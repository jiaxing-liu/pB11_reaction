#include "fusion_rates.h"
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using State=std::array<double,13>; // six species, five integrated event counts, neutrons, nuclear J/m^3
constexpr double kT=50*1.602176634e-16;
constexpr double charges[6]={1,1,1,2,2,5};
constexpr double nucleons[6]={1,2,3,3,4,11};
State derivative(const State &y,unsigned mask) {
    fusion_thermal_rates_v1 rates{};
    if (fusion_c_thermal_rates(kT,y.data(),mask,PB11_REACTIVITY_FAST,&rates))
        throw std::runtime_error("closed box out of domain or negative RK stage");
    State result{};
    for (int s=0;s<6;++s) result[s]=rates.particles.net_source[s];
    for (int c=0;c<5;++c) result[6+c]=rates.event_rate_m3_s[c];
    result[11]=rates.particles.neutron_birth;
    result[12]=rates.total_nuclear_power_W_m3;
    return result;
}
State integrate(State state,unsigned mask,double duration,int steps) {
    const double dt=duration/steps;
    for (int n=0;n<steps;++n) {
        State k1=derivative(state,mask),stage{};
        for (int s=0;s<13;++s) stage[s]=state[s]+dt*k1[s]/2;
        State k2=derivative(stage,mask);
        for (int s=0;s<13;++s) stage[s]=state[s]+dt*k2[s]/2;
        State k3=derivative(stage,mask);
        for (int s=0;s<13;++s) stage[s]=state[s]+dt*k3[s];
        State k4=derivative(stage,mask);
        for (int s=0;s<13;++s) state[s]+=dt*(k1[s]+2*k2[s]+2*k3[s]+k4[s])/6;
    }
    return state;
}
bool near(double a,double b,double scale,double tol=1e-11) {
    return std::abs(a-b)<=scale*tol;
}
void require(bool x,const char *message) {if(!x) throw std::runtime_error(message);}
void invariants(const State &initial,const State &final) {
    double z0=0,z1=0,a0=0,a1=final[11], energy=0;
    for(int s=0;s<6;++s) {
        require(final[s]>=0 && std::isfinite(final[s]),"nonnegative final state");
        z0+=initial[s]*charges[s]; z1+=final[s]*charges[s];
        a0+=initial[s]*nucleons[s]; a1+=final[s]*nucleons[s];
    }
    for(int c=0;c<5;++c) {
        double q;
        require(fusion_c_channel_q(c,&q)==0,"Q query");
        energy+=final[6+c]*q;
    }
    require(near(z0,z1,z0),"charge closure");
    require(near(a0,a1,a0),"nucleon closure including escaping neutrons");
    require(near(energy,final[12],energy),"event-weighted nuclear release closure");
}
}
int main() {
    try {
        constexpr int a[]={0,1,1,1,1},b[]={5,1,1,2,3};
        for(int ch=0;ch<5;++ch) {
            State y{}; y[a[ch]]=1e20; y[b[ch]]=1e20;
            double sv;
            require(fusion_c_thermal_reactivity(ch,kT,PB11_REACTIVITY_FAST,&sv)==0,"rate query");
            const double duration=.5/(1e20*sv);
            const State end=integrate(y,1u<<ch,duration,80);
            // Equal unlike reactants: dn/dt=-sv*n^2.
            // DD: event=sv*n^2/2 and depletion=2*event, same analytic law.
            const double n_expected=1e20/1.5;
            require(near(end[a[ch]],n_expected,1e20,1e-9),"analytic fuel depletion");
            const double events=(1e20-n_expected)/(a[ch]==b[ch]?2:1);
            require(near(end[6+ch],events,1e20,1e-9),"analytic event yield");
            invariants(y,end);
        }
        State dd{}; dd[1]=1e20;
        double st,sn;
        fusion_c_thermal_reactivity(1,kT,1,&st); fusion_c_thermal_reactivity(2,kT,1,&sn);
        const double duration=.5/(1e20*(st+sn));
        const State coarse=integrate(dd,30u,duration,80);
        const State medium=integrate(dd,30u,duration,160);
        const State fine=integrate(dd,30u,duration,320);
        invariants(dd,coarse); invariants(dd,medium); invariants(dd,fine);
        double error_coarse=0,error_fine=0;
        for(int s=0;s<12;++s) {
            error_coarse=std::max(error_coarse,std::abs(coarse[s]-medium[s])/1e20);
            error_fine=std::max(error_fine,std::abs(medium[s]-fine[s])/1e20);
        }
        require(error_fine<1e-6 && error_fine<error_coarse/8,"DD secondary network RK4 refinement");
        require(fine[9]>0 && fine[10]>0 && fine[4]>0,"both DD-born secondary channels burn");
        std::cout<<"DD secondary network scaled refinement errors "<<error_coarse<<" -> "<<error_fine<<'\n';
        std::cout<<"PASS: analytic channel burn, DD secondary network, particle and nuclear-Q ledgers\n";
    } catch(const std::exception &error) {std::cerr<<"FAIL: "<<error.what()<<'\n'; return 1;}
}
