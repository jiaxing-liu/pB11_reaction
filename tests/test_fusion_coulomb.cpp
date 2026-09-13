#include "fusion_coulomb.h"
#include "fusion_kinetics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
constexpr double keV=1.602176634e-16,ma=6.6446573450e-27;
constexpr double eps0=8.8541878188e-12,e=1.602176634e-19;
void require(bool value,const char *what) {if(!value)throw std::runtime_error(what);}
bool near(double a,double b,double relative=1e-11) {
    return std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=relative*std::max(std::abs(a),std::abs(b));
}
fusion_coulomb_energy_v1 coeff(double energy,const fusion_maxwellian_bath_v1 &bath) {
    fusion_coulomb_energy_v1 r{};
    require(fusion_c_coulomb_energy(energy,ma,2,&bath,&r)==0,"coefficient status");return r;
}
}
int main() {
    try {
        fusion_maxwellian_bath_v1 bath{1e20,3.3435837768e-27,1,keV,15};
        const double pi=std::acos(-1.),coulomb=e*e/(4*pi*eps0);
        for(double y:{.1,.5,1.,3.,10.}) {
            const double velocity=y*std::sqrt(2*bath.kT_J/bath.mass_kg);
            const double energy=ma*velocity*velocity/2,x=y*y;
            const double psi=std::erf(std::sqrt(x))-2*std::sqrt(x)*std::exp(-x)/std::sqrt(pi);
            const double psi_prime=2*std::sqrt(x)*std::exp(-x)/std::sqrt(pi);
            // Independent NRL2019 p31 relaxation-frequency route (CGS -> SI).
            const double nu0=4*pi*4*coulomb*coulomb*bath.coulomb_log*bath.density_m3/
                (ma*ma*velocity*velocity*velocity);
            const double nu_parallel=psi/x*nu0;
            const double nu_energy=2*((ma/bath.mass_kg)*psi-psi_prime)*nu0;
            const auto r=coeff(energy,bath);
            require(near(r.diffusion_J2_s,.5*ma*ma*std::pow(velocity,4)*nu_parallel),"parallel variance factor two and SI normalization");
            require(near(r.mean_energy_rate_J_s,-energy*nu_energy),"NRL energy relaxation rate");
            const double h=energy*1e-4;
            const auto plus=coeff(energy+h,bath),minus=coeff(energy-h,bath);
            const double from_flux=(plus.diffusion_J2_s-minus.diffusion_J2_s)/(2*h)-
                (1/bath.kT_J-1/(2*energy))*r.diffusion_J2_s;
            require(near(from_flux,r.mean_energy_rate_J_s,2e-7),"FP flux coefficient and Ito mean drift agree");
        }
        const auto zero=coeff(0,bath),tiny=coeff(1e-12*keV,bath);
        require(zero.diffusion_J2_s==0 && zero.mean_energy_rate_J_s>0,"zero-energy diffusion/heating limits");
        require(near(tiny.diffusion_J2_s/(1e-12*keV),2*zero.mean_energy_rate_J_s/3,1e-10),"low-speed series avoids cancellation");
        const auto ordinary=coeff(20*keV,bath);
        bath.density_m3*=2; const auto doubled=coeff(20*keV,bath); bath.density_m3/=2;
        require(near(doubled.diffusion_J2_s,2*ordinary.diffusion_J2_s),"bath density scaling");
        require(near(doubled.mean_energy_rate_J_s,2*ordinary.mean_energy_rate_J_s),"bath drift density scaling");
        fusion_coulomb_energy_v1 r{1,1};
        require(fusion_c_coulomb_energy(-keV,ma,2,&bath,&r)==PB11_STATUS_OUT_OF_RANGE && r.diffusion_J2_s==0 && r.mean_energy_rate_J_s==0,"negative-energy error clears output");
        require(fusion_c_coulomb_energy(keV,ma,2,nullptr,&r)==PB11_STATUS_INVALID_ARGUMENT,"null bath");
        require(fusion_c_coulomb_energy(keV,ma,2,&bath,nullptr)==PB11_STATUS_NULL_OUTPUT,"null coefficients");
        require(fusion_c_coulomb_energy(keV,ma,0,&bath,&r)==0 && r.diffusion_J2_s==0,"neutral test particle");

        // Coupled operator first moment: a narrow shell at 20 keV must change
        // mean energy at the independently computed single-particle NRL rate.
        for (int n:{101,201,401}) {
            std::vector<double> edges(n+1),old(n,0),trial(n),d(n-1),birth(n,0),escape(n,0);
            for (int i=0;i<=n;++i) edges[i]=40*keV*i/n;
            old[n/2]=1e12;
            for (int i=0;i<n-1;++i) d[i]=coeff(edges[i+1],bath).diffusion_J2_s;
            double heat=0; fusion_kinetic_ledger_v1 ledger{};
            constexpr double dt=1e-7;
            require(fusion_c_energy_fp_trial(n,1,dt,edges.data(),old.data(),&bath.kT_J,d.data(),birth.data(),escape.data(),0,trial.data(),&heat,&ledger)==0,"Coulomb FP trial status");
            const double mean_rate=(ledger.final_energy_J_m3-ledger.initial_energy_J_m3)/(dt*1e12);
            const double error=std::abs(mean_rate/ordinary.mean_energy_rate_J_s-1);
            std::cout<<"Coulomb first-moment relative error at "<<n<<" cells: "<<error<<'\n';
            require(error<.003,"discrete energy drift agrees with NRL mean rate");
            require(heat>0 && near(ledger.initial_energy_J_m3-ledger.final_energy_J_m3,heat,1e-9),"collision energy deposited exactly once");
        }
        std::cout<<"PASS: Coulomb normalization, limits, drift-diffusion identity and kinetic energy exchange\n";
    } catch(const std::exception &error) {std::cerr<<"FAIL: "<<error.what()<<'\n'; return 1;}
}
