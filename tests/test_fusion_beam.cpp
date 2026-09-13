#include "fusion_beam.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
constexpr double kev=1.602176634e-16,amu=1.66053906892e-27;
const double pi=std::acos(-1.);
using GK=boost::math::quadrature::gauss_kronrod<double,61>;
void require(bool condition,const char *label) {if(!condition) throw std::runtime_error(label);}
bool near(double a,double b,double r=1e-9) {
    return std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=r*std::max(std::abs(a),std::abs(b));
}
fusion_beam_window_v1 beam(int ch,double ma,double mb,double e,double t) {
    fusion_beam_window_v1 r{};
    const int status=fusion_c_beam_maxwellian_window(ch,ma,mb,e,t,&r);
    if(status) {std::cerr<<"beam failure channel="<<ch<<" E/keV="<<e/kev<<" T/keV="<<t/kev<<" status="<<status<<'\n';}
    require(status==0,"beam integral failed");
    require(std::abs(r.resolved_pair_probability+r.unresolved_pair_probability-1)<1e-10,"relative probability normalization");
    const double total=r.projectile_energy_reactivity_J_m3_s+r.target_energy_reactivity_J_m3_s;
    require(near(total,r.relative_energy_reactivity_J_m3_s+r.cm_energy_reactivity_J_m3_s,2e-9),"reactant/CM/relative energy identity");
    return r;
}

// Independent integration over target speed and its angle to the projectile.
// Does NOT use the analytic relative-speed distribution or Langevin function.
double direct_velocity_integral(int ch,double ma,double mb,double ea,double temp,int moment) {
    const double v=std::sqrt(2*ea/ma),u=std::sqrt(2*temp/mb),mu=ma*mb/(ma+mb);
    double emin=0,emax=0;
    require(fusion_c_cross_section_domain(ch,&emin,&emax)==0,"reference domain");
    const auto radial=[&](double y) {
        const double vb=u*y;
        const auto angular=[&](double costheta) {
            const double w2=v*v+vb*vb-2*v*vb*costheta,er=.5*mu*w2;
            if(er<emin || er>emax) return 0.;
            double sigma=0;
            require(fusion_c_cross_section(ch,er,&sigma)==0,"reference cross section");
            double weight=1;
            if(moment==1) weight=temp*y*y/kev;
            if(moment==2) weight=(ma*ma*v*v+mb*mb*vb*vb+2*ma*mb*v*vb*costheta)/(2*(ma+mb)*kev);
            return .5*(sigma/1e-28)*std::sqrt(w2)*weight;
        };
        return 4/std::sqrt(pi)*y*y*std::exp(-y*y)*GK::integrate(angular,-1.,1.,12,1e-10);
    };
    double sum=0;
    for(int i=0;i<20;++i) sum+=GK::integrate(radial,.5*i,.5*(i+1),12,1e-9);
    return sum*1e-28*(moment?kev:1.);
}

void cold_and_domain_limits() {
    for(int ch=0;ch<5;++ch) {
        const double ma=(ch==0?1:2)*amu,mb=(ch==0?11:ch<3?2:3)*amu;
        const double mu=ma*mb/(ma+mb),er=(ch==0?148:30)*kev,ea=er*ma/mu;
        const auto r=beam(ch,ma,mb,ea,0);
        double sigma=0; require(fusion_c_cross_section(ch,er,&sigma)==0,"cold cross section");
        require(near(r.resolved_reactivity_m3_s,sigma*std::sqrt(2*ea/ma),2e-13),"cold target sigma*v with relative energy");
        require(r.target_energy_reactivity_J_m3_s==0 && r.domain_incomplete==0,"cold target energy and complete domain");
        require(near(r.relative_energy_reactivity_J_m3_s/r.resolved_reactivity_m3_s,er),"cold relative energy moment");
        const auto zero=beam(ch,ma,mb,0,0);
        require(zero.resolved_reactivity_m3_s==0 && zero.domain_incomplete==0,"zero relative speed special case");
        double lo=0,hi=0; require(fusion_c_cross_section_domain(ch,&lo,&hi)==0,"cold high bound");
        const auto outside=beam(ch,ma,mb,2*hi*ma/mu,0);
        require(outside.resolved_reactivity_m3_s==0 && outside.unresolved_pair_probability==1 &&
                outside.domain_incomplete==1,"cold out-of-window result is explicitly unresolved");
    }
    // All relative speeds below DD's data window: analytic Maxwell speed mean.
    const double temp=1e-6*kev;
    const auto below=beam(FUSION_DD_TP,2*amu,2*amu,0,temp);
    require(below.resolved_reactivity_m3_s==0 && near(below.unresolved_pair_probability,1),"low-energy data gap visible");
    require(near(below.unresolved_relative_speed_m_s,2*std::sqrt(2*temp/(2*amu))/std::sqrt(pi),2e-12),"zero-beam analytic unresolved speed moment");
    // Numerical domain diagnostic, not a claim of nonrelativistic accuracy at
    // this artificial 20 MeV beam energy. All speeds lie above the DT fit.
    const double ea=20000*kev,t=0.001*kev,ma=2*amu,mb=3*amu;
    const auto above=beam(FUSION_DT_ALPHAN,ma,mb,ea,t);
    const double u=std::sqrt(2*t/mb),v=std::sqrt(2*ea/ma),s=v/u;
    const double exact=u*(std::exp(-s*s)/std::sqrt(pi)+(s+1/(2*s))*std::erf(s));
    require(above.resolved_reactivity_m3_s==0 && near(above.unresolved_pair_probability,1),"high-energy data gap visible");
    require(near(above.unresolved_relative_speed_m_s,exact,2e-11),"drifting Maxwell analytic speed moment");
}

void intermediate_speed_diagnostics() {
    const double ma=2*amu,mb=2*amu,temp=1e-6*kev,u=std::sqrt(2*temp/mb);
    for(double speed_ratio:{1e-6,.1,1.,3.,10.}) {
        const double ea=.5*ma*u*u*speed_ratio*speed_ratio;
        const auto r=beam(FUSION_DD_TP,ma,mb,ea,temp);
        const double expected=u*(std::exp(-speed_ratio*speed_ratio)/std::sqrt(pi)+
            (speed_ratio+1/(2*speed_ratio))*std::erf(speed_ratio));
        require(r.resolved_reactivity_m3_s==0 && near(r.unresolved_pair_probability,1),"full missing window for speed diagnostic");
        require(near(r.unresolved_relative_speed_m_s,expected,3e-12),"analytic mean speed across small/intermediate/large drift");
    }
    const double ma_dt=2*amu,mb_dt=3*amu,ea=50*kev;
    const auto cold=beam(FUSION_DT_ALPHAN,ma_dt,mb_dt,ea,0);
    double previous=1.;
    for(double t:{.1,.01,.001,.0001}) {
        const auto warm=beam(FUSION_DT_ALPHAN,ma_dt,mb_dt,ea,t*kev);
        const double error=std::abs(warm.resolved_reactivity_m3_s/cold.resolved_reactivity_m3_s-1);
        std::cout<<"DT cold-limit T/keV "<<t<<" relative error "<<error<<'\n';
        require(error<previous,"positive-temperature integration approaches cold target");
        previous=error;
    }
    require(previous<1e-4,"cold-limit convergence magnitude");
}

void independent_angular_check() {
    for(int ch:{FUSION_DT_ALPHAN,FUSION_PB11_3ALPHA}) {
        const double ma=(ch==0?1:2)*amu,mb=(ch==0?11:3)*amu;
        const double ea=(ch==0?100:20)*kev,temp=10*kev;
        const auto r=beam(ch,ma,mb,ea,temp);
        const double rate=direct_velocity_integral(ch,ma,mb,ea,temp,0);
        const double target=direct_velocity_integral(ch,ma,mb,ea,temp,1);
        const double cm=direct_velocity_integral(ch,ma,mb,ea,temp,2);
        std::cout<<"channel "<<ch<<" independent target-angle rate error "<<std::abs(r.resolved_reactivity_m3_s/rate-1)
                 <<", target-energy error "<<std::abs(r.target_energy_reactivity_J_m3_s/target-1)<<'\n';
        require(near(r.resolved_reactivity_m3_s,rate,2e-7),"analytic angular reduction agrees with direct velocities");
        require(near(r.target_energy_reactivity_J_m3_s,target,2e-7),"reaction-conditioned target energy agrees independently");
        require(near(r.cm_energy_reactivity_J_m3_s,cm,2e-7),"reaction-conditioned CM energy agrees independently");
        require(r.domain_incomplete==1,"Maxwell tail always leaves finite data domain incomplete");
        require(r.target_energy_reactivity_J_m3_s/r.resolved_reactivity_m3_s>1.6*temp,"reaction-selected target energy differs from unweighted 3T/2");
    }
}

double relative_thermal_reference(int ch,double mu,double temp) {
    double emin=0,emax=0; require(fusion_c_cross_section_domain(ch,&emin,&emax)==0,"thermal domain");
    const double lo=emin/temp,hi=emax/temp;
    std::vector<double> cuts{lo,hi};
    for(double e:{.5,1.,2.,4.,8.,16.,32.,64.}) if(e>lo && e<hi) cuts.push_back(e);
    for(double e:{101.,148.,195.,400.,530.,668.,900.,1211.,2340.,3294.,5700.})
        if(e*kev/temp>lo && e*kev/temp<hi) cuts.push_back(e*kev/temp);
    std::sort(cuts.begin(),cuts.end());
    double total=0;
    const auto kernel=[&](double y) {
        double sigma=0; require(fusion_c_cross_section(ch,temp*y,&sigma)==0,"thermal cross section");
        return y*(sigma/1e-28)*std::exp(-y);
    };
    for(std::size_t i=1;i<cuts.size();++i) total+=GK::integrate(kernel,cuts[i-1],cuts[i],15,1e-11);
    return std::sqrt(8*temp/(pi*mu))*total*1e-28;
}

void thermal_average_check() {
    for(int ch:{FUSION_DT_ALPHAN,FUSION_PB11_3ALPHA}) {
        const double ma=(ch==0?1:2)*amu,mb=(ch==0?11:3)*amu,temp=(ch==0?50:10)*kev;
        const double mu=ma*mb/(ma+mb);
        const auto kernel=[&](double y) {
            const auto r=beam(ch,ma,mb,y*temp,temp);
            return 2/std::sqrt(pi)*std::sqrt(y)*std::exp(-y)*r.resolved_reactivity_m3_s/1e-24;
        };
        double total=0;
        for(const auto interval:std::vector<std::pair<double,double>>{{0,1},{1,2},{2,4},{4,8},{8,16},{16,32},{32,64},{64,80}})
            total+=GK::integrate(kernel,interval.first,interval.second,12,2e-9);
        total*=1e-24;
        const double reference=relative_thermal_reference(ch,mu,temp);
        std::cout<<"channel "<<ch<<" double-Maxwellian versus relative-Maxwellian error "<<std::abs(total/reference-1)<<'\n';
        require(near(total,reference,2e-7),"thermal averaging of beam kernel recovers relative Maxwellian rate");
        const auto stationary=beam(ch,ma,mb,0,temp);
        require(near(stationary.target_energy_reactivity_J_m3_s,
                     mb/mu*stationary.relative_energy_reactivity_J_m3_s,2e-10),"stationary projectile energy mapping");
    }
}

void unequal_temperature_energy_check() {
    const int ch=FUSION_DT_ALPHAN;
    const double ma=2*amu,mb=3*amu,ta=5*kev,tb=12*kev;
    fusion_beam_window_v1 result{};
    require(fusion_c_thermal_pair_maxwellian_window(ch,ma,mb,ta,tb,&result)==0,"unequal-temperature pair status");
    double totals[3]={};
    for(int moment=0;moment<3;++moment) {
        const auto kernel=[&](double y) {
            const auto r=beam(ch,ma,mb,y*ta,tb);
            const double value=moment==0?r.resolved_reactivity_m3_s/1e-24:
                (moment==1?r.projectile_energy_reactivity_J_m3_s:r.target_energy_reactivity_J_m3_s)/(kev*1e-24);
            return 2/std::sqrt(pi)*std::sqrt(y)*std::exp(-y)*value;
        };
        for(const auto interval:std::vector<std::pair<double,double>>{{0,1},{1,2},{2,4},{4,8},{8,16},{16,32},{32,64},{64,80}})
            totals[moment]+=GK::integrate(kernel,interval.first,interval.second,12,2e-9);
        totals[moment]*=1e-24*(moment?kev:1.);
    }
    std::cout<<"unequal-temperature direct average rate error "<<std::abs(result.resolved_reactivity_m3_s/totals[0]-1)
             <<", a-energy error "<<std::abs(result.projectile_energy_reactivity_J_m3_s/totals[1]-1)
             <<", b-energy error "<<std::abs(result.target_energy_reactivity_J_m3_s/totals[2]-1)<<'\n';
    require(near(result.resolved_reactivity_m3_s,totals[0],2e-7),"unequal-temperature relative effective temperature");
    require(near(result.projectile_energy_reactivity_J_m3_s,totals[1],2e-7),"correlated reacting a energy");
    require(near(result.target_energy_reactivity_J_m3_s,totals[2],2e-7),"correlated reacting b energy");
    fusion_beam_window_v1 swapped{};
    require(fusion_c_thermal_pair_maxwellian_window(ch,mb,ma,tb,ta,&swapped)==0,"swapped pair status");
    require(near(swapped.resolved_reactivity_m3_s,result.resolved_reactivity_m3_s,2e-11),"species exchange leaves pair rate invariant");
    require(near(swapped.target_energy_reactivity_J_m3_s,result.projectile_energy_reactivity_J_m3_s,2e-11),"species exchange swaps energy debit");
    fusion_beam_window_v1 equal{};
    require(fusion_c_thermal_pair_maxwellian_window(ch,ma,mb,ta,ta,&equal)==0,"equal-temperature pair status");
    require(near(equal.cm_energy_reactivity_J_m3_s,1.5*ta*equal.resolved_reactivity_m3_s,2e-12),"equal-temperature CM is independent of reaction selection");
    fusion_beam_window_v1 stationary{};
    require(fusion_c_thermal_pair_maxwellian_window(ch,ma,mb,0,tb,&stationary)==0,"stationary first reactant status");
    const auto reference=beam(ch,ma,mb,0,tb);
    require(stationary.projectile_energy_reactivity_J_m3_s==0 &&
            near(stationary.target_energy_reactivity_J_m3_s,reference.target_energy_reactivity_J_m3_s,2e-11),"one cold reactant limit");
    require(fusion_c_thermal_pair_maxwellian_window(ch,ma,mb,-ta,tb,&equal)==PB11_STATUS_OUT_OF_RANGE &&
            equal.resolved_reactivity_m3_s==0,"invalid pair temperature clears output");
}

void invalid_inputs() {
    fusion_beam_window_v1 out{}; out.resolved_reactivity_m3_s=1;out.domain_incomplete=1;
    require(fusion_c_beam_maxwellian_window(0,amu,11*amu,-kev,kev,&out)==PB11_STATUS_OUT_OF_RANGE,
            "negative beam energy rejected");
    require(out.resolved_reactivity_m3_s==0 && out.domain_incomplete==0,"error clears output");
    require(fusion_c_beam_maxwellian_window(5,amu,11*amu,kev,kev,&out)==PB11_STATUS_INVALID_ARGUMENT,"invalid channel rejected");
    require(fusion_c_beam_maxwellian_window(0,0,11*amu,kev,kev,&out)==PB11_STATUS_OUT_OF_RANGE,"zero mass rejected");
    require(fusion_c_beam_maxwellian_window(0,amu,11*amu,kev,std::numeric_limits<double>::quiet_NaN(),&out)==PB11_STATUS_INVALID_ARGUMENT,"NaN rejected");
    require(fusion_c_beam_maxwellian_window(0,amu,11*amu,kev,kev,nullptr)==PB11_STATUS_NULL_OUTPUT,"null output rejected");
}
}
int main() {
    try {cold_and_domain_limits();intermediate_speed_diagnostics();independent_angular_check();thermal_average_check();unequal_temperature_energy_check();invalid_inputs();}
    catch(const std::exception &error) {std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}
    std::cout<<"PASS: beam kinematics, explicit data gaps, reaction-conditioned energy and Maxwellian limits\n";
}
