#include "fusion_rate_model.h"
#include "fusion_nuclear_data.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
constexpr double kev=1.602176634e-16;const double pi=std::acos(-1.);
using GK=boost::math::quadrature::gauss_kronrod<double,61>;
void require(bool x,const char*s){if(!x)throw std::runtime_error(s);}
double direct_velocity_integral(int ch,double ma,double mb,double ea,double temp,int moment) {
    const double v=std::sqrt(2*ea/ma),u=std::sqrt(2*temp/mb),mu=ma*mb/(ma+mb);
    double emin=0,emax=0;
    require(fusion_c_cross_section_domain(ch,&emin,&emax)==0,"reference domain");
    const auto radial=[&](double y) {
        const double vb=u*y;
        const auto angular=[&](double costheta) {
            const double w2=v*v+vb*vb-2*v*vb*costheta,er=.5*mu*w2;
            double sigma=0;
            require(fusion_c_cross_section_model(ch,FUSION_ENDPOINT_S,FUSION_PB_LOW_TB,er,&sigma)==0,"reference cross section");
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


int main(){
std::cout<<std::setprecision(17)<<"channel,beam_keV,target_keV,rate_relerr,target_energy_relerr,cm_energy_relerr\n";
for(int ch:{0,1,3,4}){
 fusion_nuclear_channel_v1 c{};fusion_nuclear_mass_v1 a{},b{};
 require(!fusion_c_nuclear_channel(ch,&c),"channel");
 require(!fusion_c_nuclear_mass(c.reactant_ids[0],&a),"mass a");
 require(!fusion_c_nuclear_mass(c.reactant_ids[1],&b),"mass b");
 double ea=(ch==1?.5:ch==3?10000:100)*kev,t=(ch==1?.01:ch==3?100:3)*kev;
 fusion_rate_model_v1 r{};require(!fusion_c_beam_maxwellian_model(ch,1,0,a.mass_kg,b.mass_kg,ea,t,&r),"model");
 double values[]={r.total.resolved_reactivity_m3_s,r.total.target_energy_reactivity_J_m3_s,r.total.cm_energy_reactivity_J_m3_s};
 std::cout<<ch<<','<<ea/kev<<','<<t/kev;
 for(int m=0;m<3;++m){double ref=direct_velocity_integral(ch,a.mass_kg,b.mass_kg,ea,t,m);double error=std::abs(values[m]-ref)/ref;std::cout<<','<<error;require(error<2e-7,"direct velocity reference");}
 std::cout<<'\n';
}
}
