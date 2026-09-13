#include "fusion_laboratory.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
using R=long double;
constexpr R c=299792458.L,c2=c*c;
bool put(R x,double &y) {
    if(!std::isfinite(x) || std::abs(x)>std::numeric_limits<double>::max()) return false;
    y=static_cast<double>(x); return x==0 || y!=0;
}
R kinetic(R rest,const R *p) {
    R pc2=0; for(int j=0;j<3;++j) pc2+=p[j]*p[j]*c2;
    return pc2/(std::sqrt(rest*rest+pc2)+rest);
}
bool write(double mass,const R *p,fusion_particle_four_vector_v1 &o) {
    o.mass_kg=mass;
    for(int j=0;j<3;++j) if(!put(p[j],o.momentum_kg_m_s[j])) return false;
    // Energy is computed from the actual rounded momentum we return.
    R rounded[3]; for(int j=0;j<3;++j) rounded[j]=o.momentum_kg_m_s[j];
    return put(kinetic(static_cast<R>(mass)*c2,rounded),o.kinetic_energy_J);
}
}
extern "C" int fusion_c_two_body_cm(double ma,double mb,double a,const double *d,
    fusion_particle_four_vector_v1 *out) {
    if(!out) return PB11_STATUS_NULL_OUTPUT;
    std::fill(out,out+2,fusion_particle_four_vector_v1{});
    if(!d || !std::isfinite(ma) || !std::isfinite(mb) || !std::isfinite(a)) return PB11_STATUS_INVALID_ARGUMENT;
    if(ma<=0 || mb<=0 || a<0) return PB11_STATUS_OUT_OF_RANGE;
    R norm2=0; for(int j=0;j<3;++j) {if(!std::isfinite(d[j])) return PB11_STATUS_INVALID_ARGUMENT;norm2+=static_cast<R>(d[j])*d[j];}
    if(std::abs(norm2-1)>2e-12L) return PB11_STATUS_OUT_OF_RANGE;
    const R ra=static_cast<R>(ma)*c2,rb=static_cast<R>(mb)*c2,A=a;
    const R t=A*(A+2*rb)/(2*(ra+rb+A));
    const R magnitude=std::sqrt(t*(t+2*ra))/c, norm=std::sqrt(norm2);
    R p[3],q[3]; for(int j=0;j<3;++j) {p[j]=magnitude*d[j]/norm;q[j]=-p[j];}
    fusion_particle_four_vector_v1 result[2]{};
    if(!write(ma,p,result[0]) || !write(mb,q,result[1])) return PB11_STATUS_NUMERICAL_FAILURE;
    if(std::abs(static_cast<R>(result[0].kinetic_energy_J)+result[1].kinetic_energy_J-A)>1e-12L*A)
        return PB11_STATUS_NUMERICAL_FAILURE;
    std::copy(result,result+2,out);return PB11_STATUS_OK;
}
extern "C" int fusion_c_boost_particles(int n,const double *v,
    const fusion_particle_four_vector_v1 *in,fusion_particle_four_vector_v1 *out,
    fusion_boost_ledger_v1 *ledger) {
    if(ledger) *ledger={};
    if(n<1) return PB11_STATUS_INVALID_ARGUMENT;
    if(out) std::fill(out,out+n,fusion_particle_four_vector_v1{});
    if(!out || !ledger) return PB11_STATUS_NULL_OUTPUT;
    if(!v || !in) return PB11_STATUS_INVALID_ARGUMENT;
    R v2=0; for(int j=0;j<3;++j) {if(!std::isfinite(v[j])) return PB11_STATUS_INVALID_ARGUMENT;v2+=static_cast<R>(v[j])*v[j];}
    const R b2=v2/c2;
    if(b2>=1) return PB11_STATUS_OUT_OF_RANGE;
    const R gamma=1/std::sqrt(1-b2),gm1=gamma*gamma*b2/(gamma+1);
    // (gamma-1)/v^2 = gamma^2/((gamma+1)c^2), including V=0.
    const R boost_factor=gamma*gamma/((gamma+1)*c2);
    try {
        std::vector<fusion_particle_four_vector_v1> result(n);
        R initial_t=0,rest_sum=0,final_t=0,initial_p[3]={},final_p[3]={};
        for(int i=0;i<n;++i) {
            const auto &x=in[i];
            if(!std::isfinite(x.mass_kg) || !std::isfinite(x.kinetic_energy_J)) return PB11_STATUS_INVALID_ARGUMENT;
            if(x.mass_kg<=0 || x.kinetic_energy_J<0) return PB11_STATUS_OUT_OF_RANGE;
            R p[3],dot=0; for(int j=0;j<3;++j) {
                if(!std::isfinite(x.momentum_kg_m_s[j])) return PB11_STATUS_INVALID_ARGUMENT;
                p[j]=x.momentum_kg_m_s[j];dot+=v[j]*p[j];
            }
            const R rest=static_cast<R>(x.mass_kg)*c2,t=x.kinetic_energy_J,on_shell=kinetic(rest,p);
            if(std::abs(on_shell-t)>1e-10L*std::max(on_shell,t)) return PB11_STATUS_OUT_OF_RANGE;
            const R factor=boost_factor*dot+gamma*(rest+t)/c2;
            R shifted[3];for(int j=0;j<3;++j) shifted[j]=p[j]+factor*v[j];
            if(!write(x.mass_kg,shifted,result[i])) return PB11_STATUS_NUMERICAL_FAILURE;
            initial_t+=t;rest_sum+=rest;final_t+=result[i].kinetic_energy_J;
            for(int j=0;j<3;++j) {initial_p[j]+=p[j];final_p[j]+=result[i].momentum_kg_m_s[j];}
        }
        R dot=0;for(int j=0;j<3;++j) dot+=v[j]*initial_p[j];
        const R expected=initial_t+gm1*(rest_sum+initial_t)+gamma*dot;
        const R er=final_t-expected, factor=boost_factor*dot+gamma*(rest_sum+initial_t)/c2;
        R pr2=0,pscale=0;for(int j=0;j<3;++j) {
            const R predicted=initial_p[j]+factor*v[j],delta=final_p[j]-predicted;
            pr2+=delta*delta;pscale+=std::abs(initial_p[j])+std::abs(factor*v[j])+std::abs(final_p[j]);
        }
        const R escale=initial_t+gm1*(rest_sum+initial_t)+std::abs(gamma*dot)+final_t;
        const R pr=std::sqrt(pr2);
        if(std::abs(er)>1e-10L*escale || pr>1e-12L*pscale) return PB11_STATUS_NUMERICAL_FAILURE;
        fusion_boost_ledger_v1 l{};
        if(!put(expected,l.expected_kinetic_energy_J) || !put(final_t,l.output_kinetic_energy_J) ||
           !put(er,l.energy_residual_J) || !put(pr,l.momentum_residual_kg_m_s)) return PB11_STATUS_NUMERICAL_FAILURE;
        std::copy(result.begin(),result.end(),out);*ledger=l;return PB11_STATUS_OK;
    } catch(...) {return PB11_STATUS_EXCEPTION;}
}
