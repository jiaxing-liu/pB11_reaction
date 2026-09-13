#include "fusion_products.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {
using Real=long double;
constexpr Real light_speed=299792458.L;
bool put(Real value,double &out) {
    if(!std::isfinite(value) || std::abs(value)>std::numeric_limits<double>::max()) return false;
    out=static_cast<double>(value);
    return value==0 || out!=0;
}
bool nonnegative(double x) {return std::isfinite(x) && x>=0;}
}
extern "C" int fusion_c_map_birth_packets(int n,const double *edges,int packets,
    const double *energy,const double *rate,double *birth,fusion_birth_mapping_v1 *out) {
    if(out) *out={};
    if(n<1 || packets<0) return PB11_STATUS_INVALID_ARGUMENT;
    if(birth) std::fill(birth,birth+n,0.);
    if(!birth || !out) return PB11_STATUS_NULL_OUTPUT;
    if(!edges || (packets && (!energy || !rate))) return PB11_STATUS_INVALID_ARGUMENT;
    for(int i=0;i<=n;++i) {
        if(!nonnegative(edges[i])) return PB11_STATUS_INVALID_ARGUMENT;
        if(i && edges[i]<=edges[i-1]) return PB11_STATUS_OUT_OF_RANGE;
    }
    for(int i=0;i<packets;++i)
        if(!nonnegative(energy[i]) || !nonnegative(rate[i])) return PB11_STATUS_INVALID_ARGUMENT;
    try {
        std::vector<Real> center(n),accumulated(n,0);
        for(int i=0;i<n;++i) center[i]=(static_cast<Real>(edges[i])+edges[i+1])/2;
        Real total_n=0,total_e=0,below_n=0,below_e=0,above_n=0,above_e=0;
        for(int p=0;p<packets;++p) {
            const Real e=energy[p],w=rate[p];
            total_n+=w;total_e+=w*e;
            if(e<center.front()) {below_n+=w;below_e+=w*e;}
            else if(e>center.back()) {above_n+=w;above_e+=w*e;}
            else {
                const auto it=std::lower_bound(center.begin(),center.end(),e);
                const auto right=static_cast<std::size_t>(it-center.begin());
                if(*it==e) accumulated[right]+=w;
                else {
                    const auto left=right-1;
                    const Real weight_right=(e-center[left])/(center[right]-center[left]);
                    accumulated[right]+=w*weight_right;
                    accumulated[left]+=w*(1-weight_right);
                }
            }
        }
        std::vector<double> result(n);
        Real mapped_n=0,mapped_e=0;
        for(int i=0;i<n;++i) {
            if(!std::isfinite(accumulated[i]) || accumulated[i]>std::numeric_limits<double>::max())
                return PB11_STATUS_NUMERICAL_FAILURE;
            result[i]=static_cast<double>(accumulated[i]);
            mapped_n+=result[i];mapped_e+=center[i]*result[i];
        }
        const Real nr=total_n-mapped_n-below_n-above_n,er=total_e-mapped_e-below_e-above_e;
        if(std::abs(nr)>1e-12L*(total_n+mapped_n+below_n+above_n) ||
           std::abs(er)>1e-12L*(total_e+mapped_e+below_e+above_e)) return PB11_STATUS_NUMERICAL_FAILURE;
        fusion_birth_mapping_v1 ledger{};
        if(!put(total_n,ledger.input_number_m3_s) || !put(total_e,ledger.input_energy_W_m3) ||
           !put(mapped_n,ledger.mapped_number_m3_s) || !put(mapped_e,ledger.mapped_energy_W_m3) ||
           !put(below_n,ledger.below_number_m3_s) || !put(below_e,ledger.below_energy_W_m3) ||
           !put(above_n,ledger.above_number_m3_s) || !put(above_e,ledger.above_energy_W_m3) ||
           !put(nr,ledger.number_residual_m3_s) || !put(er,ledger.energy_residual_W_m3))
            return PB11_STATUS_NUMERICAL_FAILURE;
        std::copy(result.begin(),result.end(),birth);*out=ledger;
        return PB11_STATUS_OK;
    } catch(...) {return PB11_STATUS_EXCEPTION;}
}
extern "C" int fusion_c_three_equal_sequential_cm(double mass_in,double available_in,
    double intermediate_in,double cosine_in,fusion_three_body_cm_v1 *out) {
    if(!out) return PB11_STATUS_NULL_OUTPUT;
    *out={};
    if(!std::isfinite(mass_in) || !std::isfinite(available_in) ||
       !std::isfinite(intermediate_in) || !std::isfinite(cosine_in)) return PB11_STATUS_INVALID_ARGUMENT;
    if(mass_in<=0 || available_in<0 || intermediate_in<0 ||
       intermediate_in>available_in || std::abs(cosine_in)>1) return PB11_STATUS_OUT_OF_RANGE;
    const Real m=mass_in,c2=light_speed*light_speed,rest=m*c2;
    const Real a=available_in,q=intermediate_in,cosine=cosine_in;
    const Real intermediate_rest=2*rest+q,parent_rest=3*rest+a,first_release=a-q;
    const Real t0=first_release*(first_release+2*intermediate_rest)/(2*parent_rest);
    const Real recoil=first_release*(first_release+2*rest)/(2*parent_rest);
    const Real primary_pc=std::sqrt(t0*(t0+2*rest));
    const Real secondary_pc=std::sqrt((q/2)*(q/2+2*rest));
    const Real gamma=1+recoil/intermediate_rest;
        const Real px=secondary_pc*std::sqrt((1-cosine)*(1+cosine))/light_speed;
    const Real pz0=primary_pc/light_speed,pzshift=gamma*secondary_pc*cosine/light_speed;
    const Real momenta_x[3]={0,px,-px};
    const Real momenta_z[3]={pz0,-pz0/2-pzshift,-pz0/2+pzshift};
    fusion_three_body_cm_v1 result{};
    for(int i=0;i<3;++i) {
        // Recover kinetic energy from momentum without subtracting rest
        // energy or two almost equal boosted secondary energies. This stays
        // nonnegative even when a secondary is exactly at rest in the CM.
        const Real pcx=momenta_x[i]*light_speed,pcz=momenta_z[i]*light_speed;
        const Real pc2=pcx*pcx+pcz*pcz;
        const Real kinetic=pc2/(std::sqrt(rest*rest+pc2)+rest);
        if(!put(kinetic,result.kinetic_energy_J[i]) ||
           !put(momenta_x[i],result.momentum_x_kg_m_s[i]) ||
           !put(momenta_z[i],result.momentum_z_kg_m_s[i])) return PB11_STATUS_NUMERICAL_FAILURE;
    }
    const Real er=static_cast<Real>(result.kinetic_energy_J[0])+result.kinetic_energy_J[1]+
        result.kinetic_energy_J[2]-a;
    const Real pr=std::abs(static_cast<Real>(result.momentum_x_kg_m_s[0])+result.momentum_x_kg_m_s[1]+
        result.momentum_x_kg_m_s[2])+std::abs(static_cast<Real>(result.momentum_z_kg_m_s[0])+
        result.momentum_z_kg_m_s[1]+result.momentum_z_kg_m_s[2]);
    if(!put(er,result.energy_residual_J) || !put(pr,result.momentum_residual_kg_m_s) ||
       std::abs(er)>1e-12L*a || pr>1e-12L*(std::abs(pz0)+2*std::abs(px)+2*std::abs(pzshift)))
        return PB11_STATUS_NUMERICAL_FAILURE;
    *out=result;return PB11_STATUS_OK;
}
