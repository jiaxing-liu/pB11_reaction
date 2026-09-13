#include "fusion_beam.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using Real=long double;
constexpr double barn=1e-28, kev=1.602176634e-16;
constexpr double tolerance=1e-11;
constexpr double gaussian_extent=40.;
using Integrator=boost::math::quadrature::gauss_kronrod<double,61>;

bool put(Real value,double &out) {
    if(!std::isfinite(value) || std::abs(value)>std::numeric_limits<double>::max()) return false;
    out=static_cast<double>(value);
    return true; // IEEE underflow is allowed; completeness is a separate flag.
}

// E[1-cos(theta)] for directions weighted by exp(a*cos(theta)), a>=0.
// Avoid cancellation as a grows and the relative velocity aligns with beam.
Real one_minus_langevin(Real a) {
    if(a<1e-3L) {
        const Real a2=a*a;
        return 1-a*(1.L/3-a2/45+2*a2*a2/945-a2*a2*a2/4725);
    }
    if(a>50) return 1/a;
    return 1/a-2/std::expm1(2*a);
}

struct Kernel {
    int channel;
    Real ma,mb,mu,ea,temp,v,u,energy_scale;
    double s;
    double probability(double x) const {
        const Real a=4*static_cast<Real>(x)*s;
        const Real factor=a==0?1:-std::expm1(-a)/a;
        const Real delta=static_cast<Real>(x)-s;
        return static_cast<double>(4*x*static_cast<Real>(x)/std::sqrt(std::acos(-1.L))*
                                  std::exp(-delta*delta)*factor);
    }
    Real relative_energy(double x) const {const Real w=u*x;return mu*w*w/2;}
    double reaction_integrand(double x,int moment) const {
        const double p=probability(x);
        if(p==0) return 0;
        const Real w=u*x,er=mu*w*w/2;
        double sigma=0;
        const int status=fusion_c_cross_section(channel,static_cast<double>(er),&sigma);
        if(status!=PB11_STATUS_OK) throw status;
        Real weight=1;
        const Real angular=one_minus_langevin(2*static_cast<Real>(x)*s);
        if(moment==1) {
            const Real delta=v-w;
            weight=mb*(delta*delta+2*v*w*angular)/(2*energy_scale);
        } else if(moment==2) weight=er/energy_scale;
        else if(moment==3) {
            const Real mb_fraction=mb/(ma+mb),delta=v-mb_fraction*w;
            weight=(ma+mb)*(delta*delta+2*v*mb_fraction*w*angular)/(2*energy_scale);
        }
        const Real value=p*x*(sigma/barn)*weight;
        if(!std::isfinite(value) || value>std::numeric_limits<double>::max())
            throw PB11_STATUS_NUMERICAL_FAILURE;
        return static_cast<double>(value);
    }
};

void add_cut(std::vector<double>& cuts,double x,double lo,double hi) {
    if(std::isfinite(x) && x>lo && x<hi) cuts.push_back(x);
}
std::vector<double> make_cuts(const Kernel &k,double lo,double hi,bool nuclear) {
    std::vector<double> cuts{lo,hi};
    for(double delta:{-40.,-16.,-8.,-4.,-2.,-1.,0.,1.,2.,4.,8.,16.,40.})
        add_cut(cuts,k.s+delta,lo,hi);
    if(nuclear) {
        // Existing fit boundaries and narrow/broad resonance scales. Gaussian
        // knots above keep the cold-target peak resolved independently.
        if(k.channel==FUSION_PB11_3ALPHA) {
            for(double e:{101.,148.,195.,400.,668.,1211.,2340.,3294.,5700.})
                add_cut(cuts,static_cast<double>(std::sqrt(2*e*kev/k.mu)/k.u),lo,hi);
        } else if(k.channel==FUSION_DT_ALPHAN)
            add_cut(cuts,static_cast<double>(std::sqrt(2*530*kev/k.mu)/k.u),lo,hi);
        else if(k.channel==FUSION_DHE3_ALPHAP)
            add_cut(cuts,static_cast<double>(std::sqrt(2*900*kev/k.mu)/k.u),lo,hi);
    }
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    return cuts;
}

template<class Function>
double integrate(const Kernel &k,double lo,double hi,bool nuclear,
                 const Function &f,double *error=nullptr) {
    if(error) *error=0;
    if(hi<=lo) return 0;
    const auto cuts=make_cuts(k,lo,hi,nuclear);
    Real total=0,errors=0;
    for(std::size_t i=1;i<cuts.size();++i) {
        double err=0;
        const double value=Integrator::integrate(f,cuts[i-1],cuts[i],15,tolerance,&err);
        if(!std::isfinite(value) || !std::isfinite(err) || value<0)
            throw PB11_STATUS_NUMERICAL_FAILURE;
        total+=value; errors+=err;
    }
    double result=0;
    if(!put(total,result) || (error && !put(errors,*error))) throw PB11_STATUS_NUMERICAL_FAILURE;
    return result;
}
}

extern "C" int fusion_c_beam_maxwellian_window(int ch,double ma_in,double mb_in,
    double ea_in,double temp_in,fusion_beam_window_v1 *out) {
    if(!out) return PB11_STATUS_NULL_OUTPUT;
    *out={};
    if(!std::isfinite(ma_in) || !std::isfinite(mb_in) ||
       !std::isfinite(ea_in) || !std::isfinite(temp_in)) return PB11_STATUS_INVALID_ARGUMENT;
    if(ma_in<=0 || mb_in<=0 || ea_in<0 || temp_in<0) return PB11_STATUS_OUT_OF_RANGE;
    double emin=0,emax=0;
    const int domain_status=fusion_c_cross_section_domain(ch,&emin,&emax);
    if(domain_status) return domain_status;
    try {
        const Real ma=ma_in,mb=mb_in,ea=ea_in,temp=temp_in,mu=ma/(1+ma/mb);
        const Real v=std::sqrt(2*ea/ma);
        fusion_beam_window_v1 result{};
        if(temp==0) {
            const Real er=mu*v*v/2;
            double er_double=0,sigma=0;
            if(!put(er,er_double)) return PB11_STATUS_NUMERICAL_FAILURE;
            const int status=fusion_c_cross_section(ch,er_double,&sigma);
            if(status==PB11_STATUS_OUT_OF_RANGE) {
                result.domain_incomplete=1;
                result.unresolved_pair_probability=1;
                if(!put(v,result.unresolved_relative_speed_m_s)) return PB11_STATUS_NUMERICAL_FAILURE;
            } else if(status) return status;
            else {
                result.resolved_pair_probability=1;
                if(!put(v*sigma,result.resolved_reactivity_m3_s) ||
                   !put(ea*v*sigma,result.projectile_energy_reactivity_J_m3_s) ||
                   !put(er*v*sigma,result.relative_energy_reactivity_J_m3_s) ||
                   !put((ea-er)*v*sigma,result.cm_energy_reactivity_J_m3_s))
                    return PB11_STATUS_NUMERICAL_FAILURE;
            }
            *out=result; return PB11_STATUS_OK;
        }
        const Real u=std::sqrt(2*temp/mb);
        const double s=static_cast<double>(v/u);
        // Beyond this ratio subtracting neighboring speeds loses useful
        // quadrature resolution. The exact cold target remains available.
        if(!std::isfinite(s) || s>1e6) return PB11_STATUS_NUMERICAL_FAILURE;
        Kernel k{ch,ma,mb,mu,ea,temp,v,u,std::max({ea,temp,static_cast<Real>(emax)}),s};
        const double lower=std::max(0.,s-gaussian_extent),upper=s+gaussian_extent;
        const double data_lower=static_cast<double>(std::sqrt(2*emin/mu)/u);
        const double data_upper=static_cast<double>(std::sqrt(2*emax/mu)/u);
        const double lo=std::clamp(data_lower,lower,upper),hi=std::clamp(data_upper,lower,upper);
        const auto probability=[&](double x){return k.probability(x);};
        const auto speed=[&](double x){return x*k.probability(x);};
        result.domain_incomplete=1;
        result.resolved_pair_probability=integrate(k,lo,hi,false,probability);
        result.unresolved_pair_probability=integrate(k,lower,lo,false,probability)+
                                           integrate(k,hi,upper,false,probability);
        const double unresolved_speed=integrate(k,lower,lo,false,speed)+
                                      integrate(k,hi,upper,false,speed);
        if(!put(u*unresolved_speed,result.unresolved_relative_speed_m_s))
            return PB11_STATUS_NUMERICAL_FAILURE;
        double rate_error=0;
        const double rate=integrate(k,lo,hi,true,[&](double x){return k.reaction_integrand(x,0);},&rate_error);
        const Real scale=u*barn;
        if(!put(scale*rate,result.resolved_reactivity_m3_s) ||
           !put(scale*rate_error,result.quadrature_error_m3_s) ||
           !put(ea*scale*rate,result.projectile_energy_reactivity_J_m3_s))
            return PB11_STATUS_NUMERICAL_FAILURE;
        double *moments[]={&result.target_energy_reactivity_J_m3_s,
                           &result.relative_energy_reactivity_J_m3_s,
                           &result.cm_energy_reactivity_J_m3_s};
        for(int moment=1;moment<=3;++moment) {
            const double integral=integrate(k,lo,hi,true,[&](double x){return k.reaction_integrand(x,moment);});
            if(!put(scale*k.energy_scale*integral,*moments[moment-1])) return PB11_STATUS_NUMERICAL_FAILURE;
        }
        const Real residual=static_cast<Real>(result.projectile_energy_reactivity_J_m3_s)+
            result.target_energy_reactivity_J_m3_s-result.relative_energy_reactivity_J_m3_s-
            result.cm_energy_reactivity_J_m3_s;
        const Real escale=static_cast<Real>(result.projectile_energy_reactivity_J_m3_s)+
            result.target_energy_reactivity_J_m3_s+result.relative_energy_reactivity_J_m3_s+
            result.cm_energy_reactivity_J_m3_s;
        if(!put(residual,result.energy_identity_error_J_m3_s) || std::abs(residual)>1e-9L*escale ||
           std::abs(result.resolved_pair_probability+result.unresolved_pair_probability-1)>1e-9 ||
           (rate>0 && rate_error>1e-8*rate)) return PB11_STATUS_NUMERICAL_FAILURE;
        *out=result;
        return PB11_STATUS_OK;
    } catch(int status) {return status;}
      catch(...) {return PB11_STATUS_EXCEPTION;}
}

extern "C" int fusion_c_thermal_pair_maxwellian_window(int ch,double ma_in,
    double mb_in,double ta_in,double tb_in,fusion_beam_window_v1 *out) {
    if(!out) return PB11_STATUS_NULL_OUTPUT;
    *out={};
    if(!std::isfinite(ma_in) || !std::isfinite(mb_in) ||
       !std::isfinite(ta_in) || !std::isfinite(tb_in)) return PB11_STATUS_INVALID_ARGUMENT;
    if(ma_in<=0 || mb_in<=0 || ta_in<0 || tb_in<0) return PB11_STATUS_OUT_OF_RANGE;
    const Real ma=ma_in,mb=mb_in,ta=ta_in,tb=tb_in,m=ma+mb,mu=ma/(1+ma/mb);
    // A zero-speed test particle with this effective target temperature has
    // exactly the required relative-velocity distribution. This is an
    // integration substitution, not the physical background temperature.
    double effective_target=0;
    if(!put(tb+(mb/ma)*ta,effective_target) ||
       (effective_target==0 && (ta>0 || tb>0))) return PB11_STATUS_NUMERICAL_FAILURE;
    fusion_beam_window_v1 result{};
    const int status=fusion_c_beam_maxwellian_window(ch,ma_in,mb_in,0,effective_target,&result);
    if(status) return status;
    if(ta==0 && tb==0) {*out=result;return PB11_STATUS_OK;}
    const Real tr=(mb*ta+ma*tb)/m;
    const Real relative_moment=result.relative_energy_reactivity_J_m3_s;
    const Real random_cm=1.5L*(ta/tr)*tb*result.resolved_reactivity_m3_s;
    const Real a=ma/m*random_cm+mb/m*(ta/tr)*(ta/tr)*relative_moment;
    const Real b=mb/m*random_cm+ma/m*(tb/tr)*(tb/tr)*relative_moment;
    const Real cm=random_cm+mu/m*((ta-tb)/tr)*((ta-tb)/tr)*relative_moment;
    if(!put(a,result.projectile_energy_reactivity_J_m3_s) ||
       !put(b,result.target_energy_reactivity_J_m3_s) ||
       !put(cm,result.cm_energy_reactivity_J_m3_s)) return PB11_STATUS_NUMERICAL_FAILURE;
    const Real residual=static_cast<Real>(result.projectile_energy_reactivity_J_m3_s)+
        result.target_energy_reactivity_J_m3_s-result.relative_energy_reactivity_J_m3_s-
        result.cm_energy_reactivity_J_m3_s;
    if(!put(residual,result.energy_identity_error_J_m3_s) ||
       std::abs(residual)>1e-12L*(a+b+relative_moment+cm)) return PB11_STATUS_NUMERICAL_FAILURE;
    *out=result;
    return PB11_STATUS_OK;
}
