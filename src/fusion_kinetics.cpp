#include "fusion_kinetics.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {
using Real=long double;
bool finite(Real value) {return std::isfinite(value);}
bool output_value(Real value,double &out) {
    if (!finite(value) || std::abs(value)>std::numeric_limits<double>::max()) return false;
    out=static_cast<double>(value);
    return value==0 || out!=0;
}
void bernoulli_pair(Real w,Real &positive,Real &negative) {
    // B(w)=w/(exp(w)-1); B(-w)=B(w)+w. No overflowing exp(+large w).
    if (std::abs(w)<1e-5L) {
        const Real w2=w*w;
        positive=1-w/2+w2/12-w2*w2/720;
        negative=positive+w;
    } else if (w>0) {
        const Real e=std::exp(-w),den=-std::expm1(-w);
        positive=w*e/den; negative=w/den;
    } else {
        const Real e=std::exp(w),den=-std::expm1(w);
        positive=-w/den; negative=-w*e/den;
    }
}
bool nonnegative(double value) {return std::isfinite(value) && value>=0;}
}

extern "C" int fusion_c_energy_fp_trial(int n,int nb,double dt,
    const double *edges,const double *old,const double *temperatures,
    const double *diffusion,const double *birth,const double *escape,
    double thermalization,double *trial,double *heat,
    fusion_kinetic_ledger_v1 *ledger) {
    if (ledger) *ledger={};
    if (n<1 || nb<0) return PB11_STATUS_INVALID_ARGUMENT;
    if (trial) std::fill(trial,trial+n,0.);
    if (heat) std::fill(heat,heat+nb,0.);
    if (!trial || !ledger || (nb && !heat)) return PB11_STATUS_NULL_OUTPUT;
    if (!edges || !old || !birth || !escape || (nb && !temperatures) ||
        (nb && n>1 && !diffusion)) return PB11_STATUS_INVALID_ARGUMENT;
    if (!std::isfinite(dt) || !std::isfinite(thermalization)) return PB11_STATUS_INVALID_ARGUMENT;
    if (dt<=0 || thermalization<0) return PB11_STATUS_OUT_OF_RANGE;
    for (int i=0;i<=n;++i) {
        if (!nonnegative(edges[i])) return PB11_STATUS_INVALID_ARGUMENT;
        if (i && edges[i]<=edges[i-1]) return PB11_STATUS_OUT_OF_RANGE;
    }
    for (int i=0;i<n;++i)
        if (!nonnegative(old[i]) || !nonnegative(birth[i]) || !nonnegative(escape[i]))
            return PB11_STATUS_INVALID_ARGUMENT;
    for (int b=0;b<nb;++b) {
        if (!std::isfinite(temperatures[b])) return PB11_STATUS_INVALID_ARGUMENT;
        if (temperatures[b]<=0) return PB11_STATUS_OUT_OF_RANGE;
    }
    try {
        const std::size_t nf=static_cast<std::size_t>(n-1);
        if (nb && nf>std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(nb))
            return PB11_STATUS_INVALID_ARGUMENT;
        const std::size_t count=nf*static_cast<std::size_t>(nb);
        std::vector<Real> energy(n),width(n),lower(n,0),diag(n,1),upper(n,0),rhs(n),state(n);
        std::vector<Real> left(count),right(count);
        for (int i=0;i<n;++i) {
            energy[i]=(static_cast<Real>(edges[i])+edges[i+1])/2;
            width[i]=static_cast<Real>(edges[i+1])-edges[i];
            diag[i]+=static_cast<Real>(dt)*(escape[i]+(i==0?static_cast<Real>(thermalization):0));
            rhs[i]=static_cast<Real>(old[i])+static_cast<Real>(dt)*birth[i];
        }
        for (int b=0;b<nb;++b) for (int f=0;f<n-1;++f) {
            const auto index=static_cast<std::size_t>(b)*nf+f;
            if (!nonnegative(diffusion[index])) return PB11_STATUS_INVALID_ARGUMENT;
            const Real distance=energy[f+1]-energy[f];
            const Real w=distance/temperatures[b]-.5L*std::log(energy[f+1]/energy[f]);
            Real bp,bm; bernoulli_pair(w,bp,bm);
            const Real coefficient=static_cast<Real>(diffusion[index])/distance;
            left[index]=coefficient*bp/width[f];
            right[index]=coefficient*bm/width[f+1];
            if (!finite(left[index]) || !finite(right[index])) return PB11_STATUS_NUMERICAL_FAILURE;
            diag[f]+=dt*left[index]; upper[f]-=dt*right[index];
            lower[f+1]-=dt*left[index]; diag[f+1]+=dt*right[index];
        }
        // Implicit Euler M-matrix. Reject loss of a positive pivot, never clip state.
        for (int i=1;i<n;++i) {
            if (!finite(diag[i-1]) || diag[i-1]<=0) return PB11_STATUS_NUMERICAL_FAILURE;
            const Real factor=lower[i]/diag[i-1];
            diag[i]-=factor*upper[i-1]; rhs[i]-=factor*rhs[i-1];
        }
        for (int i=n-1;i>=0;--i) {
            if (!finite(diag[i]) || diag[i]<=0) return PB11_STATUS_NUMERICAL_FAILURE;
            state[i]=(rhs[i]-(i<n-1?upper[i]*state[i+1]:0))/diag[i];
            if (!finite(state[i]) || state[i]<0) return PB11_STATUS_NUMERICAL_FAILURE;
        }
        std::vector<double> output(n),heat_output(nb);
        for (int i=0;i<n;++i) {
            // IEEE underflow of a positive tail is allowed; the total particle
            // and energy residuals below still bound any rounding loss. This
            // is not a population floor or a repair of negative values.
            if (state[i]>std::numeric_limits<double>::max())
                return PB11_STATUS_NUMERICAL_FAILURE;
            output[i]=static_cast<double>(state[i]);
            state[i]=output[i]; // Ledger describes the actual returned precision.
        }
        Real net_heat=0,heat_scale=0;
        for (int b=0;b<nb;++b) {
            Real bath_heat=0;
            for (int f=0;f<n-1;++f) {
                const auto index=static_cast<std::size_t>(b)*nf+f;
                const Real flux=left[index]*state[f]-right[index]*state[f+1];
                bath_heat-=dt*(energy[f+1]-energy[f])*flux;
            }
            if (!output_value(bath_heat,heat_output[b])) return PB11_STATUS_NUMERICAL_FAILURE;
            net_heat+=heat_output[b]; heat_scale+=std::abs(heat_output[b]);
        }
        Real n0=0,n1=0,e0=0,e1=0,nborn=0,eborn=0,nesc=0,eesc=0;
        for (int i=0;i<n;++i) {
            n0+=old[i]; n1+=state[i]; e0+=energy[i]*old[i]; e1+=energy[i]*state[i];
            const Real born=static_cast<Real>(dt)*birth[i];
            const Real escaped=static_cast<Real>(dt)*escape[i]*state[i];
            nborn+=born; eborn+=born*energy[i]; nesc+=escaped; eesc+=escaped*energy[i];
        }
        const Real nth=static_cast<Real>(dt)*thermalization*state[0],eth=nth*energy[0];
        const Real nr=n1-n0-nborn+nesc+nth;
        const Real er=e1-e0-eborn+eesc+eth+net_heat;
        const Real nscale=n0+nborn+nesc+nth+n1;
        const Real escale=e0+e1+eborn+eesc+eth+heat_scale;
        if (!finite(nr) || !finite(er) || std::abs(nr)>1e-10L*nscale || std::abs(er)>1e-10L*escale)
            return PB11_STATUS_NUMERICAL_FAILURE;
        fusion_kinetic_ledger_v1 result{};
        if (!output_value(n0,result.initial_number_m3) || !output_value(n1,result.final_number_m3) ||
            !output_value(e0,result.initial_energy_J_m3) || !output_value(e1,result.final_energy_J_m3) ||
            !output_value(nborn,result.born_number_m3) || !output_value(eborn,result.born_energy_J_m3) ||
            !output_value(nesc,result.escaped_number_m3) || !output_value(eesc,result.escaped_energy_J_m3) ||
            !output_value(nth,result.thermalized_number_m3) || !output_value(eth,result.thermalized_energy_J_m3) ||
            !output_value(nr,result.particle_balance_error_m3) || !output_value(er,result.energy_balance_error_J_m3))
            return PB11_STATUS_NUMERICAL_FAILURE;
        std::copy(output.begin(),output.end(),trial);
        if (nb) std::copy(heat_output.begin(),heat_output.end(),heat);
        *ledger=result;
        return PB11_STATUS_OK;
    } catch (...) {
        return PB11_STATUS_EXCEPTION;
    }
}
