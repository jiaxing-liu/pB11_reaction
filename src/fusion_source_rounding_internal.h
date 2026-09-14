#ifndef PB11_FUSION_SOURCE_ROUNDING_INTERNAL_H
#define PB11_FUSION_SOURCE_ROUNDING_INTERNAL_H
#include <cmath>
#include <limits>
namespace fusion_detail { namespace source_rounding {
using R=long double;
// The default protects every positive mapped floor packet. Optional accounting
// admits only naturally rounded zero rates whose packet AND rate are below
// half a binary64 subnormal quantum. Caller must bound aggregate missing N/E.
// Reconstruction error is measured on source scale, never old inventory.
inline bool source_roundoff_accumulate(R lost,R center,R&number,R&energy){
 if(!std::isfinite(lost)||lost<0||!std::isfinite(center)||center<=0||
    !std::isfinite(number)||number<0||!std::isfinite(energy)||energy<0)return false;
 // Upward bounds cover long-double addition/product rounding as well.
 const R up=std::numeric_limits<R>::infinity();
 const R candidateN=lost==0?number:std::nextafter(number+lost,up);
 const R candidateE=lost==0?energy:std::nextafter(energy+std::nextafter(lost*center,up),up);
 const R half_quantum=R(std::numeric_limits<double>::denorm_min())/2;
 if(!std::isfinite(candidateN)||!std::isfinite(candidateE)||
    candidateN>=half_quantum||candidateE>=half_quantum)return false;
 number=candidateN;energy=candidateE;return true;
}
inline bool floor_source_rate(R amount,double dt,double& rate,R*unrepresented=nullptr){
 if(unrepresented)*unrepresented=0;
 if(amount<0||!std::isfinite(amount)||!std::isfinite(dt)||dt<=0)return false;
 R lost=0;
 const R exact_rate=amount/R(dt);
 if(!std::isfinite(exact_rate)||std::abs(exact_rate)>std::numeric_limits<double>::max())return false;
 rate=double(exact_rate);
 if(amount>0&&rate==0){
  const R half_quantum=R(std::numeric_limits<double>::denorm_min())/2;
  if(!unrepresented||amount>=half_quantum||exact_rate>=half_quantum)return false;
  lost=amount;
 }
 const R recovered=R(rate)*dt;
 const R allowance=2*std::numeric_limits<double>::epsilon()*std::abs(amount)+
    R(std::numeric_limits<double>::denorm_min())*dt;
 if(!std::isfinite(recovered)||std::abs(recovered-amount)>allowance)return false;
 if(unrepresented)*unrepresented=lost;
 return true;
}
}}
#endif
