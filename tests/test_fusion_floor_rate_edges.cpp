#include "fusion_birth_floor_precise_internal.h"
#include "fusion_source_rounding_internal.h"
#include "fusion_kinetic_rounding_internal.h"
#include <cassert>
#include <cmath>
#include <limits>

int main() {
 using R=long double;
 using namespace fusion_detail::source_rounding;
 const R quantum=std::numeric_limits<double>::denorm_min();
 double rate=0; R lost=0;
 // Amount underflows while the rate is representable: retain the rate.
 const R amount=1.260606261366177475635e-327L;
 assert(double(amount)==0);
 assert(floor_source_rate(amount,6.249999999999312e-6,rate,&lost));
 assert(rate>0 && lost==0);
 // Aggregate first: individually unrepresentable rates can become visible.
 const R small=quantum*.3L;
 assert(floor_source_rate(small,1,rate,&lost) && rate==0 && lost==small);
 assert(floor_source_rate(small+small,1,rate,&lost) && rate>0 && lost==0);
 // Repeated unrepresentable amounts must not evade the aggregate budget.
 R missingN=0,missingE=0;
 assert(source_roundoff_accumulate(small,1,missingN,missingE));
 const R savedN=missingN,savedE=missingE;
 assert(!source_roundoff_accumulate(small,1,missingN,missingE));
 assert(missingN==savedN && missingE==savedE);
 assert(!floor_source_rate(1,0,rate,&lost));
 assert(!floor_source_rate(std::numeric_limits<R>::infinity(),1,rate,&lost));
 // Exact binary cases test inclusive borrowing limits and positive remainder.
 std::array<R,6> n{},e{};fusion_detail::precise_floor_ledger out{};
 n[0]=1;n[4]=1;e[0]=e[4]=0;
 assert(fusion_detail::project_floor_precise(.125L,1,1,1,.25,n,e,out)==0);
 assert(out.borrowed==.25L && out.correction[0]==-.125L);
 const double lower=std::nextafter(.25,0.);
 assert(fusion_detail::project_floor_precise(.125L,1,1,1,lower,n,e,out)==PB11_STATUS_OUT_OF_RANGE);
 assert(out.number[0]==0 && out.number[4]==0 && out.borrowed==0);
 assert(fusion_detail::project_floor_precise(.125L,1,.25L,1,1,n,e,out)==PB11_STATUS_OUT_OF_RANGE);
 // An empty floor must leave the temperature gate inactive.
 n={};e={};
 assert(fusion_detail::project_floor_precise(1,1,1,0,0,n,e,out)==0);
 // A trace population activates it even when its public amount rounds to zero.
 n[4]=amount;
 assert(fusion_detail::project_floor_precise(1,1,1,0,1,n,e,out)==PB11_STATUS_OUT_OF_RANGE);
 assert(out.number[4]==0);
 // Mixed ordinary and trace populations retain each species' physical energy.
 n[0]=2;e[0]=.0625L;e[4]=amount*.015625L;
 assert(fusion_detail::project_floor_precise(.125L,1,1,1,.5,n,e,out)==0);
 for(int i:{0,4}) {
  assert(out.number[i]==n[i]);
  assert(out.physical_energy[i]==e[i]);
  const R residual=std::abs(out.mapped_energy[i]+out.correction[i]-e[i]);
  const R arithmetic_bound=4*std::numeric_limits<R>::epsilon()*
   (std::abs(out.mapped_energy[i])+std::abs(out.correction[i]));
  assert(residual<=arithmetic_bound);
 }
 // Invalid input after success clears the complete output candidate.
 e[4]=-1;
 assert(fusion_detail::project_floor_precise(.125L,1,1,1,.5,n,e,out)==PB11_STATUS_INVALID_ARGUMENT);
 assert(out.number[0]==0 && out.number[4]==0 && out.borrowed==0);
}
