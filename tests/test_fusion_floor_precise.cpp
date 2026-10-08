#include "fusion_birth_floor_precise_internal.h"
#include "fusion_source_rounding_internal.h"
#include <cassert>
#include <limits>
int main(){
 std::array<long double,6> n{},e{};n[4]=1.260606261366177475635e-327L;e[4]=5.426627245485451796026e-354L;
 fusion_detail::precise_floor_ledger out{};
 assert(fusion_detail::project_floor_precise(1e-17L,1e-15L,1e4L,.1,.1,n,e,out)==0);
 assert(out.number[4]>0 && double(out.number[4])==0);
 double rate=0;assert(fusion_detail::source_rounding::floor_source_rate(out.number[4],6.249999999999312e-6,rate));assert(rate>0);
 assert(out.correction[4]<0 && out.physical_energy[4]==e[4]);
 assert(fusion_detail::project_floor_precise(1e-17L,1e-15L,1e4L,0,.1,n,e,out)==PB11_STATUS_OUT_OF_RANGE);
 assert(out.number[4]==0);
 assert(fusion_detail::project_floor_precise(1e-17L,1e-15L,1e4L,.1,0,n,e,out)==PB11_STATUS_OUT_OF_RANGE);
 n[4]=2;e[4]=1e-17L;assert(fusion_detail::project_floor_precise(1e-17L,1e-15L,1e4L,.1,.1,n,e,out)==0);
 assert(out.mapped_energy[4]+out.correction[4]==out.physical_energy[4]);
 e[4]=3e-17L;assert(fusion_detail::project_floor_precise(1e-17L,1e-15L,1e4L,.1,.1,n,e,out)==PB11_STATUS_OUT_OF_RANGE);
}
