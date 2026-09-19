#ifndef FUSION_BOUNDARY_ENERGY_HPP
#define FUSION_BOUNDARY_ENERGY_HPP
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
namespace fusion_boundary_energy {
// Binary64 N is authoritative for monoenergetic energy-domain outflow.
// Only an exact product below DBL_MIN uses this canonical representation.
// Integer significands avoid ABI-dependent long-double double rounding and FTZ.
inline bool subnormal_product(double n,double edge,double& rounded) {
 static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559 &&
               std::numeric_limits<double>::digits==53,"binary64 required");
 if(!std::isfinite(n)||!std::isfinite(edge)||n<=0||edge<=0)return false;
 uint64_t a,b;std::memcpy(&a,&n,8);std::memcpy(&b,&edge,8);
 auto unpack=[](uint64_t bits,uint64_t& sig,int& exponent){
  const unsigned e=unsigned((bits>>52)&2047);sig=bits&UINT64_C(0xfffffffffffff);
  exponent=e?int(e)-1075:-1074;if(e)sig|=UINT64_C(0x10000000000000);
 };
 uint64_t sa,sb;int ea,eb;unpack(a,sa,ea);unpack(b,sb,eb);
 using U=boost::multiprecision::uint128_t;
 const U product=U(sa)*sb;const int exponent=ea+eb;
 if(int(boost::multiprecision::msb(product))+exponent>=-1022)return false;
 const int shift=exponent+1074;U q=0;
 if(shift>=0)q=product<<shift;
 else if(-shift<128){
  const unsigned k=unsigned(-shift);q=product>>k;
  const U remainder=product-(q<<k),half=U(1)<<(k-1);
  if(remainder>half||(remainder==half&&bool(q&1)))++q;
 }
 const uint64_t bits=q.convert_to<uint64_t>();std::memcpy(&rounded,&bits,8);
 return true;
}
inline void canonicalize(double n,double edge,double& u){double rounded;
 if(subnormal_product(n,edge,rounded))u=rounded;
}
inline long double accounting(double n,double edge,double u){double rounded;
 return subnormal_product(n,edge,rounded)?static_cast<long double>(n)*edge:u;
}
}
#endif
