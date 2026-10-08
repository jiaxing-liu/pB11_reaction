// Exercise the private numerical kernel without adding a public ABI.
#include "../src/fusion_thermal_birth.cpp"
#include <iostream>

static void check(const Mapper& m,int id,long double expected_n,long double expected_e){
 long double N=m.bn[id]+m.an[id],U=m.be[id]+m.ae[id];
 for(int j=0;j<m.n;++j){auto x=m.birth[id*m.n+j];if(!std::isfinite(x)||x<0)throw std::runtime_error("invalid birth");N+=x;U+=x*m.center[j];}
 if(std::abs(N-expected_n)>1e-14L*expected_n || std::abs(U-expected_e)>1e-14L*expected_e)throw std::runtime_error("delta number/energy closure");
 for(int other=0;other<7;++other)if(other!=id){if(m.number[other]!=0||m.energy[other]!=0||m.bn[other]!=0||m.an[other]!=0)throw std::runtime_error("species leakage");for(int j=0;j<m.n;++j)if(m.birth[other*m.n+j]!=0)throw std::runtime_error("birth leakage");}
}
int main(){try{
 unsigned long long seed=431;unsigned legacy_endpoint_inversions=0,cases=0;
 auto random=[&](){seed=seed*6364136223846793005ULL+1;return static_cast<long double>(seed>>11)/9007199254740992.L;};
 for(int n:{1,8,508}){
  std::vector<double> edges(n+1);edges[0]=0;for(int j=1;j<=n;++j)edges[j]=double(std::pow(10.L,-30.L+21.L*j/n));
  for(int k=0;k<3000;++k){
   const long double mass=6.6446573357e-27L,K=k%101==0?0:std::pow(10.L,-28+16*random()),w=k%103==0?0:std::pow(10.L,-240+480*random());
   if(K>0&&K*K/K>K)++legacy_endpoint_inversions;
   Mapper m(n,edges.data());m.isotropic(k%7,mass,K,Mapper::Boost(0),w);check(m,k%7,w,w*K);++cases;
  }
  for(long double beta:{0.L,1e-12L,.01L,.2L}){
   Mapper m(n,edges.data());Mapper::Boost boost(beta);m.isotropic(4,6.6446573357e-27L,0,boost,3.L);check(m,4,3.L,3.L*boost.d*6.6446573357e-27L*c2);++cases;
  }
 }
 if(!legacy_endpoint_inversions)throw std::runtime_error("probe did not exercise rounding regression");
 std::cout<<"PASS zero-width spectra cases="<<cases<<" legacy_endpoint_inversions="<<legacy_endpoint_inversions<<"\n";
 return 0;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}}
