#include "fusion_coupled_thermal.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

int main(int argc,char**argv) {
 if(argc!=4)return 2;
 std::ifstream file(argv[1]);std::vector<double> e;double value;
 while(file>>value)e.push_back(value);
 const int n=int(e.size())-1,repetitions=std::atoi(argv[2]);
 if(n!=508||repetitions<1)return 2;
 constexpr double kev=1.602176634e-16,MeV=1.602176634e-13;
 double ne=7e19,N[6]{},Z2[6]={1,1,1,4,4,25};
 fusion_inert_ion_v1 carbon{.001*ne,12*1.66053906892e-27,36};
 N[4]=1e-5*ne;N[1]=N[3]=(ne-6*carbon.density_m3-2*N[4])/3;
 double pool=carbon.density_m3;for(double x:N)pool+=x;
 const double Ue=1.5*ne*.5*kev,Ui=1.5*pool*.1*kev;
 fusion_coupled_thermal_options_v1 o{};
 o.birth={10*MeV,40,.09184*MeV,.001*MeV,.76,0,.051,1,1,0,0,13,0,16,12,8,8};
 o.max_source_rate_error=o.max_source_debit_error=1e-5;
 o.handoff_max_L1=o.handoff_max_mean_error=.001;o.handoff_enabled=1;
 std::vector<double> s(6*n),t(6*n),sn(6*n),tn(6*n),zero(6*n),logs(48,10.);
 for(int species : {0,4}) {
  const double desired=(species==0?14700:3600)*kev;
  int j=0;while(j+1<n && (e[j]+e[j+1])/2<desired)++j;
  s[species*n+j]=1e11;
 }
 double Nnew[6]{};fusion_coupled_thermal_v1 out{};
 auto trial=[&](){return fusion_c_coupled_thermal_trial(1.5625e-6,&o,n,e.data(),N,Ue,Ui,ne,Z2,
     1,&carbon,logs.data(),s.data(),t.data(),zero.data(),zero.data(),Nnew,sn.data(),tn.data(),&out);};
 int st=trial();if(st){std::fprintf(stderr,"status %d\n",st);return st;}
 const auto start=std::chrono::steady_clock::now();
 for(int i=0;i<repetitions;++i){st=trial();if(st)return st;}
 const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 auto *destination=std::fopen(argv[3],"wb");if(!destination)return 3;
 std::fwrite(Nnew,sizeof(double),6,destination);
 std::fwrite(sn.data(),sizeof(double),sn.size(),destination);
 std::fwrite(tn.data(),sizeof(double),tn.size(),destination);
 // The ABI prefix contains doubles only; omit all possible tail padding.
 std::fwrite(&out,1,offsetof(fusion_coupled_thermal_v1,handoff_projected),destination);
 std::fwrite(out.handoff_projected,sizeof(int),6,destination);std::fclose(destination);
 std::printf("{\"seconds\":%.17g,\"trials\":%d,\"seconds_per_trial\":%.17g}\n",seconds,repetitions,seconds/repetitions);
}
