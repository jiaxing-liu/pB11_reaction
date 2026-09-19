#include "fusion_coupled_sources.h"
extern "C" void packet_layout_fixture(int n,double*v,fusion_birth_packets_v1*p){
 for(int s=0;s<2;++s)for(int ch=0;ch<5;++ch){p->events_m3[s][ch]=100*(s+1)+ch+1;
  for(int sp=0;sp<7;++sp){
   p->below_number_m3[s][ch][sp]=1000*(s+1)+100*(ch+1)+sp+1;
   p->below_energy_J_m3[s][ch][sp]=2000*(s+1)+100*(ch+1)+sp+1;
   p->above_number_m3[s][ch][sp]=3000*(s+1)+100*(ch+1)+sp+1;
   p->above_energy_J_m3[s][ch][sp]=4000*(s+1)+100*(ch+1)+sp+1;
   for(int j=0;j<n;++j)v[((s*5+ch)*7+sp)*n+j]=10000*(s+1)+1000*(ch+1)+100*(sp+1)+j+1;
  }
 }
}
