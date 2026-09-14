#ifdef NDEBUG
#undef NDEBUG
#endif
#include "fusion_coupled_thermal.h"
#include <cassert>
#include <cstring>
#include <limits>
#include <cstdio>
int main(){fusion_source_ledger_v1 l{};double inert[6]{},num[6]{};fusion_thermal_increment_v1 old{},out{},zero{};
 l.thermal_consumed_number_m3[0]=2;l.handed_off_number_m3[0]=3;
 l.heat_to_bath_J_m3[0]=5;l.heat_to_bath_J_m3[1]=8;inert[2]=2;
 assert(fusion_c_coupled_thermal_increment(&l,inert,&old)==0);
 assert(fusion_c_coupled_numerical_increment(&l,inert,num,&out)==0&&std::memcmp(&old,&out,sizeof out)==0);
 num[4]=-3;assert(fusion_c_coupled_numerical_increment(&l,inert,num,&out)==0);
 assert(out.ion_energy_J_m3==7&&out.electron_energy_J_m3==5&&out.thermal_number_m3[0]==1);
 num[0]=4;assert(fusion_c_coupled_numerical_increment(&l,inert,num,&out)==0&&out.ion_energy_J_m3==11);
 num[0]=std::numeric_limits<double>::quiet_NaN();assert(fusion_c_coupled_numerical_increment(&l,inert,num,&out)!=0&&std::memcmp(&out,&zero,sizeof out)==0);
 assert(fusion_c_coupled_numerical_increment(&l,inert,nullptr,&out)!=0&&std::memcmp(&out,&zero,sizeof out)==0);
 assert(fusion_c_coupled_numerical_increment(&l,inert,num,nullptr)==PB11_STATUS_NULL_OUTPUT);
 l={};for(double& x:inert)x=0;for(double&x:num)x=0;
 l.handed_off_energy_J_m3[0]=1e16;inert[0]=1;num[0]=-1e16;
 assert(fusion_c_coupled_thermal_increment(&l,inert,&old)==0&&old.ion_energy_J_m3==1e16);
 assert(fusion_c_coupled_numerical_increment(&l,inert,num,&out)==0&&out.ion_energy_J_m3==1);
 // Separate physical and numerical finite amounts may cancel before rounding.
 l={};inert[0]=0;num[0]=std::numeric_limits<double>::max();num[1]=num[0];
 assert(fusion_c_coupled_numerical_increment(&l,inert,num,&out)==PB11_STATUS_NUMERICAL_FAILURE&&std::memcmp(&out,&zero,sizeof out)==0);
 puts("PASS numerical fluid increment: zero parity, signed six-species sum, cancellation before rounding, atomic failures");}
