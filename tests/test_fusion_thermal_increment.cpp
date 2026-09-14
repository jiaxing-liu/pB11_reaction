#include "fusion_coupled_thermal.h"
#include <limits>
#include <iostream>
#include <stdexcept>
#include <cmath>
void need(bool p,const char*m){if(!p)throw std::runtime_error(m);}
int main(){try{
 static_assert(sizeof(fusion_thermal_increment_v1)==64,"ABI");
 fusion_source_ledger_v1 l{};double inert[6]{};fusion_thermal_increment_v1 d{};
 auto run=[&](){return fusion_c_coupled_thermal_increment(&l,inert,&d);};
 auto cleared=[&](){for(double x:d.thermal_number_m3)need(x==0,"clear N");need(d.electron_energy_J_m3==0&&d.ion_energy_J_m3==0,"clear energies");};
 need(run()==0,"zero ledger");cleared();
 l.thermal_consumed_number_m3[0]=7;l.handed_off_number_m3[0]=2;
 l.handed_off_number_m3[4]=9;
 l.thermal_consumed_energy_J_m3[0]=11;l.handed_off_energy_J_m3[4]=13;
 l.heat_to_bath_J_m3[0]=3;l.heat_to_bath_J_m3[35]=-1;
 l.heat_to_bath_J_m3[1]=-4;l.heat_to_bath_J_m3[41]=7;
 inert[0]=5;inert[5]=-2;
 // FAST birth/escape are intentionally not immediate thermal sources.
 l.external_born_number_m3[1]=123;l.external_born_energy_J_m3[1]=456;
 l.escaped_number_m3[2]=23;l.escaped_energy_J_m3[2]=56;
 need(run()==0,"signed ledger");need(d.thermal_number_m3[0]==-5&&d.thermal_number_m3[4]==9,"net particle increments");
 need(d.electron_energy_J_m3==2&&d.ion_energy_J_m3==8,"signed handoff/bath/consumption energy counted once");
 l={};for(double&x:inert)x=0;
 l.heat_to_bath_J_m3[0]=1e-30;l.thermal_consumed_number_m3[1]=1e-30;
 need(run()==0,"weak signal");const double background=1e20;
 need(background+d.electron_energy_J_m3==background,"inventory subtraction would lose heat");
 need(d.electron_energy_J_m3==1e-30&&d.thermal_number_m3[1]==-1e-30,"weak amounts preserved");
 l.escaped_energy_J_m3[5]=-1;need(run()==PB11_STATUS_INVALID_ARGUMENT,"negative ledger rejected");cleared();l.escaped_energy_J_m3[5]=0;
 inert[2]=std::numeric_limits<double>::quiet_NaN();need(run()==PB11_STATUS_INVALID_ARGUMENT,"NaN inert rejected");cleared();inert[2]=0;
 l.heat_to_bath_J_m3[4]=std::numeric_limits<double>::infinity();need(run()==PB11_STATUS_INVALID_ARGUMENT,"nonfinite signed heat rejected");cleared();l.heat_to_bath_J_m3[4]=0;
 l.heat_to_bath_J_m3[0]=l.heat_to_bath_J_m3[7]=std::numeric_limits<double>::max();need(run()==PB11_STATUS_NUMERICAL_FAILURE,"overflow rejected");cleared();
 need(fusion_c_coupled_thermal_increment(nullptr,inert,&d)==PB11_STATUS_INVALID_ARGUMENT,"null ledger");cleared();
 need(fusion_c_coupled_thermal_increment(&l,nullptr,&d)==PB11_STATUS_INVALID_ARGUMENT,"null inert");cleared();
 need(fusion_c_coupled_thermal_increment(&l,inert,nullptr)==PB11_STATUS_NULL_OUTPUT,"null output");
 std::cout<<"Thermal increment signed accounting and weak-signal tests passed\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
