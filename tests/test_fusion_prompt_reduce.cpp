#include "fusion_prompt_reduce.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <iostream>
void ck(bool b,const char*m){if(!b)throw std::runtime_error(m);}
int main(){try{
 fusion_prompt_totals_v1 out{},zero{};
 ck(fusion_c_prompt_reduce(0,nullptr,&out)==0&&!out.number_fraction_defined&&!out.energy_fraction_defined,"empty undefined");
 fusion_prompt_weight_v1 rows[]={{2,3,0,2},{3,2,0,1},{5,4,0,0}};
 ck(fusion_c_prompt_reduce(3,rows,&out)==0,"mixed status");
 ck(out.total_number==10&&out.total_energy_J==32&&out.number[0]==5&&out.energy_J[0]==20,"derived energy accounting");
 ck(out.number_fraction_low==.2&&out.number_fraction_high==.7&&out.energy_fraction_low==6./32&&out.energy_fraction_high==26./32,"weighted bounds");
 auto saved=out;ck(fusion_c_prompt_reduce(3,rows,&out)==0&&std::memcmp(&out,&saved,sizeof(out))==0,"repeat");
 for(int outcome:{0,1,2}){fusion_prompt_weight_v1 r{4,2,0,outcome};ck(fusion_c_prompt_reduce(1,&r,&out)==0,"pure category");ck(out.number_fraction_low==(outcome==2?1.:0.)&&out.number_fraction_high==(outcome==1?0.:1.),"pure fractions");}
 fusion_prompt_weight_v1 r{1,0,0,1};ck(fusion_c_prompt_reduce(1,&r,&out)==0&&out.number_fraction_defined&&!out.energy_fraction_defined,"cold energy undefined");
 r={0,3,0,0};ck(fusion_c_prompt_reduce(1,&r,&out)==0&&!out.number_fraction_defined,"zero weight");
 rows[2].api_status=1234;ck(fusion_c_prompt_reduce(3,rows,&out)==1234&&std::memcmp(&out,&zero,sizeof(out))==0,"error not unresolved");
 r={0,0,1235,0};ck(fusion_c_prompt_reduce(1,&r,&out)==1235,"zero weight error explicit");
 r={-1,1,0,2};ck(fusion_c_prompt_reduce(1,&r,&out)==PB11_STATUS_INVALID_ARGUMENT,"negative");
 r={1,1,0,3};ck(fusion_c_prompt_reduce(1,&r,&out)==PB11_STATUS_INVALID_ARGUMENT,"invalid outcome");
 r={1,std::numeric_limits<double>::infinity(),0,1};ck(fusion_c_prompt_reduce(1,&r,&out)==PB11_STATUS_INVALID_ARGUMENT,"nonfinite");
 r={std::numeric_limits<double>::max(),2,0,2};ck(fusion_c_prompt_reduce(1,&r,&out)==PB11_STATUS_NUMERICAL_FAILURE&&out.total_number==0,"overflow atomic");
 r={std::numeric_limits<double>::denorm_min(),.25,0,2};ck(fusion_c_prompt_reduce(1,&r,&out)==PB11_STATUS_NUMERICAL_FAILURE,"positive energy underflow explicit");
 r={std::numeric_limits<double>::denorm_min(),1,0,2};ck(fusion_c_prompt_reduce(1,&r,&out)==0&&out.total_number==r.number_weight&&out.number_fraction_low==1,"subnormal retained");
 std::cout<<"prompt reduction counts/energy/bounds/undefined/error/representation tests pass\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
