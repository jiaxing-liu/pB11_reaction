#include "fusion_source_state.h"
#include <vector>
#include <fstream>
#include <stdexcept>
int main(int argc,char**argv){if(argc!=2)return 1;auto ok=[](int s){if(s)throw std::runtime_error("state failure");};double edges[3]={0,2,4},s[12]{},t[12]{};for(int i=0;i<6;++i)s[2*i+1]=20*(i+1);fusion_source_state_v1*p=nullptr;ok(fusion_c_source_state_create(2,edges,s,t,3.25,0x81726354,&p));uint64_t ticket=0;ok(fusion_c_source_state_begin(p,.5,&ticket));fusion_source_ledger_v1 l{};for(int i=0;i<6;++i){s[2*i]=5*(i+1);s[2*i+1]=15*(i+1);l.heat_to_bath_J_m3[7*i]=10*(i+1);}ok(fusion_c_source_state_stage(p,ticket,s,t,&l));ok(fusion_c_source_state_commit(p,ticket));size_t n=0,w=0;ok(fusion_c_source_state_pack_size(p,&n));std::vector<unsigned char>b(n);ok(fusion_c_source_state_pack(p,b.data(),n,&w));std::ofstream f(argv[1],std::ios::binary);f.write(reinterpret_cast<const char*>(b.data()),b.size());fusion_c_source_state_destroy(p);return !f;}
