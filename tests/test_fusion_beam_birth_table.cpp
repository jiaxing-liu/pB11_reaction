#include "fusion_beam_birth_table.h"
#include "fusion_nuclear_data.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
namespace {
constexpr double kev=1.602176634e-16,mev=1.602176634e-13;
constexpr int cells=64;
void check(bool p,const char*msg){if(!p)throw std::runtime_error(msg);}
template<class C,class F>void fields(C&c,F f){f(c.reactivity_m3_s);for(auto&v:c.reactant_energy_moment_J_m3_s)f(v);for(auto&v:c.below_number_m3_s)f(v);for(auto&v:c.below_energy_J_m3_s)f(v);for(auto&v:c.above_number_m3_s)f(v);for(auto&v:c.above_energy_J_m3_s)f(v);}
fusion_beam_birth_options_v1 options(){fusion_beam_birth_options_v1 o{};o.relative_max_J=2.5*mev;o.angular_max_exponent=40;o.ground_state_q_J=91.84*kev;o.cutoff_J=.001*mev;o.l1_fraction=.76;o.narrow_peak_fraction=.051;o.continuum_peak_scale=1;o.continuation=1;o.broad_mode=13;o.relative_order=16;o.angular_order=8;o.nq=4;o.ncos=4;return o;}
std::vector<double> edges(){std::vector<double> e(cells+1);for(int j=0;j<=cells;++j)e[j]=25*mev*j/cells;return e;}
fusion_birth_coefficients_v1 coefficients(const fusion_beam_birth_v1&b){fusion_birth_coefficients_v1 c{};const auto&r=b.spectrum;c.reactivity_m3_s=r.reactivity_m3_s;for(int j=0;j<2;++j)c.reactant_energy_moment_J_m3_s[j]=r.reactant_energy_moment_J_m3_s[j];for(int j=0;j<7;++j){c.below_number_m3_s[j]=r.below_number_m3_s[j];c.below_energy_J_m3_s[j]=r.below_energy_J_m3_s[j];c.above_number_m3_s[j]=r.above_number_m3_s[j];c.above_energy_J_m3_s[j]=r.above_energy_J_m3_s[j];}return c;}
void exact(const fusion_birth_coefficients_v1&a,const fusion_birth_coefficients_v1&b){std::vector<double>x,y;fields(a,[&](double v){x.push_back(v);});fields(b,[&](double v){y.push_back(v);});check(x==y,"exact knot coefficients");}
long double relative(long double d,long double r){return r?d/r:(d?std::numeric_limits<long double>::infinity():0);}
void verify(int ch,const std::vector<double>&e,const std::vector<double>&g,const fusion_birth_coefficients_v1&c){
 fusion_nuclear_channel_v1 reaction{};check(fusion_c_nuclear_channel(ch,&reaction)==0,"channel");long double energy=0;
 fields(c,[&](double v){check(std::isfinite(v)&&v>=0,"finite positive coefficients");});
 for(int id=0;id<7;++id){long double N=c.below_number_m3_s[id]+(long double)c.above_number_m3_s[id],E=c.below_energy_J_m3_s[id]+(long double)c.above_energy_J_m3_s[id];for(int j=0;j<cells;++j){double v=g[id*cells+j];check(std::isfinite(v)&&v>=0,"finite positive source");N+=v;E+=v*((long double)e[j]+e[j+1])/2;}int m=0;for(int j=0;j<reaction.product_count;++j){int p=reaction.product_ids[j];if(p==FUSION_MASS_NEUTRON)p=6;if(p==id)++m;}check(relative(std::abs(N-m*(long double)c.reactivity_m3_s),m*(long double)c.reactivity_m3_s)<1e-10L,"particle closure");energy+=E;}
 long double expected=reaction.q_J*(long double)c.reactivity_m3_s+c.reactant_energy_moment_J_m3_s[0]+(long double)c.reactant_energy_moment_J_m3_s[1];check(relative(std::abs(energy-expected),expected)<1e-10L,"energy closure");
}
void error_check(const std::vector<double>&e,const std::vector<double>&g,const fusion_birth_coefficients_v1&c,const std::vector<double>&ref,const fusion_birth_coefficients_v1&r,const fusion_birth_table_control_v1&q){
 check(relative(std::abs((long double)c.reactivity_m3_s-r.reactivity_m3_s),r.reactivity_m3_s)<=q.max_rate_error,"off-sample rate");
 for(int j=0;j<2;++j)check(relative(std::abs((long double)c.reactant_energy_moment_J_m3_s[j]-r.reactant_energy_moment_J_m3_s[j]),r.reactant_energy_moment_J_m3_s[j])<=q.max_debit_error,"off-sample debit");
 for(int id=0;id<7;++id){long double dN=std::abs((long double)c.below_number_m3_s[id]-r.below_number_m3_s[id])+std::abs((long double)c.above_number_m3_s[id]-r.above_number_m3_s[id]);long double dE=std::abs((long double)c.below_energy_J_m3_s[id]-r.below_energy_J_m3_s[id])+std::abs((long double)c.above_energy_J_m3_s[id]-r.above_energy_J_m3_s[id]);long double N=r.below_number_m3_s[id]+(long double)r.above_number_m3_s[id],E=r.below_energy_J_m3_s[id]+(long double)r.above_energy_J_m3_s[id];for(int j=0;j<cells;++j){int k=id*cells+j;long double center=((long double)e[j]+e[j+1])/2,d=std::abs((long double)g[k]-ref[k]);dN+=d;dE+=d*center;N+=ref[k];E+=ref[k]*center;}check(relative(dN,N)<=q.max_number_L1,"off-sample complete number shape");check(relative(dE,E)<=q.max_energy_L1,"off-sample complete energy shape");}
}
void run(int ch,int slot,double ep,double lo,double hi){
 auto e=edges();auto o=options();fusion_birth_table_control_v1 q{.002,.002,.01,.01,1e-4,1e-4,512,4096,16};fusion_beam_birth_table_v1*raw=nullptr;
 int st=fusion_c_beam_birth_table_create(ch,slot,ep,lo,hi,&o,&q,cells,e.data(),&raw);if(st)std::cerr<<"create channel="<<ch<<" status="<<st<<"\n";check(st==0&&raw,"create");std::unique_ptr<fusion_beam_birth_table_v1,decltype(&fusion_c_beam_birth_table_destroy)>t(raw,fusion_c_beam_birth_table_destroy);
 fusion_beam_birth_table_info_v1 info{};check(fusion_c_beam_birth_table_info(t.get(),&info)==0,"info");check(info.channel==ch&&info.projectile_slot==slot&&info.projectile_energy_J==ep&&info.cells==cells,"fixed model metadata");check(info.knots>=2&&info.direct_evaluations>=5&&info.stored_spectral_entries<7ULL*cells*info.knots,"sparse storage");check(info.spectral_entries_evaluated>=info.stored_spectral_entries,"entry counters");
 std::vector<double> g(7*cells),ref(7*cells);for(double f:{0.,.071,.193,.347,.583,.739,.917,1.}){double T=f==0?lo:f==1?hi:std::exp((1-f)*std::log(lo)+f*std::log(hi));fusion_birth_coefficients_v1 c{};fusion_beam_birth_v1 direct{};check(fusion_c_beam_birth_table_evaluate(t.get(),T,cells,g.data(),&c)==0,"evaluate");check(fusion_c_beam_birth_grid(ch,slot,ep,T,&o,cells,e.data(),ref.data(),&direct)==0,"independent direct sample");auto rc=coefficients(direct);if(f==0||f==1){check(g==ref,"exact endpoint grid");exact(c,rc);}error_check(e,g,c,ref,rc,q);verify(ch,e,g,c);}
 auto original=e;auto original_o=o;e[1]*=.5;o.angular_order=4;fusion_birth_coefficients_v1 c{};check(fusion_c_beam_birth_table_evaluate(t.get(),lo,cells,g.data(),&c)==0,"copied grid/options");fusion_beam_birth_v1 direct{};check(fusion_c_beam_birth_grid(ch,slot,ep,lo,&original_o,cells,original.data(),ref.data(),&direct)==0,"original direct");check(g==ref,"table independent of caller mutation");
 for(double T:{std::nextafter(lo,0.),std::nextafter(hi,std::numeric_limits<double>::infinity()),std::numeric_limits<double>::quiet_NaN()}){std::fill(g.begin(),g.end(),7);fields(c,[](double&v){v=7;});check(fusion_c_beam_birth_table_evaluate(t.get(),T,cells,g.data(),&c)!=0,"bad temperature rejected");check(std::all_of(g.begin(),g.end(),[](double v){return v==0;}),"failed output grid clear");fields(c,[](double v){check(v==0,"failed coefficient clear");});}
 q.max_knots=2;q.max_evaluations=5;q.max_depth=0;q.max_rate_error=q.max_debit_error=q.max_number_L1=q.max_energy_L1=0;raw=reinterpret_cast<fusion_beam_birth_table_v1*>(1);check(fusion_c_beam_birth_table_create(ch,slot,ep,lo,hi,&original_o,&q,cells,original.data(),&raw)!=0&&raw==nullptr,"insufficient construction budget rejects atomically");
}
void subnormal_endpoint(){
 constexpr int n=800;std::vector<double> e(n+1),g(7*n),ref(7*n);
 for(int j=1;j<=n/3;++j)e[j]=1e-10*kev*std::pow(1e12,double(j-1)/(n/3-1));
 for(int j=n/3+1;j<=n;++j)e[j]=(100+24900.*(j-n/3)/(n-n/3))*kev;
 auto o=options();o.angular_order=16;o.nq=o.ncos=8;
 fusion_birth_table_control_v1 q{.01,.01,.01,.01,1e-5,1e-5,2,5,0};
 double ep=2.1399859254691007e-13,T=1.4401072590598843e-17;
 fusion_beam_birth_table_v1*raw=nullptr;
 check(fusion_c_beam_birth_table_create(3,1,ep,T,T*1.001,&o,&q,n,e.data(),&raw)==0,"subnormal table create");
 std::unique_ptr<fusion_beam_birth_table_v1,decltype(&fusion_c_beam_birth_table_destroy)>t(raw,fusion_c_beam_birth_table_destroy);
 fusion_birth_coefficients_v1 c{};fusion_beam_birth_v1 d{};
 check(fusion_c_beam_birth_table_evaluate(t.get(),T,n,g.data(),&c)==0,"subnormal table eval");
 check(fusion_c_beam_birth_grid(3,1,ep,T,&o,n,e.data(),ref.data(),&d)==0,"subnormal direct source");
 check(std::any_of(ref.begin(),ref.end(),[](double v){return std::fpclassify(v)==FP_SUBNORMAL;}),"fixture contains subnormal source entries");
 check(g==ref,"all subnormal grid values retained exactly");exact(c,coefficients(d));
}
}
int main(){try{run(3,1,2.1399859254691007e-13,.05*kev,.2*kev);run(4,0,1.5*mev,.08*kev,.12*kev);run(0,0,1.*mev,.099*kev,.101*kev);run(1,0,.5*mev,.08*kev,.12*kev);run(2,0,.5*mev,.08*kev,.12*kev);run(3,0,0.,5*kev,5.1*kev);subnormal_endpoint();fusion_c_beam_birth_table_destroy(nullptr);std::cout<<"PASS: full beam tables, exact knots, off-sample rate/debit/shape, budgets, copied inputs and clearing\n";return 0;}catch(const std::exception&e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}}
