#include "fusion_products.h"
#include "fusion_kinetics.h"
#include "fusion_coulomb.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
constexpr double kev=1.602176634e-16,mev=1000*kev,alpha=6.644657345e-27,c=299792458.;
void require(bool condition,const char *label) {if(!condition) throw std::runtime_error(label);}
bool near(double a,double b,double tolerance=1e-12) {
    return std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=tolerance*std::max(std::abs(a),std::abs(b));
}
void mapping_and_spill() {
    // Centers 1,3,7 keV. Input at 2 -> equal first/second weights.
    const double edges[]={0,2*kev,4*kev,10*kev};
    const double energies[]={2*kev,5*kev,0,9*kev};
    const double weights[]={4e10,6e10,2e10,3e10};
    double output[3]={};fusion_birth_mapping_v1 ledger{};
    require(fusion_c_map_birth_packets(3,edges,4,energies,weights,output,&ledger)==0,"mapping status");
    require(near(output[0],2e10) && near(output[1],5e10) && near(output[2],3e10),"adjacent-center conservative weights");
    require(near(ledger.below_number_m3_s,2e10) && ledger.below_energy_W_m3==0,"below hull carries number and actual energy");
    require(near(ledger.above_number_m3_s,3e10) && near(ledger.above_energy_W_m3,27e10*kev),"above hull is not clipped to final center");
    require(near(ledger.input_number_m3_s,ledger.mapped_number_m3_s+ledger.below_number_m3_s+ledger.above_number_m3_s),"input number independently conserved");
    require(near(ledger.input_energy_W_m3,ledger.mapped_energy_W_m3+ledger.below_energy_W_m3+ledger.above_energy_W_m3),"input energy independently conserved");
    const double one_edges[]={0,2*kev},one_energy[]={kev};double one_birth[1];
    require(fusion_c_map_birth_packets(1,one_edges,1,one_energy,weights,one_birth,&ledger)==0 &&
            one_birth[0]==weights[0],"single-center exact packet");
    require(fusion_c_map_birth_packets(3,edges,0,nullptr,nullptr,output,&ledger)==0 &&
            ledger.input_number_m3_s==0 && output[0]==0 && output[2]==0,"empty input is valid");
    const double negative[]={-kev};
    require(fusion_c_map_birth_packets(3,edges,1,negative,weights,output,&ledger)==PB11_STATUS_INVALID_ARGUMENT &&
            ledger.input_number_m3_s==0 && output[0]==0,"invalid packet clears outputs");
    const double huge[]={std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};
    const double exact[]={kev,kev};
    require(fusion_c_map_birth_packets(3,edges,2,exact,huge,output,&ledger)==PB11_STATUS_NUMERICAL_FAILURE &&
            output[0]==0 && ledger.input_number_m3_s==0,"unrepresentable population rejected");
}

void sequential_invariants() {
    const long double rest=static_cast<long double>(alpha)*c*c;
    const double available=8.68*mev;
    for(double fraction:{0.,.01,.1,.25-1e-12,.25,.25+1e-12,.5,.9,1.})
        for(double cosine:{-1.,-.5,0.,.5,1.}) {
            fusion_three_body_cm_v1 r{};
            const double q=fraction*available;
            require(fusion_c_three_equal_sequential_cm(alpha,available,q,cosine,&r)==0,"sequential kinematic status");
            long double energy=0,px=0,pz=0;
            for(int i=0;i<3;++i) {
                require(r.kinetic_energy_J[i]>=0,"nonnegative product kinetic energy");
                energy+=r.kinetic_energy_J[i];px+=r.momentum_x_kg_m_s[i];pz+=r.momentum_z_kg_m_s[i];
                // Check E^2-c^2 p^2=m^2 c^4 using the independently returned
                // energy and momentum arrays, scaled to the kinetic invariant.
                const long double t=r.kinetic_energy_J[i];
                const long double x=r.momentum_x_kg_m_s[i],z=r.momentum_z_kg_m_s[i];
                const long double momentum_invariant=(x*x+z*z)*c*c;
                require(std::abs(t*(t+2*rest)-momentum_invariant)<=1e-11L*std::max(momentum_invariant,rest*available),"product relativistic mass shell");
            }
            require(std::abs(energy-available)<1e-12*available,"three products conserve CM kinetic energy");
            require(std::abs(px)+std::abs(pz)<1e-12*std::sqrt(2*alpha*available),"three products conserve momentum");
            const long double pair_t=static_cast<long double>(r.kinetic_energy_J[1])+r.kinetic_energy_J[2];
            const long double pair_x=static_cast<long double>(r.momentum_x_kg_m_s[1])+r.momentum_x_kg_m_s[2];
            const long double pair_z=static_cast<long double>(r.momentum_z_kg_m_s[1])+r.momentum_z_kg_m_s[2];
            const long double invariant=(2*rest+pair_t)*(2*rest+pair_t)-(pair_x*pair_x+pair_z*pair_z)*c*c;
            const long double reconstructed_q=(invariant-4*rest*rest)/(std::sqrt(invariant)+2*rest);
            require(std::abs(reconstructed_q-q)<1e-10*available,"intermediate invariant rest energy matches supplied q");
            fusion_three_body_cm_v1 reflected{};
            require(fusion_c_three_equal_sequential_cm(alpha,available,q,-cosine,&reflected)==0,"reflected kinematic status");
            require(near(r.kinetic_energy_J[1],reflected.kinetic_energy_J[2]),"angle reflection exchanges secondary energies");
        }
    fusion_three_body_cm_v1 zero{};
    require(fusion_c_three_equal_sequential_cm(alpha,0,0,0,&zero)==0 && zero.kinetic_energy_J[0]==0 && zero.kinetic_energy_J[2]==0,"all-cold endpoint");
    require(fusion_c_three_equal_sequential_cm(alpha,available,2*available,0,&zero)==PB11_STATUS_OUT_OF_RANGE && zero.kinetic_energy_J[0]==0,"forbidden intermediate state rejected");
    require(fusion_c_three_equal_sequential_cm(alpha,available,0,1.001,&zero)==PB11_STATUS_OUT_OF_RANGE,"invalid cosine rejected");
    std::cout<<"PASS: 45 sequential events, mass shells, intermediate invariant and angle reflection\n";
}

void birth_collision_inventory() {
    // Synthetic prescribed event source, not a validated pB spectrum. Its
    // purpose is to test the full map -> kinetic trial -> explicit ledger path.
    constexpr int n=240;
    constexpr double event_rate=1e10,available=8.68*mev,dt=.002;
    std::vector<double> edges(n+1),packet_e,packet_rate,birth(n),old(n,0),trial(n),escape(n,.3),diffusion(2*(n-1));
    for(int i=0;i<=n;++i) edges[i]=.01*kev*std::pow(1e6,double(i)/n);
    for(double cosine:{-.8,-.4,0.,.4,.8}) {
        fusion_three_body_cm_v1 event{};
        require(fusion_c_three_equal_sequential_cm(alpha,available,.3*available,cosine,&event)==0,"event generation");
        for(double e:event.kinetic_energy_J) {packet_e.push_back(e);packet_rate.push_back(event_rate/5);}
    }
    fusion_birth_mapping_v1 mapped{};
    require(fusion_c_map_birth_packets(n,edges.data(),static_cast<int>(packet_e.size()),packet_e.data(),packet_rate.data(),birth.data(),&mapped)==0,"prescribed event birth mapping");
    require(mapped.below_number_m3_s==0 && mapped.above_number_m3_s==0,"complete birth grid contains synthetic source");
    require(near(mapped.mapped_number_m3_s,3*event_rate) && near(mapped.mapped_energy_W_m3,event_rate*available),"three products per event and energy preserved before evolution");
    const double temperatures[]={5*kev,5*kev};
    const fusion_maxwellian_bath_v1 baths[]={{1e20,9.1093837139e-31,1,temperatures[0],15},
                                           {1e20,3.3435837768e-27,1,temperatures[1],15}};
    for(int b=0;b<2;++b) for(int i=0;i<n-1;++i) {
        fusion_coulomb_energy_v1 r{};
        require(fusion_c_coulomb_energy(edges[i+1],alpha,2,&baths[b],&r)==0,"birth test Coulomb coefficients");
        diffusion[b*(n-1)+i]=r.diffusion_J2_s;
    }
    long double accumulated_heat=0,escaped_n=0,escaped_e=0,thermal_n=0,thermal_e=0;
    constexpr int steps=100;
    for(int step=0;step<steps;++step) {
        double heat[2];fusion_kinetic_ledger_v1 ledger{};
        require(fusion_c_energy_fp_trial(n,2,dt,edges.data(),old.data(),temperatures,diffusion.data(),birth.data(),escape.data(),1000,trial.data(),heat,&ledger)==0,"birth-collision trial");
        accumulated_heat+=heat[0]+heat[1];escaped_n+=ledger.escaped_number_m3;escaped_e+=ledger.escaped_energy_J_m3;
        thermal_n+=ledger.thermalized_number_m3;thermal_e+=ledger.thermalized_energy_J_m3;
        old.swap(trial); // Explicit caller acceptance, no hidden state.
    }
    long double final_n=0,final_e=0;
    for(int i=0;i<n;++i) {final_n+=old[i];final_e+=old[i]*(static_cast<long double>(edges[i])+edges[i+1])/2;}
    const long double born_n=3*event_rate*dt*steps,born_e=event_rate*available*dt*steps;
    const double nerror=static_cast<double>((final_n+escaped_n+thermal_n-born_n)/born_n);
    const double eerror=static_cast<double>((final_e+escaped_e+thermal_e+accumulated_heat-born_e)/born_e);
    std::cout<<"Prescribed birth -> collisions 100-step number residual "<<nerror<<", energy residual "<<eerror<<'\n';
    require(std::abs(nerror)<1e-10 && std::abs(eerror)<1e-10,"independent cumulative particle/energy budget");
    require(accumulated_heat>0 && escaped_n>0,"source transfers energy to baths and explicit escape ledger");
}
}
int main() {
    try {mapping_and_spill();sequential_invariants();birth_collision_inventory();}
    catch(const std::exception &error) {std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}
    std::cout<<"PASS: conservative product birth mapping and sequential kinematics\n";
}
