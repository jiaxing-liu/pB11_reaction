#include "fusion_pb_birth.h"
#include <iostream>
#include <iomanip>
#include <vector>
int main(){int n;double A,q;std::cout<<std::setprecision(17);while(std::cin>>n>>A>>q){std::vector<double> e(n+1),b(n);for(auto&v:e)std::cin>>v;fusion_pb_birth_v1 r{};int s=fusion_c_pb_alpha0_grid(A,q,n,e.data(),b.data(),&r);std::cout<<s<<' '<<r.primary_alpha0_energy_J<<' '<<r.secondary_min_J<<' '<<r.secondary_max_J<<' '<<r.below_number<<' '<<r.below_energy_J<<' '<<r.above_number<<' '<<r.above_energy_J;for(auto v:b)std::cout<<' '<<v;std::cout<<'\n';}}
