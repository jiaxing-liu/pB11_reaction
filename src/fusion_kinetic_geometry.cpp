#include "fusion_kinetic_geometry.h"
#include "fusion_radial_transport.h"
#include "fusion_energy_work.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace {
constexpr int bad=PB11_STATUS_INVALID_ARGUMENT,num=PB11_STATUS_NUMERICAL_FAILURE;
bool store(long double x,double& out){
 if(!std::isfinite(x)||std::abs(x)>std::numeric_limits<double>::max())return false;
 out=static_cast<double>(x);return true;
}
}
extern "C" int fusion_c_kinetic_geometry_trial(int zones,int cells,double dt,
 const double*edges,const double*va,const double*vb,const double*advection,const double*compression,
 const double*k,const double*bs,const double*bt,const double*os,const double*ot,
 double*ts,double*tt,fusion_transport_ledger_v1*ledger){
 if(zones<1||zones>300||cells<1||cells>100000||12LL*cells*(2LL*zones+1)>50000000)return bad;
 const size_t nc=6*size_t(cells),nz=size_t(zones),nf=nz+1,total=nc*nz;
 if(ts)std::fill(ts,ts+total,0.);
 if(tt)std::fill(tt,tt+total,0.);
 if(ledger)std::fill(ledger,ledger+nz,fusion_transport_ledger_v1{});
 if(!edges||!va||!vb||!advection||!compression||!k||!bs||!bt||!os||!ot||!ts||!tt||!ledger)return bad;
 if(!std::isfinite(dt)||dt<=0)return bad;
 for(int i=0;i<=cells;++i)
  if(!std::isfinite(edges[i])||edges[i]<0||(i&&edges[i]<=edges[i-1]))return bad;
 for(size_t z=0;z<nz;++z)if(!std::isfinite(compression[z]))return bad;
 try {
  // Radial kernel layout is component-major with contiguous radial cells.
  // S/T each use the same explicitly supplied species/energy conductance.
  std::vector<double> old(2*total),boundary(4*nc),conductance(2*nc*nf);
  for(size_t component=0;component<nc;++component){
   for(size_t z=0;z<nz;++z){old[component*nz+z]=os[z*nc+component];old[(nc+component)*nz+z]=ot[z*nc+component];}
   for(size_t face=0;face<nf;++face){conductance[component*nf+face]=k[component*nf+face];
    conductance[(nc+component)*nf+face]=k[component*nf+face];}
   for(size_t b=0;b<2;++b){boundary[2*component+b]=bs[2*component+b];boundary[2*(nc+component)+b]=bt[2*component+b];}
  }
  std::vector<double> radial(2*total),faces(2*nc*nf);
  std::vector<fusion_radial_ledger_v1> radial_ledgers(2*nc);
  int status=fusion_c_radial_transport_trial(zones,int(2*nc),dt,va,vb,advection,
   conductance.data(),old.data(),boundary.data(),radial.data(),faces.data(),radial_ledgers.data());
  if(status!=PB11_STATUS_OK)return status;
  std::vector<double> candidates_s(total),candidates_t(total),population(cells),worked(cells),work_faces(cells+1);
  std::vector<fusion_transport_ledger_v1> result(nz);
  for(size_t z=0;z<nz;++z)for(size_t sp=0;sp<6;++sp){
   // Accumulate both kinetic components before converting extensive moments.
   long double spatial_n=0,spatial_u=0,work=0,lower_n=0,lower_u=0,upper_n=0,upper_u=0;
   for(size_t part=0;part<2;++part){
    bool nonzero=false;
    for(int cell=0;cell<cells;++cell){
     size_t component=part*nc+sp*size_t(cells)+size_t(cell);
     population[cell]=radial[component*nz+z];nonzero|=population[cell]!=0;
     const long double exchange=static_cast<long double>(faces[component*nf+z])-faces[component*nf+z+1];
     spatial_n+=exchange;spatial_u+=exchange*(static_cast<long double>(edges[cell])+edges[cell+1])/2;
    }
    fusion_energy_work_ledger_v1 energy{};
    if(nonzero&&compression[z]!=0){
     status=fusion_c_energy_work_trial(cells,dt,compression[z],edges,population.data(),worked.data(),work_faces.data(),&energy);
     if(status!=PB11_STATUS_OK)return status;
    }else worked=population; // Validated zero operator/population, no inferred source.
    auto& candidate=part==0?candidates_s:candidates_t;
    std::copy(worked.begin(),worked.end(),candidate.begin()+z*nc+sp*size_t(cells));
    work+=static_cast<long double>(energy.work_on_particles_J_m3)*vb[z];
    lower_n+=static_cast<long double>(energy.lower_number_m3)*vb[z];
    lower_u+=static_cast<long double>(energy.lower_energy_J_m3)*vb[z];
    upper_n+=static_cast<long double>(energy.upper_number_m3)*vb[z];
    upper_u+=static_cast<long double>(energy.upper_energy_J_m3)*vb[z];
   }
   auto& r=result[z];
   if(!store(spatial_n,r.spatial_number[sp])||!store(spatial_u,r.spatial_energy_J[sp])||
      !store(work,r.work_J[sp])||!store(lower_n,r.lower_number[sp])||!store(lower_u,r.lower_energy_J[sp])||
      !store(upper_n,r.upper_number[sp])||!store(upper_u,r.upper_energy_J[sp]))return num;
  }
  std::copy(candidates_s.begin(),candidates_s.end(),ts);std::copy(candidates_t.begin(),candidates_t.end(),tt);
  std::copy(result.begin(),result.end(),ledger);return PB11_STATUS_OK;
 }catch(...){return PB11_STATUS_EXCEPTION;}
}
