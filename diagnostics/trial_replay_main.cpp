#include "fusion_capture_trial.h"
#include <cmath>
#include <iomanip>
#include <iostream>
int main(int argc,char **argv) {
 if(argc!=2){std::cerr<<"usage: fusion_trial_replay snapshot.pbtrial\n";return 2;}
 fusion_capture_replay_report_v1 report{};
 int status=fusion_capture_replay_file_v1(argv[1],128ULL*1024*1024,&report);
 std::cout<<"{\"diagnostic_status\":"<<status<<",\"replay_performed\":"<<(status==0?"true":"false")
  <<",\"physics_accepted\":false";
 if(status==0){std::cout<<",\"original_status\":"<<report.original_status<<",\"replay_status\":"<<report.replay_status
   <<",\"status_matches\":"<<(report.status_matches?"true":"false")<<",\"entry_point\":"<<report.entry_point
   <<",\"host_zone\":"<<report.host_zone<<",\"cells\":"<<report.cells<<",\"host_time_s\":";
  if(std::isfinite(report.host_time_s))std::cout<<std::setprecision(17)<<report.host_time_s;else std::cout<<"null";
  std::cout<<",\"dt_s\":";
  if(std::isfinite(report.dt_s))std::cout<<std::setprecision(17)<<report.dt_s;else std::cout<<"null";
  std::cout<<",\"kernel_identity\":\""<<report.kernel_identity<<"\",\"output_sha256\":\""<<report.output_sha256<<"\"";
 }
 std::cout<<"}\n";
 return status?1:(report.status_matches?0:2);
}
