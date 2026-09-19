#define main original_coupled_fast_main
#include "test_fusion_coupled_fast.cpp"
#undef main
int main(){try{test_source_packets();std::cout<<"PASS source packet parity, shape, moments and failure atomicity\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
