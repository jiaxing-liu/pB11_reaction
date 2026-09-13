program consumer
use iso_c_binding, only: c_double, c_int
use fusion_rate_model_fortran
implicit none
type(fusion_rate_model_v1) :: r
integer(c_int) :: status
call fusion_thermal_pair_maxwellian_model(3_c_int, FUSION_ENDPOINT_S, &
 FUSION_PB_LOW_TB,3.3435837768e-27_c_double,5.0073567512e-27_c_double, &
 3*1.602176634e-16_c_double,3*1.602176634e-16_c_double,r,status)
if(status/=0_c_int .or. r%total%resolved_reactivity_m3_s<=0.0_c_double) stop 1
end program
