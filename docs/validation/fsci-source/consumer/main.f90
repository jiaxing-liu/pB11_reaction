program consumer
 use iso_c_binding
 use fusion_spectrum_fortran, only: fusion_alpha_spectrum_v1
 use fusion_fsci_fortran
 implicit none
 real(c_double), parameter :: mev=1.602176634e-13_c_double
 real(c_double) :: edges(21),birth(20)
 type(fusion_alpha_spectrum_v1) :: result
 integer(c_int) :: rc
 integer :: i
 do i=1,21
  edges(i)=real(i-1,c_double)*.5_c_double*mev
 enddo
 call fusion_alpha_spectrum_model_grid(2_c_int,1_c_int,8.84_c_double*mev, &
  .001_c_double*mev,.76_c_double,4.2_c_double,16_c_int,16_c_int,20_c_int,edges,birth,result,rc)
 if(rc/=0)stop 1
 if(abs(result%mapped_number+result%below_number+result%above_number-3)>1e-10_c_double)stop 2
end program
