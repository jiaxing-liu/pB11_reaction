module fusion_fast_moments_fortran
 use iso_c_binding, only: c_double,c_int,c_ptr,c_loc
 implicit none
 private
 type, bind(C), public :: fusion_fast_moments_v1
  real(c_double) :: number_m3(6),energy_J_m3(6)
  real(c_double) :: charge_number_m3,charge_squared_number_m3,pressure_Pa
 end type
 public :: fusion_fast_moments
 interface
  integer(c_int) function c_moments(cells,edges,s,t,out) bind(C,name="fusion_c_fast_moments")
   import c_int,c_ptr
   integer(c_int), value :: cells
   type(c_ptr), value :: edges,s,t,out
  end function
 end interface
contains
 subroutine fusion_fast_moments(edges,s,t,out,status)
  real(c_double), intent(in),target,contiguous :: edges(:),s(:),t(:)
  type(fusion_fast_moments_v1), intent(out),target :: out
  integer(c_int), intent(out) :: status
  integer(c_int) :: cells
  out%number_m3=0;out%energy_J_m3=0;out%charge_number_m3=0
  out%charge_squared_number_m3=0;out%pressure_Pa=0;status=2
  if(size(edges)<2.or.size(edges)>100001) return
  cells=int(size(edges)-1,c_int)
  if(size(s)/=6*cells.or.size(t)/=6*cells) return
  status=c_moments(cells,c_loc(edges(1)),c_loc(s(1)),c_loc(t(1)),c_loc(out))
 end subroutine
end module
