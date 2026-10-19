! Independent read-only geometry diagnostic; no host state or indexing.
module fusion_flux_volume_fortran
 use,intrinsic::iso_c_binding,only:c_int,c_double
 implicit none
 private
 public::fusion_flux_volume_diagnostic_v1
 public::fusion_c_flux_core_enclosed_volume_diagnostic_v1
 type,bind(C)::fusion_flux_volume_diagnostic_v1
  integer(c_int)::version,angular_nodes,boundary_orientation,reserved
  real(c_double)::signed_volume_m3,absolute_volume_m3,min_abs_boundary_jacobian_m3
 end type
 interface
  function fusion_c_flux_core_enclosed_volume_diagnostic_v1(n,m,x,a,d,f,c,k,q,out) &
    bind(C,name='fusion_c_flux_core_enclosed_volume_diagnostic_v1')result(status)
   import c_int,c_double,fusion_flux_volume_diagnostic_v1
   integer(c_int),value::n,m,k
   real(c_double),intent(in)::x(*),a(*),d(*),f(*),c(*)
   real(c_double),value::q
   type(fusion_flux_volume_diagnostic_v1),intent(out)::out
   integer(c_int)::status
  end function
 end interface
end module
