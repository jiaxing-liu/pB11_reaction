module fusion_birth_floor_fortran
 use iso_c_binding
 implicit none
 private
 public::fusion_birth_floor_options_v1,fusion_birth_floor_ledger_v1,fusion_birth_floor_project
 type,bind(C)::fusion_birth_floor_options_v1
  real(c_double)::first_center_J,ion_kT_J,ion_energy_J_m3
  real(c_double)::max_center_over_ion_kT,max_ion_energy_fraction
 end type
 type,bind(C)::fusion_birth_floor_ledger_v1
  real(c_double)::born_number_m3(6),born_energy_J_m3(6)
  real(c_double)::mapped_number_m3(6),mapped_energy_J_m3(6)
  real(c_double)::ion_energy_correction_J_m3(6),energy_residual_J_m3(6)
  real(c_double)::remaining_ion_energy_J_m3
 end type
 interface
  function fusion_birth_floor_project(options,number,energy,ledger) &
    bind(C,name='fusion_c_birth_floor_project') result(status)
   import c_int,c_double,fusion_birth_floor_options_v1,fusion_birth_floor_ledger_v1
   type(fusion_birth_floor_options_v1),intent(in)::options
   real(c_double),intent(in)::number(6),energy(6)
   type(fusion_birth_floor_ledger_v1),intent(out)::ledger
   integer(c_int)::status
  end function
 end interface
end module
