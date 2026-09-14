module fusion_coupled_sources_floor_fortran
  !! Direct ISO_C_BINDING declaration for the provisional coupled birth-floor ABI.
  !!
  !! The C entry point accepts borrowed pointers.  In particular, the optional
  !! thermal-table array, beam-entry array, and inert-ion array are represented
  !! as c_ptr,value and are passed as c_null_ptr when absent.  No Fortran array
  !! descriptor or hidden optional-argument flag is part of this interface.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr
  use fusion_coupled_thermal_fortran, only : &
       fusion_coupled_thermal_options_v1, fusion_inert_ion_v1, &
       fusion_coupled_thermal_v1, fusion_handoff_diagnostics_v1, &
       FUSION_COUPLED_THERMAL_SPECIES, FUSION_COUPLED_THERMAL_CHANNELS, &
       FUSION_COUPLED_THERMAL_MAX_CELLS, FUSION_COUPLED_THERMAL_MAX_INERT, &
       FUSION_COUPLED_THERMAL_BASE_BATHS, PB11_STATUS_OK, &
       PB11_STATUS_NULL_OUTPUT, PB11_STATUS_INVALID_ARGUMENT, &
       PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_NUMERICAL_FAILURE, &
       PB11_STATUS_EXCEPTION, PB11_STATUS_UNKNOWN_METHOD
  use fusion_coupled_fast_fortran, only : fusion_fast_target_options_v1, &
       FUSION_COUPLED_FAST_CHANNELS
  use fusion_coupled_sources_fortran, only : &
       fusion_beam_table_entry_v1, fusion_beam_table_usage_v1
  use fusion_birth_floor_fortran, only : fusion_birth_floor_ledger_v1
  implicit none
  private

  ! This is the exact C layout from fusion_coupled_sources.h: two doubles,
  ! with no integer flags or Fortran-only bookkeeping.
  type, bind(C), public :: fusion_coupled_floor_limits_v1
     real(c_double) :: max_center_over_ion_kT
     real(c_double) :: max_ion_energy_fraction
  end type fusion_coupled_floor_limits_v1

  ! Re-export the ABI types used by the direct call so a fixture can obtain all
  ! declarations from this one module while their layouts remain owned by the
  ! existing external bindings.
  public :: fusion_coupled_thermal_options_v1, fusion_fast_target_options_v1
  public :: fusion_beam_table_entry_v1, fusion_beam_table_usage_v1
  public :: fusion_inert_ion_v1, fusion_coupled_thermal_v1
  public :: fusion_handoff_diagnostics_v1, fusion_birth_floor_ledger_v1
  public :: FUSION_COUPLED_THERMAL_SPECIES, FUSION_COUPLED_THERMAL_CHANNELS
  public :: FUSION_COUPLED_THERMAL_MAX_CELLS, FUSION_COUPLED_THERMAL_MAX_INERT
  public :: FUSION_COUPLED_THERMAL_BASE_BATHS, FUSION_COUPLED_FAST_CHANNELS
  public :: PB11_STATUS_OK, PB11_STATUS_NULL_OUTPUT, PB11_STATUS_INVALID_ARGUMENT
  public :: PB11_STATUS_OUT_OF_RANGE, PB11_STATUS_NUMERICAL_FAILURE
  public :: PB11_STATUS_EXCEPTION, PB11_STATUS_UNKNOWN_METHOD

  public :: c_fusion_c_coupled_sources_floor_trial

  interface
     function c_fusion_c_coupled_sources_floor_trial( &
          dt_s, options, fast_options, thermal_tables, beam_table_count, &
          beam_tables, effective_charge, cells, edges, thermal_number, &
          electron_energy, ion_energy, electron_density, &
          thermal_charge_squared, inert_count, inert, coulomb_logs, old_s, &
          old_t, external_birth, escape, new_thermal_number, new_s, new_t, &
          out, diagnostics, usage, floor_limits, floor_ledger) &
          bind(C, name="fusion_c_coupled_sources_floor_trial") result(status)
       import :: c_double, c_int, c_ptr
       real(c_double), value :: dt_s
       type(c_ptr), value :: options
       type(c_ptr), value :: fast_options
       type(c_ptr), value :: thermal_tables
       integer(c_int), value :: beam_table_count
       type(c_ptr), value :: beam_tables
       integer(c_int), value :: effective_charge
       integer(c_int), value :: cells
       type(c_ptr), value :: edges
       type(c_ptr), value :: thermal_number
       real(c_double), value :: electron_energy
       real(c_double), value :: ion_energy
       real(c_double), value :: electron_density
       type(c_ptr), value :: thermal_charge_squared
       integer(c_int), value :: inert_count
       type(c_ptr), value :: inert
       type(c_ptr), value :: coulomb_logs
       type(c_ptr), value :: old_s
       type(c_ptr), value :: old_t
       type(c_ptr), value :: external_birth
       type(c_ptr), value :: escape
       type(c_ptr), value :: new_thermal_number
       type(c_ptr), value :: new_s
       type(c_ptr), value :: new_t
       type(c_ptr), value :: out
       type(c_ptr), value :: diagnostics
       type(c_ptr), value :: usage
       type(c_ptr), value :: floor_limits
       type(c_ptr), value :: floor_ledger
       integer(c_int) :: status
     end function c_fusion_c_coupled_sources_floor_trial
  end interface

end module fusion_coupled_sources_floor_fortran
