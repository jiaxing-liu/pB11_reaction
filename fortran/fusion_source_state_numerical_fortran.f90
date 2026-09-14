module fusion_source_state_numerical_fortran
  !! Direct ISO_C_BINDING declarations for the signed numerical source-state
  !! account.  The existing source-state and source-volume modules retain the
  !! owner lifecycle wrappers and the BIND(C) ledger layouts; this small layer
  !! only exposes the four new C entry points without copying those layouts.
  !!
  !! Every opaque state, array, ledger, and output scalar is a c_ptr passed by
  !! value.  The ticket is the C uint64_t value, while snapshot epoch is a
  !! uint64_t * and is consequently also represented by a c_ptr.  Callers
  !! should pass c_loc of TARGET actuals, matching the C ABI exactly.
  use, intrinsic :: iso_c_binding, only : c_double, c_int, c_int64_t, c_ptr
  use fusion_source_state_fortran, only : fusion_source_ledger_v1
  use fusion_source_volume_fortran, only : fusion_transport_ledger_v1
  implicit none
  private

  public :: fusion_source_ledger_v1
  public :: fusion_transport_ledger_v1
  public :: c_fusion_source_state_stage_numerical
  public :: c_fusion_source_state_stage_volume_numerical
  public :: c_fusion_source_state_snapshot_numerical
  public :: c_fusion_source_state_snapshot_volume_numerical
  public :: c_fusion_c_coupled_numerical_increment

  interface
     function c_fusion_source_state_stage_numerical( &
          state, ticket, trial_s, trial_t, source, inert, numerical) &
          bind(C, name="fusion_c_source_state_stage_numerical") result(status)
       import :: c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       type(c_ptr), value :: source
       type(c_ptr), value :: inert
       type(c_ptr), value :: numerical
       integer(c_int) :: status
     end function c_fusion_source_state_stage_numerical

     function c_fusion_source_state_stage_volume_numerical( &
          state, ticket, trial_s, trial_t, trial_volume, source_volume, &
          source, inert, transport, numerical) bind(C, &
          name="fusion_c_source_state_stage_volume_numerical") result(status)
       import :: c_double, c_int, c_int64_t, c_ptr
       type(c_ptr), value :: state
       integer(c_int64_t), value :: ticket
       type(c_ptr), value :: trial_s
       type(c_ptr), value :: trial_t
       real(c_double), value :: trial_volume
       real(c_double), value :: source_volume
       type(c_ptr), value :: source
       type(c_ptr), value :: inert
       type(c_ptr), value :: transport
       type(c_ptr), value :: numerical
       integer(c_int) :: status
     end function c_fusion_source_state_stage_volume_numerical

     function c_fusion_source_state_snapshot_numerical( &
          state, accepted_s, accepted_t, cumulative, inert, numerical, &
          accepted_time, epoch) bind(C, &
          name="fusion_c_source_state_snapshot_numerical") result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: state
       type(c_ptr), value :: accepted_s
       type(c_ptr), value :: accepted_t
       type(c_ptr), value :: cumulative
       type(c_ptr), value :: inert
       type(c_ptr), value :: numerical
       type(c_ptr), value :: accepted_time
       type(c_ptr), value :: epoch
       integer(c_int) :: status
     end function c_fusion_source_state_snapshot_numerical

     function c_fusion_source_state_snapshot_volume_numerical( &
          state, accepted_s, accepted_t, cumulative, inert, transport, &
          reference_volume, accepted_volume, accepted_time, epoch, numerical) &
          bind(C, name="fusion_c_source_state_snapshot_volume_numerical") &
          result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: state
       type(c_ptr), value :: accepted_s
       type(c_ptr), value :: accepted_t
       type(c_ptr), value :: cumulative
       type(c_ptr), value :: inert
       type(c_ptr), value :: transport
       type(c_ptr), value :: reference_volume
       type(c_ptr), value :: accepted_volume
       type(c_ptr), value :: accepted_time
       type(c_ptr), value :: epoch
       type(c_ptr), value :: numerical
       integer(c_int) :: status
     end function c_fusion_source_state_snapshot_volume_numerical

     function c_fusion_c_coupled_numerical_increment( &
          ledger, inert_heat, numerical, out) bind(C, &
          name="fusion_c_coupled_numerical_increment") result(status)
       import :: c_int, c_ptr
       type(c_ptr), value :: ledger
       type(c_ptr), value :: inert_heat
       type(c_ptr), value :: numerical
       type(c_ptr), value :: out
       integer(c_int) :: status
     end function c_fusion_c_coupled_numerical_increment
  end interface

end module fusion_source_state_numerical_fortran
