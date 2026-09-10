!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module empty_mod

  use constants_mod, only : i_def
  use lfric_mpi_mod, only : global_mpi
  use log_mod,       only : log_scratch_space, log_event, log_level_error

  implicit none

  private

  integer(i_def) :: local_rank = -1
  integer(i_def), parameter :: temporary_unit = 3

contains

end module empty_mod
