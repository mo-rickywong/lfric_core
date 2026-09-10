module driver_log_mod

use constants_mod,        only: i_def, l_def
use convert_to_upper_mod, only: convert_to_upper
use config_mod,           only: config_type
use lfric_mpi_mod,        only: lfric_comm_type
use log_mod,              only: log_event,          &
                                log_set_level,      &
                                log_scratch_space,  &
                                initialise_logging, &
                                finalise_logging,   &
                                log_level_always,   &
                                log_level_error,    &
                                log_level_warning,  &
                                log_level_info,     &
                                log_level_debug,    &
                                log_level_trace

use logging_config_mod, only: key_from_run_log_level, &
                              RUN_log_level_error,    &
                              RUN_log_level_info,     &
                              RUN_log_level_debug,    &
                              RUN_log_level_trace,    &
                              RUN_log_level_warning


implicit none

public :: init_logger, final_logger
private

contains

!> @brief Initialises the logging system from a namelist.
!>
!> @param[in] config       Application namelist configuration object
!> @param[in] communicator MPI communicator to use for logging.
!> @param[in] program_name Identifies the running program.
!>
subroutine init_logger(config, communicator, program_name)

  implicit none

  type(config_type),     intent(in) :: config
  type(lfric_comm_type), intent(in) :: communicator
  character(len=*),      intent(in) :: program_name

  integer(i_def) :: log_level
  integer(i_def) :: run_log_level
  logical(l_def) :: log_to_rank_zero_only

  run_log_level         = config%logging%run_log_level()
  log_to_rank_zero_only = config%logging%log_to_rank_zero_only()

  call initialise_logging( communicator%get_comm_mpi_val(), program_name, &
                           log_to_rank_zero_only=log_to_rank_zero_only)

  select case (run_log_level)
  case( RUN_log_level_error )
    log_level = log_level_error
  case( RUN_log_level_warning )
    log_level = log_level_warning
  case( RUN_log_level_info )
    log_level = log_level_info
  case( RUN_log_level_debug )
    log_level = log_level_debug
  case( RUN_log_level_trace )
    log_level = log_level_trace
  case default
    call log_event( "Invalid option for run_log_level", log_level_error )
  end select

  call log_set_level( log_level )

  write(log_scratch_space,'(A)')                              &
      'Runtime message logging severity set to log level: '// &
      convert_to_upper(key_from_run_log_level(run_log_level))
  call log_event( log_scratch_space, log_level_always )

end subroutine init_logger

subroutine final_logger(program_name)

  implicit none

  character(len=*), intent(in) :: program_name

  ! Final logging before infrastructure is destroyed
  call log_event( program_name//' application completed.', log_level_always )

  ! Finalise the logging system
  call finalise_logging()

end subroutine final_logger

end module driver_log_mod
