!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> A module to perform linear time interpolation of lists or evaluate
!> functions at the current time.
module sci_time_interp_mod

    use constants_mod,        only: i_def, r_def, PI
    use log_mod,              only: log_event,      &
                                    log_level_info, &
                                    log_level_error

    implicit none

    public :: time_interpolate_list, &
              sinusoidal_function,   &
              diurnal_function

contains

!> @brief Interpolate a list of value at the current time
!> @details This procedure uses the current time to interpolate
!!          a given list of values and times. If the current time falls
!!          outside of the range of the provided times list, then the
!!          max/min value is used instead.
!> @param[inout] profile_now The resulting interpolated value at the
!!                           current time
!> @param[in] times      List of times, corresponding to values in prof_vary,
!!                       these times should be specifed in the units defined
!!                       by time_units.
!> @param[in] prof_vary  List of values to be interpolated
!> @param[in] time_now   Current time
subroutine time_interpolate_list( profile_now, times, &
                                  prof_vary, time_now )

  implicit none

  real(r_def),               intent(in)  :: time_now
  real(r_def), dimension(:), intent(in)  :: times
  real(r_def), dimension(:), intent(in)  :: prof_vary
  real(r_def),               intent(out) :: profile_now

  integer(i_def) :: time_index
  integer(i_def) :: n_times
  real(r_def)    :: tau

  n_times = size( times )

  if ( size(prof_vary) /= n_times ) then
    call log_event( 'prof_vary not same size as times', log_level_error )
  end if

  ! Perform linear interpolation
  if ( n_times == 1 ) then
    profile_now = prof_vary(1)

  else
    if ( time_now <= times(1) ) then
      profile_now = prof_vary(1)

    else if ( time_now >= times( n_times ) ) then
      profile_now = prof_vary(n_times)

    else
      do time_index = 1, n_times - 1

        if ( times( time_index + 1 ) > time_now ) then

          tau = ( time_now - times( time_index ) ) /              &
                ( times( time_index + 1 ) - times( time_index ) )

          profile_now = ( 1 - tau )                   &
                * prof_vary( ( time_index - 1) + 1 )  &
                + tau * prof_vary( time_index + 1 )
          exit
        end if

      end do
    end if
  end if

end subroutine time_interpolate_list

!> @brief Evaluate the value of a sine function at a certain time
!> @details Evaluates a sinusoidal function, defined by its
!!          amplitude, period and phase, at a particular time.
!> @param[out] profile_now The function value at time_now.
!> @param[in]  time_now    The time at which the function is evaluated.
!> @param[in]  amplitude   The amplitude of the function
!> @param[in]  period      The period of the function
!> @param[in]  phase       The phase of the function
subroutine sinusoidal_function( profile_now, time_now, &
                                amplitude, period, phase)

  implicit none

  real(r_def), intent(out) :: profile_now
  real(r_def), intent(in)  :: time_now
  real(r_def), intent(in)  :: amplitude
  real(r_def), intent(in)  :: period

  real(r_def), optional, intent(in) :: phase

  if ( present( phase ) ) then
    profile_now = amplitude * sin( 2.0_r_def * PI / period * time_now + phase )
  else
    profile_now = amplitude * sin( 2.0_r_def * PI / period * time_now )
  end if

end subroutine sinusoidal_function

!> @brief Evaluate the value of a diurnal function at a certain time
!> @details Evaluates a diurnal function at a particular time.
!> @param[out] profile_now       The function value at time_now.
!> @param[in]  time_now          The time at which the function is evaluated.
!> @param[in]  amplitude         The amplitude of the function
!> @param[in]  time_of_max_value The time, from the start of the run, at which
!!                               the function peaks.
!> @param[in]  length_of_day     The duration of daylight
subroutine diurnal_function( profile_now, time_now,        &
                             amplitude, time_of_max_value, &
                             length_of_day )

  implicit none

  real(r_def), intent(out) :: profile_now
  real(r_def), intent(in) :: time_now
  real(r_def), intent(in) :: amplitude
  real(r_def), intent(in) :: time_of_max_value
  real(r_def), intent(in) :: length_of_day

  real(r_def) :: xfact

  xfact  = cos( PI * ( time_of_max_value - time_now ) / &
                length_of_day )

  if ( xfact <= 0.0_r_def ) xfact = 0.0_r_def

  profile_now = amplitude * xfact ** 1.3_r_def

end subroutine diurnal_function

end module sci_time_interp_mod
