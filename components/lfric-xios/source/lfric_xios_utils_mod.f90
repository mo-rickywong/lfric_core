!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Collection of small utility procedures for lfric-xios
!>
module lfric_xios_utils_mod

  use constants_mod,            only: i_def, r_def, str_def, str_long
  use file_mod,                 only: FILE_OP_OPEN, FILE_MODE_READ
  use lfric_ncdf_dims_mod,      only: lfric_ncdf_dims_type
  use lfric_ncdf_field_mod,     only: lfric_ncdf_field_type
  use lfric_ncdf_file_mod,      only: lfric_ncdf_file_type
  use lfric_mpi_mod,            only: global_mpi
  use lfric_xios_constants_mod, only: lx_year, lx_month, lx_day, lx_second
  use lfric_xios_field_mod,     only: lfric_xios_field_type
  use log_mod,                  only: log_event, log_scratch_space, &
                                      LOG_LEVEL_ERROR, LOG_LEVEL_INFO, &
                                      LOG_LEVEL_TRACE
  use mesh_mod,                 only: mesh_type
  use xios,                     only: xios_date, xios_duration,        &
                                      xios_get_time_origin,            &
                                      xios_get_year_length_in_seconds, &
                                      xios_date_convert_to_seconds,    &
                                      operator(<), operator(+)

  implicit none

  private
  public :: parse_date_as_xios, seconds_from_date, &
            set_prime_io_mesh, prime_io_mesh_is,   &
            read_time_data, duration_from_enum

  integer(i_def), private, allocatable :: prime_io_mesh_ids(:)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Interpret a string as an XIOS date object.
  !> @description Reccomended to use ISO 8601 Meteorological Format,
  !>              yyyy-mm-ddThh:mm:ssZ
  !> @param [in] date_str The string representation of the date
  !> @result     date_obj The xios_date object represented by the input string
  !>
  function parse_date_as_xios( date_str ) result( date_obj )

    implicit none

    character(*), intent(in) :: date_str
    type(xios_date)          :: date_obj

    integer(i_def) :: year, month, day
    integer(i_def) :: hours, minutes, seconds
    integer(i_def) :: offset_hrs, offset_mins, offset_sign

    integer(i_def) :: date_index
    integer(i_def) :: time_index

    character(str_def) :: str_date
    character(str_def) :: str_time
    character(str_def) :: str_zone

    logical ::  time_zone_out_of_range

    date_index   = index(date_str,'-')
    time_index   = index(date_str,':')

    str_date = trim(adjustl(date_str(:date_index+5)))
    str_time = trim(adjustl(date_str(time_index-2:time_index+5)))
    str_zone = trim(adjustl(date_str(time_index:)))

    if (scan(str_zone,"+-") == 0) then
      str_zone = 'Z'
    else
      str_zone = trim(adjustl(date_str(time_index+6:)))
    end if

    ! Read the date
    date_index = index(str_date,'-')
    read(str_date(           1:date_index-1),*) year
    read(str_date(date_index+1:date_index+2),*) month
    read(str_date(date_index+4:date_index+5),*) day

    ! Read the time
    read(str_time(1:2),*) hours
    read(str_time(4:5),*) minutes
    read(str_time(7:8),*) seconds

    ! Read the time-zone
    if (scan(str_zone,'Z') == 0) then
      ! Time zone is non-Z format
      read(str_zone(2:3),*) offset_hrs
      read(str_zone(5:6),*) offset_mins
      if (index(str_zone,'-') == 0) then
        offset_sign = 1
      else
        offset_sign = -1
      end if
    else
      ! Time zone is UTC format
      offset_hrs  = 0
      offset_mins = 0
      offset_sign = 1
    end if

    ! Range check on time zone
    time_zone_out_of_range = .false.
    if (offset_mins > 0) then

      select case (offset_hrs*offset_sign)
      case(14:,:-12)
        time_zone_out_of_range = .true.
      end select

    else if (offset_mins == 0) then

      select case (offset_hrs*offset_sign)
      case(15:,:-13)
        time_zone_out_of_range = .true.
      end select

    end if

    if (time_zone_out_of_range) then
      call log_event('Time-zone out-of-range (-12:00 to 14:00)', &
                     log_level_error)
    end if

    date_obj = xios_date( year, month, day, hours, minutes, seconds )

  end function parse_date_as_xios

  !> @brief  Wrapper around xios_date_convert_to_seconds to nullify XIOS bug.
  !>
  !> @param[in] date            The input xios_date object
  !> @result    date_in_seconds The resulting seconds converted from "date"
  function seconds_from_date(date) result(date_in_seconds)

    implicit none

    type(xios_date), intent(in) :: date

    type(xios_date)        :: time_origin
    integer(i_def)         :: year_diff
    real(r_def)            :: date_in_seconds

    integer(i_def), parameter :: length360d = 31104000

    call xios_get_time_origin(time_origin)

    ! Get time in seconds from XIOS dates - due to a bug in XIOS, non-360day
    ! calendars do not return the correct values around the time origin so a
    ! workaround is implemented below. Also due to a bug in XIOS calendar types
    ! cannot be identified except by the number of seconds per year
    if ( date < time_origin .and. &
          xios_get_year_length_in_seconds(date%year) /= length360d ) &
      then
      year_diff = date%year - time_origin%year
      date_in_seconds = real(xios_date_convert_to_seconds(date), r_def) + &
        ( real(xios_get_year_length_in_seconds(date%year), r_def) * &
          real(year_diff, r_def) )
    else
      date_in_seconds = real(xios_date_convert_to_seconds(date), r_def)
    end if

  end function seconds_from_date

  !> @brief Registers a mesh to be used as the primary I/O mesh
  !> @param[in] mesh  The mesh object to be registered
  subroutine set_prime_io_mesh( mesh )

    implicit none

    type(mesh_type), intent(in) :: mesh

    integer(i_def), allocatable :: mesh_id_list(:)

    ! Set up array of ints to hold mesh ids, bring in previous mesh IDs if
    ! already present and add the new mesh ID to the new array
    if (allocated(prime_io_mesh_ids)) then
      allocate(mesh_id_list(size(prime_io_mesh_ids) + 1))
      mesh_id_list(1:size(prime_io_mesh_ids)) = prime_io_mesh_ids
      mesh_id_list(size(prime_io_mesh_ids)+1) = mesh%get_id()
    else
      allocate(mesh_id_list(1))
      mesh_id_list(1) = mesh%get_id()
    end if

    ! Make the new array the main array
    call move_alloc(mesh_id_list, prime_io_mesh_ids)

  end subroutine set_prime_io_mesh

  function prime_io_mesh_is( mesh ) result( mesh_is_prime_io )

    implicit none

    type(mesh_type), intent(in) :: mesh

    logical :: mesh_is_prime_io

    if (.not. allocated(prime_io_mesh_ids) ) then
      mesh_is_prime_io = .false.
    else
      mesh_is_prime_io = (any(prime_io_mesh_ids == mesh%get_id()))
    end if

  end function prime_io_mesh_is


  !>  @brief  Read time data into time axis using NetCDF tools
  !!
  !> @param[in] file_path  The path to the file containing the data controlled
  !!                       by this object
  function read_time_data(file_path) result(time_data)

    implicit none

    character(len=*),      intent(in)   :: file_path

    type(xios_date), allocatable :: time_data(:)

    ! Local variables for XIOS interface
    integer(i_def)                :: n_t, t, n_t_buf(1)
    real(r_def), allocatable      :: input_data(:)
    character(str_def)            :: time_units, ref_date_str, var_id, dim_id
    character(str_def)            :: time_meta_buf(2)
    character(str_long)           :: unit_attr
    character(str_def), parameter :: valid_units(4) = &
                                  (/'seconds', 'days   ', 'hours  ', 'months '/)
    type(xios_duration)           :: ref_time, mean_time, month_duration
    type(xios_date)               :: ref_date
    integer(i_def),  parameter :: len_date = 18 ! Length of CF date is 18 characters
    integer(i_def),  parameter :: len_delim = 7 ! Delimiter between unit and
                                                    ! date is 7 characters

    ! NetCDF reading variables
    type(lfric_ncdf_file_type)  :: file_ncdf
    type(lfric_ncdf_dims_type)  :: time_dim
    type(lfric_ncdf_field_type) :: time_var

    if (global_mpi%get_comm_rank() == 0) then
      file_ncdf = lfric_ncdf_file_type( trim(file_path)//".nc", &
                                        open_mode=FILE_OP_OPEN, &
                                        io_mode=FILE_MODE_READ )

      ! Some JULES surface ancils have non-CF-compliant time representation
      ! so we need to account for that for the time being
      if (file_ncdf%contains_var("time")) then
        var_id = "time"
        dim_id = "time"
      else if (file_ncdf%contains_var("month_number")) then
        var_id = "month_number"
        dim_id = "month_number"
      else
        call log_event( "Invalid representation of time in file ["// &
                        trim(file_path)//"]", log_level_error)
      end if

      ! Get size of time axis from file
      time_dim = lfric_ncdf_dims_type(trim(dim_id), file_ncdf)
      n_t = time_dim%get_size()
      allocate( input_data( n_t ) )
      allocate( time_data( n_t ) )

      ! Read the time data from the ancil file
      time_var = lfric_ncdf_field_type(trim(var_id), file_ncdf)
      call time_var%read_data(input_data)

      ! Read time units and reference date from file
      if (trim(var_id) == "time") then
        unit_attr = time_var%get_char_attribute("units")
        time_units = unit_attr( 1 : len(trim(unit_attr))-len_date-len_delim )
        if ( .not. any( valid_units == trim(adjustl(time_units)) ) ) then
          write( log_scratch_space,'(A,A)' ) "Invalid units of ["//trim(time_units)// &
                                            "] for time axis in file: "// trim(file_path)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        ref_date_str = unit_attr( len(trim(unit_attr))-len_date : len(trim(unit_attr)) )
      else if (trim(var_id) == "month_number") then
        ! Non CF files are given the bare minimum treatment
        time_units = "months"
        ref_date_str = "1970-01-01 00:00:00"
      end if

      call file_ncdf%close_file()

    end if ! comm_rank == 0

    ! Broadcast data across MPI ranks
    n_t_buf(1) = n_t
    time_meta_buf = [ref_date_str, time_units]
    call global_mpi%broadcast(n_t_buf, 1, 0)
    if (global_mpi%get_comm_rank() /= 0) allocate(input_data(n_t_buf(1)))
    call global_mpi%broadcast(input_data, n_t_buf(1), 0)
    call global_mpi%broadcast(time_meta_buf, size(time_meta_buf, 1)*str_def, 0)

    ref_date_str = time_meta_buf(1)
    time_units = time_meta_buf(2)
    ref_date = parse_date_as_xios(trim(adjustl(ref_date_str)))

    ! Convert input time data to xios_date type
    do t = 1, size(input_data)
      ref_time = xios_duration(0, 0, 0, 0, 0, 0)
      if ( trim(adjustl(time_units)) == "seconds" ) then
        ref_time%second = input_data(t)
      else if ( trim(adjustl(time_units)) == "hours" ) then
        ref_time%hour = input_data(t)
      else if ( trim(adjustl(time_units)) == "days" ) then
        ref_time%day = input_data(t)
      else if ( trim(adjustl(time_units)) == "months" ) then
        ! Offset months backwards to account for monthly mean
        ref_time%month = input_data(t) - 1
      end if
      time_data(t) = ref_date + ref_time
    end do

    ! Correct "months" data to be monthly mean - centred on middle of month.
    ! This can't be done above as the conversion to seconds can only be done
    ! with xios_date objects, not xios_durations
    if ( time_units == "months" ) then
      mean_time = xios_duration(0, 0, 0, 0, 0, 0)
      month_duration = xios_duration(0, 1, 0, 0, 0, 0)
      do t = 1, n_t_buf(1)
        mean_time%second = ( &
          seconds_from_date(time_data(t)+month_duration) - &
          seconds_from_date(time_data(t)) ) / 2
        time_data(t) = time_data(t) + mean_time
      end do
    end if

  end function read_time_data


  !> @brief  Construct an XIOS duration object from an integer time enum
  !!
  !> @param[in] time_enum  The integer time enum to be converted to an XIOS duration object
  function  duration_from_enum(time_enum) result(duration)

    implicit none

    integer(i_def), intent(in) :: time_enum
    type(xios_duration)        :: duration

    duration = xios_duration(0, 0, 0, 0, 0, 0)

    if (mod(time_enum, lx_year) == 0) then
      duration%year = time_enum / lx_year
    else if (mod(time_enum, lx_month) == 0) then
      duration%month = time_enum / lx_month
    else if (mod(time_enum, lx_day) == 0) then
      duration%day = time_enum / lx_day
    else if (mod(time_enum, lx_second) == 0) then
      duration%second = time_enum / lx_second
    else
      call log_event( "Unable to construct XIOS duration from time enum", &
                      log_level_error )
    end if

  end function duration_from_enum

end module lfric_xios_utils_mod
