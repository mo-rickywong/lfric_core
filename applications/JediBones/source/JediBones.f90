!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp JediBones program

!> @brief Main program used to illustrate how to write LFRic applications.

!> @details Calls init, step and finalise routines from a driver module

program JediBones

  use cli_mod,                 only: parse_command_line
  use constants_mod,           only: precision_real
  use driver_collections_mod,  only: init_collections, final_collections
  use driver_comm_mod,         only: init_comm, final_comm
  use driver_config_mod,       only: init_config, final_config
  use driver_log_mod,          only: init_logger, final_logger
  use driver_modeldb_mod,      only: modeldb_type
  use driver_time_mod,         only: init_time, final_time
  use lfric_mpi_mod,           only: global_mpi
  use log_mod,                 only: log_event,       &
                                     log_level_trace, &
                                     log_scratch_space

  use JediBones_mod, only: JediBones_S1_namelists,&
                           JediBones_S2_namelists

  use JediBones_driver_mod, only: initialise, step, finalise

  use season_cast_mod, only: roll_cast

  implicit none

  ! The technical and scientific state
  type(modeldb_type) :: modeldb_s1
  type(modeldb_type) :: modeldb_s2

  character(*), parameter :: S1_Title = "DustyBones"
  character(*), parameter :: S2_Title = "FleshyBones"

  character(*), parameter   :: program_name = "JediBones"
  character(:), allocatable :: S1_nml_file
  character(:), allocatable :: S2_nml_file

  call parse_command_line( S1_nml_file, filename2=S2_nml_file )

  call modeldb_S1%config%initialise(S1_Title)
  call modeldb_S2%config%initialise(S2_Title)

  modeldb_S1%mpi => global_mpi
  modeldb_S2%mpi => global_mpi

  ! Both modeldbs point to a global_mpi
  ! So just use the on in season 1.
  ! May need to speak to Mike about this
  call init_comm( "JediBones", modeldb_s1 )
! call init_comm( "JediBones", modeldb_s2 ) ???????

  ! Read in the 2nd file first as settings in module scope configs will be the last
  ! config read in.
  call init_config( S2_nml_file, JediBones_S2_namelists, config=modeldb_s2%config )
  call init_config( S1_nml_file, JediBones_S1_namelists, config=modeldb_s1%config )

  call init_logger( modeldb_s1%mpi%get_comm(), program_name )

  write(log_scratch_space,'(A)')                          &
      'Application built with '// trim(precision_real) // &
      '-bit real numbers.'
  call log_event( log_scratch_space, log_level_trace )

  call roll_cast( modeldb_s1%config, modeldb_s2%config )

  call init_collections()

  call init_time( modeldb_s1 )
! call init_time( modeldb_s2 ) ?????????

  deallocate( S1_nml_file, S2_nml_file )

  ! Create the depository field collection and place it in modeldb
  call modeldb_s1%fields%add_empty_field_collection("depository")
  call modeldb_s1%io_contexts%initialise(program_name, 100)

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise( program_name, modeldb_s1 )

  do while (modeldb_s1%clock%tick())
    call step( program_name, modeldb_s1 )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name, modeldb_s1 )

  call final_time( modeldb_s1 )
  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb_s1 )

end program JediBones
