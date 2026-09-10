!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------
module check_global_mesh_mod

  use constants_mod,   only: i_def, str_def, str_max_filename
  use config_mod,      only: config_type
  use global_mesh_mod, only: global_mesh_type
  use log_mod,         only: log_event,         &
                             log_scratch_space, &
                             log_level_error

  use global_mesh_collection_mod, only: global_mesh_collection

  !------------------------------
  ! Configuration modules
  !------------------------------
  use base_mesh_config_mod, only: key_from_geometry,       &
                                  key_from_topology,       &
                                  geometry_spherical,      &
                                  geometry_planar,         &
                                  TOPOLOGY_FULLY_PERIODIC, &
                                  topology_non_periodic
  implicit none

  private
  public :: check_global_mesh

contains

!> @brief Basic validation that global meshes are suitable
!!        for the specified configuration.
!> @param[in]  config      Configuration object.
!> @param[in]  mesh_names  Global meshes held in application
!!                         global mesh collection object.
subroutine check_global_mesh( config, mesh_names )

  implicit none

  type(config_type),  intent(in) :: config
  character(str_def), intent(in) :: mesh_names(:)

  integer(i_def) :: topology
  integer(i_def) :: geometry

  logical :: valid_geometry
  logical :: valid_topology

  type(global_mesh_type), pointer :: global_mesh

  integer(i_def) :: i

  geometry = config%base_mesh%geometry()
  topology = config%base_mesh%topology()

  do i=1, size(mesh_names)

    global_mesh => global_mesh_collection%get_global_mesh(mesh_names(i))

    ! Check mesh has valid domain geometry
    !=====================================
    valid_geometry = .false.
    select case ( geometry )

    case ( geometry_spherical )
      if ( global_mesh%is_geometry_spherical() ) valid_geometry = .true.

    case ( geometry_planar )
      if ( global_mesh%is_geometry_planar() ) valid_geometry = .true.

    end select

    if ( .not. valid_geometry ) then
      write(log_scratch_space, '(A)')        &
          'Mesh (' // trim(mesh_names(i)) // &
          ') in file is not valid as a ' //  &
          trim(key_from_geometry(geometry)) // ' domain geometry'
      call log_event(log_scratch_space, log_level_error )
    end if


    ! Check mesh has valid domain toplogy
    !=====================================
    valid_topology = .false.
    select case ( topology )

    case ( TOPOLOGY_FULLY_PERIODIC )
      if ( global_mesh%is_topology_periodic() ) valid_topology = .true.

    case ( topology_non_periodic )
      if ( global_mesh%is_topology_non_periodic() ) valid_topology = .true.

    end select

    if ( .not. valid_topology ) then
      write(log_scratch_space, '(A)')           &
          'Mesh (' // trim(mesh_names(i)) //    &
          ') in file does not have a valid ' // &
          trim(key_from_topology(topology)) // ' topology'
      call log_event(log_scratch_space, log_level_error )
    end if

  end do

nullify(global_mesh)

end subroutine check_global_mesh

end module check_global_mesh_mod
