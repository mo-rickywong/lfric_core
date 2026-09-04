!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used
!-----------------------------------------------------------------------------

!> @brief Module to access mesh enumerations.
module sci_mesh_enums_mod

  use constants_mod, only: i_def
  use log_mod,       only: log_event, log_level_error
  use mesh_mod,      only: mesh_type

  use base_mesh_config_mod, only:                                  &
          config_geometry_spherical    => geometry_spherical,      &
          config_geometry_planar       => geometry_planar,         &
          config_topology_periodic     => topology_fully_periodic, &
          config_topology_non_periodic => topology_non_periodic

  implicit none

  private

  public :: geometry_spherical, geometry_planar
  public :: topology_periodic, topology_non_periodic

  public :: get_mesh_geometry, get_mesh_topology

  ! These will get switched to something hardcoded for the science components
  ! at a later date to break dependence on base_mesh_config_mod
  integer(i_def), parameter :: geometry_spherical      = config_geometry_spherical    ! 157
  integer(i_def), parameter :: geometry_planar         = config_geometry_planar       ! 358
  integer(i_def), parameter :: topology_periodic       = config_topology_periodic     ! 492
  integer(i_def), parameter :: topology_non_periodic   = config_topology_non_periodic ! 157

contains

!---------------------------------------------------------------------------
!> @brief  Returns mesh geometry enumeration
!> @param[in] mesh   Mesh object to query
!> @return geometry  Geometry enumeration
!>
function get_mesh_geometry(mesh) result(geometry)

  implicit none

  type(mesh_type), intent(in) :: mesh

  integer(i_def) :: geometry

   if (mesh%is_geometry_spherical()) then
     geometry = geometry_spherical
   else if (mesh%is_geometry_planar()) then
     geometry = geometry_planar
   else
     call log_event('Unsupported mesh geometry', log_level_error)
   end if

 end function get_mesh_geometry



!---------------------------------------------------------------------------
!> @brief  Returns mesh topology enumeration
!> @param[in] mesh   Mesh object to query
!> @return topology  Topology enumeration
!>
function get_mesh_topology(mesh) result(topology)

  implicit none

  type(mesh_type), intent(in) :: mesh

  integer(i_def) :: topology

   if (mesh%is_topology_periodic()) then
     topology = topology_periodic
   else if (mesh%is_topology_non_periodic()) then
     topology = topology_non_periodic
   else
     call log_event('Unsupported mesh topology', log_level_error)
   end if

 end function get_mesh_topology

end module sci_mesh_enums_mod
