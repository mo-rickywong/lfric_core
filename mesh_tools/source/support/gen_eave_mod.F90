!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Module to derive eave mesh(es) from a given gencube_ps_type
!>          generation stragtegy.
!> @details Type implements the ugrid_generator_type interface to
!>          construct eave meshes derived from a given cubed-sphere mesh.
!>
!>      +---+                       +---------+
!>      | 5 |                       |         |
!>  +---+---+---+---+               |  +---+  |
!>  | 1 | 2 | 3 | 4 | --->      6 X |  | ? |  |
!>  +---+---+---+---+               |  +---+  |
!>      | 6 |                       |         |
!>      +---+                       +---------+
!>
!>  6 x eave meshes are generated where the centre cells of each are linked
!>  (via intergrid map) to the common cells on the corresponding panel of the
!>  cubed-sphere parent.
!>
!>  Each eave mesh can be described as a grid topology with non-periodic domain
!>  boundaries. The mappings and node/cell centre coords are linked to the
!>  corresponding cells on the cubed-sphere. The eave mesh extends past the
!>  parent panel by an addition number of cells (eave depth). The cell
!>  size in the eaves is kept consistent to that of the parent cubed-sphere
!>  mesh.
!-------------------------------------------------------------------------------
module gen_eave_mod
!-------------------------------------------------------------------------------

  use constants_mod,                  only: r_def, i_def, str_def, l_def,     &
                                            str_long, str_longlong,           &
                                            PI, radians_to_degrees,           &
                                            degrees_to_radians, rmdi, imdi,   &
                                            emdi
  use coord_transform_mod,            only: ll2xyz, xyz2ll
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use log_mod,                        only: log_event, log_scratch_space, &
                                            LOG_LEVEL_ERROR, LOG_LEVEL_INFO
  use reference_element_mod,          only: W, S, E, N, SWB, SEB, NWB, NEB
  use rotation_mod,                   only: rotate_mesh_coords, &
                                            TRUE_NORTH_POLE_LL, &
                                            TRUE_NULL_ISLAND_LL

  use ugrid_generator_mod,            only: ugrid_generator_type

  use mesh_config_mod, only: coord_sys_ll,       &
                             coord_sys_xyz,      &
                             key_from_geometry,  &
                             key_from_topology,  &
                             key_from_coord_sys, &
                             geometry_spherical, &
                             topology_periodic


  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Mesh Vertex directions: local aliases for reference_element_mod values
  integer(i_def), parameter :: NW = NWB
  integer(i_def), parameter :: NE = NEB
  integer(i_def), parameter :: SE = SEB
  integer(i_def), parameter :: SW = SWB

  ! For a cubesphere these panels are the 6 faces of the domain cube
  integer(i_def), parameter :: NPANELS = 6
  integer(i_def), parameter :: PANEL_ROTATIONS(NPANELS) = (/ 0, 0, 1, 1, -1, 0 /)

  ! Prefix for error messages
  character(*),       parameter :: PREFIX = "[Cubed-Sphere Mesh] "

  ! flag to print out mesh data for debugging purposes
  logical(l_def),     parameter :: DEBUG = .false.

  !-------------------------------------------------------------------------------

  type, extends(ugrid_generator_type), public :: gen_eave_type

    private

    character(str_def) :: mesh_name
    integer(i_def)     :: geometry  = geometry_spherical
    integer(i_def)     :: topology  = topology_periodic
    integer(i_def)     :: coord_sys = emdi
    character(str_def) :: coord_units(2)

    integer(i_def)     :: edge_cells
    real(r_def)        :: domain_extents(2,4,NPANELS)
    integer(i_def)     :: npanels   = NPANELS
    real(r_def)        :: north_pole(2)
    real(r_def)        :: null_island(2)
    real(r_def)        :: equatorial_latitude = 0.0_r_def

    character(str_longlong) :: constructor_inputs

    integer(i_def)     :: nsmooth
    real(r_def)        :: stretch_factor = 1.0_r_def
    logical(l_def)     :: rotate_mesh    = .false.
    logical(l_def)     :: periodicity(2) = .false.

    ! Connectivity and coordinates
    integer(i_def), allocatable :: cell_next(:,:,:)
    integer(i_def), allocatable :: nodes_on_cell(:,:,:)

    integer(i_def), allocatable :: edges_on_cell(:,:,:)
    integer(i_def), allocatable :: nodes_on_edge(:,:,:)
    real(r_def),    allocatable :: node_coords(:,:,:)
    real(r_def),    allocatable :: cell_coords(:,:,:)


    ! Intergrid maps
    ! Only 1 map to the parent cubed-sphere
    integer(i_def)     :: nmaps = 1
    character(str_def) :: parent_mesh_name
    integer(i_def)     :: parent_edge_cells

    integer(i_def)     :: eave_depth = imdi
    integer(i_def)     :: eave_edge_cells
    integer(i_def)     :: cpp
    integer(i_def)     :: npp
    integer(i_def)     :: epp

    integer(i_def), allocatable :: cell_maps(:,:)
    integer(i_def)     :: void_id

    type(global_mesh_map_collection_type), allocatable :: global_mesh_maps

    integer(i_def) :: max_num_faces_per_node

    logical :: generated = .false.

  contains
    procedure :: generate
    procedure :: get_number_of_panels
    procedure :: get_metadata
    procedure :: get_dimensions
    procedure :: get_coordinates
    procedure :: get_connectivity
    procedure :: get_global_mesh_maps
    procedure :: calc_cell_maps
    procedure :: compare_eave
    procedure :: is_generated
    procedure :: populate_ugrid_2d

    procedure :: clear
    final     :: gen_eave_final

  end type gen_eave_type

!-------------------------------------------------------------------------------
  interface gen_eave_type
    module procedure gen_eave_constructor
  end interface gen_eave_type
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> @brief   Constructor for gen_cs_eave_type.
!> @details Accepts mesh dimension for initialisation and validation.
!>
!> @param[in] mesh_name          Name of this mesh topology
!> @param[in] edge_cells         Number of cells per panel edge of the cubed-sphere.
!>                               Each panel will contain edge_cells*edge_cells faces.
!> @param[in] nsmooth            Number of smoothing passes to be performed on mesh nodes.
!>                               Each panel will contain edge_cells*edge_cells faces.
!> @param[in] coord_sys          Coordinate system to position nodes.
!> @param[in, optional] rotate_mesh
!>                               Logical to indicate rotation of the resulting mesh
!> @param[in, optional] target_north_pole
!>                               If rotating the cubed-sphere, then these are the target
!>                               coordinates [longitude, latitude] (degrees)
!>                               to move the reference north pole.
!> @param[in, optional] target_null_island
!>                               If rotating the cubed-sphere, then these are the target
!>                               coordinates [longitude, latitude] (degrees)
!>                               to move the reference null island.
!> @param[in, optional] target_mesh_names
!>                               Names of meshes to map to.
!> @param[in, optional] target_edge_cells
!>                               Number of cells per panel edge of the meshes to map to.
!> @param[in, optional] stretch_factor
!>                               Attracts points to the North (< 1) or South (> 1) to give
!>                               a variable resolution mesh
!>
!> @return    self               Instance of gen_cs_eave_type
!-------------------------------------------------------------------------------
 function gen_eave_constructor( gen_cs, eave_depth ) &
                                    result( self )

  use gencube_ps_mod,  only: gencube_ps_type
  use mesh_config_mod, only: coord_sys_from_key, &
                             geometry_from_key,  &
                             topology_non_periodic

  implicit none

  class(gencube_ps_type), intent(in) :: gen_cs
  integer(i_def),         intent(in) :: eave_depth

  type( gen_eave_type ) :: self

  character(str_def) :: geometry_str
  character(str_def)     :: coord_sys_str

  call gen_cs%get_metadata( mesh_name    = self%parent_mesh_name, &
                            edge_cells_x = self%edge_cells,       &
                            north_pole   = self%north_pole,       &
                            null_island  = self%null_island,      &
                            coord_sys    = coord_sys_str,         &
                            geometry     = geometry_str,          &
                            void_cell    = self%void_id )

  self%equatorial_latitude = gen_cs%get_equatorial_latitude()

  self%rotate_mesh = gen_cs%get_rotate_mesh()
  self%north_pole  = self%north_pole*degrees_to_radians
  self%null_island = self%null_island*degrees_to_radians

  self%mesh_name  = trim(self%parent_mesh_name)//'_eave'

  self%nsmooth    = 0
  self%nmaps      = 1
  self%coord_sys  = coord_sys_from_key(coord_sys_str)
  self%topology   = topology_non_periodic
  self%geometry   = geometry_from_key(geometry_str)

  self%eave_depth = eave_depth
  self%eave_edge_cells = self%edge_cells + 2* self%eave_depth
  self%cpp = self%eave_edge_cells**2
  self%npp = (self%eave_edge_cells+1)**2
  self%epp = 2*(self%eave_edge_cells**2 + self%eave_edge_cells)

  ! There are a maximum of 4 faces around a node in this type of mesh
  self%max_num_faces_per_node = 4
  self%nmaps = 1

  write(self%constructor_inputs,'(A,I0)')  &
      'cs_strategy=<gencube_ps_type,"'//   &
      trim(self%parent_mesh_name)//'">;'// &
      'eave_depth=', self%eave_depth

  return
end function gen_eave_constructor

!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the set of cells to which it is adjacent.
!>          (Private Routine)
!> @details Allocates and populates the instance's cell_next(:,:) array
!>          with the id of each cell to which the index cell is adjacent.
!>
!> @param[in]   gen_cube   Generator strategy for a cubed-sphere
!> @param[out]  cell_next  A rank 2 (4,ncells)-sized array containing the
!>                         adjacency map.
!-------------------------------------------------------------------------------
subroutine calc_adjacency(self, cell_next)

  implicit none

  class(gen_eave_type),    intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: cell_next(:,:,:)

  integer(i_def) :: edge_cells, ncells, cpp, i
  integer(i_def) :: cell, astat, panel_number

  integer(i_def), allocatable :: panel_edge_cells_west(:)
  integer(i_def), allocatable :: panel_edge_cells_south(:)
  integer(i_def), allocatable :: panel_edge_cells_east(:)
  integer(i_def), allocatable :: panel_edge_cells_north(:)

  integer(i_def), allocatable :: panel_next(:,:)
  integer(i_def), allocatable :: panel_edge_cells(:,:,:)

  integer(i_def) :: eave_depth
  integer(i_def) :: eave_edge_cells

  edge_cells = self%edge_cells
  eave_depth = self%eave_depth

  eave_edge_cells = edge_cells + 2*eave_depth
  cpp             = eave_edge_cells**2
  ncells          = cpp*NPANELS



  allocate(cell_next(4, cpp, NPANELS), stat=astat)
  if (astat /= 0)                                               &
      call log_event( PREFIX//"Failure to allocate cell_next.", &
                      LOG_LEVEL_ERROR )

  allocate(panel_edge_cells_west(edge_cells))
  allocate(panel_edge_cells_south(edge_cells))
  allocate(panel_edge_cells_east(edge_cells))
  allocate(panel_edge_cells_north(edge_cells))

  cell_next = 0

  ! Panels are arranged and numbered as indicated
  !==============================================
  !      +---+
  !      | 5 |
  !  +---+---+---+---+
  !  | 1 | 2 | 3 | 4 |
  !  +---+---+---+---+
  !      | 6 |
  !      +---+
  !
  ! Cells in each panel are numbered from NW panel corner
  ! i.e. 3x3 panel number #1 is:
  !
  ! +-----+
  ! |1|2|3|
  ! |-+-+-|
  ! |4|5|6|
  ! |-+-+-|
  ! |7|8|9|
  ! +-----+

  allocate(panel_next(4, NPANELS))
  allocate(panel_edge_cells (eave_edge_cells,4,NPANELS))

  ! Ordering : W,S,E,N
  !  panel_next(:,1) = [4,6,2,5]
  !  panel_next(:,2) = [1,6,3,5]
  !  panel_next(:,3) = [2,6,4,5]
  !  panel_next(:,4) = [3,6,1,5]
  !  panel_next(:,5) = [1,2,3,4]
  !  panel_next(:,6) = [1,4,3,2]

  call get_panel_edge_cell_ids( eave_edge_cells, panel_edge_cells )

  cell_next = 0_i_def

  ! Default settings
  panel_number = 1
  do cell=1, cpp
    ! Default: W, S, E, N
    cell_next(:, cell, panel_number) = (/ cell - 1,                 & ! W
                                          cell + eave_edge_cells,   & ! S
                                          cell + 1,                 & ! E
                                          cell - eave_edge_cells /)   ! N
  end do

  ! Panel I
  !===========================================================================

  ! Since each panel will start with local ids of 1 use
  ! panel 1 to mark the connectivy on the edge of the eave
  ! mesh domains.
  do i=1, eave_edge_cells
    ! Top edge
    cell = panel_edge_cells(i, N, panel_number )
    cell_next(N, cell, panel_number) = self%VOID_ID

    ! Right edge
    cell = panel_edge_cells(i, E, panel_number )
    cell_next(E, cell, panel_number) = self%VOID_ID

    ! Bottom edge
    cell = panel_edge_cells(i, S, panel_number )
    cell_next(S, cell, panel_number) = self%VOID_ID

    ! Left edge
    cell = panel_edge_cells(i, W, panel_number )
    cell_next(W, cell, panel_number) = self%VOID_ID
  end do


  return
end subroutine calc_adjacency

!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the four vertices which comprise it.
!>          (Private Routine)
!> @details Allocates and populates the instance's verts_on_cell(:,:) array with
!>          the vertices which form each cell.
!>
!> @param[in]   gen_cube       Generator strategy for a cubed-sphere
!> @param[out]  verts_on_cell  A rank 2 (4,ncells)-sized integer array of vertices
!>                             which constitute each cell.
!-------------------------------------------------------------------------------
subroutine calc_face_to_node(self, nodes_on_cell)

  implicit none

  class(gen_eave_type),    intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: nodes_on_cell(:,:,:)

  integer(i_def) :: edge_cells, ncells, cpp
  integer(i_def) :: cell,  nxf, astat
  integer(i_def) :: cell_id
  integer(i_def) :: node_id
  integer(i_def) :: col
  integer(i_def) :: row
  integer(i_def) :: loop
  integer(i_def) :: panel_id
  integer(i_def) :: eave_depth,  eave_edge_cells



  edge_cells = self%edge_cells
  eave_depth = self%eave_depth

  eave_edge_cells = edge_cells + 2*eave_depth
  cpp        = eave_edge_cells**2
  ncells     = cpp*NPANELS

  allocate(nodes_on_cell(4, cpp, NPANELS), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate nodes_on_cell.", &
                      LOG_LEVEL_ERROR )

  nodes_on_cell = 0
  cell = 1
  nxf = 1

  ! NW node of every cell in panels 1
  panel_id = 1
  node_id = 0
  cell_id = 0
  do row=1, eave_edge_cells

    do col=1, eave_edge_cells
      cell_id = cell_id + 1
      node_id = node_id + 1
      nodes_on_cell(NW, cell_id, panel_id) = node_id
      if (col == eave_edge_cells) then
        node_id = node_id + 1
        nodes_on_cell(NE, cell_id, panel_id) = node_id
      end if
    end do

    if (row == eave_edge_cells) then
      cell_id = cell_id - eave_edge_cells
      do col=1, eave_edge_cells
        cell_id = cell_id + 1
        node_id = node_id + 1
        nodes_on_cell(SW, cell_id, panel_id) = node_id
        if (col == eave_edge_cells) then
          node_id = node_id + 1
          nodes_on_cell(SE, cell_id, panel_id) = node_id
        end if
      end do
    end if

  end do


  ! Fill in the rest now the nodes have been numbered
  do loop=1,2
    panel_id = 1
    node_id = 0
    cell_id = 0
    do row=1, eave_edge_cells
      do col=1, eave_edge_cells
        cell_id = cell_id + 1
        node_id = node_id + 1

        if (self%cell_next(E, cell_id, panel_id) /= self%void_id) then

          if (nodes_on_cell(NE, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( NW, self%cell_next(E, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(NE, cell_id, panel_id) =                        &
                       nodes_on_cell( NW,                                   &
                                      self%cell_next(E, cell_id, panel_id), &
                                      panel_id )
            end if
          end if

          if (nodes_on_cell(SE, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( SW, self%cell_next(E, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(SE, cell_id, panel_id) =                        &
                       nodes_on_cell( SW,                                   &
                                      self%cell_next(E, cell_id, panel_id), &
                                      panel_id )
            end if
          end if

        end if

        if (self%cell_next(S, cell_id, panel_id) /= self%void_id) then
          if (nodes_on_cell(SE, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( NE, self%cell_next(S, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(SE, cell_id, panel_id) =                        &
                       nodes_on_cell( NE,                                   &
                                      self%cell_next(S, cell_id, panel_id), &
                                      panel_id )
            end if
          end if

          if (nodes_on_cell(SW, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( NW, self%cell_next(S, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(SW, cell_id, panel_id) =                        &
                       nodes_on_cell( NW,                                   &
                                      self%cell_next(S, cell_id, panel_id), &
                                      panel_id )
            end if
          end if
        end if


        if (self%cell_next(N, cell_id, panel_id) /= self%void_id) then

          if (nodes_on_cell(NW, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( SW, self%cell_next(N, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(NW, cell_id, panel_id) =                        &
                       nodes_on_cell( SW,                                   &
                                      self%cell_next(N, cell_id, panel_id), &
                                      panel_id )
            end if
          end if

          if (nodes_on_cell(NE, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( SE, self%cell_next(N, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(NE, cell_id, panel_id) =                        &
                       nodes_on_cell( SE,                                   &
                                      self%cell_next(N, cell_id, panel_id), &
                                      panel_id )
            end if
          end if
        end if



        if (self%cell_next(W, cell_id, panel_id) /= self%void_id) then

          if (nodes_on_cell(NW, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( NE, self%cell_next(W, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then

              nodes_on_cell(NW, cell_id, panel_id) =                        &
                       nodes_on_cell( NE,                                   &
                                      self%cell_next(W, cell_id, panel_id), &
                                      panel_id )

            end if
          end if

          if (nodes_on_cell(SW, cell_id, panel_id) == 0_i_def) then
            if (nodes_on_cell( SE, self%cell_next(W, cell_id, panel_id), &
                               panel_id ) /= 0_i_def ) then
              nodes_on_cell(SW, cell_id, panel_id) =                        &
                       nodes_on_cell( SE,                                   &
                                      self%cell_next(W, cell_id, panel_id), &
                                      panel_id )
            end if
          end if
        end if

      end do ! columns
    end do ! rows
  end do ! loops

  return
end subroutine calc_face_to_node

!-------------------------------------------------------------------------------
!> @brief   Calculates the edges which are found on each cell and the
!>          pair of vertices which are found on each edge.
!>          (Private Routine)
!> @details Allocates and populates both the edges_on_cell and
!>          verts_on_edge arrays for the instance.
!>
!> @param[in]   gen_cube       Generator strategy for a cubed-sphere
!> @param[out]  edges_on_cell  A rank-2 (4,ncells)-sized integer array of
!>                             the edges found on each cell.
!> @param[out]  verts_on_edge  A rank-2 (2,2*ncells)-sized integer array
!>                             of the vertices found on each edge.
!-------------------------------------------------------------------------------
subroutine calc_edges( self, edges_on_cell, nodes_on_edge)

  implicit none

  class(gen_eave_type),    intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: edges_on_cell(:,:,:)
  integer(i_def), allocatable, intent(out) :: nodes_on_edge(:,:,:)

  integer(i_def),parameter :: edges_per_cell = 4
  integer(i_def) :: edge_cells, cpp
  integer(i_def) :: cell,  nxf, astat !, edge

  integer(i_def) :: n_edges
  integer(i_def) :: row
  integer(i_def) :: col
  integer(i_def) :: panel_id
  integer(i_def) :: eave_edge_cells

  edge_cells      = self%edge_cells
  eave_edge_cells = self%eave_edge_cells

  cpp        = eave_edge_cells**2
  n_edges    = 2*(eave_edge_cells**2 + eave_edge_cells)

  allocate(edges_on_cell(edges_per_cell, cpp, NPANELS), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate edges_on_cell.", &
                      LOG_LEVEL_ERROR )

  allocate(nodes_on_edge(2, n_edges, NPANELS), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate nodes_on_edge.", &
                      LOG_LEVEL_ERROR )

  edges_on_cell = self%VOID_ID
  nodes_on_edge = self%VOID_ID
  cell = 1
  nxf = 1
  panel_id = 1

  ! Top row of panel
  do cell=1, eave_edge_cells

    edges_on_cell(N, cell, panel_id) = nxf
    edges_on_cell(W, cell, panel_id) = nxf+1
    edges_on_cell(S, cell, panel_id) = nxf+2

    nodes_on_edge(1, nxf, panel_id)   = self%nodes_on_cell(NW, cell, panel_id)
    nodes_on_edge(2, nxf, panel_id)   = self%nodes_on_cell(NE, cell, panel_id)
    nodes_on_edge(1, nxf+1, panel_id) = self%nodes_on_cell(SW, cell, panel_id)
    nodes_on_edge(2, nxf+1, panel_id) = self%nodes_on_cell(NW, cell, panel_id)
    nodes_on_edge(1, nxf+2, panel_id) = self%nodes_on_cell(SE, cell, panel_id)
    nodes_on_edge(2, nxf+2, panel_id) = self%nodes_on_cell(SW, cell, panel_id)

    if ( cell == eave_edge_cells ) then
      edges_on_cell(E, cell, panel_id)  = nxf+3
      nodes_on_edge(1, nxf+3, panel_id) = self%nodes_on_cell(SE, cell, panel_id)
      nodes_on_edge(2, nxf+3, panel_id) = self%nodes_on_cell(NE, cell, panel_id)
      nxf = nxf + 4
    else
      nxf = nxf + 3
    end if

    ! Push W edge to W neighbour
    if ( self%cell_next(W, cell, panel_id) /= self%VOID_ID ) then
      edges_on_cell(E, self%cell_next(W, cell, panel_id), panel_id) = &
           edges_on_cell(W, cell, panel_id)
    end if
  end do

  cell = eave_edge_cells
  ! Remainder of panel
  do row=2, eave_edge_cells
    do col=1, eave_edge_cells

      cell = cell+1
      edges_on_cell(W, cell, panel_id)  = nxf
      edges_on_cell(S, cell, panel_id)  = nxf+1
      nodes_on_edge(1, nxf, panel_id)   = self%nodes_on_cell(SW, cell, panel_id)
      nodes_on_edge(2, nxf, panel_id)   = self%nodes_on_cell(NW, cell, panel_id)
      nodes_on_edge(1, nxf+1, panel_id) = self%nodes_on_cell(SE, cell, panel_id)
      nodes_on_edge(2, nxf+1, panel_id) = self%nodes_on_cell(SW, cell, panel_id)

      if ( col == eave_edge_cells ) then
        edges_on_cell(E, cell, panel_id) = nxf+2
        nodes_on_edge(1, nxf+2, panel_id) = self%nodes_on_cell(SE, cell, panel_id)
        nodes_on_edge(2, nxf+2, panel_id) = self%nodes_on_cell(NE, cell, panel_id)
        nxf = nxf + 3
      else
        nxf = nxf + 2
      end if

      ! Copy N edge from N cell
      edges_on_cell(N, cell, panel_id) = edges_on_cell(S, self%cell_next(N, cell, panel_id), panel_id)

      ! Push W edge to W neighbour
      if ( self%cell_next(W, cell, panel_id) /= self%VOID_ID ) then
        edges_on_cell(E, self%cell_next(W, cell, panel_id), panel_id) = edges_on_cell(W, cell, panel_id)
      end if

    end do ! columns
  end do ! rows

  return
end subroutine calc_edges

!-------------------------------------------------------------------------------
!> @brief   Calculates the coordinates of vertices in the mesh.
!> @details Assigns an (x,y) lat-long coordinate to each mesh
!>          vertex according to its Cartesian position in the mesh.
!>
!> @param[in]   gen_cube       Generator strategy for a cubed-sphere
!> @param[out]  vert_coords    A rank 2 (2,ncells)-sized real array of long and
!>                             lat coordinates (degrees) respectively for
!>                             each vertex.
!> @param[out]  coord_units_x  Units of x-coordinate.
!> @param[out]  coord_units_y  Units of y-coordinate.
!-------------------------------------------------------------------------------
subroutine calc_coords( self,        &
                        node_coords, &
                        coord_units )


  implicit none

  class(gen_eave_type), intent(inout)  :: self

  real(r_def), allocatable, intent(out) :: node_coords(:,:,:)
  character(str_def), intent(out) :: coord_units(2)

  integer(i_def) ::  x, y, astat, cpp, node !, node0

  integer(i_def) :: n_cells
  integer(i_def) :: n_nodes
  integer(i_def) :: panel_id
  integer(i_def) :: eave_edge_cells
  integer(i_def) :: edge_cells
  real(r_def)    :: eave_panel_size

  real(r_def)    :: min_x,max_x
  real(r_def)    :: min_y,max_y

  real(r_def)    :: lat, long
  real(r_def)    :: x0, y0, z0
  real(r_def)    :: xs, ys, zs

  real(r_def)    :: dlambda
  real(r_def)    :: lambda1, lambda2
  real(r_def)    :: t1, t2

  real(r_def), parameter :: pio4 = PI*0.25_r_def ! 45 degrees
  real(r_def), parameter :: pio2 = PI*0.50_r_def ! 90 degrees

  integer(i_def) :: sw_domain_cell
  integer(i_def) :: se_domain_cell
  integer(i_def) :: ne_domain_cell
  integer(i_def) :: nw_domain_cell
  integer(i_def) :: sw_domain_node
  integer(i_def) :: se_domain_node
  integer(i_def) :: ne_domain_node
  integer(i_def) :: nw_domain_node

  edge_cells      = self%edge_cells
  eave_edge_cells = self%eave_edge_cells
  cpp             = eave_edge_cells**2
  n_cells         = cpp
  n_nodes         = (eave_edge_cells+1)**2


  if (allocated(node_coords)) deallocate(node_coords)

  allocate(node_coords(2, n_nodes, NPANELS), stat=astat)

  if (astat /= 0)                                                 &
      call log_event( PREFIX//"Failure to allocate node_coords.", &
                      LOG_LEVEL_ERROR )

  node_coords = 0.0_r_def
  dlambda     = pio2/edge_cells  ! dlamba in radians

  eave_panel_size = eave_edge_cells*dlambda

  node = 1

  ! Panels I-VI
  do y=1, eave_edge_cells + 1
    lambda2 = (y-1)*dlambda - eave_panel_size*0.5_r_def
    t2 = tan(lambda2)
    do x=1, eave_edge_cells + 1
      lambda1 = (x-1)*dlambda - eave_panel_size*0.5_r_def
      t1 = tan(lambda1)

      ! Panel I
      xs = 1.0_r_def/sqrt(1.0_r_def + t1*t1 + t2*t2)
      ys = xs*t1
      zs = xs*t2
      call xyz2ll(xs, ys, zs, long, lat)
      node_coords(1, node, 1) = long
      node_coords(2, node, 1) = -lat


      ! Panel II
      x0 = -ys
      y0 =  xs
      z0 =  zs
      call xyz2ll(x0, y0, z0, long, lat)
      node_coords(1, node, 2) = long
      node_coords(2, node, 2) = -lat


      ! Panel III
      x0 = -xs
      y0 = -ys
      z0 =  zs
      call xyz2ll(x0, y0, z0, long, lat)
      node_coords(1, node, 3) = long
      node_coords(2, node, 3) = -lat


      ! Panel IV
      x0 =  ys
      y0 = -xs
      z0 =  z0
      call xyz2ll(x0, y0, z0, long, lat)
      node_coords(1, node, 4) = long
      node_coords(2, node, 4) = -lat


      ! Panel V
      x0 = -ys
      y0 =  zs
      z0 = -xs
      call xyz2ll(x0, y0, z0, long, lat)
      node_coords(1, node, 5) = long
      node_coords(2, node, 5) = -lat


      ! Panel VI
      x0 = -ys
      y0 = -zs
      z0 =  xs
      call xyz2ll(x0, y0, z0, long, lat)
      node_coords(1, node, 6) = long
      node_coords(2, node, 6) = -lat

      node = node + 1

    end do
  end do

  ! Output units from xyz2ll are in radians
  coord_units(:) = 'radians'

  sw_domain_cell = self%eave_edge_cells*(self%eave_edge_cells - 1) + 1
  se_domain_cell = self%eave_edge_cells**2
  ne_domain_cell = self%eave_edge_cells
  nw_domain_cell = 1

  do panel_id=1, NPANELS
    sw_domain_node = self%nodes_on_cell( SW, sw_domain_cell, panel_id )
    se_domain_node = self%nodes_on_cell( SE, se_domain_cell, panel_id )
    ne_domain_node = self%nodes_on_cell( NE, ne_domain_cell, panel_id )
    nw_domain_node = self%nodes_on_cell( NW, nw_domain_cell, panel_id )

    min_x = minval(self%node_coords(1, :, panel_id))
    max_x = maxval(self%node_coords(1, :, panel_id))
    min_y = minval(self%node_coords(2, :, panel_id))
    max_y = maxval(self%node_coords(2, :, panel_id))

    self%domain_extents(:, SW, panel_id) = [min_x,min_y] ! self%node_coords(:, sw_domain_node, panel_id )
    self%domain_extents(:, SE, panel_id) = [max_x,min_y] ! self%node_coords(:, se_domain_node, panel_id )
    self%domain_extents(:, NE, panel_id) = [max_x,max_y] ! self%node_coords(:, ne_domain_node, panel_id )
    self%domain_extents(:, NW, panel_id) = [max_x,min_y] ! self%node_coords(:, nw_domain_node, panel_id )
  end do

  return
end subroutine calc_coords


!-------------------------------------------------------------------------------
!> @brief Populates the arguments with the dimensions defining
!>        the mesh.
!>
!> @param[in]   self                   The gen_eave_type instance reference.
!> @param[out]  num_nodes              The number of nodes on the mesh.
!> @param[out]  num_edges              The number of edges on the mesh.
!> @param[out]  num_faces              The number of faces on the mesh.
!> @param[out]  num_nodes_per_face     The number of nodes around each face.
!> @param[out]  num_edges_per_face     The number of edges around each face.
!> @param[out]  num_nodes_per_edge     The number of nodes around each edge.
!> @param[out]  max_num_faces_per_node The maximum number of faces surrounding
!>                                     each node.
!-------------------------------------------------------------------------------
subroutine get_dimensions(self, num_nodes, num_edges, num_faces,           &
                                num_nodes_per_face, num_edges_per_face,    &
                                num_nodes_per_edge, max_num_faces_per_node )

  implicit none

  class(gen_eave_type), intent(in) :: self

  integer(i_def), intent(out) :: num_nodes
  integer(i_def), intent(out) :: num_edges
  integer(i_def), intent(out) :: num_faces
  integer(i_def), intent(out) :: num_nodes_per_face
  integer(i_def), intent(out) :: num_edges_per_face
  integer(i_def), intent(out) :: num_nodes_per_edge
  integer(i_def), intent(out) :: max_num_faces_per_node

  integer(i_def) :: edge_cells, cpp, ncells

  edge_cells = self%edge_cells
  cpp        = edge_cells*edge_cells
  ncells     = cpp*NPANELS

  num_faces = ncells
  num_nodes = ncells + 2
  num_edges = ncells * 2

  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2

  max_num_faces_per_node = self%max_num_faces_per_node

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!> @brief   Populates the argument array with the coordinates of
!>          the mesh's vertices.
!> @details Exposes the instance's vert_coords array to the caller.
!>
!> @param[out]  node_coordinates  The argument to receive the vert_coords data.
!> @param[out]  cell_coordinates  Cell centre coordinates
!> @param[out]  domain_extents    Principal coordiantes describing domain.
!> @param[out]  coord_units_x     Units of x-coordinate.
!> @param[out]  coord_units_y     Units of y-coordinate.
!-------------------------------------------------------------------------------
subroutine get_coordinates(self, node_coordinates, &
                                 cell_coordinates, &
                                 domain_extents,   &
                                 coord_units_x,    &
                                 coord_units_y)

  implicit none

  class(gen_eave_type), intent(in)  :: self

  real(r_def), allocatable, intent(out) :: node_coordinates(:,:)
  real(r_def), allocatable, intent(out) :: cell_coordinates(:,:)
  real(r_def), allocatable, intent(out) :: domain_extents(:,:)

  character(str_def),     intent(out) :: coord_units_x
  character(str_def),     intent(out) :: coord_units_y

  if (allocated(node_coordinates)) deallocate(node_coordinates)
  if (allocated(cell_coordinates)) deallocate(cell_coordinates)
  if (allocated(domain_extents))   deallocate(domain_extents)

  allocate(node_coordinates, source=self%node_coords(:,:,1))
  allocate(cell_coordinates, source=self%cell_coords(:,:,1))
  allocate(domain_extents,   source=self%domain_extents(:,:,1))

  coord_units_x    = self%coord_units(1)
  coord_units_y    = self%coord_units(2)

  return
end subroutine get_coordinates

!-------------------------------------------------------------------------------
!> @brief   Populates the argument arrays with the corresponding mesh
!>          connectivity information.
!> @details Implements the connectivity-providing interface required
!>          by the ugrid writer.
!>
!>  @param[out]  face_node_connectivity  Face-node connectivity.
!>  @param[out]  face_edge_connectivity  Face-edge connectivity.
!>  @param[out]  face_face_connectivity  Face-face connectivity.
!>  @param[out]  edge_node_connectivity  Edge-node connectivity.
!-------------------------------------------------------------------------------
subroutine get_connectivity(self, face_node_connectivity, &
                                  face_edge_connectivity, &
                                  face_face_connectivity, &
                                  edge_node_connectivity )
  implicit none

  class(gen_eave_type), intent(in) :: self

  integer(i_def), allocatable, intent(out) :: face_node_connectivity(:,:)
  integer(i_def), allocatable, intent(out) :: face_edge_connectivity(:,:)
  integer(i_def), allocatable, intent(out) :: face_face_connectivity(:,:)
  integer(i_def), allocatable, intent(out) :: edge_node_connectivity(:,:)

  if (allocated(face_node_connectivity)) deallocate(face_node_connectivity)
  if (allocated(face_edge_connectivity)) deallocate(face_edge_connectivity)
  if (allocated(face_face_connectivity)) deallocate(face_face_connectivity)
  if (allocated(edge_node_connectivity)) deallocate(edge_node_connectivity)

  allocate(face_node_connectivity, source=self%nodes_on_cell(:,:,1) )
  allocate(face_edge_connectivity, source=self%edges_on_cell(:,:,1) )
  allocate(face_face_connectivity, source=self%cell_next(:,:,1)     )
  allocate(edge_node_connectivity, source=self%nodes_on_edge(:,:,1) )

  return
end subroutine get_connectivity

!-------------------------------------------------------------------------------
!> @brief  Gets the global mesh map collection which uses
!>         this mesh as the source mesh
!>
!> @return global_mesh_maps global_mesh_map_collection_type
!-------------------------------------------------------------------------------
function get_global_mesh_maps(self) result (global_mesh_maps)

  implicit none

  class(gen_eave_type), target, intent(in) :: self

  type(global_mesh_map_collection_type), pointer  :: global_mesh_maps

  global_mesh_maps => self%global_mesh_maps

  return
end function get_global_mesh_maps

!-------------------------------------------------------------------------------
!> @brief   Generates the mesh and connectivity.
!> @details Calls each of the instance methods which calculate the
!>          specified mesh and populate the arrays.
!>
!> @param[in,out]  self  The gen_eave_type instance reference.
!-------------------------------------------------------------------------------
subroutine generate(self)


  use gencube_ps_mod, only: stretch_mesh, calc_cell_centres
  use rotation_mod,   only: rotate_mesh_coords

  implicit none

  class(gen_eave_type), intent(inout) :: self

  integer(i_def) :: panel_id
!  integer(i_def) :: shift_value

  integer(i_def), allocatable :: nodes_on_cell(:,:)
  real(r_def),    allocatable :: node_coords(:,:)
  real(r_def),    allocatable :: cell_coords(:,:)

  call calc_adjacency(self, self%cell_next)
  call calc_face_to_node(self, self%nodes_on_cell)
  call calc_edges(self, self%edges_on_cell, self%nodes_on_edge)

  do panel_id=2,6
    self%cell_next(:,:,panel_id)     = self%cell_next(:,:,1)
    self%nodes_on_cell(:,:,panel_id) = self%nodes_on_cell(:,:,1)
    self%edges_on_cell(:,:,panel_id) = self%edges_on_cell(:,:,1)
    self%nodes_on_edge(:,:,panel_id) = self%nodes_on_edge(:,:,1)
  end do

  call self%calc_cell_maps()

  ! Co-ord output from calc_coords in radians
  call calc_coords(self, self%node_coords, &
                         self%coord_units )

  call orient_lfric(self, panel_rotations)

  allocate( nodes_on_cell(4,self%cpp) )
  allocate( node_coords(2,self%npp) )

  ! Apply stretch transform
  if (self%equatorial_latitude /= 0.0_r_def) then
    do panel_id=1, npanels
      call stretch_mesh( self%equatorial_latitude, self%node_coords(:,:, panel_id) )
    end do
  end if

  ! Apply rotation transform
  if (self%rotate_mesh) then
    do panel_id=1, npanels
      node_coords(:,:) = self%node_coords(:,:,panel_id)
      call rotate_mesh_coords( self%north_pole, node_coords )
      self%node_coords(:,:,panel_id) = node_coords(:,:)
    end do
  end if

  ! Calculate cell centres
  if ( allocated(self%cell_coords) ) deallocate(self%cell_coords)
  allocate( self%cell_coords(2,self%cpp,npanels) )

  do panel_id=1, npanels
    nodes_on_cell(:,:) = self%nodes_on_cell(:,:,panel_id)
    node_coords(:,:)   = self%node_coords(:,:,panel_id)
    call calc_cell_centres( nodes_on_cell, node_coords, cell_coords )
    self%cell_coords(:,:,panel_id) = cell_coords(:,:)
  end do

  self%generated = .true.

  return
end subroutine generate

!-------------------------------------------------------------------------------
!> @brief   Reorients the cubed-sphere to be compatible with
!>          the orientation assumed by the LFRic infrastructure.
!>          (Private Routine)
!> @details Performs circular shifts on appropriate panels
!>
!> @param[in,out] gen_cube             Generator strategy for a cubed-sphere
!> @param[in]     panel_rotation_array Specifies the amount of rotation for
!>                                     each panel
!-------------------------------------------------------------------------------
subroutine orient_lfric(self, panel_rotation_array)

  implicit none

  class(gen_eave_type), intent(inout) :: self
  integer(i_def),       intent(in)    :: panel_rotation_array(:)

  integer(i_def) :: panel_id
  integer(i_def) :: shift_value

  ! Apply panel orientation to match Cubed-Sphere strategy
  do panel_id=1, size(panel_rotation_array)

    shift_value = panel_rotation_array(panel_id)
    self%domain_extents(:,:,panel_id) = cshift( self%domain_extents(:, :, panel_id), &
                                                shift_value, 2 )
    self%nodes_on_cell(:,:,panel_id)  = cshift( self%nodes_on_cell(:,:,panel_id),    &
                                                shift_value, 1 )
    self%cell_next(:,:,panel_id)      = cshift( self%cell_next(:,:,panel_id),        &
                                                shift_value, 1 )
    self%edges_on_cell(:,:,panel_id)  = cshift( self%edges_on_cell(:,:,panel_id),    &
                                                shift_value, 1 )
  end do

  return
end subroutine orient_lfric


!-----------------------------------------------------------------------------
!> @brief Returns the number of panels in the mesh topology.
!> @description Panels are a subset of cells in the mesh domain which may
!>              exhibit common properties.
!> @return answer Number of panels resulting from this generation strategy.
!-----------------------------------------------------------------------------
function get_number_of_panels( self ) result( answer )

  implicit none

  class(gen_eave_type), intent(in) :: self

  integer(i_def) :: answer

  answer = self%npanels

end function get_number_of_panels


!-----------------------------------------------------------------------------
!> @brief Returns mesh metadata information.
!>
!> @param[out]  mesh_name          Optional, Name of mesh instance to generate
!> @param[out]  geometry           Optional, Mesh domain surface type.
!> @param[out]  topology           Optional, Mesh boundary/connectivity type
!> @param[out]  coord_sys          Optional, Coordinate system to position nodes.
!> @param[out]  periodic_x         Optional, Domain periodicity in x-axis
!> @param[out]  periodic_y         Optional, Domain periodicity in y-axis
!> @param[out]  edge_cells_x       Optional, Number of panel edge cells (x-axis).
!> @param[out]  edge_cells_y       Optional, Number of panel edge cells (y-axis).
!> @param[out]  constructor_inputs Optional, Inputs used to create this mesh from
!>                                           the this ugrid_generator_type
!> @param[out]  nmaps              Optional, Number of maps to create with this mesh
!>                                           as source mesh.
!> @param[out]  rim_depth          Optional, Depth of LBC mesh rim (in cells)
!> @param[out]  void_id          Optional, Cell ID for null connectivity.
!> @param[out]  target_mesh_names  Optional, Mesh names of the target meshes that
!>                                           this mesh has maps for.
!> @param[out]  maps_edge_cells_x  Optional, Number of panel edge cells (x-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out]  maps_edge_cells_y  Optional, Number of panel edge cells (y-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out]  north_pole         Optional, [Longitude, Latitude] of north pole
!>                                           used for domain orientation (degrees)
!> @param[out]  null_island        Optional, [Longitude, Latitude] of null
!>                                           island used for domain orientation (degrees)
!-----------------------------------------------------------------------------
subroutine get_metadata( self,               &
                         mesh_name,          &
                         geometry,           &
                         topology,           &
                         coord_sys,          &
                         periodic_xy,        &
                         edge_cells_x,       &
                         edge_cells_y,       &
                         constructor_inputs, &
                         nmaps,              &
                         rim_depth,          &
                         eave_depth,         &
                         void_cell,          &
                         target_mesh_names,  &
                         maps_edge_cells_x,  &
                         maps_edge_cells_y,  &
                         north_pole,         &
                         null_island,        &
                         equatorial_latitude )

  implicit none

  class(gen_eave_type), intent(in)  :: self

  character(str_def), optional, intent(out) :: mesh_name
  character(str_def), optional, intent(out) :: geometry
  character(str_def), optional, intent(out) :: topology
  character(str_def), optional, intent(out) :: coord_sys
  logical(l_def),     optional, intent(out) :: periodic_xy(2)
  integer(i_def),     optional, intent(out) :: edge_cells_x
  integer(i_def),     optional, intent(out) :: edge_cells_y
  integer(i_def),     optional, intent(out) :: nmaps
  integer(i_def),     optional, intent(out) :: rim_depth
  integer(i_def),     optional, intent(out) :: eave_depth
  integer(i_def),     optional, intent(out) :: void_cell

  character(str_longlong), optional, intent(out) :: constructor_inputs

  character(str_def), allocatable, optional, intent(out) :: target_mesh_names(:)
  integer(i_def),     allocatable, optional, intent(out) :: maps_edge_cells_x(:)
  integer(i_def),     allocatable, optional, intent(out) :: maps_edge_cells_y(:)

  real(r_def), optional, intent(out) :: north_pole(2)
  real(r_def), optional, intent(out) :: null_island(2)
  real(r_def), optional, intent(out) :: equatorial_latitude

  if (present(mesh_name))    mesh_name      = trim(self%mesh_name)
  if (present(geometry))     geometry       = key_from_geometry(self%geometry)
  if (present(topology))     topology       = key_from_topology(self%topology)
  if (present(coord_sys))    coord_sys      = key_from_coord_sys(self%coord_sys)
  if (present(periodic_xy))  periodic_xy    = .true.
  if (present(edge_cells_x)) edge_cells_x   = self%edge_cells
  if (present(edge_cells_y)) edge_cells_y   = self%edge_cells
  if (present(nmaps))        nmaps          = self%nmaps
  if (present(rim_depth))    rim_depth      = imdi
  if (present(eave_depth))   eave_depth     = self%eave_depth
  if (present(void_cell))    void_cell      = self%VOID_ID

  if (present(north_pole))     north_pole(:)  = radians_to_degrees * self%north_pole(:)
  if (present(null_island))    null_island(:) = radians_to_degrees * self%null_island(:)
  if (present(equatorial_latitude)) &
                          equatorial_latitude = radians_to_degrees * self%equatorial_latitude

  if (present(constructor_inputs)) constructor_inputs = trim(self%constructor_inputs)

  if (self%nmaps > 0) then
    if (present(target_mesh_names)) target_mesh_names = self%parent_mesh_name
    if (present(maps_edge_cells_x)) maps_edge_cells_x = 1
    if (present(maps_edge_cells_x)) maps_edge_cells_y = 1
  end if

  return
end subroutine get_metadata

!-------------------------------------------------------------------------------
!> @brief Subroutine to manually deallocate any memory used by the object.
!-------------------------------------------------------------------------------
subroutine clear(self)

  implicit none

  class (gen_eave_type), intent(inout) :: self

  if (allocated(self%cell_next))     deallocate( self%cell_next     )
  if (allocated(self%nodes_on_cell)) deallocate( self%nodes_on_cell )
  if (allocated(self%edges_on_cell)) deallocate( self%edges_on_cell )
  if (allocated(self%nodes_on_edge)) deallocate( self%nodes_on_edge )
  if (allocated(self%node_coords))   deallocate( self%node_coords   )
  if (allocated(self%cell_coords))   deallocate( self%cell_coords   )

  if (allocated(self%global_mesh_maps)) then
    call self%global_mesh_maps%clear()
    deallocate( self%global_mesh_maps )
  end if


  return
end subroutine clear

!-----------------------------------------------------------------------------
!> @brief Finaliser for the cubedsphere ugrid mesh generator (gen_eave_type)
!-----------------------------------------------------------------------------
subroutine gen_eave_final(self)

  implicit none

  type(gen_eave_type), intent(inout) :: self

  call self%clear()

  return
end subroutine gen_eave_final


subroutine get_panel_edge_cell_ids( edge_cells, panel_edge_cells )

  implicit none

!
!          North
!      o---->>>----o         Cell ids on panel edges are
!      |           |         listed in direction shown:
! West Y   Panel   Y East
!      |           |         panel_edge_cells( cell_ids,
!      o---->>>----o                           side of the panel,
!          South                               panel number )
!

  integer(i_def), intent(in)  :: edge_cells
  integer(i_def), intent(out) :: panel_edge_cells(edge_cells,4,NPANELS)

  integer(i_def) :: i, cpp, panel

  cpp = edge_cells*edge_cells

  do panel=1, NPANELS
    ! Panel edge ordering W,S,E,N
    do i=1, edge_cells
      panel_edge_cells(i,W,panel) = (panel-1)*cpp + edge_cells*(i-1) + 1
      panel_edge_cells(i,S,panel) = (panel-1)*cpp + (cpp-edge_cells) + i
      panel_edge_cells(i,E,panel) = (panel-1)*cpp + edge_cells*(i)
      panel_edge_cells(i,N,panel) = (panel-1)*cpp + i
    end do
  end do
end subroutine get_panel_edge_cell_ids

subroutine calc_cell_maps(self)

  implicit none

  class(gen_eave_type), intent(inout) :: self

  integer(i_def) :: panel_edge_cells(2)
  integer(i_def) :: full_edge_cells(2)
  integer(i_def) :: flap_depth(2)

  integer(i_def), allocatable :: tmp_panel_ids(:)
  integer(i_def), allocatable :: cs_panel_ids(:,:)


  integer(i_def) :: cpp
  integer(i_def) :: full_cpp


  integer(i_def) :: start_id
  integer(i_def) :: start_row, end_row
  integer(i_def) :: start_index, end_index

  integer(i_def) :: panel_row, panel_id, i

  panel_edge_cells(:) = self%edge_cells
  flap_depth(:)       = self%eave_depth

  cpp = panel_edge_cells(1)*panel_edge_cells(2)

  allocate( tmp_panel_ids(cpp) )
  allocate( cs_panel_ids( panel_edge_cells(1), &
                          panel_edge_cells(2)) )

  full_edge_cells(1)  = panel_edge_cells(1) + 2* flap_depth(1)
  full_edge_cells(2)  = panel_edge_cells(2) + 2* flap_depth(2)


  full_cpp = full_edge_cells(1)*full_edge_cells(2)

  if ( allocated(self%cell_maps)) deallocate(self%cell_maps)
  allocate( self%cell_maps(full_cpp, npanels) )

  ! Set array to be void cell as this will blank out the map everywhere
   self%cell_maps(:,:) = self%void_id

  ! Get the starting index that would correspond to the real cubedsphere panel
  ! rows down  = self%flap_depth(2) + 1
  ! coloums in = self%flap_depth(1) + 1
  ! So start with local id:
  start_row = flap_depth(2) + 1
  end_row   = start_row + panel_edge_cells(2) - 1


  !  start_column = self%flap_depth(1) + 1
  !  end_column   = self%flap_depth(1) + panel_edge_cells(1)
  do panel_id=1, npanels

    ! Generate expected cubedsphere ids for the
    ! given panel.
    start_id = (panel_id-1)*cpp + 1

    do i=1, cpp
      ! Initialise panel ids to map
      tmp_panel_ids(i) = start_id
      start_id = start_id + 1
    end do

    cs_panel_ids = reshape( tmp_panel_ids, [ panel_edge_cells(1), &
                                             panel_edge_cells(2)] )

    ! Create eave map to cs_panel_ids
    panel_row = 0
    do i=start_row, end_row
      panel_row   = panel_row + 1
      start_index = (i-1)*full_edge_cells(1) + flap_depth(1) + 1
      end_index   = start_index + panel_edge_cells(1) - 1
      self%cell_maps(start_index:end_index,panel_id) = cs_panel_ids(:, panel_row)
    end do

  end do

end subroutine calc_cell_maps

subroutine compare_eave(self, cs)

  use coord_transform_mod, only: rebase_longitude_range
  use gencube_ps_mod,      only: gencube_ps_type

  implicit none

  class(gen_eave_type), intent(in) :: self

  class(gencube_ps_type), intent(in) :: cs

  integer(i_def) :: panel_id
  integer(i_def) :: nodes_per_cell

  real(r_def),    allocatable :: node_coordinates(:,:)
  real(r_def),    allocatable :: cell_coordinates(:,:)
  real(r_def),    allocatable :: domain_extents(:,:)
  integer(i_def), allocatable :: face_node_connectivity(:,:)
  integer(i_def), allocatable :: face_edge_connectivity(:,:)
  integer(i_def), allocatable :: face_face_connectivity(:,:)
  integer(i_def), allocatable :: edge_node_connectivity(:,:)
  character(str_def) :: coord_units_x
  character(str_def) :: coord_units_y

  character(str_def) :: cmessage

  integer(i_def) :: node_index
  integer(i_def) :: cs_cell_id
  integer(i_def) :: cs_node_id
  integer(i_def) :: eave_cell_id
  integer(i_def) :: eave_node_id
  logical :: out_of_tol

  real(r_def), parameter :: tol = 1.0e-10_r_def
  real(r_def) :: cs_node_coords(2)
  real(r_def) :: eave_node_coords(2)
  real(r_def) :: diff_lat
  real(r_def) :: diff_lon

  nodes_per_cell = 4

  call cs%get_connectivity( face_node_connectivity, &
                            face_edge_connectivity, &
                            face_face_connectivity, &
                            edge_node_connectivity )

  call cs%get_coordinates(  node_coordinates, &
                            cell_coordinates, &
                            domain_extents,   &
                            coord_units_x,    &
                            coord_units_y     )

  do panel_id=1, 6

    out_of_tol=.false.

    write(cmessage,'(A,E9.2)') 'All nodes on panel are within tolerance ', tol
    write(log_scratch_space, '(A,I0,A)') 'PANEL ', panel_id,': '//trim(cmessage)

    do eave_cell_id=1, self%cpp
      cs_cell_id = self%cell_maps(eave_cell_id, panel_id)

      if (cs_cell_id /= self%void_id) then
        do node_index=1, nodes_per_cell
          cs_node_id   = face_node_connectivity( node_index, cs_cell_id )
          eave_node_id = self%nodes_on_cell( node_index, eave_cell_id, panel_id )

          cs_node_coords(:)   = node_coordinates(:, cs_node_id)
          eave_node_coords(:) = self%node_coords(:, eave_node_id, panel_id)


          diff_lon = cs_node_coords(1) - eave_node_coords(1)*radians_to_degrees
          diff_lat = cs_node_coords(2) - eave_node_coords(2)*radians_to_degrees
          if ( ( ABS(diff_lon) > tol ) .or. &
               ( ABS(diff_lat) > tol ) ) then

            out_of_tol=.true.
            write(cmessage, '(A,I0,A,E5.1,A,2F20.6)') 'PANEL ', panel_id,': '// &
                'Node coordinates (degrees) out of tolerance cf (',tol,'): ',   &
                diff_lon, diff_lat
            call log_event(cmessage, log_level_error)

          end if

        end do ! Nodes per cell
      end if ! If not a void cell

    end do ! eave cells

  end do ! panels

end subroutine compare_eave


function is_generated(self) result(answer)

  implicit none
  class(gen_eave_type), intent(in) :: self

  logical :: answer

  answer = self%generated

  return
end function is_generated


subroutine populate_ugrid_2d( self, ugrid_2d, panel_id )

  use ugrid_2d_mod,    only: ugrid_2d_type
  use mesh_config_mod, only: key_from_coord_sys

  implicit none

  class(gen_eave_type), intent(inout) :: self
  class(ugrid_2d_type), intent(inout) :: ugrid_2d

  integer(i_def), intent(in) :: panel_id
  character(str_def) :: name

  integer(i_def), allocatable :: cell_map(:,:,:)
  real(r_def) :: factor

  character(str_def) :: coord_units(2)

  call ugrid_2d%set_dimensions(        &
               num_nodes = self%npp,   &
               num_edges = self%epp,   &
               num_faces = self%cpp,   &
               num_nodes_per_face = 4, &
               num_edges_per_face = 4, &
               num_nodes_per_edge = 4, &
               max_num_faces_per_node = 4 )

  if (self%coord_units(1) == 'radians') then
    factor = radians_to_degrees
  else
    factor = 1.0_r_def
  end if
  coord_units = 'degrees'

  call ugrid_2d%set_coords(                                     &
           node_coords = self%node_coords(:,:,panel_id)*factor, &
           face_coords = self%cell_coords(:,:,panel_id)*factor, &
           north_pole  = self%north_pole*factor,                &
           null_island = self%null_island*factor,               &
           coord_sys   = key_from_coord_sys(self%coord_sys),    &
           units_xy    = coord_units )

  call ugrid_2d%set_connectivity(                              &
           nodes_on_faces = self%nodes_on_cell(:,:,panel_id),  &
           edges_on_faces = self%edges_on_cell(:,:,panel_id),  &
           faces_on_faces = self%cell_next(:,:, panel_id),     &
           nodes_on_edges = self%nodes_on_edge(:,:, panel_id), &
           void_cell      = self%void_id )

  write(name,'(A,I0)') trim(self%parent_mesh_name)//'-eave-', panel_id

  call ugrid_2d%set_metadata(                                             &
           mesh_name          = name,                                     &
           geometry           = key_from_geometry(self%geometry),         &
           topology           = key_from_topology(self%topology),         &
           npanels            = 1,                                        &
           domain_extents     = self%domain_extents(:,:,panel_id)*factor, &
           periodic_xy        = self%periodicity,                         &
           edge_cells_x       = self%eave_edge_cells,                     &
           edge_cells_y       = self%eave_edge_cells,                     &
           eave_depth         = self%eave_depth,                          &
           constructor_inputs = self%constructor_inputs,                  &
           nmaps              = 1,                                        &
           target_mesh_names  = [self%parent_mesh_name] )

  allocate(cell_map(1,1,self%cpp))
  cell_map = reshape( self%cell_maps(:,panel_id), [1,1,self%cpp] )

  if ( .not. allocated(self%global_mesh_maps) ) then
    allocate(self%global_mesh_maps)
  else
    call self%global_mesh_maps%clear()
  end if

  call self%global_mesh_maps%add_global_mesh_map(1, 2, cell_map)
  call ugrid_2d%set_mesh_maps( self%global_mesh_maps )

end subroutine populate_ugrid_2d

end module gen_eave_mod
