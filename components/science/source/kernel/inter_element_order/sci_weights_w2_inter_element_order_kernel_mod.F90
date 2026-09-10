!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes weights for the W2 inter element-order mapping.
!> @details Compute the weights for the mapping higher to lower-order as the
!!          basis functions of the higher-order space evaluated at the dofs
!!          of the lower-order space inside the coarse mesh cells. Inverts this
!!          matrix to get the weights for the mapping lower to higher-order.
!!          The segment centre quadrature rule gives the fine dof locations.

module sci_weights_w2_inter_element_order_kernel_mod

use argument_mod,          only: arg_type, func_type,       &
                                 GH_FIELD, GH_REAL,         &
                                 GH_READ, GH_WRITE,         &
                                 GH_COARSE, GH_FINE,        &
                                 GH_BASIS,                  &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 ANY_DISCONTINUOUS_SPACE_2, &
                                 ANY_DISCONTINUOUS_SPACE_3, &
                                 CELL_COLUMN,               &
                                 GH_QUADRATURE_XYoZ
use constants_mod,         only: i_def, r_def, imdi
use kernel_mod,            only: kernel_type
use matrix_invert_mod,     only: matrix_invert_plu
use reference_element_mod, only: W, S, E, N, B, T

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>
! Note: W2 spaces are currently declared as ANY_DISCONTINUOUS_SPACE due to
! #4052.
type, public, extends(kernel_type) :: weights_w2_inter_element_order_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/ &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1, mesh_arg=GH_FINE  ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3, mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3, mesh_arg=GH_COARSE)  &
       /)
  type(func_type) :: meta_funcs(1) = (/func_type(ANY_DISCONTINUOUS_SPACE_2, GH_BASIS)/)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: weights_w2_inter_element_order_code
end type weights_w2_inter_element_order_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: weights_w2_inter_element_order_code

contains

  !> @brief Computes weights for the scalar inter element-order mapping
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in]     dummy_fine               An unused fine field supplied
  !!                                         to ensure that appropriate fine
  !!                                         field properties are provided.
  !> @param[in]     dummy_coarse             An unused coarse field supplied
  !!                                         to ensure that appropriate coarse
  !!                                         field properties are provided.
  !> @param[in,out] weights_fe_to_fv         The coarse grid 2D field that
  !!                                         contains the weights to be
  !!                                         computed here for mapping high
  !!                                         order to low order
  !> @param[in,out] weights_fv_to_fe         The coarse grid 2D field that
  !!                                         contains the weights to be
  !!                                         computed here for mapping low
  !!                                         order to high order
  !> @param[in]     ndf_fine                 Num of DoFs per cell on the fine
  !!                                         grid
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     ndf_coarse               Num of DoFs per cell on the coarse
  !!                                         grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     basis_coarse             Basis functions on the coarse grid
  !> @param[in]     undf_weights             Total num of DoFs on the weights
  !!                                         field for this mesh partition
  !> @param[in]     map_weights              DoFmap of cells on the weights
  !!                                         fields
  !> @param[in]     nqp_xy                   Number of quadrature points in the
  !!                                         horizontal
  !> @param[in]     nqp_z                    Number of quadrature points in the
  !!                                         vertical
  !> @param[in]     weights_xy               Weights for the horizontal
  !!                                         quadrature points
  !> @param[in]     weights_z                Weights for the vertical
  !!                                         quadrature points
  subroutine weights_w2_inter_element_order_code(nlayers,                 &
                                                 cell_map,                &
                                                 ncell_fine_per_coarse_x, &
                                                 ncell_fine_per_coarse_y, &
                                                 ncell_fine,              &
                                                 dummy_fine,              &
                                                 dummy_coarse,            &
                                                 weights_fe_to_fv,        &
                                                 weights_fv_to_fe,        &
                                                 ndf_fine,                &
                                                 undf_fine,               &
                                                 map_fine,                &
                                                 ndf_coarse,              &
                                                 undf_coarse,             &
                                                 map_coarse,              &
                                                 basis_coarse,            &
                                                 undf_weights,            &
                                                 map_weights,             &
                                                 nqp_xy,                  &
                                                 nqp_z,                   &
                                                 weights_xy,              &
                                                 weights_z)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf_fine, ndf_coarse
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse
    integer(kind=i_def), intent(in) :: undf_weights
    integer(kind=i_def), intent(in) :: nqp_xy, nqp_z

    real(kind=r_def),    intent(in) :: basis_coarse(3, ndf_coarse, nqp_xy, nqp_z)
    real(kind=r_def),    intent(in) :: weights_xy(nqp_xy)
    real(kind=r_def),    intent(in) :: weights_z(nqp_z)

    ! Fields
    real(kind=r_def),    intent(inout) :: weights_fe_to_fv(undf_weights)
    real(kind=r_def),    intent(inout) :: weights_fv_to_fe(undf_weights)
    real(kind=r_def),    intent(in)    :: dummy_fine(undf_fine)
    real(kind=r_def),    intent(in)    :: dummy_coarse(undf_coarse)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ncell_fine, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf_coarse)
    integer(kind=i_def), intent(in) :: map_weights(ndf_coarse)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, &
                                                ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: df_c, df_f, x_idx, y_idx
    integer(kind=i_def) :: j, k, x, y, qp_xy, qp_z, component
    integer(kind=i_def) :: nqp_x

    integer(kind=i_def) :: map_idx(ndf_fine, ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), allocatable :: map_quad(:,:)

    real(kind=r_def) :: weights_matrix(ndf_coarse, ndf_coarse)
    real(kind=r_def) :: inv_weights_matrix(ndf_coarse, ndf_coarse)

    real(kind=r_def) :: coeff

    !---------------------------------------------------------------------------
    ! Explanation
    !---------------------------------------------------------------------------
    ! This kernel takes the basis functions of a higher-order W2 space,
    ! evaluated using the segment centre quadrature rule on a lattice.
    ! An appropriate choice of lattice dimensions results in all of the dofs
    ! of (element_order_h + 1)^2 embedded, lowest-order, W2 cells lying on the
    ! nodes of the lattice.
    !
    ! Let W be the wind field. For the ith dof of the higher-order cell d_i,
    ! let b_i be the corresponding basis function and n_i be the normal
    ! direction of the dof. Then by basis function decomposition:
    !
    !   W(x) = sum_i b_i(x) * (n_i \dot W(d_i))
    !
    ! Letting x be dof of a lowest-order cell, we get that the weights for the
    ! mapping of the wind field evaluated at the higher-order dofs to the
    ! wind field evaluated at the lowest-order dofs are the basis functions
    ! evaluated at the dofs of the lowest-order cells.
    !
    ! As a result, from the basis functions evaluated on the lattice, we can
    ! extract the values at the dofs of the lowest-order cells, ignoring
    ! unused lattice evaluations.
    !
    ! The reference image for element_order_h=1 is:
    !
    !   Top-down view:                           Side view:
    !   1    2    3    4    5                    1    2    3    4    5
    !             N                                         T
    !   -----v---------v-----    5   |            -----^---------^-----    3
    !   |         |         |        |            |         |         |
    !   >    *    >    *    >    4   |            |         |         |
    !   |         |         |        |            |         |         |
    ! W -----v---------v----- E  3  y_idx       W >    *    >    *    > E  2
    !   |         |         |        |            |         |         |
    !   >    *    >    *    >    2   |            |         |         |
    !   |         |         |        |            |         |         |
    !   -----v---------v-----    1   v            -----^---------^-----    1
    !             S                                         B
    !   --------x_idx------->
    !
    ! The arrows represent lowest-order dofs in that direction, and the stars
    ! represent lowest-order dofs normal to the plane. Each integer coordinate
    ! of the lattice represents a point at which the basis functions have been
    ! evaluated.
    !
    ! We have 4 embedded cells per higher-order cell, so we need a lattice of
    ! dimension 5x5x3 to evaluate the basis functions at every lowest-order
    ! dof.
    !
    ! We call the higher-order function space the "coarse space", since it is
    ! stored on a coarse mesh, and likewise the lowest-order function space is
    ! called the "fine space".

    !---------------------------------------------------------------------------
    ! Define maps
    !---------------------------------------------------------------------------
    ! This map takes a dof in the a lowest-order cell on the fine mesh, and the
    ! coordinates of the cell in the local embedding. It returns an index in
    ! the range 1 and ndf_coarse.
    !
    ! Dofs which are not used are set to imdi, which is large enough to trigger
    ! a segmentation fault when weights_matrix is indexed with it.
    j = 1
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do df_f = 1, ndf_fine
          if (df_f == E .and. x_idx /= ncell_fine_per_coarse_x) then
            ! East, internal dof. Will be counted by the West dof of the next
            ! cell in the x-direction.
            map_idx(df_f, x_idx, y_idx) = imdi
          else if (df_f == N .and. y_idx /= 1) then
            ! North, internal dof. Will be counted by the South dof of the
            ! next cell in the y-direction.
            map_idx(df_f, x_idx, y_idx) = imdi
          else
            map_idx(df_f, x_idx, y_idx) = j
            j = j + 1
          end if
        end do
      end do
    end do

    ! Map from (x,y) coordinates on the lattice to the xy quadrature points, as
    ! in the quadrature_xyoz initialisation routine.
    nqp_x = int(sqrt(real(nqp_xy, r_def)), i_def)
    allocate(map_quad(nqp_x, nqp_x))
    k = 1
    do x = 1, nqp_x
      do y = 1, nqp_x
        map_quad(x, y) = k
        k = k + 1
      end do
    end do

    !---------------------------------------------------------------------------
    ! Extract the basis functions evaluated at the fine dofs
    !---------------------------------------------------------------------------
    ! Now we use the maps to extract the basis function evaluations which
    ! correspond to the dofs of the lowest-order cells.
    !
    ! We may assume the cells are square, since element_order_h modifies both
    ! directions.
    !
    ! For East-West or North-South oriented dofs, we only need to evaluate the
    ! basis functions at the second quadrature point in the z-direction since
    ! this corresponds to the second layer (see diagram above). This kernel
    ! will only be called for element_order_v = 0 so this can be assumed.
    ! Likewise, for Bottom-Top oriented dofs, we only need to evaluate the
    ! basis functions at the first and third quadrature points in the
    ! z-direction.
    weights_matrix(:,:) = 0.0_r_def
    coeff = 1.0_r_def/real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, r_def)
    do x = 1, nqp_x
      do y = 1, nqp_x
        if (mod(x+1,2)==0 .and. mod(y+1,2)==0) then
          ! Odd coordinates are cell corners, no W2 dofs so do nothing
          cycle
        end if
        qp_xy = map_quad(x, y)

        if (mod(x+1,2)==0 .and. mod(y,2)==0) then
          ! Odd x so an East or West face,  even y so at the centre of a face.
          ! Therefore East-West oriented dof
          component = 1
          qp_z = 2
          if (x == nqp_x) then
            ! On East boundary, use East dof
            df_f = E
            x_idx = (nqp_x-1)/2
            y_idx = y/2
          else
            ! Not on East boundary, store on West dof to avoid overcounting
            df_f = W
            x_idx = (x+1)/2
            y_idx = y/2
          end if
          j = map_idx(df_f, x_idx, y_idx)

        else if (mod(x,2)==0 .and. mod(y+1,2)==0) then
          ! Odd y so a North or South face, even x so at the centre of a face.
          ! Therefore North-South oriented dof
          component = 2
          qp_z = 2
          if (y == nqp_x) then
            ! On North boundary, use North dof
            df_f = N
            x_idx = x/2
            y_idx = 1
          else
            ! Not on North boundary, store on South dof to avoid overcounting
            df_f = S
            x_idx = x/2
            y_idx = (nqp_x-y)/2
          end if
          j = map_idx(df_f, x_idx, y_idx)

        else if (mod(x,2)==0 .and. mod(y,2)==0) then
          ! Even x and y so at the centre of a cell. Therefore we can evaluate
          ! the Bottom dofs here and the Top dofs in the usual loop
          component = 3
          x_idx = x/2
          y_idx = y/2

          df_f = B
          qp_z = 1
          j = map_idx(df_f, x_idx, y_idx)
          do df_c = 1, ndf_coarse
            weights_matrix(df_c, j) = basis_coarse(component,df_c,qp_xy,qp_z)*coeff
          end do

          df_f = T
          qp_z = 3
          j = map_idx(df_f, x_idx, y_idx)
        end if

        ! Loop over the higher-order dofs to fill in the rows of the matrix
        ! with the correct component of the basis functions, evaluated at the
        ! quadrature points we specified above
        do df_c = 1, ndf_coarse
            weights_matrix(df_c, j) = basis_coarse(component,df_c,qp_xy,qp_z)*coeff
        end do
      end do
    end do

    !---------------------------------------------------------------------------
    ! Invert the weights matrix on each coarse cell
    !---------------------------------------------------------------------------
    ! The resultant values represent the weights for the mapping from the
    ! lowest to the higher-order space. The matrix does not typically have an
    ! LU decomposition, so use matrix inversion with pivoting
    call matrix_invert_plu(weights_matrix, inv_weights_matrix, ndf_coarse)

    !---------------------------------------------------------------------------
    ! Flatten the weights matrices
    !---------------------------------------------------------------------------
    do df_c = 1, ndf_coarse
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          do df_f = 1, ndf_fine
            ! Only evaluate East and North dofs at the East and North edges of
            ! the coarse cell to avoid overcounting
            if (df_f == E .and. x_idx /= ncell_fine_per_coarse_x) cycle
            if (df_f == N .and. y_idx /= 1) cycle

            j = map_idx(df_f, x_idx, y_idx)

            weights_fe_to_fv(map_weights(df_c) + j - 1) = &
              weights_matrix(df_c, j)

            ! Take the transpose of the inverted matrix
            weights_fv_to_fe(map_weights(df_c) + j - 1) = &
              inv_weights_matrix(j, df_c)
          end do
        end do
      end do
    end do

    !---------------------------------------------------------------------------
    ! Tidy up
    !---------------------------------------------------------------------------
    if (allocated(map_quad)) deallocate(map_quad)

  end subroutine weights_w2_inter_element_order_code

end module sci_weights_w2_inter_element_order_kernel_mod
