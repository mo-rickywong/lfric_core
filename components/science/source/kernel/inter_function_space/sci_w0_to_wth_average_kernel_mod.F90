!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the average of nearest W0 space to Wtheta space.
!> @details Kernel to average a W0 lower-level field to the horizontal parts
!!          of a Wtheta field. The method is valid for the bottom DoFs, horizontal
!!          and lower level of the lowest-order finite elements on a cubed-sphere mesh.

module sci_w0_to_wth_average_kernel_mod

  use argument_mod,      only : arg_type,                   &
                                GH_FIELD, GH_REAL,          &
                                GH_READWRITE, GH_READ,      &
                                CELL_COLUMN
  use constants_mod,     only : i_def, r_def
  use fs_continuity_mod, only : W0, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: w0_to_wth_average_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                     &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W0)      &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: w0_to_wth_average_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: w0_to_wth_average_code

  contains

  !> @brief Computes a 1-2-1 Filter from W0 to Wtheta space.
  !> @param[in]     nlayers         Number of layers
  !> @param[in,out] field_wth       Output field from Filter on Wtheta space
  !> @param[in]     field_w0        Input field for filter on W0 spaceal
  !> @param[in]     ndf_wtheta      Number of degrees of freedom per cell for Wtheta
  !> @param[in]     undf_wtheta     Number of unique degrees of freedom for Wtheta
  !> @param[in]     map_wtheta      Dofmap for the cell at the base of the column for Wtheta
  !> @param[in]     ndf_w0          Number of degrees of freedom per cell for W0
  !> @param[in]     undf_w0         Number of unique degrees of freedom for W0
  !> @param[in]     map_w0          Dofmap for the cell at the base of the column for W0
  subroutine w0_to_wth_average_code(nlayers,                               &
                                    field_wth, field_w0,                   &
                                    ndf_wtheta, undf_wtheta, map_wtheta,   &
                                    ndf_w0, undf_w0, map_w0)
    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta
    integer(kind=i_def), intent(in) :: ndf_w0, undf_w0

    integer(kind=i_def), intent(in), dimension(ndf_wtheta) :: map_wtheta
    integer(kind=i_def), intent(in), dimension(ndf_w0)     :: map_w0

    real(kind=r_def), intent(inout), dimension(undf_wtheta) :: field_wth
    real(kind=r_def), intent(in),    dimension(undf_w0)     :: field_w0

    ! Internal variables
    integer(kind=i_def) :: df, k

    do k = 0, nlayers
      field_wth(map_wtheta(1) + k) = 0.0_r_def
      do df = 1,4 ! Loop at the Bottom
        field_wth(map_wtheta(1) + k) = field_wth(map_wtheta(1) + k) + &
                    field_w0(map_w0(df) + k)/4.0_r_def
      end do
    end do

  end subroutine w0_to_wth_average_code

end module sci_w0_to_wth_average_kernel_mod
