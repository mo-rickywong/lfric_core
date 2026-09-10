!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief An algorithm for performing the coupling. In a real model
!>       the coupled data could be processed with kernel calls, but
!>       here it is a simple test, so there is no processing
module coupled_alg_mod

  use constants_mod,                  only: i_def, r_def, str_longlong
  use coupler_exchange_2d_mod,        only: coupler_exchange_2d_type
  use coupling_mod,                   only: coupling_type, &
                                            get_coupling_from_collection
  use driver_modeldb_mod,             only: modeldb_type
  use field_collection_iterator_mod,  only: field_collection_iterator_type
  use field_collection_mod,           only: field_collection_type
  use field_mod,                      only: field_type
  use field_parent_mod,               only: field_parent_type
  use log_mod,                        only: log_event,         &
                                            log_scratch_space, &
                                            log_level_info,    &
                                            log_level_error
  use sci_field_minmax_alg_mod,       only: get_field_minmax
  implicit none

  private

  public :: coupled_alg

contains

  !> @details An algorithm for performing the coupling
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine coupled_alg(modeldb)

    implicit none

    ! Prognostic fields
    type(modeldb_type), intent(inout)   :: modeldb

    character(str_longlong), pointer       :: cpl_component_name

#ifdef MCT
    type( field_type ),pointer             :: field_1
    type( field_type ),pointer             :: field_2
    type(coupling_type), pointer           :: coupling_ptr
    integer(i_def), pointer                :: local_index(:)
    type( field_collection_type ), pointer :: cpl_snd_2d
    type( field_collection_type ), pointer :: cpl_rcv_2d
    type( field_collection_iterator_type)  :: iter
    class( field_parent_type ), pointer    :: field
    type( field_type ), pointer            :: field_ptr
    type(coupler_exchange_2d_type)         :: coupler_exchange_2d
    integer(i_def)                         :: ierror
    real(r_def)                            :: field_min, field_max
#endif

    call log_event( "coupled: Running algorithm", log_level_info )

#ifdef MCT

    ! Extract the version that was actaully placed in the collection
    coupling_ptr => get_coupling_from_collection(modeldb%values, "coupling" )
    local_index => coupling_ptr%get_local_index()

    call modeldb%values%get_value("cpl_name", cpl_component_name)
    if (trim(cpl_component_name) == "lfric_o") then
      ! Load in data to be coupled with other component.
      ! For now, just use fixed data
      cpl_snd_2d => modeldb%fields%get_field_collection("cpl_snd_2d")
      call cpl_snd_2d%get_field("field_1",  field_1)
      call get_field_minmax(field_1, field_min, field_max)
      write(log_scratch_space, &
        '("Sent field (field_1) min= ",f8.3,", max=",f8.3)' ) &
        field_min, field_max
      call log_event( log_scratch_space, log_level_info )

      call iter%initialise(cpl_snd_2d)
      do
        if ( .not. iter%has_next() ) exit
        field => iter%next()
        select type(field)
        type is (field_type)
          field_ptr => field
          call coupler_exchange_2d%initialise(field_ptr, local_index)
          call coupler_exchange_2d%set_time(modeldb%clock)
          call coupler_exchange_2d%copy_from_lfric(ierror)
          call coupler_exchange_2d%clear()
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM: coupled_alg: field ", &
                trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, log_level_error )
        end select
      end do
    end if

    if (trim(cpl_component_name) == "lfric_i") then
      cpl_rcv_2d => modeldb%fields%get_field_collection("cpl_rcv_2d")
      call iter%initialise(cpl_rcv_2d)
      do
        if ( .not. iter%has_next() ) exit
        field => iter%next()
        select type(field)
        type is (field_type)
          field_ptr => field
          call coupler_exchange_2d%initialise(field_ptr, local_index)
          call coupler_exchange_2d%set_time(modeldb%clock)
          ! Receive field from coupler to LFRic
          call coupler_exchange_2d%copy_to_lfric(ierror)
          call coupler_exchange_2d%clear()
        class default
          write(log_scratch_space, '(2A)' ) "PROBLEM: coupled_alg: field ", &
                trim(field%get_name())//" is NOT field_type"
          call log_event( log_scratch_space, log_level_error )
        end select
      end do
      call cpl_rcv_2d%get_field("field_2",  field_2)
      call get_field_minmax(field_2, field_min, field_max)
      write(log_scratch_space, &
        '("Received field (field_2) min= ",f8.3,", max=",f8.3)' ) &
        field_min, field_max
      call log_event( log_scratch_space, log_level_info )
    end if
#endif

    call log_event( "coupled: finished algorithm", log_level_info )

  end subroutine coupled_alg

end module coupled_alg_mod
