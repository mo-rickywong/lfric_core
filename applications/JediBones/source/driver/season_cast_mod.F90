!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Initialisation functionality for the JediBones miniapp

!> @details Handles init of prognostic fields and through the call to
!>          runtime_contants the coordinate fields and fem operators

module season_cast_mod

  use config_mod, only: config_type
  use log_mod,    only: log_scratch_space

  implicit none

  contains

  !> @details Initialises everything needed to run the JediBones miniapp
  !> @param[in,out] modeldb  The structure that holds model state
  !> @param[in]     mesh     Representation of the mesh the code will run on
  !> @param[in,out] chi      The co-ordinate field
  !> @param[in,out] panel_id 2d field giving the id for cubed sphere panels

  subroutine roll_cast(Season1, Season2)

    implicit none

    type(config_type), intent(in) :: Season1
    type(config_type), intent(in) :: Season2


    print*, '========================================================'
    print*, '== SEASON 1 :' // trim(Season1%name())
    print*, '========================================================'

    if (Season1%namelist_exists('master')) then
      print*, ''
      print*, '--=: Jedi Master :=--'
      print*, 'Name: '//trim(Season1%master%name())
      write(log_scratch_space,'(I0)') Season1%master%age()
      print*, 'Age: '// trim(log_scratch_space)
      write(log_scratch_space,'(F4.2)') Season1%master%height()
      print*, 'Height: '//trim(log_scratch_space)
      write(log_scratch_space,'(I0)') Season1%master%street_cred()
      print*, 'Reputation: '//trim(log_scratch_space)
      write(log_scratch_space,'(L1)') Season1%master%ownseries()
      print*, 'Own Series: '//trim(log_scratch_space)
    end if

    if (Season1%namelist_exists('sith')) then
      print*, ''
      print*, '--=: Sith Lord :=--'
      print*, 'Name: '//trim(Season1%sith%name())
      write(log_scratch_space,'(I0)') Season1%sith%age()
      print*, 'Age: '// trim(log_scratch_space)
      write(log_scratch_space,'(F4.2)') Season1%sith%height()
      print*, 'Height: ', trim(log_scratch_space)
      write(log_scratch_space,'(I0)') Season1%sith%street_cred()
      print*, 'Reputation: ', trim(log_scratch_space)
      write(log_scratch_space,'(L1)') Season1%sith%ownseries()
      print*, 'Own Series: ', trim(log_scratch_space)
    end if

    if (Season1%namelist_exists('donkey')) then
      print*, ''
      print*, '--=: Transport Vechicle :=--'
      print*, 'Name: '//trim(Season1%donkey%name())
      write(log_scratch_space,'(I0)') Season1%donkey%age()
      print*, 'Age: '//trim(log_scratch_space)
      print*, 'Fuel: '//trim(Season1%donkey%fuel())
      write(log_scratch_space,'(F7.3)') Season1%donkey%weight()
      print*, 'Weight: '// trim(adjustl(log_scratch_space))
      write(log_scratch_space,'(F7.3)') Season1%donkey%capacity()
      print*, 'Capacity: ', trim(adjustl(log_scratch_space))
    end if

    if (Season1%namelist_exists('pleb')) then
      print*, ''
      print*, '--=: The Sidkick :=--'
      print*, 'Name: '//trim(Season1%pleb%name())
      write(log_scratch_space,'(I0)') Season1%pleb%age()
      print*, 'Age: '//trim(log_scratch_space)
      print*, 'Race: '//trim(Season1%pleb%race())
      print*, 'Expendable: ', Season1%pleb%for_fodder()
    end if

    if (Season1%namelist_exists('weapon')) then
      print*, ''
      print*, '--=: The Weapon :=--'
      print*, 'Name: '//trim(Season1%weapon%name())
      print*, 'One-Handed: ', Season1%weapon%onehanded()
      print*, 'Defensive: ', Season1%weapon%defensive()
      write(log_scratch_space,'(I0)') Season1%weapon%strike_damage()
      print*, 'Strike damage: '// trim(log_scratch_space)
      print*, 'Special move: ', Season1%weapon%special_move()
    end if

    print*, ''
    print*, '========================================================'
    print*, '== SEASON 2 :' // trim(Season2%name())
    print*, '========================================================'

    if (Season2%namelist_exists('master')) then
      print*, ''
      print*, '--=: Jedi Master :=--'
      print*, 'Name: '//trim(Season2%master%name())
      write(log_scratch_space,'(I0)') Season2%master%age()
      print*, 'Age: '// trim(log_scratch_space)
      write(log_scratch_space,'(F4.2)') Season2%master%height()
      print*, 'Height: ', trim(log_scratch_space)
      write(log_scratch_space,'(I0)') Season2%master%street_cred()
      print*, 'Reputation: ', trim(log_scratch_space)
      write(log_scratch_space,'(L1)') Season2%master%ownseries()
      print*, 'Own Series: ', trim(log_scratch_space)
    end if

    if (Season2%namelist_exists('sith')) then
      print*, ''
      print*, '--=: Sith Lord :=--'
      print*, 'Name: '//trim(Season2%sith%name())
      write(log_scratch_space,'(I0)') Season2%sith%age()
      print*, 'Age: '// trim(log_scratch_space)
      write(log_scratch_space,'(F4.2)') Season2%sith%height()
      print*, 'Height: ', trim(log_scratch_space)
      write(log_scratch_space,'(I0)') Season2%sith%street_cred()
      print*, 'Reputation: ', trim(log_scratch_space)
      write(log_scratch_space,'(L1)') Season2%sith%ownseries()
      print*, 'Own Series: ', trim(adjustl(log_scratch_space))
    end if

    if (Season2%namelist_exists('donkey')) then
      print*, ''
      print*, '--=: Transport Vechicle :=--'
      print*, 'Name: '//trim(Season2%donkey%name())
      write(log_scratch_space,'(I0)') Season2%donkey%age()
      print*, 'Age: '// trim(log_scratch_space)
      print*, 'Fuel: '//trim(Season2%donkey%fuel())
      write(log_scratch_space,'(F7.3)') Season2%donkey%weight()
      print*, 'Weight: '// trim(log_scratch_space)
      write(log_scratch_space,'(F7.3)') Season2%donkey%capacity()
      print*, 'Capacity: ', trim(adjustl(log_scratch_space))
    end if

    if (Season2%namelist_exists('pleb')) then
      print*, ''
      print*, '--=: The Sidkick :=--'
      print*, 'Name: '//trim(Season2%pleb%name())
      write(log_scratch_space,'(I0)') Season2%pleb%age()
      print*, 'Age: ', trim(log_scratch_space)
      print*, 'Race: '//trim(Season2%pleb%race())
      print*, 'Expendable: ', Season2%pleb%for_fodder()
    end if

    if (Season2%namelist_exists('weapon')) then
      print*, ''
      print*, '--= The Weapon =--'
      print*, 'Name: '//trim(Season2%weapon%name())
      print*, 'One-Handed: ', Season2%weapon%onehanded()
      print*, 'Defensive: ', Season2%weapon%defensive()
      write(log_scratch_space,'(I0)') Season2%weapon%strike_damage()
      print*, 'Strike damage: '// trim(log_scratch_space)
      print*, 'Special move: ', Season2%weapon%special_move()
    end if

    print*, ''
    print*, ''
    print*, ''
    print*, '====================================================='
    print*, '======================================='
    print*, '========================'
    print*, '============'
    print*, '====='

  end subroutine roll_cast

end module season_cast_mod
