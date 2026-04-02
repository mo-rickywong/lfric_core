!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> JediBones knows what configuration it needs.
!>
module JediBones_mod

  implicit none

  private

  character(*), public, parameter ::                   &
      JediBones_S1_namelists(8) = [ 'base_mesh     ', &
                                     'extrusion     ', &
                                     'finite_element', &
                                     'partitioning  ', &
                                     'planet        ', &
                                     'master        ', &
                                     'sith          ', &
                                     'weapon        ' ]

  character(*), public, parameter ::                   &
      JediBones_S2_namelists(8) = [ 'base_mesh     ', &
                                     'extrusion     ', &
                                     'finite_element', &
                                     'partitioning  ', &
                                     'planet        ', &
                                     'sith          ', &
                                     'donkey        ', &
                                     'weapon        ' ]


end module JediBones_mod
