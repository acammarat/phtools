!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! phonchar. Copyright (C) 2022 Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Program to calculate the atomic character
! of phonon eigenvectors obtained from PHONOPY
! ( https://phonopy.github.io/phonopy )
!
! If used for production, you should cite
! Phys. Rev. B 103, 035406 (2021)
! https://doi.org/10.1103/PhysRevB.103.035406
! where the formulation is reported in section V "Atomic character of the phonon modes" 
! of the supplemental material.
!
! Only the real part of the eigenvector is considered, as the purpose
! is to analyse the contribution to the atomic motions.
! 
!
!    This file is part of phonchar.
!
!    phonchar is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    phonchar is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with phonchar.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! v 2.2
! Removed unnecessary variables
!
! v 2.1
! The line-skip for the search in the band.yaml file is now general
!
! v 2.0
! It is possible to specify up to ngmax groups of atoms
! If the number of atomic groups is 2, the weight is calculated
! as in v 1.0.
!
! v 1.0
! It is possible to specify only two groups of atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Format of the input file
!
! int                     number of atomic groups
! int                     number of atoms in group 1
! int                     label of the first atom in the group
! ...
! int                     label of the last atom in the group
! int                     number of atoms in group 2
! int                     label of the first atom in the group
! ...
! int                     label of the last atom in the group
! ...                     number of atoms in group 3 - if any, and similar
!                         input structure as above
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module integer_to_string
contains
  function i2a(i) result(out)
    character(:), allocatable :: out
    integer(4), intent(in) :: i
    character(range(i)+2) :: x

    write(x,'(i0)') i

    out = trim(x)
    
  end function i2a
end module integer_to_string

module var
  ! global parameters
  integer, parameter :: natmax = 1000
  character(3), parameter :: version = '2.2'
  character(8), parameter :: progname = 'PHONCHAR'
  ! from input files
  integer, save :: nat, nqp, neig, nag
  integer, save :: natg(natmax)
  integer, save, allocatable :: ag(:,:)
  real(8), save, allocatable :: freq(:,:), eig(:,:,:)
end module var

program phonchar

  call init

  call read_eigen
  call calc_char

  call deallocate_all
  call credits

  stop
end program phonchar
