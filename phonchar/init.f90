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

subroutine init
  use var, only: nag, ag, natmax, natg
  implicit none
  integer :: i, j
  character(256) :: infile, in2file
  logical :: file_exists

  call show_logo

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '' ) ) then
     write(*,'(a)') ' Syntax: phonchar <setting file> <band.yaml file>'
     write(*,*)
     stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  call get_command_argument(2,in2file)
  if ( (in2file == '-h') .or. (in2file == '' ) ) then
     write(*,'(a)') ' Syntax: phonchar <setting file> <band.yaml file>'
     write(*,*)
     stop
  end if
  inquire(file=infile,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input yaml file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  open(unit=20,file=infile,action='READ')

  write(*,'(2a)') ' Reading settings from file: ', trim(infile)
  read(20,*) nag  ! number of atomic groups
  if ( nag < 2 ) then
     write(*,'(a59)') ' ERROR: the number of atomic groups must be greater than 1.'
     write(*,*)
     stop
  end if
  allocate ( ag(nag,natmax), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for ag'

  do i = 1, nag
    read(20,*) natg(i) ! number of atoms in i-th group
    if ( natg(i) > natmax ) then
       write(*,'(a)') ' ERROR: the number of atoms in the groups exceeds natmax.'
       write(*,*)
       stop
    end if

    do j = 1, natg(i)
      read(20,*) ag(i,j)    ! atom j in group i
    end do
  end do

  close(20)

  return
end subroutine init
