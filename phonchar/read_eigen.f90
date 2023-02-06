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

subroutine read_eigen
  use var, only: nat, nqp, neig, eig, freq
  use integer_to_string
                 
  implicit none
  integer(4) :: i, j, k, ix
  integer(8) :: l
  character(200) :: dum, in2file

  call get_command_argument(2,in2file)
  open(unit=10,file=in2file, action='READ')

  write(*,'(2a)') ' Reading eigenvectors and frequencies from file: ', trim(in2file)
  write(*,'(a)', advance='no') ' Checking if input file contains eigenvectors: '
  do
    read(10,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, no eigenvectors found.'
      write(*,*)
      stop 
    else if ( dum == 'eigenvector:' ) then
      write(*,'(a)') 'yes'
      exit
    end if
  end do
  call fseek(10, 0, 0, i)
  l=ftell(10)
  do
    read(10,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, nqpoints string not found.'
      write(*,*)
      stop 
    else if ( dum == 'nqpoint:' ) then
      backspace(10)
      read(10,*) dum, nqp
      exit
    end if
  end do

  write(*,'(2a)') ' Number of q-points: ', i2a(nqp)
  call fseek(10, 0, 0, i)
  l=ftell(10)
  do
    read(10,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(a)') 'reached end of file, natom string not found.'
      write(*,*)
      stop 
    else if ( dum == 'natom:' ) then
      backspace(10)
      read(10,*) dum, nat
      write(*,'(2a)') ' Number of atoms: ', i2a(nat)
      exit
    end if
  end do

  neig=nat*3

  allocate ( eig(nqp,neig,neig), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for eig'
  allocate ( freq(nqp,neig), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for freq'

  do
    read(10,*,iostat=i) dum
    if ( i<0 ) then
      write(*,'(*(a))') 'reached end of file, the input file ', trim(in2file), ' is not complete.'
      write(*,*)
      stop 
    else if ( dum == 'phonon:' ) then
      exit
    end if
  end do
  
  do i = 1, nqp
    do k = 1,3
      read(10,*)
    end do

    do j = 1, neig
      read(10,*)
      read(10,*) dum, freq(i,j)
      read(10,*)

      do k = 1, nat
        read(10,*)
        do ix = 1, 3
          read(10,*) dum, dum, eig(i,j,(k-1)*3+ix)
        end do
      end do

    end do
    read(10,*)
  end do

  close(10)

  write(*,*) 'Reading input file done.'

  return
end subroutine read_eigen
