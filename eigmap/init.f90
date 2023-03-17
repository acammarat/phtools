!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! eigmap. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Program to calculate the map between different eigenvector
! via the scalar product
!
! If used for production, you should cite
! Phys. Rev. B XX, XXXXX (XXXX)
! https://doi.org/10.1103/xxx
! where the formulation is reported in section II "Eigenvector map based on
! atomic displacements" of the Supplemental Material.
!
!    This file is part of eigmap.
!
!    eigmap is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    eigmap is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with phonchar.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init
  use String_mod, only: String_type
  use functions, only: i2a
  use var, only: qmatndref, freqndref, qmatndcomp, freqndcomp, &
                 atmap, nq, atoms_ref

  implicit none
  integer :: i, j, k, natmap
  character(256) :: infile
  type(String_type) :: filename
  logical :: file_exists, fl_chkunique
  character(1) :: comment

  call show_logo

  call get_command_argument(1,infile)
  if ( (infile == '-h') .or. (infile == '') ) then
    write(*,'(a)') ' Syntax: eigmap <setting file>'
    write(*,*)
    stop
  end if

  inquire(file=infile,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: setting file ',trim(infile),' not found.'
    write(*,*)
    stop
  end if

  open(unit=10,file=infile,action='READ')
  write(*,'(*(a))') ' Reading settings from file: ', trim(infile)

  read(10,'(a)') qmatndref
  filename%value = trim(qmatndref)
  filename%Parts = filename%split(filename%value, delim = " ")
  qmatndref = trim(filename%Parts(1)%record)
  inquire(file=qmatndref,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(qmatndref),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  reference eigenvector file: ', trim(qmatndref)

  read(10,'(a)') freqndref
  filename%value = trim(freqndref)
  filename%Parts = filename%split(filename%value, delim = " ")
  freqndref = trim(filename%Parts(1)%record)
  inquire(file=freqndref,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(freqndref),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  reference frequency file: ', trim(freqndref)

  read(10,'(a)') qmatndcomp
  filename%value = trim(qmatndcomp)
  filename%Parts = filename%split(filename%value, delim = " ")
  qmatndcomp = trim(filename%Parts(1)%record)
  inquire(file=qmatndcomp,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(qmatndcomp),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  eigenvector file to compare: ', trim(qmatndcomp)

  read(10,'(a)') freqndcomp
  filename%value = trim(freqndcomp)
  filename%Parts = filename%split(filename%value, delim = " ")
  freqndcomp = trim(filename%Parts(1)%record)
  inquire(file=freqndcomp,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(freqndcomp),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  frequency file to compare: ', trim(freqndcomp)

  open(unit=20,file=freqndref,action='READ')
  read(20,*)
  read(20,*) comment, comment , i , comment, i, comment, i, comment, nq
  close(20)

  atoms_ref = nq/3
  allocate ( atmap(atoms_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for atmap'

  ! atom map is initialized as the identity
  do i = 1, atoms_ref
    atmap(i) = i
  end do

  read(10,*,iostat=i) natmap
  if ( natmap > 0 ) then
    write(*,'(a)') '  atom mapping:'
    do i = 1, natmap
      read(10,*) j, k
      write(*,'(*(a))') '    ',i2a(j),' -> ',i2a(k)
      atmap(j) = k
    end do
  end if

  !check that the mapping is unique
  fl_chkunique = .false.
  do i = 1, atoms_ref-1
    do j = i+1, atoms_ref
      if ( atmap(i) == atmap(j) ) then
        write(*,'(*(a))') '  non unique map: ',i2a(i), ' -> ',i2a(atmap(i)),' ; ',i2a(j), ' -> ',i2a(atmap(j))
        fl_chkunique = .true.
      end if
    end do
  end do

  if ( fl_chkunique ) then
    write(*,'(a)') ' ERROR: map is not unique.'
    write(*,*)
    stop
  end if

  close(10)

  write(*,'(a)') ' done.'

  return
end subroutine init
