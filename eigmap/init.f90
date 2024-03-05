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
  use var, only: qmatnd_ref, freqnd_ref, qmatnd_comp, freqnd_comp, atmap, &
                 nq_ref, atoms_UC_ref, natoms_UC_ref, &
                 pos_eq_UC_ref, qgm_ref, mass_UC_ref, &
                 side_eq_UC_ref, &
                 nq_comp, atoms_UC_comp, natoms_UC_comp, &
                 pos_eq_UC_comp, qgm_comp, mass_UC_comp, &
                 side_eq_UC_comp, &
                 fu_mul

  implicit none
  integer :: i, j, k, natmap, natom_types_ref, natom_types_comp
  real(8) :: r_tmp, x_tmp, y_tmp, z_tmp
  real(8), allocatable :: mass_pertype_ref(:), mass_pertype_comp(:)
  character(2), allocatable :: at_tmp(:), at_pertype_ref(:), at_pertype_comp(:)
  character(256) :: infile, word
  character(256) :: pos_ref, pos_comp
  logical :: file_exists, fl_chkunique

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

  ! reference configuration
  read(10,*) pos_ref
  inquire(file=pos_ref,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(pos_ref),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  POSCAR reference file: ', trim(pos_ref)

  read(10,*) natom_types_ref
  write(*,'(*(a))') '  Number of atom types: ',i2a(natom_types_ref)

  allocate ( at_pertype_ref(natom_types_ref), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for at_pertype'
  allocate ( mass_pertype_ref(natom_types_ref), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for mass_pertype'
  do i = 1, natom_types_ref
     read(10,*) at_pertype_ref(i), mass_pertype_ref(i) ! atom symbol, uma
  end do

  open(unit=15,file=pos_ref,action='READ') 
  read(15,*)
  read(15,*) r_tmp
  if ( abs(r_tmp - 1.d0) > tiny(1.d0) ) then
    write(*,'(a)') '  ERROR: the scale factor in the reference POSCAR file must be 1.0.'
    write(*,*)
    stop
  end if
  do i = 1, 3
    read(15,*) side_eq_UC_ref(i,:) ! Ang
  end do
  
  allocate ( at_tmp(natom_types_ref), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for at_tmp'
  read(15,*) at_tmp(:) ! atom symbols
  do i = 1, natom_types_ref
    if ( at_pertype_ref(i) /= at_tmp(i) ) then
      write(*,'(a)') '  ERROR: atomic types in the POSCAR file do not match'
      write(*,'(a)') '         those in the input file.'
      write(*,*)
      stop
    end if
  end do

  deallocate ( at_tmp, stat = i )
  if ( i/=0 ) stop 'Deallocation failed for at_tmp'

  allocate ( natoms_UC_ref(natom_types_ref), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for natoms_uc_ref'
  read(15,*) natoms_UC_ref(:) ! number of atoms for each type

  atoms_UC_ref = 0
  do i = 1, natom_types_ref
     atoms_UC_ref = atoms_UC_ref + natoms_UC_ref(i) ! total number of atoms
  end do
  nq_ref = 3 * atoms_UC_ref ! number of degrees of freedom

  allocate ( mass_UC_ref(atoms_UC_ref), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for mass_uc_ref'

  k = 0
  do i = 1, natom_types_ref
     do j = 1, natoms_UC_ref(i)
        k = k + 1
        mass_UC_ref(k) = mass_pertype_ref(i) ! uma
     end do
  end do
  
  allocate ( pos_eq_UC_ref(atoms_UC_ref,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for pos_eq_uc_ref'

  read(15,*) word
  do i = 1, atoms_UC_ref
     read(15,*) pos_eq_UC_ref(i,:) ! Ang or adim
     if ( word == 'Direct' ) then
        x_tmp = side_eq_UC_ref(1,1)*pos_eq_UC_ref(i,1) + side_eq_UC_ref(2,1)*pos_eq_UC_ref(i,2) + side_eq_UC_ref(3,1)*pos_eq_UC_ref(i,3)
        y_tmp = side_eq_UC_ref(1,2)*pos_eq_UC_ref(i,1) + side_eq_UC_ref(2,2)*pos_eq_UC_ref(i,2) + side_eq_UC_ref(3,2)*pos_eq_UC_ref(i,3)
        z_tmp = side_eq_UC_ref(1,3)*pos_eq_UC_ref(i,1) + side_eq_UC_ref(2,3)*pos_eq_UC_ref(i,2) + side_eq_UC_ref(3,3)*pos_eq_UC_ref(i,3)
        pos_eq_UC_ref(i,1) = x_tmp
        pos_eq_UC_ref(i,2) = y_tmp
        pos_eq_UC_ref(i,3) = z_tmp
     end if
  end do
  close(15)

  read(10,*) qmatnd_ref
  inquire(file=qmatnd_ref,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(qmatnd_ref),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  reference eigenvector file: ', trim(qmatnd_ref)

  read(10,*) freqnd_ref
  inquire(file=freqnd_ref,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(freqnd_ref),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  reference frequency file: ', trim(freqnd_ref)

  read(10,*) qgm_ref(:)

  ! comparison configuration
  read(10,*) pos_comp
  inquire(file=pos_comp,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(pos_comp),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  POSCAR comparison file: ', trim(pos_comp)

  read(10,*) natom_types_comp
  write(*,'(*(a))') '  Number of atom types: ',i2a(natom_types_comp)

  allocate ( at_pertype_comp(natom_types_comp), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for at_pertype_comp'
  allocate ( mass_pertype_comp(natom_types_comp), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for mass_pertype_comp'
  do i = 1, natom_types_comp
     read(10,*) at_pertype_comp(i), mass_pertype_comp(i) ! atom symbol, uma
  end do

  open(unit=15,file=pos_comp,action='READ') 
  read(15,*)
  read(15,*) r_tmp
  if ( abs(r_tmp - 1.d0) > tiny(1.d0) ) then
    write(*,'(a)') '  ERROR: the scale factor in the comparison POSCAR file must be 1.0.'
    write(*,*)
    stop
  end if
  do i = 1, 3
    read(15,*) side_eq_UC_comp(i,:) ! Ang
  end do
  
  allocate ( at_tmp(natom_types_comp), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for at_tmp'
  read(15,*) at_tmp(:) ! atom symbols
  do i = 1, natom_types_comp
    if ( at_pertype_comp(i) /= at_tmp(i) ) then
      write(*,'(a)') '  ERROR: atomic types in the POSCAR file do not match'
      write(*,'(a)') '         those in the input file.'
      write(*,*)
      stop
    end if
  end do

  deallocate ( at_tmp, stat = i )
  if ( i/=0 ) stop 'Deallocation failed for at_tmp'

  allocate ( natoms_UC_comp(natom_types_comp), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for natoms_uc_comp'
  read(15,*) natoms_UC_comp(:) ! number of atoms for each type

  atoms_UC_comp = 0
  do i = 1, natom_types_comp
     atoms_UC_comp = atoms_UC_comp + natoms_UC_comp(i) ! total number of atoms
  end do
  nq_comp = 3 * atoms_UC_comp ! number of degrees of freedom

  allocate ( mass_UC_comp(atoms_UC_comp), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for mass_uc_comp'

  k = 0
  do i = 1, natom_types_comp
     do j = 1, natoms_UC_comp(i)
        k = k + 1
        mass_UC_comp(k) = mass_pertype_comp(i) ! uma
     end do
  end do
  
  allocate ( pos_eq_UC_comp(atoms_UC_comp,3), stat = i )
  if ( i /= 0 ) stop 'Allocation failed for pos_eq_uc_comp'

  read(15,*) word
  do i = 1, atoms_UC_comp
     read(15,*) pos_eq_UC_comp(i,:) ! Ang or adim
     if ( word == 'Direct' ) then
        x_tmp = side_eq_UC_comp(1,1)*pos_eq_UC_comp(i,1) + side_eq_UC_comp(2,1)*pos_eq_UC_comp(i,2) + side_eq_UC_comp(3,1)*pos_eq_UC_comp(i,3)
        y_tmp = side_eq_UC_comp(1,2)*pos_eq_UC_comp(i,1) + side_eq_UC_comp(2,2)*pos_eq_UC_comp(i,2) + side_eq_UC_comp(3,2)*pos_eq_UC_comp(i,3)
        z_tmp = side_eq_UC_comp(1,3)*pos_eq_UC_comp(i,1) + side_eq_UC_comp(2,3)*pos_eq_UC_comp(i,2) + side_eq_UC_comp(3,3)*pos_eq_UC_comp(i,3)
        pos_eq_UC_comp(i,1) = x_tmp
        pos_eq_UC_comp(i,2) = y_tmp
        pos_eq_UC_comp(i,3) = z_tmp
     end if
  end do
  close(15)

  if ( atoms_UC_comp > atoms_UC_ref ) then
    if ( mod(atoms_UC_comp, atoms_UC_ref) /= 0 ) then
      write(*,'(a)') ' ERROR: formula unit in comparison structure is not an integer multiple of the reference one.'
      write(*,*)
      stop
    end if
    fu_mul = atoms_UC_comp / atoms_UC_ref
  else if ( atoms_UC_comp < atoms_UC_ref ) then
    if ( mod(atoms_UC_ref, atoms_UC_comp) /= 0 ) then
      write(*,'(a)') ' ERROR: formula unit in reference structure is not an integer multiple of the comparison one.'
      write(*,*)
      stop
    end if
    fu_mul = atoms_UC_ref / atoms_UC_comp
  else
    fu_mul = 1
  end if

  read(10,*) qmatnd_comp
  inquire(file=qmatnd_comp,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(qmatnd_comp),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  comparison eigenvector file: ', trim(qmatnd_comp)

  read(10,*) freqnd_comp
  inquire(file=freqnd_comp,exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(*(a))') ' ERROR: input file ',trim(freqnd_comp),' not found.'
    write(*,*)
    stop
  end if
  write(*,'(*(a))') '  comparison frequency file: ', trim(freqnd_comp)

  read(10,*) qgm_comp(:)

  allocate ( atmap(atoms_UC_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for atmap'

  ! atom map is initialized as the identity
  do i = 1, atoms_UC_ref
    atmap(i) = i
  end do

  read(10,*,iostat=i) natmap
  if ( natmap /= 0 .and. fu_mul /= 1 ) then
    write(*,'(a)') ' ERROR: atom mapping allowed only if the two structures have the same number of atoms'
    write(*,*)
    stop
  end if
  if ( natmap > 0 ) then
    write(*,'(a)') '  atom mapping:'
    do i = 1, natmap
      read(10,*) j, k
      write(*,'(*(a))') '    ',i2a(j),' -> ',i2a(k)
      atmap(j) = k
    end do
  else if ( natmap < 0 ) then
    write(*,'(a)') ' ERROR: the number of atoms to map cannot be negative.'
    write(*,*)
    stop
  end if

  !check that the mapping is unique
  fl_chkunique = .false.
  do i = 1, atoms_UC_ref-1
    do j = i+1, atoms_UC_ref
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
  write(*,*)

  return
end subroutine init
