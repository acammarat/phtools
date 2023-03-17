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

subroutine read_eigen
  use var, only: qmatndref, freqndref, qmatndcomp, freqndcomp, &
                 npunique, nq, eig_ref, eig_comp, freq_ref, freq_comp, &
                 atmap
  use functions, only: i2a

  implicit none
  integer :: i, j, k, npG, npH, npS, nqcomp
  integer :: i1
  real(8) :: x, y, z 
  real(8), allocatable :: tmp_e(:), tmp1_e(:)
  complex(8), allocatable :: tmp_eig(:)
  character(1) :: comment
 
  write(*,'(a)') ' Reading eigenvectors and frequencies:'

  open(unit=20,file=freqndref,action='READ')
  read(20,*)
  read(20,*) comment, comment , npG , comment, npH, comment, npS, comment, nq

  open(unit=21,file=freqndcomp,action='READ')
  read(21,*)
  read(21,*) comment, comment, i, comment, j, comment, k, comment, nqcomp

  if ( (i /= npG) .or. (j /=npH) .or. (k /= npS) .or. (nqcomp /= nq) ) then
    write(*,'(a)') ' ERROR: the q-sets or the number of modes are not consistent among the reference and the files to compare'
    write(*,*)
    stop
  end if

  ! number of "unique" q-points (i.e. Gamma + set H + set S) 
  npunique = npG + npH + npS
  write(*,'(*(a))') '  q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)
  write(*,'(*(a))') '  Number of q-points: ', i2a(npunique)

  allocate ( eig_ref(npunique,nq,nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for eig_ref'
  allocate ( eig_comp(npunique,nq,nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for eig_comp'
  allocate ( freq_ref(npunique,nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for freq_ref'
  allocate ( freq_comp(npunique,nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for freq_comp'
  allocate ( tmp_e(2*nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp_e'
  allocate ( tmp1_e(2*nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp1_e'
  allocate ( tmp_eig(nq), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp_eig'

  read(20,*)
  read(20,*)
  read(21,*)
  read(21,*)

  ! read eigenvectors and frequencies
  open(unit=10,file=qmatndref,action='READ')
  read(10,*)
  open(unit=11,file=qmatndcomp,action='READ')
  read(11,*)

  write(*,'(a)') '  Reading eigenvectors and frequencies at'
  do k = 1, npunique 
    read(20,*) ! we don't need the q-vector
    read(20,*) comment, i, x, y, z ! but we print the reduced components for check
    write(*,'(a,a,a,3(f10.5,1x),a,a)') '   ref  q',i2a(k), ' ', x, y, z
    read(21,*) ! we don't need the q-vector
    read(21,*) comment, i, x, y, z ! but we print the reduced components for check
    write(*,'(a,a,a,3(f10.5,1x),a,a)') '   comp q',i2a(k), ' ', x, y, z
    do i = 1, nq
      read(10,*) tmp_e(:)
      read(11,*) tmp1_e(:)
      read(20,*) freq_ref(k,i)
      read(21,*) freq_comp(k,i)
      i1 = 0
      do j = 1, 2*nq-1, 2
        i1 = i1 + 1
        eig_ref(k,i,i1) = dcmplx( tmp_e(j), tmp_e(j+1) ) ! adim
        eig_comp(k,i,i1) = dcmplx( tmp1_e(j), tmp1_e(j+1) ) ! adim
      end do
    end do

    ! rearranging the components of reference eigenvectors according to the mapping
    do i = 1, nq
      tmp_eig(:) = eig_comp(k,i,:)
      do j = 1, nq/3 
        eig_comp(k,i, (j-1)*3+1:(j-1)*3+3 ) = tmp_eig( (atmap(j)-1)*3+1:(atmap(j)-1)*3+3 )
      end do
    end do

  end do

  close(20)
  close(21)
  close(10)
  close(11)

  deallocate ( tmp_e, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp_e'
  deallocate ( tmp1_e, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp1_e'
  deallocate ( tmp_eig, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp_eig'

  write(*,'(a)') ' done.'

  return
end subroutine read_eigen
