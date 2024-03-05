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
! https://xxx
! where the formulation is reported
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
  use var, only: qmatnd_ref, freqnd_ref, qmatnd_comp, freqnd_comp, &
                 npunique_ref, nq_ref, eig_ref, freq_ref, &
                 vec_ref, vec_red_ref, ncells_vec_ref, &
                 npunique_comp, nq_comp, eig_comp, freq_comp, &
                 vec_comp, vec_red_comp, ncells_vec_comp, &
                 atmap, fu_mul
  use functions, only: i2a, lcm

  implicit none
  integer :: i, j, k, kk, npG, npH, npS, i1, mcm, sign_vec
  integer :: p(3)
  integer, allocatable :: q(:,:)
  real(8), allocatable :: tmp_e(:)
  complex(8), allocatable :: tmp_eig(:)
  character(1) :: comment
 
  write(*,'(a)') ' Reading eigenvectors and frequencies:'

  write(*,*)
  write(*,'(a)') '  reference structure'
  ! reference eigenvector and frequencies
  open(unit=20,file=freqnd_ref,action='READ')
  read(20,*)
  read(20,*) comment, comment , npG , comment, npH, comment, npS, comment, i
  if ( i /= nq_ref ) then 
    write(*,'(a)') ' ERROR: the number of modes in ',trim(freqnd_ref),' is not consistent'
    write(*,'(a)') '        with the number of modes calculated from the geometry file.'
    write(*,*)
    stop
  end if

  ! number of "unique" q-points (i.e. Gamma + set H + set S) 
  npunique_ref = npG + npH + npS
  write(*,'(*(a))') '  q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)
  write(*,'(*(a))') '  Number of q-points: ', i2a(npunique_ref)

  allocate ( eig_ref(npunique_ref,nq_ref,nq_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for eig_ref'
  allocate ( freq_ref(npunique_ref,nq_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for freq_ref'
  allocate ( vec_ref(npunique_ref,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for vec_ref'
  allocate ( vec_red_ref(npunique_ref,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for vec_red_ref'
  allocate ( q(npunique_ref,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for q'
  allocate ( tmp_e(2*nq_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp_e'
  allocate ( tmp_eig(nq_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp_eig'

  read(20,*)
  read(20,*)

  ! read eigenvectors and frequencies
  open(unit=10,file=qmatnd_ref,action='READ')
  read(10,*)

  write(*,'(a)') '  Reading eigenvectors and frequencies at'
  do k = 1, npunique_ref
    read(20,*) comment, vec_ref(k,:)
    read(20,*) comment, comment, vec_red_ref(k,:)
    write(*,'(a,a,a,3(1x,f10.5))') '   ref  q',i2a(k), ' ', vec_red_ref(k,:)
    do i = 1, nq_ref
      read(10,*) tmp_e(:)
      read(20,*) freq_ref(k,i)
      i1 = 0
      do j = 1, 2*nq_ref-1, 2
        i1 = i1 + 1
        eig_ref(k,i,i1) = dcmplx( tmp_e(j), tmp_e(j+1) ) ! adim
      end do
    end do

    ! extract integer denominators from reduced vector components
    do i = 1, 3
      if ( abs(vec_red_ref(k,i)) > tiny(1.d0) ) then
         sign_vec = 1
         if ( vec_red_ref(k,i) < 0.d0 ) sign_vec = -1
         call real_to_rational ( dble(sign_vec)*vec_red_ref(k,i), p(i), q(k,i) )
         p(i) = sign_vec*p(i)
      else
         p(i) = 0
         q(k,i) = 1
      end if
    end do

  end do

  close(20)
  close(10)

  deallocate ( tmp_e, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp_e'
  deallocate ( tmp_eig, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp_eig'

  ! calculate commensurate supercell
  do i = 1, 3    
     if (npunique_ref == 1) then
        mcm = q(1,i)
     else        
        mcm = 1
        do k = 1, npunique_ref-1
           do kk = k+1, npunique_ref
              mcm = lcm(mcm, lcm(q(k,i), q(kk,i)))
           end do
        end do
     end if         
     ncells_vec_ref(i) = abs(mcm)
  end do
  write(*,'(*(a))') '  commensurate supercell: ', i2a(ncells_vec_ref(1)),' x ', i2a(ncells_vec_ref(2)),' x ', i2a(ncells_vec_ref(3))

  deallocate ( q, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for q'


  write(*,*)
  write(*,'(a)') '  comparison structure'
  ! comparison eigenvector and frequencies
  open(unit=20,file=freqnd_comp,action='READ')
  read(20,*)
  read(20,*) comment, comment , npG , comment, npH, comment, npS, comment, i
  if ( i /= nq_comp ) then 
    write(*,'(a)') ' ERROR: the number of modes in ',trim(freqnd_comp),' is not consistent'
    write(*,'(a)') '        with the number of modes calculated from the geometry file.'
    write(*,*)
    stop
  end if

  ! number of "unique" q-points (i.e. Gamma + set H + set S) 
  npunique_comp = npG + npH + npS
  write(*,'(*(a))') '  q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)
  write(*,'(*(a))') '  Number of q-points: ', i2a(npunique_comp)

  allocate ( eig_comp(npunique_comp,nq_comp,nq_comp), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for eig_comp'
  allocate ( freq_comp(npunique_comp,nq_comp), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for freq_comp'
  allocate ( vec_comp(npunique_comp,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for vec_comp'
  allocate ( vec_red_comp(npunique_comp,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for vec_red_comp'
  allocate ( q(npunique_comp,3), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for q'
  allocate ( tmp_e(2*nq_comp), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp_e'
  allocate ( tmp_eig(nq_comp), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for tmp_eig'

  read(20,*)
  read(20,*)

  ! read eigenvectors and frequencies
  open(unit=10,file=qmatnd_comp,action='READ')
  read(10,*)

  write(*,'(a)') '  Reading eigenvectors and frequencies at'
  do k = 1, npunique_comp
    read(20,*) comment, vec_comp(k,:)
    read(20,*) comment, comment, vec_red_comp(k,:)
    write(*,'(a,a,a,3(1x,f10.5))') '   ref  q',i2a(k), ' ', vec_red_comp(k,:)
    do i = 1, nq_comp
      read(10,*) tmp_e(:)
      read(20,*) freq_comp(k,i)
      i1 = 0
      do j = 1, 2*nq_comp-1, 2
        i1 = i1 + 1
        eig_comp(k,i,i1) = dcmplx( tmp_e(j), tmp_e(j+1) ) ! adim
      end do
    end do

    ! extract integer denominators from reduced vector components
    do i = 1, 3
      if ( abs(vec_red_comp(k,i)) > tiny(1.d0) ) then
         sign_vec = 1
         if ( vec_red_comp(k,i) < 0.d0 ) sign_vec = -1
         call real_to_rational ( dble(sign_vec)*vec_red_comp(k,i), p(i), q(k,i) )
         p(i) = sign_vec*p(i)
      else
         p(i) = 0
         q(k,i) = 1
      end if
    end do

    ! rearranging the components of comparison eigenvectors according to the mapping
    if ( fu_mul == 1 ) then
      do i = 1, nq_comp
        tmp_eig(:) = eig_comp(k,i,:)
        do j = 1, nq_comp/3 
          eig_comp(k,i, (j-1)*3+1:(j-1)*3+3 ) = tmp_eig( (atmap(j)-1)*3+1:(atmap(j)-1)*3+3 )
        end do
      end do
    end if

  end do

  close(20)
  close(10)

  deallocate ( tmp_e, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp_e'
  deallocate ( tmp_eig, stat = i )
  if ( i /= 0 ) stop ' ERROR: Deallocation failed for tmp_eig'

  ! calculate commensurate supercell
  do i = 1, 3    
     if (npunique_comp == 1) then
        mcm = q(1,i)
     else        
        mcm = 1
        do k = 1, npunique_comp-1
           do kk = k+1, npunique_comp
              mcm = lcm(mcm, lcm(q(k,i), q(kk,i)))
           end do
        end do
     end if         
     ncells_vec_comp(i) = abs(mcm)
  end do
  write(*,'(*(a))') '  commensurate supercell: ', i2a(ncells_vec_comp(1)),' x ', i2a(ncells_vec_comp(2)),' x ', i2a(ncells_vec_comp(3))

  write(*,*)
  write(*,'(a)') ' done.'
  write(*,*)

  return
end subroutine read_eigen

! This subroutine finds the best integer p and q such that p/q is close to a given real number x
subroutine real_to_rational ( x, p, q )
  use functions, only: gcd
  implicit none
  real(8), intent(in) :: x
  integer, intent(out) :: p, q
  integer :: f!, gcd
  real(8) :: r, e, best 
  
  p = 1 
  q = 1         
  best = x * 6.d0    

  do 
     r = dble(p) / dble(q)                
     e = x - r                  
     if ( abs(e) <= best ) then 
        best = abs(e) * 0.125d0             
        f = gcd(p,q)                    
        if ( abs(e) < 0.000001d0 ) exit                
     end if
     if ( e > 0.d0 ) then 
        p = p + ceiling( e * q )    
     else if ( e < 0.d0 ) then    
        q = q + 1                       
     end if
  end do

  return        
end subroutine real_to_rational


!! This subroutine returns the least common multiplier (lcm) of a pair of integers (a,b)
!integer function lcm ( a, b )
!  implicit none
!  integer, intent(in) :: a, b
!  integer :: gcd
!  
!  lcm = a * b / gcd(a,b)
!  
!  return
!end function lcm
!
!
!! This subroutine returns the greatest common divisor (gcd) of a pair of integers (a,b)
!integer function gcd ( a, b )
!  implicit none
!  integer, intent(in) :: a, b
!  integer :: aa, bb, t
!  
!  aa = a
!  bb = b
!  do while ( bb /= 0 )
!     t = bb
!     bb = mod(aa,bb)
!     aa = t
!  end do
!  gcd = abs(aa)
!  
!  return
!end function gcd
