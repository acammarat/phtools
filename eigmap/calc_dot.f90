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

subroutine calc_dot
  use var, only: rad_to_deg, &
                 atoms_UC_ref, npunique_ref, nq_ref, eig_ref, freq_ref, &
                 vec_ref, ncells_vec_ref, side_eq_UC_ref, pos_eq_UC_ref, mass_UC_ref, &
                 atoms_UC_comp, npunique_comp, nq_comp, eig_comp, freq_comp, &
                 vec_comp, ncells_vec_comp, side_eq_UC_comp, pos_eq_UC_comp, mass_UC_comp, &
                 fu_mul!, qgm_ref, qgm_comp
                 
  use functions, only: i2a, lcm!, inv

  implicit none
  integer :: i, k1, j1, k2, j2, k, j, i1, i2, i3, ix1, n, n1, nrep, mcm
  integer :: natoms_tot_ref, natoms_tot_comp, ncells_tot_ref, ncells_tot_comp
  integer :: side_mul_ref(3), side_mul_comp(3), ncells_ref(3), ncells_comp(3)
  real(8) :: xt, yt, zt
  real(8) :: side_ref(3), side_comp(3)!, E_acu_ref(3,3), E_acu_comp(3,3), M(3,3)
  real(8), allocatable :: u_ref(:), u_comp(:), pos_eq_EC_ref(:,:,:), pos_eq_EC_comp(:,:,:)
  real(8) :: alphamin, refmod, compmod, alpha

  
  ! check if one structure is an integer multiple of the other
  do i = 1, 3
    side_ref(i) = sqrt(dot_product(side_eq_UC_ref(i,:), side_eq_UC_ref(i,:)))
    side_comp(i) = sqrt(dot_product(side_eq_UC_comp(i,:), side_eq_UC_comp(i,:)))
  end do

  side_mul_ref(:) = 1
  side_mul_comp(:) = 1
  if ( fu_mul /= 1 ) then
    if ( atoms_UC_ref >= atoms_UC_comp ) then
      side_mul_comp(:) = nint(side_ref(:)/side_comp(:))
      nrep = side_mul_comp(1)*side_mul_comp(2)*side_mul_comp(3)
    else
      side_mul_ref(:) = nint(side_comp(:)/side_ref(:))
      nrep = side_mul_ref(1)*side_mul_ref(2)*side_mul_ref(3)
    end if
    if ( nrep /= fu_mul ) then
      write(*,'(a)') ' ERROR: the cell multiplicity obtained from the ratio between corresponding lattice vectors is not the same one'
      write(*,'(a)') '        obtained from the ratio between the two formula units.'
      write(*,*)
      stop
    end if
  end if
  
  write(*,'(*(a))') ' Ratio between the formula units: ',i2a(fu_mul)

  do i = 1, 3
     mcm = lcm(ncells_vec_ref(i), ncells_vec_comp(i))
     ncells_ref(i) = abs(mcm)
     ncells_comp(i) = abs(mcm)
  end do

  ncells_ref(:) = side_mul_ref(:)*ncells_ref(:)
  ncells_comp(:) = side_mul_comp(:)*ncells_comp(:)
  write(*,'(*(a))') ' Reference supercell: ', i2a(ncells_ref(1)),' x ', i2a(ncells_ref(2)),' x ', i2a(ncells_ref(3))
  write(*,'(*(a))') ' Comparison supercell: ', i2a(ncells_comp(1)),' x ', i2a(ncells_comp(2)),' x ', i2a(ncells_comp(3))

  ncells_tot_ref = ncells_ref(1)*ncells_ref(2)*ncells_ref(3)
  natoms_tot_ref = atoms_UC_ref*ncells_tot_ref
  ncells_tot_comp = ncells_comp(1)*ncells_comp(2)*ncells_comp(3)
  natoms_tot_comp = atoms_UC_comp*ncells_tot_comp

  allocate ( u_ref(3*natoms_tot_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for u_ref'
  allocate ( pos_eq_EC_ref(atoms_UC_ref,3,ncells_tot_ref), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for pos_eq_EC_ref'
  allocate ( u_comp(3*natoms_tot_comp), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for u_comp'
  allocate ( pos_eq_EC_comp(atoms_UC_comp,3,ncells_tot_comp), stat = i )
  if ( i /= 0 ) stop ' ERROR: Allocation failed for pos_eq_EC_comp'

  ! reference supercell
  n = 0
  do i3 = 1, ncells_ref(3)
     do i2 = 1, ncells_ref(2)
        do i1 = 1, ncells_ref(1)
          n = n + 1
          xt = dble(i1-1)*side_eq_UC_ref(1,1) + dble(i2-1)*side_eq_UC_ref(2,1) + dble(i3-1)*side_eq_UC_ref(3,1)
          yt = dble(i1-1)*side_eq_UC_ref(1,2) + dble(i2-1)*side_eq_UC_ref(2,2) + dble(i3-1)*side_eq_UC_ref(3,2)
          zt = dble(i1-1)*side_eq_UC_ref(1,3) + dble(i2-1)*side_eq_UC_ref(2,3) + dble(i3-1)*side_eq_UC_ref(3,3)
          pos_eq_EC_ref(:,1,n) = pos_eq_UC_ref(:,1) + xt
          pos_eq_EC_ref(:,2,n) = pos_eq_UC_ref(:,2) + yt
          pos_eq_EC_ref(:,3,n) = pos_eq_UC_ref(:,3) + zt
        end do
     end do
  end do

  ! comparison supercell
  n = 0
  do i3 = 1, ncells_comp(3)
     do i2 = 1, ncells_comp(2)
        do i1 = 1, ncells_comp(1)
          n = n + 1
          xt = dble(i1-1)*side_eq_UC_comp(1,1) + dble(i2-1)*side_eq_UC_comp(2,1) + dble(i3-1)*side_eq_UC_comp(3,1)
          yt = dble(i1-1)*side_eq_UC_comp(1,2) + dble(i2-1)*side_eq_UC_comp(2,2) + dble(i3-1)*side_eq_UC_comp(3,2)
          zt = dble(i1-1)*side_eq_UC_comp(1,3) + dble(i2-1)*side_eq_UC_comp(2,3) + dble(i3-1)*side_eq_UC_comp(3,3)
          pos_eq_EC_comp(:,1,n) = pos_eq_UC_comp(:,1) + xt
          pos_eq_EC_comp(:,2,n) = pos_eq_UC_comp(:,2) + yt
          pos_eq_EC_comp(:,3,n) = pos_eq_UC_comp(:,3) + zt
        end do
     end do
  end do

  ! M is the rotation to apply on the acoustic modes in comp to obtain
  ! the acoustic modes in ref: E_acu_ref = M.E_acu_comp

!  do i = 1, 3
!    E_acu_ref(:,i) = dble(eig_ref(qgm_ref(1), qgm_ref(i+1), 1:3)/sqrt(mass_UC_ref(1)))
!    E_acu_comp(:,i) = dble(eig_comp(qgm_comp(1), qgm_comp(i+1), 1:3)/sqrt(mass_UC_comp(1)))
!  end do
!  M = matmul(E_acu_ref,inv(E_acu_comp))

  write(*,'(a)') ' Calculating phonon mapping:'


  open(unit=100,file='eigmap.txt',action='write')
  write(100,'(a)') '#    1  ,   2  ->   4  ,   5  ,     6     ,     7'
  write(100,'(a)') '# k_ref, j_ref -> k_cmp, j_cmp,    ang (°), freq_diff'
  
  open(unit=200,file='eigscal.txt',action='write')
  write(200,'(a)') '#    1  ,   2  ->   4  ,   5  ,     6   ,     7'
  write(200,'(a)') '# k_ref, j_ref -> k_cmp, j_cmp,    ang (°), freq_diff'

  do k1 = 1, npunique_ref
  do j1 = 1, nq_ref
    alphamin = 1.d10
  
    k = 0
    j = 0
    do k2 = 1, npunique_comp
    do j2 = 1, nq_comp
      
      ! mode contribution to unitary ref displacement
      n1 = 0
      do i = 1, atoms_UC_ref
        n = 0
        do i3 = 1, ncells_ref(3)
        do i2 = 1, ncells_ref(2)
        do i1 = 1, ncells_ref(1)
          n = n + 1
          do ix1 = 1, 3
            n1 = n1 + 1
            u_ref(n1) = dble( exp(dcmplx( 0.d0, dot_product(vec_ref(k1,:),pos_eq_EC_ref(i,:,n)) ))*eig_ref(k1,j1,(i-1)*3+ix1) )
            u_ref(n1) = u_ref(n1)/sqrt(mass_UC_ref(i))
          end do
        end do
        end do
        end do
      end do

      ! mode contribution to unitary comp displacement
      n1 = 0
      do i = 1, atoms_UC_comp
        n = 0
        do i3 = 1, ncells_comp(3)
        do i2 = 1, ncells_comp(2)
        do i1 = 1, ncells_comp(1)
          n = n + 1
          do ix1 = 1, 3
            n1 = n1 + 1
            u_comp(n1) = dble( exp(dcmplx( 0.d0, dot_product(vec_comp(k2,:),pos_eq_EC_comp(i,:,n)) ))*eig_comp(k2,j2,(i-1)*3+ix1) )
            u_comp(n1) = u_comp(n1)/sqrt(mass_UC_comp(i))
          end do
        end do
        end do
        end do
      end do

      ! rotate the atom components of u_comp according to M
!      do i = 1, natoms_tot_comp
!        u_comp((i-1)*3+1:(i-1)*3+3) = matmul(M,u_comp((i-1)*3+1:(i-1)*3+3))
!      end do

      refmod = sqrt( dot_product( u_ref(:), u_ref(:)) )
      compmod = sqrt( dot_product( u_comp(:), u_comp(:)) )
  
      ! alpha is the Euclidean angle
      alpha = abs(dot_product( u_ref(:), u_comp(:) ))
      alpha = alpha / (refmod * compmod)
      if ( alpha - 1.d0 > tiny(1.d0) ) then
        write(*,'(a,E29.18)') '  Warning: cos(alpha) = ', alpha
        alpha = 1.d0
      end if
      alpha = acos(alpha)
  
      write(200,'(i5,1x,i5,1x,a3,i4,1x,i4,1x,f8.2,1x,f12.7)') k1, j1, ' . ', k2, j2, alpha*rad_to_deg, &
           freq_comp(k2,j2)-freq_ref(k1,j1)
  
      if ( alpha < alphamin ) then
        k = k2
        j = j2
        alphamin = alpha
      end if

!      write(*,*) u_ref
!      write(*,*) u_comp
!      stop
      
    end do
    end do
  
    write(100,'(i5,1x,i5,1x,a3,i4,1x,i4,1x,f8.2,1x,f12.7)') k1, j1, ' . ', k, j, alphamin*rad_to_deg, &
         freq_comp(k,j)-freq_ref(k1,j1)
  
  end do
  end do

  close(100)
  write(*,'(a)') '  Phonon mapping written in eigmap.txt'
  close(200)
  write(*,'(a)') '  Scalar products written in eigscal.txt'


  return
end subroutine calc_dot


