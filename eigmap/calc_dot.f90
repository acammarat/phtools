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
  use var, only: npunique, nq, eig_ref, eig_comp, freq_ref, &
                 freq_comp, rad_to_deg, pi

  implicit none
  integer :: k1, j1, k2, j2, k, j, i1
  real(8) :: alphamin, eigrefmod, eigcompcmod, thetaC_re, thetaC_im
  real(8) :: rho, theta_H, phi, theta, eigref_re_mod, eigcomp_re_mod, alpha
  real(8) :: theta_kj, phi_kj, theta_H_kj
  
  write(*,'(a)',advance='no') ' Calculating phonon mapping...'

  open(unit=100,file='eigmap.txt',action='write')
  write(100,'(a)') '#    1  ,   2  ->   4  ,   5  ,     6   ,     7    ,     8   ,     9     ,    10'
  write(100,'(a)') '# k_ref, j_ref -> k_cmp, j_cmp, disp_ang, freq_diff, Eucl_ang, pseudo_ang, Herm_ang'

  open(unit=101,file='eigscal.txt',action='write')
  write(101,'(a)') '#    1  ,   2  ->   4  ,   5  ,     6   ,     7    ,     8   ,     9     ,    10'
  write(101,'(a)') '# k_ref, j_ref -> k_cmp, j_cmp, disp_ang, freq_diff, Eucl_ang, pseudo_ang, Herm_ang'

  ! The angle definitions follow K. Scharnhorst, Acta Applicandae Mathematicae 69, 95 (2001)
  ! https://doi.org/10.1023/A:1012692601098

  do k1 = 1, npunique
  do j1 = 1, nq
    alphamin = 1.d10
    theta_kj = 1.d10
    phi_kj = 1.d10
    theta_H_kj = 1.d10

    k = 0
    j = 0
    do k2 = 1, npunique
    do j2 = 1, nq

      ! We can safely take only the real part of the square root because the
      ! scalar product yields the squared norm of the eigenvector
      eigrefmod = dble(sqrt( dot_product( eig_ref(k1,j1,:), eig_ref(k1,j1,:)) ))
      eigcompcmod = dble(sqrt( dot_product( eig_comp(k2,j2,:), eig_comp(k2,j2,:)) ))

      thetaC_re = dble( dot_product( eig_ref(k1,j1,:), eig_comp(k2,j2,:) )/(eigrefmod * eigcompcmod) )
      thetaC_im = aimag( dot_product( eig_ref(k1,j1,:), eig_comp(k2,j2,:) )/(eigrefmod * eigcompcmod) )

      ! cos theta_C(a,b) = rho.exp(i.phi) eq. (4)
      rho = sqrt( thetaC_re*thetaC_re + thetaC_im*thetaC_im )

      ! Hermitian angle eq. (5)
      theta_H = acos(rho)

      ! Kasner's pseudoangle eq. (4)
      if ( abs(thetaC_re) < tiny (1.d0) ) then
        phi = sign(0.5d0,thetaC_re) * pi
      else
        phi = atan(thetaC_im / thetaC_re)
      end if

      ! Euclidean angle eq. (2)
      theta = 0.d0
      do i1 = 1, nq
        theta = theta + dble(eig_ref(k1,j1,i1))*dble(eig_comp(k2,j2,i1)) + &
                                aimag(eig_ref(k1,j1,i1))*aimag(eig_comp(k2,j2,i1))
      end do
      theta = theta / (eigrefmod * eigcompcmod)
      theta = acos(theta)

      ! Euclidean angle of real part only (e.g. for displacement comparison)
      eigref_re_mod = sqrt( dot_product( dble(eig_ref(k1,j1,:)) , dble(eig_ref(k1,j1,:)) ) )
      eigcomp_re_mod = sqrt( dot_product( dble(eig_comp(k2,j2,:)) , dble(eig_comp(k2,j2,:)) ) )
      alpha = dot_product( dble(eig_ref(k1,j1,:)) , dble(eig_comp(k2,j2,:)) )
      alpha = alpha / ( eigref_re_mod * eigcomp_re_mod )
      ! We care only about the relative direction the displacements as the motion is harmonic
      alpha = acos(abs(alpha))

      write(101,'(i4,i4,a3,i4,i4,f8.2,1x,f12.7,3(1x,f8.2))') k1, j1, ' . ', k2, j2, alpha*rad_to_deg, &
           freq_comp(k2,j2)-freq_ref(k1,j1), theta*rad_to_deg, phi*rad_to_deg, theta_H*rad_to_deg

      if ( alpha < alphamin ) then
        k = k2
        j = j2
        alphamin = alpha
        theta_kj = theta
        phi_kj = phi
        theta_H_kj = theta_H
      end if
      
    end do
    end do

    write(100,'(i4,i4,a3,i4,i4,f8.2,1x,f12.7,3(1x,f8.2))') k1, j1, ' . ', k, j, alphamin*rad_to_deg, &
         freq_comp(k,j)-freq_ref(k1,j1), theta_kj*rad_to_deg, phi_kj*rad_to_deg, theta_H_kj*rad_to_deg

  end do
  end do

  write(*,'(a)') ' done.'

  close(100)
  write(*,'(a)') ' Phonon mapping written in eigmap.txt'

  close(101)
  write(*,'(a)') ' Scalar products written in eigscal.txt'

  return
end subroutine calc_dot


