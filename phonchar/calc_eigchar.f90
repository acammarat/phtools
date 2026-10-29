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

subroutine calc_eigchar
  use io_units, only: out_eigchar
  use pars, only: version
  use in_user, only: nag, ag, natg
  use in_yaml
  use functions
  use refconf 

  implicit none
!  integer, external :: lcm
  integer :: j, k, l, m
  real(8) :: wshift
  real(8) :: s(nqp), w(nqp,nq), s1(nag,nqp,nq), maxl(nag)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phonon character based on contribution to eigenvector components
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,'(a)',advance='no') ' Calculating phonon character based on eigenvector components... '

  s1(:,:,:) = 0.d0
  s(:) = 0.d0
  do k = 1, nqp
    ! the amplitude is the same for all the modes at fixed q

    do j = 1, nq
  
      do l = 1, nag ! loop over number of atomic groups
        do m = 1, natg(l) ! loop over number of atoms in group l
          ! the contribution is evaluated as the square modulus of the eigenvector atom component
          s1(l,k,j) = s1(l,k,j) + realpart(dot_product ( eig(k,j,(ag(l,m)-1)*3+1:(ag(l,m)-1)*3+3) , eig(k,j,(ag(l,m)-1)*3+1:(ag(l,m)-1)*3+3) ))
        end do
        s(k) = s(k) + s1(l,k,j)
      end do

    end do

  end do

  if ( nag == 2 ) then
    do k = 1, nqp
      do j = 1, nq
        w(k,j) = (-s1(1,k,j)+s1(2,k,j)) / s(k)
      end do
      wshift = (maxval(w(k,:))+minval(w(k,:)))/2.d0
      w(k,:) = w(k,:) - wshift
      w(k,:) = w(k,:) / maxval(w(k,:))
    end do
  else
    do k = 1, nqp
      do j = 1, nq
        maxl(:) = s1(:,k,j)
        w(k,j) = dble(maxloc(maxl,dim=1))
      end do
    end do
  end if
  
  open(unit=out_eigchar,file='eigchar.dat')
  write(out_eigchar,'(*(a))') '# phonchar v. ', version
  write(out_eigchar,'(*(a))') '# character: eig_i(k,j).eig_i(k,j)*'
  write(out_eigchar,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(nq),' ; groups: ',i2a(nag)
  if ( nag > 2 ) then
    write(out_eigchar,'(a)') '# the weight value is equal to the group label with greatest projection'
  else
    write(out_eigchar,'(a)') '# smaller weights correspond to largest projections on group 1'
  end if

  ! unified mode index is for use with phind
  write(out_eigchar,'(a)') '# q-point, freq[THz], weight, mode index, unified mode index'
  do j = 1, nq
    do k = 1, nqp
      m = (k-1)*nq + j
      write(out_eigchar,'(a,1x,f20.10,1x,f7.2,2(1x,a))') i2a(k), freq(k,j), w(k,j), i2a(j), i2a(m)
    end do
    write(out_eigchar,*) 
    write(out_eigchar,*) 
  end do
  close(out_eigchar)

  write(*,'(a)') ' done.'

  write(*,'(a)') ' Eigenvector character written in eigchar.dat'
  write(*,*)

  return
end subroutine calc_eigchar

