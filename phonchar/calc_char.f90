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

subroutine calc_char
  use var, only: nqp, neig, eig, freq, nag, ag, natg, version
  use integer_to_string

  implicit none
  integer :: i, j, k, ix, l
  real(8) :: s(nqp), w(nqp,neig)
  real(8) :: s1(nag,nqp,neig), maxl(nag)
  
  write(*,'(a)',advance='no') ' Calculating phonon character... '

  s1(:,:,:) = 0.d0
  s(:) = 0.d0
  do i = 1, nqp
    do j = 1, neig

      do l = 1, nag
        do k = 1, natg(l)
          do ix = 1, 3
            s1(l,i,j) = s1(l,i,j) + eig(i,j,(ag(l,k)-1)*3+ix)*eig(i,j,(ag(l,k)-1)*3+ix)
          end do
        end do
        s(i) = s(i) + s1(l,i,j)
      end do


    end do
  end do

  if ( nag == 2 ) then

    do i = 1, nqp
      do j = 1, neig
        w(i,j) = (-s1(1,i,j)+s1(2,i,j)) / s(i)
      end do
    end do

  else

    do i = 1, nqp
      do j = 1, neig
        maxl(:) = s1(:,i,j)
        w(i,j) = dble(maxloc(maxl,dim=1))
      end do
    end do

  end if
  
  open(unit=30,file='phchar.dat')
  write(30,'(*(a))') '# phonchar v. ', version
  write(30,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(neig),' ; groups: ',i2a(nag)
  if ( nag > 2 ) then
    write(30,'(a)') '# the weight value is equal to the group label with greatest projection'
  else
    write(30,'(a)') '# smaller weights correspond to largest projections on group 1'
  end if

  write(30,'(a28)') '# q-point, freq[THz], weight'
  do j = 1, neig
    do i = 1, nqp
      write(30,'(a,1x,E12.6,1x,E12.6)') i2a(i), freq(i,j), w(i,j)
    end do
    write(30,*) 
    write(30,*) 
  end do
  close(30)

  write(*,'(a)') 'done.'

  write(*,'(a)') ' Output written in phchar.dat'

  return
end subroutine calc_char



