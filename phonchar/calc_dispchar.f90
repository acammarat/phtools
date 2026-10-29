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

subroutine calc_dispchar
  use io_units, only: out_dispchar
  use pars, only: version, othread_factor
  use in_user, only: nag, ag, natg, max_num_othreads
  use in_yaml
  use functions
  use refconf 
  use omp_lib

  implicit none
!  integer, external :: lcm
  integer :: i, j, k, ix, l, ncells_tot, i1, i2 ,i3, m, n, sign_vec, mcm
  integer :: ncells(3), q(3), p(3)
  real(8) :: xt, yt, zt, wshift, cpustart, cpuend, ompstart, ompend
  real(8) :: s_tmp, s(nqp), w(nqp,nq), s1(nag,nqp,nq), maxl(nag)
  real(8), allocatable :: pos_eq_EC(:,:,:)
  complex(8) :: eexp 
  character(1000) :: frac(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phonon character based on contribution to unit displacement
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  write(*,'(a)') ' Calculating phonon character based on atom contribution to u_i = m_i^(-1/2) exp(ik.r) eig_i(k,j)... '
  write(*,'(a)') '  # q-point, commensurate supercell'

  s1(:,:,:) = 0.d0
  s(:) = 0.d0
  do k = 1, nqp
    ! the amplitude is the same for all the modes at fixed q

    call cpu_time(cpustart)
    ompstart = omp_get_wtime()

    ! print the q-point as a fraction
    do i = 1, 3
       if ( abs(vec(k,i)) > tiny(1.d0) ) then
          sign_vec = 1
          if ( vec(k,i) < 0.d0 ) sign_vec = -1
          call real_to_rational ( dble(sign_vec)*vec(k,i), p(i), q(i) )
          p(i) = sign_vec*p(i)
       else
          p(i) = 0
          q(i) = 1
       end if
       call print_fraction ( p(i), q(i), frac(i) )
    end do
    write(*,'(2x,8a)',advance='no') i2a(k), ' (', trim(frac(1)), ',', trim(frac(2)), ',', trim(frac(3)), '), '

    ! supercell commensurate with vec(k,:)
    do i = 1, 3
       if ( p(i) == 0  ) then
          mcm = 1
       else
          mcm = lcm(p(i), q(i))
       end if
       ncells(i) = abs(mcm)
    end do

    write(*,'(6a)', advance='no') i2a(ncells(1)),' x ', i2a(ncells(2)),' x ', i2a(ncells(3))
    ncells_tot = ncells(1)*ncells(2)*ncells(3)

    if ( ncells_tot > othread_factor ) then ! the factor is chosen according to local benchmarks
      call omp_set_num_threads(max_num_othreads)
    else
      call omp_set_num_threads(1)
    end if

    allocate ( pos_eq_EC(atoms_UC,3,ncells_tot), stat = i )
    if ( i /= 0 ) stop 'Allocation failed for pos_eq_EC'

    ! generate supercell atomic positions according to the reference configuration
    do i = 1, atoms_UC
       n = 0
       do i3 = 1, ncells(3)
          do i2 = 1, ncells(2)
             do i1 = 1, ncells(1)
                n = n + 1
                xt = dble(i1-1)*side_UC(1,1) + dble(i2-1)*side_UC(2,1) + dble(i3-1)*side_UC(3,1)
                yt = dble(i1-1)*side_UC(1,2) + dble(i2-1)*side_UC(2,2) + dble(i3-1)*side_UC(3,2)
                zt = dble(i1-1)*side_UC(1,3) + dble(i2-1)*side_UC(2,3) + dble(i3-1)*side_UC(3,3)
                pos_eq_EC(i,1,n) = pos_eq_UC(i,1) + xt
                pos_eq_EC(i,2,n) = pos_eq_UC(i,2) + yt
                pos_eq_EC(i,3,n) = pos_eq_UC(i,3) + zt
             end do
          end do
       end do
    end do

    do j = 1, nq
  
      do l = 1, nag ! loop over number of atomic groups
        do m = 1, natg(l) ! loop over number of atoms in group l
          do ix = 1, 3 ! loop over Cartesian components

            s_tmp = 0.d0
            !$omp parallel do schedule ( dynamic ) private (eexp) reduction(+:s_tmp)
            do n = 1, ncells_tot 
              eexp = exp( dcmplx( 0.d0, dot_product( vec(k,:), pos_eq_EC(ag(l,m),:,n) ) ) )
              ! it is not necessary to consider 2 times the real part of the displacement
              ! to calculate the weight; also, better to consider the quadratic displacement as
              ! the direction is not relevant
              s_tmp = s_tmp + realpart( eig(k,j,(ag(l,m)-1)*3+ix) * eexp ) * realpart( eig(k,j,(ag(l,m)-1)*3+ix) * eexp )
            end do
            !$omp end parallel do

          end do
          s1(l,k,j) = s1(l,k,j) + s_tmp/sqrt(mass_UC(ag(l,m)))
        end do
        s1(l,k,j) = s1(l,k,j)/sqrt(dble(ncells_tot))
        s(k) = s(k) + s1(l,k,j)
      end do

    end do

    deallocate ( pos_eq_EC, stat = i )
    if ( i /= 0 ) stop 'Deallocation failed for pos_eq_EC'

    call cpu_time(cpuend)
    ompend = omp_get_wtime()
    write(*,'(*(a))') ', CPU: ', trim(f2a(cpuend-cpustart)), ', OMP: ', trim(f2a(ompend-ompstart)), ', ratio: ', trim(f2a((cpuend-cpustart)/(ompend-ompstart)))

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
  
  open(unit=out_dispchar,file='dispchar.dat')
  write(out_dispchar,'(*(a))') '# phonchar v. ', version
  write(out_dispchar,'(*(a))') '# character: atom contribution to unitary phonon displacement u_i = m_i^(-1/2) exp(ik.r) eig_i(k,j)'
  write(out_dispchar,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(nq),' ; groups: ',i2a(nag)
  if ( nag > 2 ) then
    write(out_dispchar,'(a)') '# the weight value is equal to the group label with greatest projection'
  else
    write(out_dispchar,'(a)') '# smaller weights correspond to largest projections on group 1'
  end if

  ! unified mode index is for use with phind
  write(out_dispchar,'(a)') '# q-point, freq[THz], weight, mode index, unified mode index'
  do j = 1, nq
    do k = 1, nqp
      n = (k-1)*nq + j
      write(out_dispchar,'(a,1x,f20.10,1x,f7.2,2(1x,a))') i2a(k), freq(k,j), w(k,j), i2a(j), i2a(n)
    end do
    write(out_dispchar,*) 
    write(out_dispchar,*) 
  end do
  close(out_dispchar)

  write(*,'(a)') ' Group displacement character written in dispchar.dat'
  write(*,*)

  return
end subroutine calc_dispchar
