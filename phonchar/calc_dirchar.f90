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

subroutine calc_dirchar
  use io_units, only: out_dirchar
  use pars, only: version, othread_factor, warning_string
  use in_user, only: nag, ag, natg, max_num_othreads, dir, rotax, scanang
  use in_yaml
  use functions
  use refconf 
  use omp_lib

  implicit none
  integer :: i, j, k, ix, l, ncells_tot, i1, i2 ,i3, m, n, sign_vec, mcm, scandiv
  integer :: ncells(3), q(3), p(3)
  real(8) :: xt, yt, zt, cpustart, cpuend, ompstart, ompend
  real(8) :: dpcos, mod_dir, mod_u, mod_u2, w, w2, rotang, mod_rotdir
  real(8) :: u(nag,nqp,nq,3), mass_tot(nag), s1(3), rotdir(3), rotdir_max(3)
  real(8), allocatable :: pos_eq_EC(:,:,:)
  complex(8) :: eexp
  character(len=:), allocatable :: outfile
  character(1000) :: frac(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate phonon character based on the projection of a reference direction vector
! onto the unit displacement and on the relative displacement among all the atom groups
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  write(*,'(a)') ' Calculating phonon character based on displacement projections... '
  write(*,'(a)') '  # q-point, commensurate supercell'


  ! total mass of each group
  mass_tot(:) = 0.d0
  do l = 1, nag ! loop over number of atomic groups
    do m = 1, natg(l) ! loop over number of atoms in group l
      mass_tot(l) = mass_tot(l) + mass_UC(ag(l,m))
    end do
  end do

  u(:,:,:,:) = 0.d0
  do k = 1, nqp

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

    ! the amplitude Q(k,j) is the same for all the modes at fixed q and
    ! it is assumed to be dcmplx(1.d0, 1.d0)
    do j = 1, nq
  
      do l = 1, nag ! loop over number of atomic groups
        do m = 1, natg(l) ! loop over number of atoms in group l

            s1(:) = 0.d0
            !$omp parallel do schedule ( dynamic ) private (eexp) reduction(+:s1)
            do n = 1, ncells_tot 
              eexp = exp( dcmplx( 0.d0, dot_product( vec(k,:), pos_eq_EC(ag(l,m),:,n) ) ) )
              ! it is not necessary to consider 2 times the real part of the displacement
              s1(:) = s1(:) + realpart( eig(k,j,(ag(l,m)-1)*3+1:(ag(l,m)-1)*3+3) * eexp )
            end do
            !$omp end parallel do

          u(l,k,j,:) = u(l,k,j,:) + s1(:)*sqrt(mass_UC(ag(l,m)))
        end do
        u(l,k,j,:) = u(l,k,j,:)/sqrt(dble(ncells_tot))/mass_tot(l) ! displacement of the center mass of the group [Ang]
      end do

    end do

    deallocate ( pos_eq_EC, stat = i )
    if ( i /= 0 ) stop 'Deallocation failed for pos_eq_EC'

    call cpu_time(cpuend)
    ompend = omp_get_wtime()
    write(*,'(*(a))') ', CPU: ', trim(f2a(cpuend-cpustart)), ', OMP: ', trim(f2a(ompend-ompstart)), ', ratio: ', trim(f2a((cpuend-cpustart)/(ompend-ompstart)))

  end do

  mod_dir = sqrt(dot_product(dir(:),dir(:)))
  ! dirchar_l.dat: the weight is the angle between the non-rotated direction and the
  !                displacement u(l,k,j,:) of the center mass of group l in mode (k,j)

  do l = 1, nag

    outfile='dispchar_'//i2a(l)//'.dat'
    open(unit=out_dirchar,file=outfile,action='write')
    write(out_dirchar,'(*(a))') '# phonchar v. ', version
    write(out_dirchar,'(*(a))') '# displacement character: u.dir ; group ',i2a(l)
    write(out_dirchar,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(nq)
    write(out_dirchar,'(a)') '# weight w1 is the angle (u,dir), weight w2 is the angle from acos(abs(dpcos))'
    ! unified mode index is for use with phind
    write(out_dirchar,'(a)') '# q-point, freq[THz], w1, w2, mod(u), u.dir, mode index, unified mode index'
    do j = 1, nq
      do k = 1, nqp
        mod_u = sqrt(dot_product(u(l,k,j,:),u(l,k,j,:)))
        if ( abs(mod_u) < tiny(1.d0) ) then
          w = 0.d0
        else
          dpcos = dot_product(u(l,k,j,:),dir(:))/(mod_dir*mod_u)
          if ( (abs(dpcos)-1.d0)>tiny(1.d0) ) then
            write(0,*) warning_string, ' dpcos ',dpcos
            dpcos = sign(1.d0, dpcos)
          end if
          w = rad2deg(acos(dpcos))
          w2 = rad2deg(acos(abs(dpcos)))
        end if
        n = (k-1)*nq + j
        write(out_dirchar,'(a,1x,f12.6,1x,f7.2,1x,f7.2,1x,E9.3,1x,E9.3,2(1x,a))') i2a(k), freq(k,j), w, w2, mod_u, dot_product(u(l,k,j,:),dir(:)), i2a(j), i2a(n)
      end do
      write(out_dirchar,*) 
      write(out_dirchar,*) 
    end do
    close(out_dirchar)

  end do

  write(*,'(a)') ' Displacement character u(l).dir in dispchar_l.dat files'

  do i1 = 1, nag-1
  do i2 = i1+1, nag

    outfile='dispchar_'//i2a(i1)//'_'//i2a(i2)//'.dat'
    open(unit=out_dirchar,file=outfile,action='write')
    write(out_dirchar,'(*(a))') '# phonchar v. ', version
    write(out_dirchar,'(*(a))') '# displacement character: u(l1=',i2a(i1),').u(l2=',i2a(i2),')'
    write(out_dirchar,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(nq)
    write(out_dirchar,'(a)') '# weight w1 is the angle (u,dir), weight w2 is the angle from acos(abs(dpcos))'
    ! unified mode index is for use with phind
    write(out_dirchar,'(a)') '# q-point, freq[THz], w1, w2, mod(u(l1)), mod(u(l2)), u(l1).u(l2), mode index, unified mode index'
    do j = 1, nq
      do k = 1, nqp
        mod_u = sqrt(dot_product(u(i1,k,j,:),u(i1,k,j,:)))
        mod_u2 = sqrt(dot_product(u(i2,k,j,:),u(i2,k,j,:)))
        if ( abs(mod_u)<tiny(1.d0) .or. abs(mod_u2)<tiny(1.d0) ) then
          w = 0.d0
        else
          dpcos = dot_product(u(i1,k,j,:),u(i2,k,j,:))/(mod_u*mod_u2)
          if ( (abs(dpcos)-1.d0)>tiny(1.d0) ) then
            write(0,*) warning_string, ' dpcos ',dpcos
            dpcos = sign(1.d0, dpcos)
          end if
          w = rad2deg(acos(dpcos))
          w2 = rad2deg(acos(abs(dpcos)))
        end if
        n = (k-1)*nq + j
        write(out_dirchar,'(a,1x,f12.6,1x,f7.2,1x,f7.2,1x,2(E9.3,1x),E9.3,2(1x,a))') i2a(k), freq(k,j), w, w2, mod_u, mod_u2, dot_product(u(i1,k,j,:),u(i2,k,j,:)), i2a(j), i2a(n)
      end do
      write(out_dirchar,*) 
      write(out_dirchar,*) 
    end do
    close(out_dirchar)

  end do
  end do

  write(*,'(a)') ' Displacement character u(l1).u(l2) in dispchar_l1_l2.dat files'


  ! dirchar_l1_maxproj.dat: the weight is the angle (u,R(dir)) which maximises the projection u.dir 
  if ( abs(scanang(2)-scanang(1))<tiny(1.d0) ) then
    write(*,'(a)') ' Initial and final scan angle are the same, no scan will be performed.'
  else

    scandiv = floor((scanang(2)-scanang(1))/scanang(3))
    do l = 1, nag
  
      outfile='dispchar_'//i2a(l)//'_maxproj.dat'
      open(unit=out_dirchar,file=outfile,action='write')
      write(out_dirchar,'(*(a))') '# phonchar v. ', version
      write(out_dirchar,'(*(a))') '# displacement character: u.R(dir) ; group ',i2a(l)
      write(out_dirchar,'(*(a))') '# q-points: ',i2a(nqp),' ; bands: ',i2a(nq)
      write(out_dirchar,'(a)') '# the weight is the rotation angle for R(dir) which maximises the projection u.R(dir)'
      ! unified mode index is for use with phind
      write(out_dirchar,'(a)') '# q-point, freq[THz], weight, mod(u), u.R(dir), mode index, unified mode index'
      do j = 1, nq
        do k = 1, nqp
          mod_u = sqrt(dot_product(u(l,k,j,:),u(l,k,j,:)))
          w = 1.d3
          if ( abs(mod_u) > tiny(1.d0) ) then
            do ix = 0, scandiv
              rotang = scanang(1) + scanang(3)*ix
              call rotvec(dir,rotax,rotang,rotdir)
              mod_rotdir = sqrt(dot_product(rotdir(:),rotdir(:)))
              dpcos = dot_product(u(l,k,j,:),rotdir(:))/(mod_rotdir*mod_u)
              if ( (abs(dpcos)-1.d0)>tiny(1.d0) ) then
                write(0,*) warning_string, ' dpcos ',dpcos
                dpcos = sign(1.d0, dpcos)
              end if
              if ( rad2deg(acos(dpcos)) < w ) then
                w = rad2deg(acos(dpcos))
                rotdir_max(:) = rotdir(:)
              end if
            end do
          end if
          n = (k-1)*nq + j
          write(out_dirchar,'(a,1x,f12.6,1x,f7.2,1x,E9.3,1x,E9.3,2(1x,a))') i2a(k), freq(k,j), w, mod_u, dot_product(u(l,k,j,:),rotdir_max(:)), i2a(j), i2a(n)
        end do
        write(out_dirchar,*) 
        write(out_dirchar,*) 
      end do
      close(out_dirchar)
  
    end do
    write(*,'(a)') ' Character with maximised projection u.R(dir) in dispchar_l_maxproj.dat files'

  end if



  return
end subroutine calc_dirchar

