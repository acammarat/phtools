!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!                                                                             
! qpoints. Copyright (C) Antonio Cammarata 
! https://nano.cvut.cz/researchers/antonio-cammarata 
! https://orcid.org/0000-0002-5691-0682              
!                                       
! Extracts phonon eigenvectors and eigenvalues from 
! the file qpoints.yaml generated by PHONOPY 
! ( https://phonopy.github.io/phonopy )        
!                                       
! If used for production, you should cite 
! Phys. Rev. B XX, XXXXX (XXXX)    
! https://doi.org/10.1103/xxx 
!                                                                    
!    This file is part of qpoints.                                        
!                                                                      
!    qpoints is free software: you can redistribute it and/or modify 
!    it under the terms of the GNU General Public License as published by 
!    the Free Software Foundation, either version 3 of the License, or 
!    (at your option) any later version.                            
!                                                                 
!    qpoints is distributed in the hope that it will be useful, 
!    but WITHOUT ANY WARRANTY; without even the implied warranty of 
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      
!    GNU General Public License for more details.                     
! 
!    You should have received a copy of the GNU General Public License        
!    along with phonchar.  If not, see <http://www.gnu.org/licenses/>. 
!                 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine qextract
  use var, only: npoints, side_eq_RC, q, version, nskip, skipmod
  use functions, only: i2a

  implicit none
  integer :: i, j, k, l, natoms, nq
  integer :: npG, npH, npS
  real(8), parameter :: toldelta = 1.d-4
  real(8) :: freq, er, ei, kx, ky, kz
  character(20) :: word
  character(1) :: dum
  logical :: file_exists, warnH

  inquire(file='qpoints.yaml',exist=file_exists)
  if ( .not. (file_exists) ) then
    write(*,'(a)') ' ERROR: qpoints.yaml not found.'
    write(*,*)
    stop
  end if

  npG = 0
  npH = 0
  npS = 0
  do i = 1, npoints
    if ( (abs(q(i,1))<tiny(1.d0)) .and. (abs(q(i,2))<tiny(1.d0)) .and. (abs(q(i,3))<tiny(1.d0)) ) then
      npG = 1
    else if ( (abs(abs(q(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(q(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(q(i,3))-0.5d0)<tiny(1.d0)) ) then
      npH = npH + 1
    else
      npS = npS + 1
    end if
  end do
  write(*,*)
  write(*,'(*(a))') ' q-set composition: G ',i2a(npG),', H ',i2a(npH),', S ',i2a(npS)

  do i = 1, npG
    if ( .not. ((abs(q(1,1))<tiny(1.d0)) .and. (abs(q(1,2))<tiny(1.d0)) .and. (abs(q(1,3))<tiny(1.d0))) ) then
      write(*,'(a)') '  WARNING: Gamma is in the set but is not the first point, the output cannot be used for phind.'
      write(*,'(a)') '           Use q4phi to obtain an ordered set.'
    end if
  end do

  warnH = .true.
  do i = npG+1, npG+npH
    if ( warnH ) then
      if ( (npG == 1) .and. (.not. ((abs(abs(q(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(q(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(q(i,3))-0.5d0)<tiny(1.d0)))) ) then
        write(*,'(a)') '  WARNING: H set present but not following Gamma, the output cannot be used for phind.'
        write(*,'(a)') '           Use q4phi to obtain an ordered set.'
        warnH = .false.
      else if ( (npG == 0) .and. (.not. ((abs(abs(q(i,1))-0.5d0)<tiny(1.d0)) .or. (abs(abs(q(i,2))-0.5d0)<tiny(1.d0)) .or. (abs(abs(q(i,3))-0.5d0)<tiny(1.d0)))) ) then
        write(*,'(a)') '  WARNING: H set present but not as first set, the output cannot be used for phind.'
        write(*,'(a)') '           Use q4phi to obtain an ordered set.'
        warnH = .false.
      end if
    end if
  end do

  do k = 1, npoints
    if ( (q(k,1)+0.5d0<=tiny(1.d0)) .or. (q(k,1)-0.5d0>tiny(1.d0)) .or. &
         (q(k,2)+0.5d0<=tiny(1.d0)) .or. (q(k,2)-0.5d0>tiny(1.d0)) .or. &
         (q(k,3)+0.5d0<=tiny(1.d0)) .or. (q(k,3)-0.5d0>tiny(1.d0)) ) then
      write(*,'(a,a,3(1x,f10.5))') ' WARNING: q-components out of (-1/2,1/2] range: q',i2a(k),q(k,:)
      write(*,'(a)') '          the output cannot be used for phind.'
    end if
  end do

  ! check if the provided q-point set
  ! contains any complex conjugated couple

  do k = 1 , npoints-1
    do i = k+1, npoints
      if ( ( abs(q(k,1)+q(i,1)) < toldelta ) .and. ( abs(q(k,2)+q(i,2)) < toldelta ) .and. &
           ( abs(q(k,3)+q(i,3)) < toldelta ) ) then
        write(*,'(*(a))') ' WARNING: vectors ',i2a(k),' and ',i2a(i), ' are a conjugated couple.'
        write(*,'(a)') '          the output cannot be used for phind.'
        write(*,'(a)') '          Use q4phi to remove the conjugated couples.'
      end if
    end do
  end do


  write(*,*)
  write(*,'(a)') ' Extracting eigenvectors and frequencies:'

  open(unit=15,file='qpoints.yaml',action='READ')

  read(15,*)
  read(15,*) word, natoms
  nq = natoms*3

  do i = 1, 5
    read(15,*)
  end do

  open(unit=20,file='qmatrix.nd')
  write(20,'(*(a))') '# qpoints v. ', version
  open(unit=22,file='freq.nd')
  write(22,'(*(a))') '# qpoints v. ', version
  write(22,'(*(a))') '# Gamma: ',i2a(npG),' H: ',i2a(npH),' S: ', i2a(npS),' nq: ', i2a(nq)
  write(22,'(*(a))',advance='no') '# Skip modes: ',i2a(nskip)
  do i = 1, nskip-1
    write(22,'(*(a))',advance='no') '  ',i2a(skipmod(i,1)),' ',i2a(skipmod(i,2))
  end do
  write(22,'(*(a))') '  ',i2a(skipmod(nskip,1)),' ',i2a(skipmod(nskip,2))
  write(22,'(a)') '# Frequencies [THz] without 2pi at q-point:'

  do i = 1, npoints

    kx = q(i,1)*side_eq_RC(1,1) + q(i,2)*side_eq_RC(2,1) + q(i,3)*side_eq_RC(3,1)
    ky = q(i,1)*side_eq_RC(1,2) + q(i,2)*side_eq_RC(2,2) + q(i,3)*side_eq_RC(3,2)
    kz = q(i,1)*side_eq_RC(1,3) + q(i,2)*side_eq_RC(2,3) + q(i,3)*side_eq_RC(3,3)
    write(22,'(a2,3(f20.15,1x),a9)') '# ', kx, ky, kz, ' [Ang^-1]'
    write(22,'(a2,i4,1x,3(f20.15,1x),a13)') '# ', i, q(i,:), ' [red. coor.]'

    read(15,*)
    read(15,*)
    do j = 1, nq
      read(15,*)
      read(15,*) word, freq
      write(22,'(f20.10)') freq
      read(15,*)
      do k = 1, natoms-1
        read(15,*)
        do l = 1, 3
          read(15,*) dum, dum, er, ei
          write(20,'(2(f17.14,1x))',advance='no') er, ei
        end do
      end do

      read(15,*)
      read(15,*) dum, dum, er, ei
      write(20,'(2(f17.14,1x))',advance='no') er, ei
      read(15,*) dum, dum, er, ei
      write(20,'(2(f17.14,1x))',advance='no') er, ei
      read(15,*) dum, dum, er, ei
      write(20,'(f17.14,1x,f17.14)') er, ei
    end do

    write(*,'(a3,a,a3,3(f8.5,1x))') '  q', i2a(i), ' : ', q(i,:)

    read(15,*)

  end do

  deallocate ( q, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for q'

  close(15)
  close(20)
  close(22)

end subroutine qextract