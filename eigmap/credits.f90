
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! eigmap. Copyright (C) Antonio Cammarata
! https://nano.cvut.cz/researchers/antonio-cammarata
! https://orcid.org/0000-0002-5691-0682
! 
! Program to calculate the atomic character
! of phonon eigenvectors obtained from PHONOPY
! ( https://phonopy.github.io/phonopy )
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
!    along with eigmap.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine credits
  use var, only: progname

  write(*,*)
  write(*,'(*(a))') ' Suggested reference for the acknowledgment of ',progname,' usage.'
  write(*,*)
  write(*,'(*(a))') ' The users of ',progname, ' have little formal obligations'
  write(*,'(a)') ' specified in the GNU General Public License, http://www.gnu.org/copyleft/gpl.txt .'
  write(*,'(a)') ' However, it is common practice in the scientific literature,'
  write(*,'(a)') ' to acknowledge the efforts of people that have made the research possible.'
  write(*,'(a)') ' In this spirit, please cite '
  write(*,*)
  write(*,'(a)') ' Normal Dynamics: solving Newton’s equations of motion in the reciprocal space' 
  write(*,'(a)') ' A. Cammarata et al. '
  write(*,'(a)') ' https://doi.org/10.1103/xxx'
  write(*,*)
  write(*,*)

end subroutine credits
