
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
  write(*,'(a)') ' Fine control of lattice thermal conductivity in low-dimensional materials'
  write(*,'(a)') ' A. Cammarata, T. Polcar, Phys. Rev. B 103, 035406 (2021)'
  write(*,'(a)') ' https://doi.org/10.1103/PhysRevB.103.035406'
  write(*,*)
  write(*,'(a)') ' where the formulation used to calculate the phonon atomic character'
  write(*,'(a)') ' is reported in section V "Atomic character of the phonon modes" of the Supplemental Material.'
  write(*,*)

end subroutine credits
