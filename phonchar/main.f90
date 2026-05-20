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
! If used for production, you should cite
! Phys. Rev. B 103, 035406 (2021)
! https://doi.org/10.1103/PhysRevB.103.035406
! where the original formulation is reported in section V "Atomic character of the phonon modes" 
! of the supplemental material.
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
!
! v 2.5
! - In dirchar_l.dat, dirchar_l_maxproj.dat and dirchar_l1_l2.dat, added a second weight w2 as an angle in [0,90] deg
!   to take into account that parallel and antiparallel displacement directions might be considered equivalent like
!   for example, in the case of layer sliding.
! - corrected a bug when calculating the weight for dirchar_l1_l2.dat 
!
! v 2.4
! Modulus of the displacements is written in dispchar_l.dat, dirchar_l_maxproj.dat and dirchar_l1_l2.dat.
!
! v 2.3
! - at each q-point, the weight of the two-group case is normalised to 1 and centered at 0
! - phonon atomic group character is now calculated in different ways:
!     1) eigchar.dat: the weight is the square modulus of the eigenvector atom component eig_i(k,j).eig_i(k,j)*
!     2) dirchar_l.dat: the weight is the angle between the input direction and the
!        displacement u(l,k,j,:) of the center mass of group l in mode (k,j)
!     3) dirchar_l_maxproj.dat: the weight is the angle between the center mass displacement u 
!        and the rotated direction which maximises the scalar product dir.u
!     4) dirchar_l1_l2.dat: the weight is the angle in the scalar product u(l1,k,j,:).u(l2,k,j,:)
!        among all the possible group couples (l1,l2)
!     5) dispchar.dat: the weight is the atom contribution to unitary phonon displacement
!        u_i = m_i^(-1/2) exp(ik.r) eig_i(k,j)' 
! - at each q-point, the displacement u is generated commensurate with the q-point
!
!
! v 2.2
! Removed unnecessary variables
!
! v 2.1
! The line-skip for the search in the band.yaml file is now general
!
! v 2.0
! It is possible to specify up to ngmax groups of atoms
! If the number of atomic groups is 2, the weight is calculated
! as in v 1.0.
!
! v 1.0
! It is possible to specify only two groups of atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Format of the input file
!
! int int int             1 or 0, flag to specify the character to calculate:
!                         eigenvalues, direction, displacement.
! char [char]             yaml file [POSCAR]: if yaml file is of kind qpoints.yaml
!                         the corresponding POSCAR geometry must be provided
! real real real          Cartesian (not crystallographic!) drift vector components
! real real real          rotation axis
! real real real          starting angle, final angle, step
!                         if starting angle = final angle, no scan is performed
! int                     number of atomic groups
! int                     number of atoms in group 1
! int                     label of the first atom in the group
! ...
! int                     label of the last atom in the group
! int                     number of atoms in group 2
! int                     label of the first atom in the group
! ...
! int                     label of the last atom in the group
! ...                     number of atoms in group 3 - if any, and similar
!                         input structure as above
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program phonchar
  use in_user, only: run_flag

  implicit none

  call init

  call read_eigen

  if ( run_flag(1) == 1 ) then
    call calc_eigchar
  end if

  if ( run_flag(2) == 1 ) then
    call calc_dirchar
  end if

  if ( run_flag(3) == 1 ) then
    call calc_dispchar
  end if

  call deallocate_all
  call credits

end program phonchar
