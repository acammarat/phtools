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

subroutine deallocate_all
  use var, only: freq_ref, freq_comp, eig_ref, eig_comp, atmap
  implicit none
  integer :: i
   
  deallocate ( atmap, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for atmap'

  deallocate ( eig_ref, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for eig_UC'
  deallocate ( freq_ref, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for freq_UC'

  deallocate ( eig_comp, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for eig_EC'
  deallocate ( freq_comp, stat = i )
  if ( i /= 0 ) stop 'Deallocation failed for freq_EC'
  
  return
end subroutine deallocate_all
