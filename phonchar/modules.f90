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
! where the formulation is reported in section V "Atomic character of the phonon modes" 
! of the supplemental material.
!
! Only the real part of the eigenvector is considered, as the purpose
! is to analyse the contribution to the atomic motions.
! 
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

module functions
  implicit none
contains

  integer function gcd ( a, b )
    implicit none
    integer, intent(in) :: a, b
    integer :: aa, bb, t
    
    aa = a
    bb = b
    do while ( bb /= 0 )
       t = bb
       bb = mod(aa,bb)
       aa = t
    end do
    gcd = abs(aa)
    
    return
  end function gcd

  integer function lcm ( a, b )
    implicit none
    integer, intent(in) :: a, b
    
    lcm = a * b / gcd(a,b)
    
    return
  end function lcm

  real(8) function rad2deg(radians)
    implicit none
    real(8), intent(in) :: radians  
    real(8), parameter :: pi = acos(-1.d0)
    rad2deg = radians * 180.d0 / pi
  end function rad2deg

  function i2a(i) result(out)
    character(:), allocatable :: out
    integer(4), intent(in) :: i
    character(range(i)+2) :: x

    write(x,'(i0)') i

    out = trim(x)
  end function i2a

  function f2a(f) result(out)
    character(20) :: out
    real(8), intent(in) :: f
    character(20) :: x

    write(x,'(f20.4)') f

    out = adjustl(x)
  end function f2a

  function wordcount(line) result(out)
    implicit none
    character(256) :: line
    integer :: i, len_line, nword, out
    logical :: in_word

    len_line = len_trim(line)
    nword = 0
    in_word = .false.

    ! Loop through each character
    do i = 1, len_line
        if (line(i:i) /= ' ' .and. line(i:i) /= char(9)) then
            if (.not. in_word) then
                nword = nword + 1
                in_word = .true.
            end if
        else
            in_word = .false.
        end if
    end do

    out = nword
  end function wordcount

end module functions

module pars
  ! strings
  character(3), parameter :: version = '2.5'
  character(8), parameter :: progname = 'PHONCHAR'
  character(6), parameter :: output_format = 'g22.14', screen_format = 'f12.4'
  character(7), parameter :: error_string = 'ERROR: '
  character(9), parameter :: warning_string = 'WARNING: '
  integer(8), parameter :: othread_factor = 1000000 ! chosen according to local benchmarks
  real(8), parameter :: pi = acos(-1.d0)
end module pars

module in_user
  integer, parameter :: natmax = 10000
  integer, save :: nag, max_num_othreads
  integer, save :: natg(natmax), run_flag(3)
  integer, save, allocatable :: ag(:,:)
  real(8), save :: dir(3), rotax(3), scanang(3)
end module in_user

module in_yaml
  character(256), save :: inyaml
  integer, save :: nqp, nq
  real(8), save, allocatable :: freq(:,:), vec(:,:)
  complex(8), save, allocatable :: eig(:,:,:)
end module in_yaml

module refconf
  integer, save :: atoms_UC, atom_types
  integer, save, allocatable :: natoms_UC(:)
  real(8), save :: side_UC(3,3)
  real(8), save, allocatable :: mass_UC(:), pos_eq_UC(:,:)
end module refconf

module io_units
  integer, parameter :: input             = 10
  integer, parameter :: inp_yaml          = 11
  integer, parameter :: inp_poscar        = 12
  integer, parameter :: out_dispchar      = 30
  integer, parameter :: out_eigchar       = 31
  integer, parameter :: out_dirchar       = 32
end module io_units

module masses
  integer, parameter :: n_phonopy = 92
  character(2), parameter, dimension(n_phonopy) :: at_phonopy = (/"H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U "/)
  real(8), parameter, dimension(n_phonopy) :: mass_phonopy =(/1.00794d0, 4.002602d0, 6.941d0, 9.012182d0, 10.811d0, 12.0107d0, 14.0067d0, 15.9994d0, 18.9984032d0, 20.1797d0, 22.98976928d0, 24.3050d0, 26.9815386d0, 28.0855d0, 30.973762d0, 32.065d0, 35.453d0, 39.948d0, 39.0983d0, 40.078d0, 44.955912d0, 47.867d0, 50.9415d0, 51.9961d0, 54.938045d0, 55.845d0, 58.933195d0, 58.6934d0, 63.546d0, 65.38d0, 69.723d0, 72.64d0, 74.92160d0, 78.96d0, 79.904d0, 83.798d0, 85.4678d0, 87.62d0, 88.90585d0, 91.224d0, 92.90638d0, 95.96d0, 98.d0, 101.07d0, 102.90550d0, 106.42d0, 107.8682d0, 112.411d0, 114.818d0, 118.710d0, 121.760d0, 127.60d0, 126.90447d0, 131.293d0, 132.9054519d0, 137.327d0, 138.90547d0, 140.116d0, 140.90765d0, 144.242d0, 145.d0, 150.36d0, 151.964d0, 157.25d0, 158.92535d0, 162.500d0, 164.93032d0, 167.259d0, 168.93421d0, 173.054d0, 174.9668d0, 178.49d0, 180.94788d0, 183.84d0, 186.207d0, 190.23d0, 192.217d0, 195.084d0, 196.966569d0, 200.59d0, 204.3833d0, 207.2d0, 208.98040d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 227.d0, 232.03806d0, 231.03588d0, 238.02891d0/)
end module masses

module String_mod

    use iso_fortran_env, only: IK => int32
    implicit none
    public

    type :: CharVec_type
        character(:)    , allocatable   :: record
    end type CharVec_type

    type :: String_type
        character(:)      , allocatable   :: value          !< The string value.
        type(CharVec_type), allocatable   :: Parts(:)       !< The string parts.
        integer(IK)                       :: nPart = 0_IK   !< The number of parts in the string.
    contains
        procedure, nopass :: split
    end type String_type

contains

    function split(string,delim,npart) result(Parts)

        implicit none
        character(len=*)    , intent(in)            :: string, delim

        integer(IK)         , intent(out), optional :: npart

        type(CharVec_type)  , allocatable           :: Parts(:)
        integer(IK)         , allocatable           :: PartEnd(:)
        integer(IK)         , allocatable           :: PartBegin(:)
        integer(IK)                                 :: dlmlenMinusOne
        integer(IK)                                 :: strlen, dlmlen, npartMax, ipart, ibeg, iend, i
        logical                                     :: npartIsPresent

        dlmlen = len(delim)
        strlen = len(string)
        npartIsPresent = present(npart)

        ! if dlm is empty, return the whole string split character by character

        if (dlmlen==0_IK) then
            allocate(Parts(strlen))
            do ipart = 1, strlen
                Parts(ipart)%record = string(ipart:ipart)
            end do
            if (npartIsPresent) npart = strlen
            return
        end if

        npartMax = 1_IK + strlen / dlmlen ! There can be at most strlen + 1 splits
        allocate(PartBegin(npartMax), PartEnd(npartMax)) ! This will contain the beginning and the ends of the splits.
        dlmlenMinusOne = dlmlen - 1_IK

        ibeg = 0_IK
        ipart = 1_IK
        PartBegin(ipart) = 1_IK
        loopParseString: do

            ibeg = ibeg + 1_IK
            iend = ibeg + dlmlenMinusOne

            if (strlen<iend) then ! the remaining part of the string is shorter than the delim
                PartEnd(ipart) = strlen
                exit loopParseString
            elseif ( string(ibeg:iend) == delim ) then
                PartEnd(ipart) = ibeg - 1_IK
                ipart = ipart + 1_IK
                PartBegin(ipart) = iend + 1_IK
                ibeg = iend
            end if

        end do loopParseString

        allocate(Parts(ipart))
        do i = 1, ipart
            Parts(i)%record = string(PartBegin(i):PartEnd(i))
        end do
        if (present(npart)) npart = ipart

        deallocate(PartBegin, PartEnd)

    end function split
    
end module String_mod
