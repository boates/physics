!
!  Program: pdos.f90
!  Author: Brian Boates
!
!======================================================
!
!  Analyze and average a VASP DOSCAR file from a
!  projected DOS calculation for a given set of atoms
!  and s, p, d character (shifted by Efermi)
!
!======================================================

implicit none

integer i, j, k, a, b
integer maxAtoms, maxBins
parameter (maxAtoms=512,maxBins=1000)

integer atom_i, atom_f, natom, nbins, bin
integer atom_list(maxAtoms), len_atom_list, use_atom

real*8 s, p, d, tot
real*8 emin, emax, efermi, E(maxBins)
real*8 s_tot(maxBins), p_tot(maxBins), d_tot(maxBins), tot_tot(maxBins)
real*8 atom_k, s_k, p_k, d_k, tot_k

character*128 d1, d2, d3, d4, d5, d6, d7, d8, d9
character*128 doscar, fatom

!=============================!
! RETRIEVE USER INPUT

! PROCAR filename
write(*,*) "DOSCAR file:"
read(*,*) doscar

! which atoms to consider
write(*,*) "atom_i, atom_f (indexing starts at 1, not 0)"
write(*,*) " - OR give 0,0 for all atoms"
write(*,*) " - OR give -1,-1 to specify a file with atom indices:"
read(*,*) atom_i, atom_f
if (atom_i.eq.-1.and.atom_f.eq.-1) then
  write(*,*) "Name of atom index file:"
  read(*,*) fatom
endif
!==============================

!=========================
! OPEN DOSCAR FILE
open(1,file=doscar,status='old',ERR=100)

! remove header
read(1,*) natom, d2, d3, d4 
read(1,*) d1, d2, d3, d4, d5
read(1,*) d1
read(1,*) d1
read(1,*) d1, d2
read(1,*) emax, emin, nbins, efermi, d5

! read in total DOS, only keep the energy axis
do i=1, nbins
  read(1,*) E(i), d2, d3
  E(i) = E(i) - Efermi  ! shift by Efermi
enddo
!==========================

!=========================
! ATOMS LIST
! check if all atoms should be considered
if (atom_i.eq.0.and.atom_f.eq.0) then
  atom_i = 1
  atom_f = natom
  ! use subsequent range of atoms
  len_atom_list = 0
  do i = atom_i, atom_f
    atom_list(i) = i
    len_atom_list = len_atom_list + 1
  enddo

! check if atom indices are in a file
elseif (atom_i.eq.atom_f.eq.-1) then
  open(10,file=fatom,status='old',ERR=110)
  do i = 1, natom
    read(10,*,END=10) atom_list(i)
    len_atom_list = len_atom_list + 1
  enddo
  10 continue
  close(10)

else
  ! get range of atoms
  len_atom_list = 0
  do i = atom_i, atom_f
    atom_list(i) = i
    len_atom_list = len_atom_list + 1
  enddo
endif
!===========================!

!=============================!
!  READ IN DOS FOR ALL ATOMS  
do a=1, natom

  read(1,*) ! this is the emax, emin, nbins, efermi, d5 line (should always be same)

  ! check if atom is wanted
  use_atom = 0
  do i=1, len_atom_list
    if (a.eq.atom_list(i)) then
      use_atom = 1
    endif
  enddo

  ! check that atom is wanter
  if (use_atom.eq.1) then

    ! loop over bins
    do b=1, nbins

      ! read in pdos for atom a and bin b
      read(1,*) d1, s, p, d
      tot = s + p + d
      s_tot(b) = s_tot(b) + s
      p_tot(b) = p_tot(b) + p
      d_tot(b) = d_tot(b) + d
      tot_tot(b) = tot_tot(b) + tot

    enddo
  endif
enddo
!=============================!

!=================================================
! WRITE s, p, d, AND tot DOS FOR DESIRED ATOMS
open(2,file='pdos.dat',status='unknown')

! header info
write(2,*) '# from file: DOSCAR'!, procar
write(2,*) '# using atoms: ', atom_i, 'to', atom_f
write(2,*) '# E, total, s, p, d'

! write DOS's
5 format(f12.8,x,f12.8,x,f12.8,x,f12.8,x,f12.8)
do b=1, nbins
  ! energies, tot, s, p, & d DOS's
  write(2,5) E(b), tot_tot(b), s_tot(b), p_tot(b), d_tot(b)
enddo
!=================================================

close(2)

100 continue

close(1)

110 continue

END
