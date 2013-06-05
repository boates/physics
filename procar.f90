!
!  Program: procar.f90
!  Author: Brian Boates
!
!======================================================
!
!  Analyze and average a VASP PROCAR file for specific
!  atoms, bands, k-point, and s/p/d/tot character
!
!======================================================

implicit none

integer i, j, k, a, b, c
integer maxAtoms, maxBands, maxKpts, maxBins
parameter (maxAtoms=1024,maxBands=2000,maxKpts=1024,maxBins=1000)

integer atom_i, atom_f, band_i, band_f, kpt_i, kpt_f
integer nbins, bin, natom, nkpt, nband, iband
integer atom_list(maxAtoms), band_list(maxBands), kpt_list(maxKpts)
integer len_atom_list, len_band_list, len_kpt_list
integer use_atom, use_band, use_kpt

real*8 sigma, E, w, arg1, arg2
real*8 eband, oband, emin, emax, de, energies(maxBins)
real*8 s(maxBins), p(maxBins), d(maxBins), tot(maxBins)
real*8 atom_k, s_k, p_k, d_k, tot_k

character*128 d1, d2, d3, d4, d5, d6, d7, d8, d9
character*128 procar, fatom, fband, fkpt

!=====================!
! RETRIEVE USER INPUT !
!=====================!

! PROCAR filename
write(*,*) "PROCAR file:"
read(*,*) procar

! which atoms to consider
write(*,*) "atom_i, atom_f (indexing starts at 1, not 0)"
write(*,*) " - OR give 0,0 for all atoms"
write(*,*) " - OR give -1,-1 to specify a file with atom indices:"
read(*,*) atom_i, atom_f
if (atom_i.eq.-1.and.atom_f.eq.-1) then
  write(*,*) "Name of atom index file:"
  read(*,*) fatom
endif

! which bands to consider
write(*,*) "band_i, band_f (indexing starts at 1, not 0)"
write(*,*) " - OR give 0,0 for all bands"
write(*,*) " - OR give -1,-1 to specify a file with band indices:"
read(*,*) band_i, band_f
if (band_i.eq.-1.and.band_f.eq.-1) then
  write(*,*) "Name of atom index file:"
  read(*,*) fband
endif

! which kpts to consider
write(*,*) "kpt_i, kpt_f (indexing starts at 1, not 0)"
write(*,*) " - OR give 0,0 for all kpts"
write(*,*) " - OR give -1,-1 to specify a file with kpt indices:"
read(*,*) kpt_i, kpt_f
if (kpt_i.eq.-1.and.kpt_f.eq.-1) then
  write(*,*) "Name of kpt index file:"
  read(*,*) fkpt
endif

! number of bins
write(*,*) "Number of bins (~200):"
read(*,*) nbins

! smearing parameter
write(*,*) "Amount of smearing (i.e. 0.1):"
read(*,*) sigma


!==============================!
! DETERMINE BAND EMIN AND EMAX !
! AND GET HEADER INFO          !
!==============================!

! open PROCAR file
open(1,file=procar,status='old',ERR=100)

! header info
read(1,*) 
read(1,*) d1, d2, d3, nkpt, d4, d5, d6, nband, d7, d8, d9, natom

! loop over kpts
do i = 1, nkpt
  read(1,*) ! blank line
  read(1,*) ! current kpt

  ! loop over bands
  do j = 1, nband
    read(1,*) ! blank line
    read(1,*) d1, d2, d3, d4, eband

    ! check for new emin
    if (eband.lt.emin) then
      emin = eband
    endif

    ! check for new emax
    if (eband.gt.emax) then
        emax = eband
    endif

    read(1,*) ! blank line
    read(1,*) ! ion, s, p, d, tot line

    ! loop over atoms
    do k = 1, natom
      read(1,*)
    enddo

    read(1,*) ! tot line
  enddo

  read(1,*) ! blank line
enddo

close(1)

!=================================================

! round off emin and emax
emin = float( int(emin -1.5) )
emax = float( int(emax +1.5) )

! total energy range
de = emax - emin

!=================================================

!=====================================!
! DETERMINE ATOM, BAND, AND KPT LISTS !
!=====================================!

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

!=========================
! BANDS LIST

! check if all bands should be considered
if (band_i.eq.0.and.band_f.eq.0) then
  band_i = 1
  band_f = nband
  len_band_list = 0
  do i = band_i, band_f
    band_list(i) = i
    len_band_list = len_band_list + 1
  enddo

! check if band indices are in a file
elseif (band_i.eq.-1.and.band_f.eq.-1) then
  open(11,file=fband,status='old',ERR=111)
  do i = 1, nband
    read(11,*,END=11) band_list(i)
    len_band_list = len_band_list + 1
  enddo
  11 continue
  close(11)

else
  ! get range of bands
  len_band_list = 0
  do i = band_i, band_f
    band_list(i) = i
    len_band_list = len_band_list + 1
  enddo
endif

!=========================
! KPTS LIST

! check if all kpts should be considered
if (kpt_i.eq.0.and.kpt_f.eq.0) then
  kpt_i = 1
  kpt_f = nkpt
  ! use subsequent range of kpts
  len_kpt_list = 0
  do i = kpt_i, kpt_f
    kpt_list(i) = i
    len_kpt_list = len_kpt_list + 1
  enddo

! check if kpt indices are in a file
elseif (kpt_i.eq.-1.and.kpt_f.eq.-1) then
  open(12,file=fkpt,status='old',ERR=112)
  do i = 1, nkpt
    read(12,*,END=12) kpt_list(i)
    len_kpt_list = len_kpt_list + 1
  enddo
  12 continue
  close(12)

else
  ! get range of kpts
  len_kpt_list = 0
  do i = kpt_i, kpt_f
    kpt_list(i) = i
    len_kpt_list = len_kpt_list + 1
  enddo
endif

! write(*,*) len_atom_list, len_band_list, len_kpt_list
! write(*,*) atom_list(320), band_list(304), kpt_list(4)

!=================================================

! re-open PROCAR file
open(1,file=procar,status='old',ERR=100)

! remove header
read(1,*) 
read(1,*) 

! loop over kpts
do i = 1, nkpt
  read(1,*) ! blank line
  read(1,*) ! d1, d2, d3, d4, d5, d6, d7, d8, d9, weights(i)

  ! check if current kpt should be considered
  use_kpt = 0
  do c = 1, len_kpt_list
    if (i.eq.kpt_list(c)) then
      use_kpt = 1
    endif
  enddo
!  if (use_kpt.ne.1) then
!    write(*,*) 'use_kpt = ',use_kpt
!  endif

  ! loop over bands
  do j = 1, nband
    read(1,*) ! blank line
    read(1,*) d1, iband, d3, d4, eband, d5, d6, oband
    read(1,*) ! blank line
    read(1,*) ! ion, s, p, d, tot line

    ! check if current band should be considered
    use_band = 0
    do b = 1, len_band_list
      if (iband.eq.band_list(b)) then
        use_band = 1
      endif
    enddo
!    if (use_band.ne.1) then
!      write(*,*) 'use_band = ',use_band
!    endif

    ! loop over atoms
    do k = 1, natom

      ! read s, p, d, tot for current atom
      read(1,*) atom_k, s_k, p_k, d_k, tot_k

      ! check if current atom should be considered
      use_atom = 0
      do a = 1, len_atom_list
        if (atom_k.eq.atom_list(a)) then
          use_atom = 1
        endif
      enddo
!      if (use_atom.ne.1) then
!        write(*,*) 'use_atom = ',use_atom
!      endif

      ! If the current atom, band, and kpt should be considered
      if (use_kpt.eq.1.and.use_band.eq.1.and.use_atom.eq.1) then

        ! build each DOS with smearing
        if (sigma.gt.0) then

          ! loop over energy bins
          do b = 1, nbins

            ! simplify some variable names
            E = emin + b*(de/nbins)
!!!!            E = emin + (b-1)*(de/nbins)
            w = (1.0 / len_kpt_list) *2 ! equal kpt weighting & 2 e per band
            arg1 = (E + de/nbins - eband) / (2**0.5 * sigma)
            arg2 = (E - eband) / (2**0.5 * sigma)

            s(b)   = s(b)   + w/2. *( erf(arg1) - erf(arg2) ) *s_k
            p(b)   = p(b)   + w/2. *( erf(arg1) - erf(arg2) ) *p_k
            d(b)   = d(b)   + w/2. *( erf(arg1) - erf(arg2) ) *d_k
            tot(b) = tot(b) + w/2. *( erf(arg1) - erf(arg2) ) *tot_k

          enddo

        ! if sigma = 0 (no smearing)
        elseif (sigma.eq.0) then
          bin = int( (eband-emin)*(nbins/de) +0.5 )
          s(bin)   = s(bin)   + s_k
          p(bin)   = p(bin)   + p_k
          d(bin)   = d(bin)   + d_k
          tot(bin) = tot(bin) + tot_k
        endif
      
      endif

    enddo

    read(1,*) ! tot line
  enddo

  read(1,*) ! blank line
enddo

! create energy array
do i = 1, nbins
  energies(i) = emin + i*(de/nbins)
enddo

!===========================================================

! write s, p, d, and tot DOS for desired atoms/bands/kpts to file
open(2,file='pdos.dat',status='unknown')

! header info
write(2,*) '# from file: PROCAR'!, procar
write(2,*) '# using atoms: ', atom_i, 'to', atom_f
write(2,*) '# using bands: ', band_i, 'to', band_f
write(2,*) '# using kpts: ', kpt_i, 'to', kpt_f
write(2,*) '# E, total, s, p, d'

! write DOS's
5 format(f12.8,x,f12.8,x,f12.8,x,f12.8,x,f12.8)
do i = 1, nbins

  ! energies, tot, s, p, & d DOS's
  write(2,5) energies(i), tot(i), s(i), p(i), d(i)

enddo

close(2)

goto 100

110 continue
goto 100

111 continue
close(10)
goto 100

112 continue
close(10)
close(11)
goto 100

100 continue

close(1)

END
