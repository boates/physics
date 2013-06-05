! program nn_dist_two_species.f90
!****************************************************
! Create TRAJEC.cnn and nn_dist.hist files
! for one species w.r.t. another
!****************************************************

implicit none

integer maxSteps, maxAtoms, maxBins, maxWrite, i, j, k, l
parameter (maxSteps=100000,maxAtoms=100000,maxBins=1000,maxWrite=100)
integer natom, natom1, natom2, tstep, nsteps, nn(maxAtoms)
integer nbins, nwrite, bin, index, list(maxAtoms), pbc_round(3)
real*8 r_nn(maxAtoms), sorted(maxAtoms), arg
real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
real*8 ax, ay, az, dx, dy, dz, r, binsize
real*8 hist(maxWrite,maxBins), input_value(3)
character*2 typat(maxAtoms), typat1, typat2
character*128 fxyz

! Retrieve user input
write(*,*) "TRAJEC.xyz file:"
read(*,*) fxyz
write(*,*) "ax, ay, az (angstroms):"
read(*,*) ax, ay, az
write(*,*) "Number of bins (~200):"
read(*,*) nbins
write(*,*) "Number of neighbour distributions to write (max=natom2 or 100):"
read(*,*) nwrite
write(*,*) "Atom type 1, type 2:"
read(*,*) typat1, typat2
write(*,*) "Number of species 1, species 2:"
read(*,*) natom1, natom2

! Determine the binsize
binsize = ( (ax**2+ay**2+az**2)**0.5 / 2.0 ) / nbins

! Open necessary files
open(1,file=fxyz,status='old',ERR=90)
open(2,file='TRAJEC.cnn',status='unknown')
open(3,file='nn_dist.hist',status='unknown')
read(1,*) natom
rewind(1)

! Check nwrite parameter
if (nwrite.gt.natom2) then
  write(*,*) "nwrite=",nwrite," is too large, setting nwrite=natom2"
  nwrite = natom2
endif
if (nwrite.gt.100) then
  write(*,*) "nwrite=",nwrite," is still too large, setting nwrite=100"
  nwrite = 100
endif

! Check sum of natom1 + natom2 = natom
if (natom1+natom2.ne.natom) then
  write(*,*) 'natom1 + natom2 != natom --- exiting...'
  goto 110
endif

! Write TRAJEC.cnn header
write(2,*) "# comment = "
write(2,*) "# a = ", ax
write(2,*) "# b = ", ay
write(2,*) "# c = ", az
write(2,*) "# natom = ", natom1
write(2,*) "# nneighbors = ", nwrite
write(2,*) "#"
write(2,*) "#"
write(2,*) "#"
write(2,*) "# units = angstrom"


! Loop over all timesteps in file
do i = 1, maxSteps

  ! Read natom and tstep from xyz file
  read(1,*,END=100) natom
  read(1,*) tstep

  ! Read in trajectories
  do j = 1, natom
    read(1,*) typat(j), X(j), Y(j), Z(j)
  enddo

  ! Compute the neighbour list for each atom
  do j = 1, natom

    ! If current atom is of type 1
    if (typat(j).eq.typat1) then

      do k = 1, natom

        r_nn(k) = 9999.0

        ! If "other" atom is of type 2
        if (typat(k).eq.typat2) then

          ! Calculate the distance between atom j & k
          dx = X(k) - X(j)
          dy = Y(k) - Y(j)
          dz = Z(k) - Z(j)
          input_value(1) = dx/ax
          input_value(2) = dy/ay
          input_value(3) = dz/az
          do l = 1, 3
            pbc_round(l) = int(input_value(l))
            if (abs(input_value(l)-pbc_round(l)).ge.0.5) then
              if (input_value(l).gt.0) pbc_round(l) = pbc_round(l) + 1
              if (input_value(l).lt.0) pbc_round(l) = pbc_round(l) - 1
            endif
          enddo
          dx = dx - ax*pbc_round(1)
          dy = dy - ay*pbc_round(2)
          dz = dz - az*pbc_round(3)
          r = ( dx**2 + dy**2 + dz**2 )**0.5

          ! Distance between j & k is r
          r_nn(k) = r

        endif
      enddo

      ! Sort in terms of ascending distance for ordered neighbour array
      do k = 1, nwrite
        arg = 9999.0
        do l = 1, natom
          if (r_nn(l).lt.arg) then
            arg = r_nn(l)
            index = l
          endif
        enddo
        r_nn(index) = 9999.0
        sorted(k) = arg
        nn(k) = index - 1   ! Shifted to maintain dumb original python indexing
      enddo
    
      ! Re-assign the sorted values to original r_nn array
      do l = 1, nwrite
        r_nn(l) = sorted(l)
      enddo

      ! Write neighbour histograms
      do k = 1, nwrite
        bin = int(r_nn(k)/binsize + 0.5)
        hist(k,bin) = hist(k,bin) + 1
      enddo

      ! Write line to TRAJEC.cnn file
      5 format(A2,x,f12.8,x,f12.8,x,f12.8,100(x,i3))
      write(2,5) typat(j), X(j), Y(j), Z(j), (nn(k),k=1,nwrite)

    endif
  enddo
enddo

100 continue

! Write the neighbour distance histograms
nsteps = i - 1
6 format(f12.8,100(x,f12.8))
do i = 1, nbins
  write(3,6) i*binsize, (hist(k,i)/(2*(natom-1)*nsteps),k=1,nwrite)
!  write(3,6) i*binsize, (hist(k,i)/(natom*nsteps),k=1,nwrite)
enddo

!write(*,*) nsteps, natom
110 continue

close(1)
close(2)
close(3)

90 continue

END
