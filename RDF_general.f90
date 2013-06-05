! program RDF_general.f90
!****************************************************
! Caclulation of pair correlation function g(r) for
! a general shaped cell.
!
! To calculate the distances in the reduced zone
! scheme, I am creating a 3x3x3 supercell and
! calculating all of those distances and using 
! the smallest. Yes, this is inefficient but for
! a general cell the algorithm would be difficult
! to write. I don't see anyone else writing this
! code so don't complain :P after all it is in
! fortran so it is still fast.
!****************************************************

implicit none

integer maxSteps, maxAtoms, i, j, k, l, m, n, c
parameter (maxSteps=20000,maxAtoms=1024)
integer natom, tstep, bin
real*8 ax, ay, az, bx, by, bz, cx, cy, cz, alat, r
real*8 x(maxAtoms), y(maxAtoms), z(maxAtoms)
real*8 rX(3*3*3), rY(3*3*3), rZ(3*3*3)
real*8 dX(3*3*3), dY(3*3*3), dZ(3*3*3), dR(3*3*3)
real*8 rmin, rmax, bin_size, pi, V, rho
real*8 rdf(10000), denominator
character*2 typat(maxAtoms)
character*128 xyz, pos

! Retrieve user input
write(*,*) "xyz file:"
read(*,*) xyz
write(*,*) "POSCAR file:"
read(*,*) pos
write(*,*) "Max r:"
read(*,*) rmax
write(*,*) "Bin size:"
read(*,*) bin_size

! Open necessary files
open(1,file=xyz,status='old',ERR=90)
open(3,file='RDF.dat',status='unknown')
open(2,file=pos,status='old',ERR=95)

! Read POSCAR header
read(2,*) 
read(2,*) alat
read(2,*) ax, ay, az
read(2,*) bx, by, bz
read(2,*) cx, cy, cz
read(2,*) natom
close(2)

! Scale the lattice vector components
ax = ax*alat
ay = ay*alat
az = az*alat
bx = bx*alat
by = by*alat
bz = bz*alat
cx = cx*alat
cy = cy*alat
cz = cz*alat

! Define necessary parameters
pi = 3.1415926535897931
V = ax*(by*cz-bz*cy) + ay*(bz*cx-bx*cz) + az*(bx*cy-by*cx)
rho = natom / V

! Perform the calculation
do i = 1, maxSteps

  ! Read xyz info lines
  read(1,*,END=100) natom
  read(1,*) tstep

  ! Read in atomic coordinates for current timestep
  do j = 1, natom
    read(1,*,END=100) typat(j), x(j), y(j), z(j)
  enddo

  do j = 1, natom
    do k = j, natom

      ! Replicate the k atom in 3x3x3 supercell
      c = 0
      rmin = rmax
      do l = -1, 1
        do m = -1, 1
          do n = -1, 1
            c = c + 1
            rX(c) = x(k) + l*ax + m*bx + n*cx
            rY(c) = y(k) + l*ay + m*by + n*cy
            rZ(c) = z(k) + l*az + m*bz + n*cz

            ! Calculate distances and find the smallest
            dX(c) = x(j) - rX(c)
            dY(c) = y(j) - rY(c)
            dZ(c) = z(j) - rZ(c)
            dR(c) = (dX(c)**2 + dY(c)**2 + dZ(c)**2)**0.5
            if (dR(c).lt.rmin) rmin = dR(c)

          enddo
        enddo
      enddo

      ! Add to the g(r) histogram
      bin = int( rmin/bin_size )

      rdf(bin) = rdf(bin) + 1.0

    enddo
  enddo

enddo

100 continue

tstep = i - 1

! Scale g(r) and write to file
5 format(x,f12.4,x,f12.4)
do i = 1, int(rmax/bin_size) - 1

  r = i*bin_size + bin_size/2.0

  write(3,5) r, rdf(i) / (natom-1) / rho / (4*pi*r**2*bin_size) / tstep *2

enddo

95 continue

close(1)
close(3)

90 continue

END
