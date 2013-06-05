! program unwrap_PBC_general.f90
!******************************************************************
! Author: Brian Boates
!
! Unwrap xyz file, write to unwrapped.xyz
!
! NOTE: This algorithm assumes a reasonably small stepsize
!       with respect to the lattice constants. Is partially
!       taken from F.26 from Allen-Tildesley.
!******************************************************************
implicit none

! n = max # of timesteps, m = max # of atoms
integer n, m, i, j, k
parameter (n=50000,m=512)
integer natoms(n), timesteps(n)
real*4 Rx(n,m), Ry(n,m), Rz(n,m)
real*4 newRx(n,m), newRy(n,m), newRz(n,m)
real*4 dx(n,m), dy(n,m), dz(n,m), dRx, dRy, dRz
real*4 ax, ay, az, bx, by, bz, cx, cy, cz, amag, bmag, cmag
real*4 iax, iay, iaz, ibx, iby, ibz, icx, icy, icz, det
character*2 typat(n,m)
character*128 fxyz

! Ask user for required input, which can be given via
! an input file (still echoes the requests though...)
write(*,*) "Name of xyz file with wrapped coordinates"
read(*,*) fxyz
write(*,*) "a lattice vector components (ax,ay,az)"
read(*,*) ax, ay, az
write(*,*) "b lattice vector components (bx,by,bz)"
read(*,*) bx, by, bz
write(*,*) "c lattice vector components (cx,cy,cz)"
read(*,*) cx, cy, cz

! open input and output files
open(1,file=fxyz,status='old',ERR=90)
open(2,file='unwrapped.xyz')

! Calculate vector magnitudes
amag = (ax**2 + ay**2 + az**2)**0.5
bmag = (bx**2 + by**2 + bz**2)**0.5
cmag = (cx**2 + cy**2 + cz**2)**0.5

! Create inverse lattice vectors for coord transformation
det = ax*(cz*by-cy*bz) - bx*(cz*ay-cy*az) + cx*(bz*ay-by*az)
iax = (cz*by - cy*bz) / det
iay = (cy*az - cz*ay) / det
iaz = (bz*ay - by*az) / det
ibx = (cx*bz - cz*bx) / det
iby = (cz*ax - cx*az) / det
ibz = (bx*az - bz*ax) / det
icx = (cy*bx - cx*by) / det
icy = (cx*ay - cy*ax) / det
icz = (by*ax - bx*ay) / det

!write(*,*) ax, ay, az
!write(*,*) bx, by, bz
!write(*,*) cx, cy, cz
!write(*,*) amag, bmag, cmag
!write(*,*) det
!write(*,*) iax, iay, iaz
!write(*,*) ibx, iby, ibz
!write(*,*) icx, icy, icz

! Read xyz file, convert to reduced coordinates & shift origin
do i = 1, n

  read(1,*,END=100) natoms(i)
  read(1,*) timesteps(i)

  do j=1,natoms(i)

    read(1,*) typat(i,j), Rx(i,j), Ry(i,j), Rz(i,j)

!    if (i.eq.100.and.j.eq.1) then
!       write(*,*)
!      write(*,*) Rx(i,j), Ry(i,j), Rz(i,j)
!    endif

!    newRx(i,j) = (iax*Rx(i,j) + iay*Ry(i,j) + iaz*Rz(i,j)) !+ 0.5
!    newRy(i,j) = (ibx*Rx(i,j) + iby*Ry(i,j) + ibz*Rz(i,j)) !+ 0.5
!    newRz(i,j) = (icx*Rx(i,j) + icy*Ry(i,j) + icz*Rz(i,j)) !+ 0.5

!    Rx(i,j) = newRx(i,j)
!    Ry(i,j) = newRy(i,j)
!    Rz(i,j) = newRz(i,j)

    Rx(i,j) = Rx(i,j) !+ 0.5
    Ry(i,j) = Ry(i,j) !+ 0.5
    Rz(i,j) = Rz(i,j) !+ 0.5

    if (i.eq.1) then
      write(*,*) Rx(i,j), Ry(i,j), Rz(i,j)
    endif
    if (Rx(i,j)-0.5.gt.1.0.or.Ry(i,j)-0.5.gt.1.0.or.Rz(i,j)-0.5.gt.1.0) then
      write(*,*) 'WARNING:', Rx(i,j), Ry(i,j), Rz(i,j)
    endif

  enddo
enddo

100    continue

! Initialize dx, dy, & dz arrays
do j = 1, natoms(1)

  dx(1,j) = 0.0
  dy(1,j) = 0.0
  dz(1,j) = 0.0

enddo

! Find all dx, dy, & dz values from the wrapped data
do i = 2, n

  do j = 1, natoms(i)

    dRx = Rx(i,j) - Rx(i-1,j)
    dRy = Ry(i,j) - Ry(i-1,j)
    dRz = Rz(i,j) - Rz(i-1,j)

    dx(i,j) = dx(i-1,j) + (dRx - anint(dRx))
    dy(i,j) = dy(i-1,j) + (dRy - anint(dRy))
    dz(i,j) = dz(i-1,j) + (dRz - anint(dRz))

  enddo
enddo

! Unwrap the data
do i = 1, n

  if (natoms(i).eq.0) then
    goto 80
  endif

  do j=1,natoms(i)

    ! Calculate the unwrapped coordinate
    Rx(i,j) = Rx(1,j) + dx(i,j)
    Ry(i,j) = Ry(1,j) + dy(i,j)
    Rz(i,j) = Rz(1,j) + dz(i,j)

  enddo
enddo

80     continue

! Write to file
do i = 1, n

  if (natoms(i).eq.0) then
    goto 70
  endif

  write(2,*) natoms(i)
  write(2,*) timesteps(i)

  do j = 1, natoms(i)

!    if (i.eq.100.and.j.eq.1) then
!      write(*,*) Rx(i,j), Ry(i,j), Rz(i,j)
!    endif

    ! Convert back to Cartesian coordinates (re-shift origin)
!    newRx(i,j) = ax*(Rx(i,j)-0.5) + ay*(Ry(i,j)-0.5) + az*(Rz(i,j)-0.5)
!    newRy(i,j) = bx*(Rx(i,j)-0.5) + by*(Ry(i,j)-0.5) + bz*(Rz(i,j)-0.5)
!    newRz(i,j) = cx*(Rx(i,j)-0.5) + cy*(Ry(i,j)-0.5) + cz*(Rz(i,j)-0.5)

    newRx(i,j) = ax*Rx(i,j) + ay*Ry(i,j) + az*Rz(i,j)
    newRy(i,j) = bx*Rx(i,j) + by*Ry(i,j) + bz*Rz(i,j)
    newRz(i,j) = cx*Rx(i,j) + cy*Ry(i,j) + cz*Rz(i,j)

    Rx(i,j) = newRx(i,j)
    Ry(i,j) = newRy(i,j)
    Rz(i,j) = newRz(i,j)

!    if (i.eq.100.and.j.eq.1) then
!      write(*,*) Rx(i,j), Ry(i,j), Rz(i,j)
!    endif

    ! Write unwrapped coordinates to file
    write(2,*) typat(i,j), Rx(i,j), Ry(i,j), Rz(i,j)

  enddo
enddo

70     continue

close(1)
close(2)

90     continue

END
