! program molecules.f90
!****************************************************
! Scan a xyz file and find which molecules exist
!****************************************************

implicit none

integer maxSteps, maxAtoms, i, j, k
parameter (maxSteps=100000,maxAtoms=1024)
integer natom, tstep
real*8 x(maxAtoms), y(maxAtoms), z(maxAtoms)
real*8 dx, dy, dz
character*2 typat(maxAtoms)
character*128 fxyz

! Retrieve user input
write(*,*) "Name of xyz File:"
read(*,*) fxyz

! Open necessary files
open(1,file=fxyz,status='old',ERR=90)
open(2,file='output.dat',status='unknown')


! Perform the calculation
!do i = 1, maxSteps

  ! Read in xyz snapshot
!  read(1,*,END=100) natom
!  read(1,*) tstep
!  do j = 1, natom
!    read(1,*) typat(j), x(j), y(j), z(j)
!  enddo

  

!enddo

100 continue

close(1)
close(2)

90 continue

END
