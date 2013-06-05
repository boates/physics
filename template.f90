! program template.f90
!****************************************************
! Describe code
!****************************************************

implicit none

integer maxSteps, maxAtoms, natom, i, j
parameter (maxSteps=100000,maxAtoms=1024)
real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
real*8 dx, dy, dz
character*2 typat, d1, d2, d3, d4
character*128 fin

! Retrieve user input
write(*,*) "Filename:"
read(*,*) fin

! Open necessary files
open(1,file=fin,status='old',ERR=90)
open(2,file='output.dat',status='unknown')

! Perform the calculation
do i=1,maxSteps

  do j = 1,natom

    read(1,*,END=100) X(j), Y(j), Z(j)

    write(2,*) X(j), Y(j), Z(j)

  enddo

enddo

100 continue

close(1)
close(2)

90 continue

END
