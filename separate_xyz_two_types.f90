! program separate_xyz_two_types.f90
!****************************************************
! Separate xyz file into xyz files for type1 &
! type2 atoms
!****************************************************

implicit none

integer maxSteps, i, j
parameter (maxSteps=100000)
real*8 x, y, z
integer tstep, natom, natom1, natom2
character*2 typat, type1, type2
character*128 fxyz

! Retrieve user input
write(*,*) "TRAJEC.xyz:"
read(*,*) fxyz
write(*,*) "natom1, natom2:"
read(*,*) natom1, natom2
write(*,*) "typat1, typat2:"
read(*,*) type1, type2

! Open necessary files
open(1,file=fxyz,status='old',ERR=90)
open(2,file='TRAJEC_type1.xyz',status='unknown')
open(3,file='TRAJEC_type2.xyz',status='unknown')

! Perform the calculation
do i= 1, maxSteps

  read(1,*,END=100)  natom
  read(1,*)  tstep
  write(2,*) natom1
  write(2,*) tstep
  write(3,*) natom2
  write(3,*) tstep

  do j = 1, natom
    read(1,*) typat, x, y, z
    if (typat.eq.type1) write(2,*) typat, x, y, z
    if (typat.eq.type2) write(3,*) typat, x, y, z
  enddo

enddo

100 continue

close(1)
close(2)
close(3)

90 continue

END
