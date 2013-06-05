!       program VDOS.f90
!***********************************************************
!       Calculate normalized VDOS from normalized VACF
!***********************************************************

implicit none

integer i, j, Nsteps, maxArray
parameter (maxArray=20000)
real*8 pi, f, t, X, dt, dummy
real*8 sum,vacf(maxArray),cosXvacf(maxArray,maxArray)
character*128 fin

! Define pi
pi = 2.0*asin(1.0)

! Retrieve user input
write(*,*) "VACF.dat file from VACF code"
read(*,*) fin
write(*,*) "Enter timestep in a.u."
read(*,*) dt

! Open necessary files
open(1,file=fin,status='old')
open(2,file="VDOS.dat",status='unknown')

! Get the number of steps
Nsteps = 0
100 read(1,*,END=200) 
Nsteps = Nsteps + 1
goto 100
200 continue
close(1)
open(1,file=fin,status='old')

! Read vacf from file
do i=1,Nsteps
  read(1,*) t, dummy, vacf(i)
end do
close(1)

! Convert timestep to fs
dt = dt * 2.418800d-5 * 1500.0
dt = dt / 1500.0

! Calculate the cosXvacf array
do i=0,1500
  f = i * 0.1
  do j = 1,Nsteps
    cosXvacf(i,j) = dcos(2.0*pi*f*(j-1)*dt)*vacf(j)
  enddo
enddo

! Calculate VDOS and write to file
do i=0,1500
  sum = cosXvacf(i,1) + cosXvacf(i,Nsteps/2)
  X = 4.0
  do j=2,Nsteps/2-1
    sum = sum + X*cosXvacf(i,j)
    X = 6.0 - X
  enddo
  write(2,*) i*0.1, 4.0*sum*dt!/3.0
enddo

close(2)

END
