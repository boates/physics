! gibbs-solid.f90
! Brian's modified version of Amanuel's code

implicit none

integer i, j, npts, natom
real*8 k_B,h,T,P,dos,w_A,w_S,w_E,x,v,v1,v2,hv,kBT,n,Lx,Ly,Lz,E
real*8 dosXw_A(:),dosXw_S(:),dosXw_E(:)
real*8 int_dosXw_A,int_dosXw_S,int_dosXw_E
real*8 vol,q,enth,G,sum_1,sum_2,sum_3,D
real*8 alat,ax,ay,az,bx,by,bz,cx,cy,cz
character*12 fdos, fpos
parameter (k_B=1.3806503E-23,h=6.626068E-34,q=1.60217646E-19)  
allocatable :: dosXw_A,dosXw_S,dosXw_E

! Retrieve user input
write(*,*) "name of VDOS.dat file in THZ normalized to 3N:"
read(*,*) fdos
write(*,*) "name of POSCAR file for volume calculation:"
read(*,*) fpos
write(*,*) "temperature in K:"
read(*,*) T

! Read lattice vector info and natom from POSCAR file
open(20,file=fpos,status='old')
read(20,*)
read(20,*) alat
read(20,*) ax, ay, az
read(20,*) bx, by, bz
read(20,*) cx, cy, cz
read(20,*) natom
close(20)

! Calculate the volume as the scalar triple product
vol = (ax*(by*cz-bz*cy) + ay*(bz*cx-bx*cz) + az*(bx*cy-by*cx))*alat**3

! Determine size of VDOS file
open(10,file=fdos,status='old')
npts=0
100 read(10,*,END=200) 
npts = npts + 1
goto 100
200 continue
REWIND (10)

! Allocate arrays accordingly
allocate(dosXw_A(npts),dosXw_S(npts),dosXw_E(npts))

! Do the calculation
kBT = k_B*T
do i = 1, npts
  read(10,*)  v, dos
  if (i==1) v1=v
  if (i==2) v2=v  
  hv = h*v*1.0E+12
  x= hv/KBT
  w_A = -dlog(dexp(-x/2.0d0)/(1.0d0-dexp(-x)))
  W_S = x/(dexp(x)-1.0d0)-log(1.0d0-dexp(-x))
  w_E = x/2.0d0 + x/(dexp(x)-1.0d0)
  if (v==0.0d0) then
    dosXw_A(i) = 0.0d0
    dosXw_S(i) = 0.0d0
    dosXw_E(i) = 0.0d0
  else
   dosXw_A(i) = dos*w_A
   dosXw_S(i) = dos*w_S
   dosXw_E(i) = dos*w_E
  endif
enddo 
close(10)

sum_1=dosXw_A(1)+dosXw_A(npts)
sum_2=dosXw_S(1)+dosXw_S(npts)
sum_3=dosXw_E(1)+dosXw_E(npts)
D=4.0d0
do i = 2,npts-1
   sum_1 = sum_1 + D*dosXw_A(i)
   sum_2 = sum_2 + D*dosXw_S(i)
   sum_3 = sum_3 + D*dosXw_E(i)
   D=6.0d0-D
enddo

int_dosXw_A=sum_1*(v2-v1)/3.0d0
int_dosXw_S=sum_2*(v2-v1)/3.0d0
int_dosXw_E=sum_3*(v2-v1)/3.0d0

kBT = kBT/q

! Print the results to screen
write(*,*) "Entropy per atom in k_B:", int_dosXw_S/dfloat(natom)
write(*,*) "Entropy per atom in eV/K:", int_dosXw_S/dfloat(natom)*8.617343d-05
write(*,*) "TS per atom in eV:", T*int_dosXw_S/dfloat(natom)*8.617343d-05

!open(30,file='entropy.dat')
!5 format(f12.10)
!write(30,5) int_dosXw_S/dfloat(natom)*8.617343d-05
!close(30)

END
