implicit none

integer i,j, npts,NE,natoms,natoms1
real*8 k_B,h,T,P,dos,w_A_s,w_S_s,w_A_g,w_S_g,w_E,x,v,v1,v2,hv,kBT,n,L,E,E_tot,E_ave
real*8 dosXw_A_s(:),dosXw_S_s(:),int_dosXw_A_s,int_dosXw_A_g,vol,q,enth,G,sum_1,sum_2,D
real*8 xx,pi,f,s0,dosXw_A_g(:),dosXw_S_g(:),dos_g(:),dos_s,int_dosXw_S_s,int_dosXw_S_g
real*8 a1,a2,a3,a4,Del,m,amu2Kg,y,fy,z,lambda,cut_off,dummy
character*128 f_vdos
parameter (k_B=1.3806503E-23,h=6.626068E-34,q=1.602E-19,amu2Kg=1.66053886E-27)  !k_B in J/K& h in J.s 
allocatable :: dosXw_A_s,dosXw_S_s,dosXw_A_g,dosXw_S_g,dos_g

write(*,*) 'name of VDOS.dat file in THz normalized to 3N:'
read(*,*) f_vdos
write(*,*) 'number of atoms:'
read(*,*) natoms
write(*,*) 'number of atoms of this species:'
read(*,*) natoms1
write(*,*) 'volume in angstroms^3:'
read(*,*) vol
write(*,*) 'temperature in K:'
read(*,*) T
write(*,*) 'mass in amu:'
read(*,*) m
write(*,*) 'fluidicity factor:'
read(*,*) f

open(10,file=f_vdos,status='old')
open(11,FILE ='dos_s.dat',status='unknown')
open(12,FILE ='dos_g.dat',status='unknown')

npts=0
100 read(10,*,END=200) cut_off,dummy
npts = npts + 1
goto 100
200 continue
REWIND (10)

allocate(dosXw_A_s(npts),dosXw_S_s(npts),dos_g(npts))

kBT = k_B*T 
pi = 2.0*DASIN(1.0D0)

do i = 1, npts
  read(10,*) v,dos
  if (v>cut_off) dos=0.0d0
  if (i==1) then 
   v1=v
   s0=dos
  endif
  if (i==2) v2=v
  hv = h*v*1.0E+12
  x= hv/KBT
  w_A_s = -dlog(dexp(-x/2.0d0)/(1.0d0-dexp(-x)))
  w_S_s = x/(dexp(x)-1.0d0)-log(1.0d0-dexp(-x))
  xx=pi*s0*v/(6.0d0*f*dfloat(natoms1))
  dos_g(i) = s0/(1.0d0+xx**2)
  dos_s = dos-dos_g(i)
  if (v==0.0d0) then
    dosXw_A_s(i) = 0.0d0
    dosXw_S_s(i) = 0.0d0
  else
   dosXw_A_s(i) = dos_s*w_A_s
   dosXw_S_s(i) = dos_s*w_S_s
  endif
write(11,*) v,dos_s
write(12,*) v,dos_g(i)
enddo
close(10)

!write(*,*) s0

!integrate with simpson's rule
sum_1=dosXw_A_s(1)+dosXw_A_s(npts)
sum_2=dosXw_S_s(1)+dosXw_S_s(npts)
D=4.0d0
do i = 2,npts-1
   sum_1 = sum_1 + D*dosXw_A_s(i)
   sum_2 = sum_2 + D*dosXw_S_s(i)
   D=6.0d0-D
enddo

int_dosXw_A_s=sum_1*(v2-v1)/3.0d0
int_dosXw_S_s=sum_2*(v2-v1)/3.0d0

vol=vol*(1.0E-10)**3
m=m*amu2Kg
lambda= h/(dsqrt(2.0d0*pi*m*k_B*T))

s0=s0*1.0E-12

a1=2.0d0*s0/(9.0d0*natoms1)
a2=sqrt(pi*k_B*T/m)
a3=(natoms1/vol)**(1.0d0/3.0d0)
a4=(6.0d0/pi)**(2.0d0/3.0d0)
Del=a1*a2*a3*a4
y=(f/Del)**(1.5d0)
fy=y*f
z=(1.0d0+fy+fy**2-fy**3)/((1.0d0-fy)**3)
w_S_g=2.5d0-log(natoms1/vol*lambda**3)+log(z)+fy*(3.0d0*fy-4.0d0)/((1.0d0-fy)**2)
w_S_g=w_S_g/3.0d0
w_A_g=0.5d0-w_S_g

allocate(dosXw_A_g(npts),dosXw_S_g(npts))

do i=1,npts
  dosXw_A_g(i) = dos_g(i)*w_A_g
  dosXw_S_g(i) = dos_g(i)*w_S_g
enddo

!integrate with simpson's rule
sum_1=dosXw_A_g(1)+dosXw_A_g(npts)
sum_2=dosXw_S_g(1)+dosXw_S_g(npts)
D=4.0d0
do i = 2,npts-1
   sum_1 = sum_1 + D*dosXw_A_g(i)
   sum_2 = sum_2 + D*dosXw_S_g(i)
   D=6.0d0-D
enddo

int_dosXw_A_g=sum_1*(v2-v1)/3.0d0
int_dosXw_S_g=sum_2*(v2-v1)/3.0d0

kBT = kBT/q !in eV

!vol = L**3

Write(*,*) "Entropy per atom in k_B:", (int_dosXw_S_s+int_dosXw_S_g)/dfloat(natoms)
Write(*,*) "Entropy per atom in eV/K:", (int_dosXw_S_s+int_dosXw_S_g)/dfloat(natoms)*8.6173423E-5
Write(*,*) "Entropy per atom in k_B: solid part", int_dosXw_S_s/dfloat(natoms)
Write(*,*) "Entropy per atom in k_B: gas part", int_dosXw_S_g/dfloat(natoms)
Write(*,*) "TS per atom in eV:", kBT*(int_dosXw_S_s+int_dosXw_S_g)/dfloat(natoms)

END

