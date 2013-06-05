!Calculates the universal equation to solve the fluidicity factor f
implicit none
integer N
real*8 f,D,x,k,amu2Kg,T,m,s0,L,pi,V,a1,a2,a3,a4,freq
parameter(k=1.3806503E-23,amu2Kg=1.66053886E-27) !,m=6.941d0*amu2Kg)  

OPEN (9,FILE ='dos.dat',status='old')
OPEN (10,FILE ='latt_const_P',status='old')

read(10,*) L
close(10)
L=L*1.0E-10
V=L**3

write(*,*) "enter the # of atoms"
read(*,*) N
write(*,*) "enter temp. in K"
read(*,*) T
write(*,*) "mass in amu"
read(*,*) m

read(9,*) freq, s0
close(9)

pi = 2.0*ASIN(1.0D0)

s0=s0*1.0E-12 !changes 1/THz(1/10^12*1/s) into s  
m=m*amu2Kg

a1=2.0d0*s0/(9.0d0*N)
a2=sqrt(pi*k*T/m)
a3=(N/V)**(1.0d0/3.0d0)
a4=(6.0d0/pi)**(2.0d0/3.0d0)
D=a1*a2*a3*a4

!write(*,*) "Delta:", D

do f = 0.0d0,1.0d0,0.000001d0 
 x=2.0d0*f**(7.5)*D**(-4.5)-6.0d0*f**5*D**(-3)-f**(3.5)*D**(-1.5)+6.0d0*f**(2.5)*D**(-1.5)+2.0d0*f-2.0d0
if (x>-0.00001 .and. x<0.00001) then
write(*,*) "========================================="
write(*,*) "Fluidicity factor:", f
write(*,*) "========================================="
endif
enddo

END
