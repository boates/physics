implicit none

integer i,j, npts,NE,natoms,natoms_mg,natoms_si,natoms_o
real*8 m_mg,m_si,m_o,a1_mg,a2_mg,a1_si,a2_si,a1_o,a2_o,factrl_mg,factrl_si,factrl_o
real*8 k_B,h,T,P,dos,w_A_s,w_S_s,w_A_g,w_S_g,w_E,x,v,v1,v2,hv,kBT,n,L,E,E_tot,E_ave,ds
real*8 dosXw_A_s(:),dosXw_S_s(:),int_dosXw_A_s,int_dosXw_A_g,vol,q,enth,G,sum_1,sum_2,D,df1,df2
real*8 xx,pi,f,s0,dosXw_A_g(:),dosXw_S_g(:),dos_g(:),dos_s,int_dosXw_S_s,int_dosXw_S_g,F_id3
real*8 a1,a2,a3,a4,Del,m,amu2Kg,y,fy,z,lambda,cut_off,F_id1,F_id2,F_id,dF,k_B2,factrl,U,U_withTS_e
!parameter (k_B=1.3806503E-23,h=6.626068E-34,q=1.602E-19,amu2Kg=1.66053886E-27)  !k_B in J/K& h in J.s 
parameter (k_B=1.3806503d0,h=6.626068d0,amu2Kg=1.66053886d0,k_B2=8.617385E-5)  !k_B in J/K& h in J.s 
allocatable :: dosXw_A_s,dosXw_S_s,dosXw_A_g,dosXw_S_g,dos_g

write(*,*) "enter lattice const. in angst."
read(*,*) L
write(*,*) "enter temp. in K"
read(*,*) T
write(*,*) "enter dF (id->DFT/cl)"
read(*,*) dF 
write(*,*) "enter ^=1 potential energy:KE will be added" 
read(*,*) U

U=U+1.5d0*k_B2*T

vol=L**3

pi = 2.0*DASIN(1.0D0)

m_mg=24.305
m_si=28.086d0
m_o=15.99d0
natoms_mg=27
natoms_si=27
natoms_o=81

factrl_mg=1.00d0
factrl_si=1.00d0
factrl_o=1.00d0

do i = 1, natoms_mg
 factrl_mg=factrl_mg*dfloat(i)
enddo 
do i = 1, natoms_si
 factrl_si=factrl_si*dfloat(i)
enddo 
do i = 1, natoms_o
 factrl_o=factrl_o*dfloat(i)
enddo 


!Mg
a1_mg=0.001d0*vol/dfloat(natoms_mg)
a2_mg=2.0d0*pi*M_mg*amu2Kg*k_B*T/(h**2)
a2_mg=a2_mg**(3.0d0/2.0d0)
!Si
a1_si=0.001d0*vol/dfloat(natoms_si)
a2_si=2.0d0*pi*M_si*amu2Kg*k_B*T/(h**2)
a2_si=a2_si**(3.0d0/2.0d0)
!O
a1_o=0.001d0*vol/dfloat(natoms_O)
a2_o=2.0d0*pi*M_O*amu2Kg*k_B*T/(h**2)
a2_o=a2_o**(3.0d0/2.0d0)

a1=0.001d0*vol

F_id1=-k_B2*T*(dfloat(natoms_mg)*log(a1*a2_mg)-log(factrl_mg)) 
F_id2=-k_B2*T*(dfloat(natoms_si)*log(a1*a2_si)-log(factrl_si))
F_id3=-k_B2*T*(dfloat(natoms_o)*log(a1*a2_o)-log(factrl_o))
F_id=(F_id1+F_id2+F_id3)/dfloat(natoms_mg+natoms_si+natoms_o) !per atom (exact)


write(*,*) "F_id"
write(*,*) F_id

write(*,*) "Entropy in k_B/atom"
write(*,*) (U-F_id-dF)/(T*k_B2) 



END
