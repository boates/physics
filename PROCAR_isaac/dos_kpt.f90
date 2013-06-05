PROGRAM dos_kpt 
!
!       program dos_kpt.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer i, j, k, nkpt, neig, nedos 
integer maxline, counter
parameter(maxline=10000000)      ! this just sets a max number of lines that can be read in
real(8) E, emin, emax, dE, weight, eigenvalue, sigma
real(8) dft_sum, gw_sum
real(8), allocatable :: DFT_dos(:,:)
real(8), allocatable :: GW_dos(:,:)
real(8), allocatable :: wtk(:)
real(8), allocatable :: DFT_eig_array(:,:)
real(8), allocatable :: GW_eig_array(:,:)
real, parameter :: pi = 3.14159 ! cannot be changed
character(128) fin, fin2
character(6) dummya, dummyb, dummyc, dummyd, dummye 

nedos = 2101 
sigma = 0.1

write(6,*) 'Name of eigenvalue data file '
!read(5,*) fin
fin='eig.dat'

open(1,file=fin,status='old',ERR=100)

read(1,*) dummya, dummyb, nkpt, dummyc, neig
read(1,*)

write(6,*) 'Smearing is set at: ', sigma

allocate(DFT_dos(3,nedos))
allocate(GW_dos(3,nedos))
allocate(wtk(nkpt))
allocate(DFT_eig_array(nkpt,neig))
allocate(GW_eig_array(nkpt,neig))

counter = 2

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

write(6,*) 'Name of wtk data file '
!read(5,*) fin2
fin2='wtk.dat'

open(2,file=fin2,status='old',ERR=101)

do i=1,nkpt/5   
    read(2,*,END=101) (wtk(j + (i-1)*6),j=1,6)
end do

read(2,*,END=101) (wtk(j + (i-1)*6),j=1, nkpt - (i-1)*6)

101  continue

close(2)

print *, wtk

rewind(1)               ! this puts the read head back at top of file
read(1,*)               ! skips the header
read(1,*)               ! skips the header

DFT_dos(:,:) = 0.0
GW_dos(:,:) = 0.0
DFT_eig_array(:,:) = 0.0
GW_eig_array(:,:) = 0.0

do j=1,nkpt
    do k=1,neig
        read(1,*) DFT_eig_array(j,k), GW_eig_array(j,k)
    end do
end do

!emin = min(minval(DFT_eig_array), minval(GW_eig_array)) - 2.0  ! so we can see +/- 2 eV
!emax = max(maxval(DFT_eig_array), maxval(GW_eig_array)) + 2.0

!emin = -14
!emax = 0.5
emin = -10
emax = 10
dE = (emax - emin)/nedos

do i=1,nedos

    E = (i - 1)*dE + emin  ! should I change this to middle..?
                           ! no, this is the left side of the box for
                           ! integration

    do j=1,nkpt

        weight = 2.0*wtk(j)  ! two e per state

        do k=1,neig

            eigenvalue = DFT_eig_array(j,k)
!            DFT_dos(2,i) = DFT_dos(2,i) + weight*(-(1.0/2.0)*sqrt(1.0*sigma)*sqrt(pi)*erf((eigenvalue - (E + dE))/(sqrt(1.0*sigma))) - -(1.0/2.0)*sqrt(1.0*sigma)*sqrt(pi)*erf((eigenvalue - E)/(sqrt(1.0*sigma))))
            DFT_dos(2,i) = DFT_dos(2,i) + weight*(1.0/2.0)*(erf(sqrt(2.0)*(E + dE -eigenvalue)/(2*sigma)) - erf(sqrt(2.0)*(E - eigenvalue)/(2*sigma)))
            eigenvalue = GW_eig_array(j,k)
!            GW_dos(2,i) = GW_dos(2,i) + weight*(-(1.0/2.0)*sqrt(1.0*sigma)*sqrt(pi)*erf((eigenvalue - (E + dE))/(sqrt(1.0*sigma))) - -(1.0/2.0)*sqrt(1.0*sigma)*sqrt(pi)*erf((eigenvalue - E)/(sqrt(1.0*sigma))))
            GW_dos(2,i) = GW_dos(2,i) + weight*(1.0/2.0)*(erf(sqrt(2.0)*(E + dE -eigenvalue)/(2*sigma)) - erf(sqrt(2.0)*(E - eigenvalue)/(2*sigma)))

        end do 

    end do

end do

open(3,file="dft_dos.dat")

dft_sum = 0.0

do i=1,nedos
  E = (i-1)*dE + emin + dE/2
  dft_sum = dft_sum + DFT_dos(2,i)
  write(3,*) E, DFT_dos(2,i), dft_sum
end do

close(3)

open(3,file="gw_dos.dat")

gw_sum = 0.0

do i=1,nedos
  E = (i-1)*dE + emin + dE/2
  gw_sum = gw_sum + GW_dos(2,i)
  write(3,*) E, GW_dos(2,i), gw_sum
end do

close(3)

deallocate(wtk)
deallocate(DFT_eig_array)
deallocate(GW_eig_array)
deallocate(DFT_dos)
deallocate(GW_dos)

END PROGRAM dos_kpt 
