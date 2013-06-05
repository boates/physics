!       program wigner-kirkwood.f
!***********************************************************
!       Calculate the Wigner-Kirkwood expansion for
!       the 1st order quantum correction to the ionic
!       free energies.
!
!       FORCES.fxyz must be in eV/Angstrom, T in K
!       will get the free energy correction at each
!       step in eV/atom
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=1000000,maxAtoms=1024)
        real*8 hbar, kB, mass, units_conversion, factor, temp
        parameter (hbar=6.58211899E-16,kB=8.617343E-05,mass=14.007)
        parameter (units_conversion=9.6485341E+27)
        integer natom, tstep
        real*8 T(maxSteps),Fx(maxAtoms),Fy(maxAtoms),Fz(maxAtoms)
        real*8 sum, delta_F
        character*128 fxyz, ftemp

        ! Format of FORCES.fxyz is 'Fx Fy Fz' in eV/Ang
        write(*,*) 'FORCES.fxyz file:'
        read(*,*) fxyz

        ! Format of temperature.dat is 'tstep T'
        write(*,*) 'temperature.dat file:'
        read(*,*) ftemp

        write(*,*) 'Number of atoms'
        read(*,*) natom
        if (natom.gt.maxAtoms) then
          goto 105
        endif

        ! Open the input and output files
        open(1,file=fxyz,status='old',ERR=90)
        open(2,file=ftemp,status='old',ERR=90)

        ! Read the header of the temperature.dat file
        read(2,*)

        open(11,file='WK_expansion.dat')
        write(11,*) '# tstep delta_F (eV/atom)'

        factor = (hbar**2 / (24*kB**2)) * units_conversion

        do i=1,maxSteps

          read(2,*,END=100) tstep, temp
          sum = 0.0

          do j=1,natom

            read(1,*) Fx(j), Fy(j), Fz(j)

            sum = sum + (Fx(j)**2 + Fy(j)**2 + Fz(j)**2)/3.0 /mass

          enddo
          
          delta_F = sum*(factor/temp**2)

          ! Write 1st order correction in eV/atom
          write(11,*) tstep, delta_F / natom

        enddo

 100    continue

        close(1)
        close(2)
        close(11)

        goto 90

 105    continue

        write(*,*) 'natom = ',natom,' exceeds maxAtoms = ',maxAtoms

 90     continue

        END

