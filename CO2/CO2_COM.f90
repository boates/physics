!       program CO2_COM.f90
!***********************************************************
!       Create an xyz file of CO2 molecular COM's 
!       from TRAJEC.mol
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, nC, nO, nmol, nn1(maxAtoms), nn2(maxAtoms)
        integer pbc_round(3)
        real*8 input_value(3)
        real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
        real*8 Xcom(maxAtoms), Ycom(maxAtoms), Zcom(maxAtoms)
        real*8 ax, ay, az, dx, dy, dz, r, d
        real*8 Cmass, Omass, CO2mass, A(3), B(3)
        character*2 typat, d1, d2, d3
        character*128 fmol

        ! Retrieve user input
        write(*,*) 'Name of .mol file:'
        read(*,*) fmol

        ! Open input and output files
        open(1,file=fmol,status='old',ERR=90)
        open(2,file='CO2_COM.xyz')

        ! Read cbn header info
        read(1,*,END=100) d1, d2, d3, natom
        read(1,*) d1, d2, d3, nC
        read(1,*) d1, d2, d3, nO
        read(1,*) 
        read(1,*) d1, d2, d3, ax
        read(1,*) d1, d2, d3, ay
        read(1,*) d1, d2, d3, az
        read(1,*) 
        read(1,*) 
        read(1,*) 

        ! Define halfbox distance and atomic masses (g/mol or amu)
        d = ( (ax/2)**2 + (ay/2)**2 + (az/2)**2 )**0.5 / 2.0
        Cmass = 12.0107
        Omass = 15.9994
        CO2mass = Cmass + 2*Omass

        ! Begin calculation
        do i=1,maxSteps

          ! Read in from .mol file
          do j=1,natom
            read(1,*,END=100) typat, X(j), Y(j), Z(j), nn1(j), nn2(j)
          enddo

          ! Initiate molecule counting for current timestep
          nmol = 0

          ! Loop over carbons (possible molecules)
          do j=1,nC

            ! If this carbon is part of a molecule (as defined)
            if ((nn1(j).ne.0).and.(nn2(j).ne.0)) then

              ! Add one molecule to the count
              nmol = nmol + 1

              ! Calculate COM coordinates for molecule
              dx = X(j) - X(nn1(j))
              dy = Y(j) - Y(nn1(j))
              dz = Z(j) - Z(nn1(j))
              input_value(1) = dx/ax
              input_value(2) = dy/ay
              input_value(3) = dz/az
              pbc_round(1) = int(input_value(1))
              pbc_round(2) = int(input_value(2))
              pbc_round(3) = int(input_value(3))
              if (abs(input_value(1)-pbc_round(1)).ge.0.5) then
                if (input_value(1).gt.0) pbc_round(1) = pbc_round(1) + 1
                if (input_value(1).lt.0) pbc_round(1) = pbc_round(1) - 1
              endif
              if (abs(input_value(2)-pbc_round(2)).ge.0.5) then
                if (input_value(2).gt.0) pbc_round(2) = pbc_round(2) + 1
                if (input_value(2).lt.0) pbc_round(2) = pbc_round(2) - 1
              endif
              if (abs(input_value(3)-pbc_round(3)).ge.0.5) then
                if (input_value(3).gt.0) pbc_round(3) = pbc_round(3) + 1
                if (input_value(3).lt.0) pbc_round(3) = pbc_round(3) - 1
              endif
              A(1) = dx - ax*pbc_round(1)
              A(2) = dy - ay*pbc_round(2)
              A(3) = dz - az*pbc_round(3)
              r = ( A(1)**2 + A(2)**2 + A(3)**2 )**0.5
              if (r.gt.d) then
                write(*,*) "WARNING --- PROBLEM WITH PBC"
                goto 100
              endif
              
              dx = X(j) - X(nn2(j))
              dy = Y(j) - Y(nn2(j))
              dz = Z(j) - Z(nn2(j))
              input_value(1) = dx/ax
              input_value(2) = dy/ay
              input_value(3) = dz/az
              pbc_round(1) = int(input_value(1))
              pbc_round(2) = int(input_value(2))
              pbc_round(3) = int(input_value(3))
              if (abs(input_value(1)-pbc_round(1)).ge.0.5) then
                if (input_value(1).gt.0) pbc_round(1) = pbc_round(1) + 1
                if (input_value(1).lt.0) pbc_round(1) = pbc_round(1) - 1
              endif
              if (abs(input_value(2)-pbc_round(2)).ge.0.5) then
                if (input_value(2).gt.0) pbc_round(2) = pbc_round(2) + 1
                if (input_value(2).lt.0) pbc_round(2) = pbc_round(2) - 1
              endif
              if (abs(input_value(3)-pbc_round(3)).ge.0.5) then
                if (input_value(3).gt.0) pbc_round(3) = pbc_round(3) + 1
                if (input_value(3).lt.0) pbc_round(3) = pbc_round(3) - 1
              endif
              B(1) = dx - ax*pbc_round(1)
              B(2) = dy - ay*pbc_round(2)
              B(3) = dz - az*pbc_round(3)
              r = ( B(1)**2 + B(2)**2 + B(3)**2 )**0.5
              if (r.gt.d) then
                write(*,*) "WARNING --- PROBLEM WITH PBC"
                goto 100
              endif
              
              ! Final COM result for molecule
              Xcom(nmol) = (1/CO2mass) * ( Cmass*X(j) + Omass*(2*X(j)+A(1)+B(1)) )
              Ycom(nmol) = (1/CO2mass) * ( Cmass*Y(j) + Omass*(2*Y(j)+A(2)+B(2)) )
              Zcom(nmol) = (1/CO2mass) * ( Cmass*Z(j) + Omass*(2*Z(j)+A(3)+B(3)) )

            endif

          enddo

          ! Write to CO2_COM.xyz file
          write(2,'(i2)') nmol
          write(2,'(i6)') i

          ! Wrap and write COM's to file
          do j=1,nmol

 200        if (Xcom(j).lt.0.0) then
              Xcom(j) = Xcom(j) + ax
              goto 200
            endif
 210        if (Xcom(j).gt.ax) then
              Xcom(j) = Xcom(j) - ax
              goto 210
            endif
 300        if (Ycom(j).lt.0.0) then
              Ycom(j) = Ycom(j) + ay
              goto 300
            endif
 310        if (Ycom(j).gt.ay) then
              Ycom(j) = Ycom(j) - ay
              goto 310
            endif
 400        if (Zcom(j).lt.0.0) then
              Zcom(j) = Zcom(j) + az
              goto 400
            endif
 410        if (Zcom(j).gt.az) then
              Zcom(j) = Zcom(j) - az
              goto 410
            endif

            ! Use H to denote molecular symbol for the COM's
            write(2,'(A1,x,f10.6,x,f10.6,x,f10.6)') 'H',Xcom(j),Ycom(j),Zcom(j)

          enddo

        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
