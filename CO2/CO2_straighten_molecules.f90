!       program CO2_straighten_molecules.f90
!***********************************************************
!       Straighten's CO2 molecules, for one tstep only!
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, maxBins, i, j, k, f
        parameter (maxSteps=100000,maxAtoms=1024, maxBins=1000)
        integer natom, nC, nO, nSteps, nStraight
        integer nn(maxAtoms,2), pbc_round(3), sC(maxAtoms)
        real*8 X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
        real*8 ax, ay, az, dx, dy, dz, r
        real*8 A(3), B(3), C(3), Amag, Bmag, input_value(3)
        character*8 typat(maxAtoms), d1, d2, d3
        character*532 fmol, fmt

        ! Retrieve user input
        write(*,*) 'Name of single snapshot .mol file'
        read(*,*) fmol
        write(*,*) 'Number of molecules to straighten'
        read(*,*) nStraight

        ! Open .mol and output files
        open(1,file=fmol,status='old',ERR=90)
        open(2,file='STRAIGHT.xyz')

        ! Read in .mol header
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

        ! Read in .mol information
        do j=1,natom
          read(1,*,END=100) typat(j), X(j), Y(j), Z(j), nn(j,1), nn(j,2)
        enddo

        ! If all molecules
        if (nStraight.eq.nC) then
          do i=1,nC
            sC(i) = i
          enddo
          goto 50          
        endif

        ! If not all carbon atoms
        write(*,*) 'Indicies of carbons of molecules to straighten'
        write(*,*) 'i.e. 1,3,11,26,27,32'
        read(*,*) (sC(i),i=1,nStraight)
 
 50     continue

        ! Do the requested molecular straightening
        do j=1,nStraight

          if ((nn(sC(j),1).ne.0).and.(nn(sC(j),2).ne.0)) then

            dx = X(sC(j)) - X(nn(sC(j),1))
            dy = Y(sC(j)) - Y(nn(sC(j),1))
            dz = Z(sC(j)) - Z(nn(sC(j),1))
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
            if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
              write(*,*) "Warning: Problem with PBC --- Exiting..."
              goto 100
            endif

            dx = X(sC(j)) - X(nn(sC(j),2))
            dy = Y(sC(j)) - Y(nn(sC(j),2))
            dz = Z(sC(j)) - Z(nn(sC(j),2))
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
            if (r.gt.((ax**2+ay**2+az**2)**0.5)/2.0) then
              write(*,*) "Warning: Problem with PBC --- Exiting..."
              goto 100
            endif

            ! Calculate the straightened vector
            Amag = (A(1)**2+A(2)**2+A(3)**2)**0.5
            Bmag = (B(1)**2+B(2)**2+B(3)**2)**0.5
            C(1) = X(sC(j)) + Bmag*A(1)/Amag
            C(2) = Y(sC(j)) + Bmag*A(2)/Amag
            C(3) = Z(sC(j)) + Bmag*A(3)/Amag

            ! Make sure new coordinate is wrapped
 10         if (C(1).gt.ax) then
              C(1) = C(1) - ax
              goto 10
            endif
 11         if (C(1).lt.0.0) then
               C(1) = C(1) + ax
              goto 11
            endif
 20         if (C(2).gt.ay) then
              C(2) = C(2) - ay
              goto 20
            endif
 21         if (C(2).lt.0.0) then
              C(2) = C(2) + ay
              goto 21
            endif
 30         if (C(3).gt.az) then
              C(3) = C(3) - az
              goto 30
            endif
 31         if (C(3).lt.0.0) then
              C(3) = C(3) + az
              goto 31
            endif

            ! Set coordinates to straightened/wrapped ones 
            X(nn(sC(j),2)) = C(1)
            Y(nn(sC(j),2)) = C(2)
            Z(nn(sC(j),2)) = C(3)

          endif
        enddo

        ! Write the straightened snapshot to .xyz file
        write(2,'(i2)') natom
        write(2,'(A1)') '1'
        do i=1,natom
          write(2,'(A1,x,f10.6,x,f10.6,x,f10.6)') typat(i), X(i), Y(i), Z(i)
        enddo

 100    continue

        close(1)
        close(2)

 90     continue

        END
