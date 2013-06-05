!       program SDF.f
!***********************************************************
!       Cacluates the 'Spatial Distribution Function'
!       Requires a cbn file.
!***********************************************************
        implicit none
        
        integer maxSteps, maxAtoms, i, j, k
        parameter (maxSteps=100000,maxAtoms=1024)
        integer natom, tstep, B(maxAtoms)
        real r, theta, X(maxAtoms), Y(maxAtoms), Z(maxAtoms)
        real alat(3), CX(maxAtoms), CY(maxAtoms), CZ(maxAtoms)
        real dx, dy, dz, ox, oy, oz, oox, ooy, ooz
        real dop, ntp, ndp, nip, nto, ndo, nio
        real do_oo, doop, dooo, ndrp, ndro
        character*8 typat, d1, d2, d3, d4, d5

        ! Open input and output files
        open(1,file='TRAJEC.cbn',status='old',ERR=90)
        open(2,file='circle_all.dat')
        open(3,file='circle_nearest_pole.dat')
        open(4,file='circle_nearest_origin.dat')

        ! Read in cbn file header
        read(1,*) 
        read(1,*) d1, d2, d3, alat(1)
        read(1,*) d1, d2, d3, alat(2)
        read(1,*) d1, d2, d3, alat(3)
        read(1,*) d1, d2, d3, d4, d5, natom
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 
        read(1,*) 

        ! Write output file headers
        write(2,*) '# r theta bonding'
        write(3,*) '# r theta bonding'
        write(4,*) '# r theta bonding'

        ! BEGIN ANALYSIS
        do i=1,maxSteps

          ! Read in current snapshot data from cbn file
          do j=1,natom
            read(1,*,END=100) typat, X(j), Y(j), Z(j), B(j)
          enddo

          do j=1,natom

            ! If particle is bonded
            if (B(j).ne.-1) then 

              do k=1,natom
                dx = X(k) - X(j)
                dy = Y(k) - Y(j)
                dz = Z(k) - Z(j)
                X(k) = dx - int(dx/alat(1) + 0.5)*alat(1)
                Y(k) = dy - int(dy/alat(2) + 0.5)*alat(2)
                Z(k) = dz - int(dz/alat(3) + 0.5)*alat(3)
              enddo

              ! Coords of atom that j is bonded to
              ox = X( B(j) )
              oy = Y( B(j) )
              oz = Z( B(j) )

              ! Do some center shifting of coords
              do k=1,natom
                CX(k) = X(k) - ox/2.0
                CY(k) = Y(k) - oy/2.0
                CZ(k) = Z(k) - oz/2.0
                CX(k) = CX(k) - int(CX(k)/alat(1) + 0.5)*alat(1)
                CY(k) = CY(k) - int(CY(k)/alat(2) + 0.5)*alat(2)
                CZ(k) = CZ(k) - int(CZ(k)/alat(3) + 0.5)*alat(3)
              enddo

              X(j) = CX(j)
              Y(j) = CY(j)
              Z(j) = CZ(j)
              ox = CX( B(j) )
              oy = CY( B(j) )
              oz = CZ( B(j) )

              ! Acronyms of variables from Isaac's original code
              dop = (X(j)**2 + Y(j)**2 + Z(j)**2)**0.5
              ntp = -1
              ndp = 1000000
              nip = maxSteps*natom
              nto = -1
              ndo = 1000000
              nio = maxSteps*natom

              ! Meat of the code
              do k=1,natom

                if (k.ne.j.and.k.ne.B(j)) then

                  oox = CX(k)
                  ooy = CY(k)
                  ooz = CZ(k)

                  do_oo = (oox**2 + ooy**2 + ooz**2)**0.5
                  theta = acos( (X(j)*oox + Y(j)*ooy + Z(j)*ooz)/(dop*do_oo) )

                  doop = ( (oox-X(j))**2 + (ooy-Y(j))**2 + (ooz-Z(j))**2 )**0.5
                  dooo = ( (oox-ox)**2 + (ooy-oy)**2 + (ooz-oz)**2 )**0.5

                  if (min(doop,dooo).lt.ndp) then
                    nip = k
                    ndp = min(doop,dooo)
                    ntp = theta
                    ndrp = do_oo
                  endif

                  if (do_oo.lt.ndo) then
                    nio = k
                    ndo = do_oo
                    nto = theta
                    ndro = do_oo
                  endif

                  write(2,*) do_oo, theta, B(k)

                endif
              enddo

              write(3,*) ndrp, ntp, B(nip)
              write(4,*) ndro, nto, B(nio)

            endif
          enddo
        enddo

 100    continue

        close(1)
        close(2)
        close(3)
        close(4)

 90     continue

        END
