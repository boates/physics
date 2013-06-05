!       program cbn_from_xyz_general.f90
!***********************************************************
!       Translation of Isaac's bonding_analysis.py
!       into fortran for speed
!       However, this assumes constant lattice constants
!***********************************************************
        implicit none
        
        ! n = max # of allowed timesteps
        integer n, i, j, k, l, m, g
        parameter (n=10000000)
        integer natom, tstep, Nsteps
        integer c, s, p, o, oo, molecule_number
        integer closest_to_p, closest_to_o
        real*8 r(3), alat, rX(27), rY(27), rZ(27)
        real*8 ax, ay, az, bx, by, bz, cx, cy, cz, amag, bmag, cmag
        real*4 POS_array(n,4)
        real*8 dX(27), dY(27), dZ(27), dR(27)
        real*8 sanity, distance, minimum_distance
        character*2 typat(1024)
        character*128 fxyz, fpos

        ! Get the xyz filename and alat from the user
        write(*,*) 'Name of xyz file'
        read(*,*) fxyz
        write(*,*) 'Name of POSCAR file for lattice information'
        read(*,*) fpos

        ! Open the input and output files
        open(1,file=fxyz,status='old',ERR=90)
        open(4,file=fpos,status='old',ERR=90)
        open(2,file='TRAJEC.cbn')
        open(3,file='diss.dat')

        ! Read lattice information from POSCAR header
        read(4,*)
        read(4,*) alat
        read(4,*) ax, ay, az
        read(4,*) bx, by, bz
        read(4,*) cx, cy, cz

        ! Scale the lattice vector components by alat
        ax = ax*alat
        ay = ay*alat
        az = az*alat
        bx = bx*alat
        by = by*alat
        bz = bz*alat
        cx = cx*alat
        cy = cy*alat
        cz = cz*alat
        amag = (ax**2 + ay**2 + az**2)**0.5
        bmag = (bx**2 + by**2 + bz**2)**0.5
        cmag = (cx**2 + cy**2 + cz**2)**0.5

        Nsteps = 0

        do i=1,n

          read(1,*,END=100) natom
          read(1,*) tstep
          Nsteps = Nsteps + 1

          do j=1,natom

            read(1,*) typat(j), r(1), r(2), r(3)
            POS_array(j+natom*(i-1),1) = r(1)/0.5291772
            POS_array(j+natom*(i-1),2) = r(2)/0.5291772
            POS_array(j+natom*(i-1),3) = r(3)/0.5291772
            POS_array(j+natom*(i-1),4) = -10
          
          enddo

        enddo

 100    continue

        close(1)
        
        sanity = ( (amag/2)**2 + (bmag/2)**2 + (cmag/2)**2 )**0.5

        s = 0

        do while (s.lt.Nsteps)

          molecule_number = 0

          ! counts over particles
          p = 0

          do while (p.lt.natom)

            if (POS_array(s*natom+p+1,4).eq.-10) then

              ! counts over other particles
              o = 0
              minimum_distance = sanity
              closest_to_p = p

              do while (o.lt.natom)

                ! Replicate the k atom in 3x3x3 supercell
                c = 0
                distance = sanity
                do l = -1, 1
                  do m = -1, 1
                    do g = -1, 1
                      c = c + 1
                      rX(c) = POS_array(s*natom+o+1,1) + l*ax + m*bx + g*cx
                      rY(c) = POS_array(s*natom+o+1,2) + l*ay + m*by + g*cy
                      rZ(c) = POS_array(s*natom+o+1,2) + l*az + m*bz + g*cz

                      ! Calculate distances and find the smallest
                      dX(c) = POS_array(s*natom+p+1,1) - rX(c)
                      dY(c) = POS_array(s*natom+p+1,2) - rY(c)
                      dZ(c) = POS_array(s*natom+p+1,3) - rZ(c)
                      dR(c) = (dX(c)**2 + dY(c)**2 + dZ(c)**2)**0.5
                      if (dR(c).lt.distance) distance = dR(c)

                    enddo
                  enddo
                enddo

                if (distance.gt.sanity) then
                  write(*,*) "Warning, problem with pbc"
                  goto 90
                endif

                if (distance.lt.minimum_distance) then
                  if (distance.ne.0.0) then
                    minimum_distance = distance
                    closest_to_p = o
                  endif
                endif

                o = o + 1

              enddo

              minimum_distance = sanity

              o = closest_to_p
              oo = 0

              if (POS_array(s*natom+o+1,4).eq.-10) then

                do while (oo.lt.natom)

                  ! Replicate the k atom in 3x3x3 supercell
                  c = 0
                  distance = sanity
                  do l = -1, 1
                    do m = -1, 1
                      do g = -1, 1
                        c = c + 1
                        rX(c) = POS_array(s*natom+oo+1,1) + l*ax + m*bx + g*cx
                        rY(c) = POS_array(s*natom+oo+1,2) + l*ay + m*by + g*cy
                        rZ(c) = POS_array(s*natom+oo+1,3) + l*az + m*bz + g*cz

                        ! Calculate distances and find the smallest
                        dX(c) = POS_array(s*natom+o+1,1) - rX(c)
                        dY(c) = POS_array(s*natom+o+1,2) - rY(c)
                        dZ(c) = POS_array(s*natom+o+1,3) - rZ(c)
                        dR(c) = (dX(c)**2 + dY(c)**2 + dZ(c)**2)**0.5
                        if (dR(c).lt.distance) distance = dR(c)

                      enddo
                    enddo
                  enddo

                  if (distance.gt.sanity) then
                    write(*,*) "Warning, problem with pbc"
                    goto 90
                  endif

                  if (distance.lt.minimum_distance) then
                    if (distance.ne.0.0) then
                      minimum_distance = distance
                      closest_to_o = oo
                    endif
                  endif

                  oo = oo + 1

                enddo

                if (closest_to_p.eq.o.and.closest_to_o.eq.p) then

                  molecule_number = molecule_number + 1

                  POS_array(s*natom+p+1,4) = o
                  POS_array(s*natom+o+1,4) = p

                else
                  POS_array(s*natom+p+1,4) = -1

                endif

              else
                POS_array(s*natom+p+1,4) = -1

              endif

            endif

            p = p + 1

          enddo
          
          write(3,*) 2.0*molecule_number/natom

          s = s + 1

        enddo

        close(3)

        ! Write the cbn file
        write(2,*) '# cbn file created with cbn_from_xyz.x'
        write(2,*) '# a =', amag!, ax, ay, az
        write(2,*) '# b =', bmag!, bx, by, bz
        write(2,*) '# c =', cmag!, cx, cy, cz
        write(2,*) '# number_of_particles =', natom
        write(2,*) '# number_of_neighbors =', 1
        write(2,*) '#'
        write(2,*) '#'
        write(2,*) '#'
        write(2,*) '# units = bohr'

        do i=1,Nsteps
          do j=1,natom
            write(2,*) typat(j),POS_array((i-1)*natom+j,1), &
                                POS_array((i-1)*natom+j,2), &
                                POS_array((i-1)*natom+j,3), &
                                int( POS_array((i-1)*natom+j,4) )
          enddo
        enddo

        close(2)

 90     continue

        END
