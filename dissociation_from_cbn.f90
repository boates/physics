!       program dissociation_from_cbn.f
!***********************************************************
!       Calculate the percentage of dissociated vs. molecular
!       atoms over time from a cbn file.
!***********************************************************
        implicit none
        
        ! n = max # of timesteps, m = max # of atoms
        integer maxAtoms, maxSteps, i, j, k
        parameter (maxSteps=100000,maxAtoms=10000)
        integer natom, n
        integer dis, mol
        real*8 f_dis, f_mol, x, y, z
        character*8 typat, d1, d2, d3, d4, d5
        character*128 fin

        ! Retrieve user input
        write(*,*) 'Name of cbn file:'
        read(*,*) fin

        ! Open cbn and output files
        open(1,file=fin,status='old',ERR=90)
        open(2,file='dissociation.dat')
        write(2,*) '# tstep, f_dis, f_mol'

        ! Read in cbn header info
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) d1, d2, d3, natom
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)

        do i=1,maxSteps

            ! Initialize counting variables as zero
            dis = 0
            mol = 0

          do j=1,natom

            ! Read in data for each atom at the current timestep
            read(1,*,END=100) typat, x, y, z, n

            ! Check if dissociated or not
            if (n.lt.0) then
              dis = dis + 1
            endif
            if (n.ge.0) then
              mol = mol + 1
            endif

          enddo

          ! Check for atom conservation
          if ((1.0*dis + 1.0*mol)/natom.ne.1.0) then
            write(*,*) 'dis + mol != natom'
            write(*,*) 'atoms not conserved, exiting...'
            goto 100
          endif

          ! Caclulate fraction values
          f_dis = (1.0*dis) / natom
          f_mol = (1.0*mol) / natom
          
          ! Write to file
          write(2,*) i, f_dis, f_mol

        enddo

 100    continue

        close(1)

        close(2)

        goto 110

 90     continue

        write(*,*) 'No file:',fin

 110    continue

        END
