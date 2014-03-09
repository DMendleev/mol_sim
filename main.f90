!--------------------------------------------------------------------------!
!                            SAW - Homopolymer                             !
!                        Lamella - Block Copolymer                         !
!--------------------------------------------------------------------------!

PROGRAM main

   USE sumup
   USE variables
   USE moves
   USE util
   USE util_random
   
   IMPLICIT NONE
   SAVE

   LOGICAL             :: linit, didimove, initrandom
   INTEGER             :: i, j, nstep, nbead,length, avlength
   REAL                :: t, energy, engmove, engmovetot, aveng, rm
   REAL                :: rg2, rg2_para, rg2_perp, rend, rend_para, rend_perp
   REAL                :: rg2tot, rg2tot_para, rg2tot_perp, rendtot, rendtot_para, rendtot_perp
   INTEGER             :: iounit, rx, ry, rz, check_condition, c, naccept, sample
   INTEGER             :: vec(3)

   open(unit=11,action='read',file='input.dat',status='old')
   open(unit=13,action='write',file='results.dat')
   open(unit=42,action='write',file='movie.xyz')

   iounit = 11

   read(iounit,*) 
   read(iounit,*) nseed
   read(iounit,*)
   read(iounit,*) nstep  
   read(iounit,*)
   read(iounit,*) t
   read(iounit,*)
   read(iounit,*) nbead
   read(iounit,*)
   read(iounit,*) linit
   read(iounit,*) 
   read(iounit,*) initrandom
   read(iounit,*)
   read(iounit,*) pmreptation,pmpivot
   read(iounit,*) 
   read(iounit,*) mono, diblock, triblock
   read(iounit,*)
   read(iounit,*) sample  
   close(iounit)
  
    
   L  = nbead
   naccept = 0

   dir(:,1)  = (/1,0,0/) 
   dir(:,2)  = (/-1,0,0/)
   dir(:,3)  = (/0,1,0/)
   dir(:,4)  = (/0,-1,0/)
   dir(:,5)  = (/0,0,1/)
   dir(:,6)  = (/0,0,-1/)
   
   call ranset(nseed)

   if (linit) then
      call initia(nbead,initrandom)
      call totenergy(nbead,energy)
   else
      open(unit=12,action='read',file='config_final.dat',status='old')
      allocate(x(nbead),y(nbead),z(nbead),jij(nbead,nbead))
          do i=1,nbead
          read(12,*) x(i),y(i),z(i)
          enddo
      read(12,*) energy
   endif


   close(12)
   close(21)
   close(78)

   engmovetot = energy
   call writexyz(42,nbead,0)
    
   do i = 1,nstep
      rm = random()
      if (rm <= pmreptation) then
         call slither(nbead, engmovetot, t, didimove)
         if (didimove .eqv. .true.) then
             naccept = naccept + 1
             energy  = energy + engmovetot
!            call writexyz(42,nbead,i) 
!            write(99,*) i,energy/i,length/i
         endif
      else
         call pivot(nbead,engmovetot,t,didimove)
         if (didimove .eqv. .true.) then
             naccept = naccept + 1
             energy  = energy + engmovetot
!            call writexyz(42,nbead,i)
!            write(99,*) i,energy/i,length/i
         endif
      endif

      if(mod(i,sample)==0)then
         call properties(nbead, rg2, rg2_para, rg2_perp, rend, rend_para, rend_perp)      
         rg2tot       = rg2tot + rg2
         rg2tot_para  = rg2tot_para + rg2_para
         rg2tot_perp  = rg2tot_perp + rg2_perp
         rendtot      = rendtot + rend
         rendtot_para = rendtot_para + rend_para
         rendtot_perp = rendtot_perp + rend_perp
      end if
   
      check_condition = 0
      do j=1,nbead-1
        vec(1) = x(j+1)- x(j)
        vec(2) = y(j+1)- y(j)
        vec(3) = z(j+1)- z(j)
        c = dot_product(vec,vec)
        if(c.ne.1)then
           check_condition = 1
           write(*,*)i
           stop
        endif
      end do

   end do
   
   open(unit=77,action='write',file='config_final.dat') 
   do i=1,nbead
      write(77,*) x(i),y(i),z(i)
   enddo
   write(77,*) engmovetot
   
!   call PROPERTIES(nbead, rg2, rend)
   rg2tot       = rg2tot / (nstep)*sample
   rg2tot_para  = rg2tot_para / (nstep)*sample
   rg2tot_perp  = rg2tot_perp / (nstep)*sample
   rendtot      = rendtot / (nstep)*sample
   rendtot_para = rendtot_para / (nstep)*sample
   rendtot_perp = rendtot_perp / (nstep)*sample
   energy       = energy/(nstep)

   write(13,*)'-------------------------------------------------------------'
   write(13,*)'                           RESULTS                           '
   write(13,*)'-------------------------------------------------------------'
   write(13,*) '|  1. Seed                              = ', nseed
   write(13,*) '|  2. Temperature [J/kb]                = ', t
   write(13,*) '|  3. Number of Steps                   = ', nstep
   write(13,*) '|  4. Energy                            = ', energy
   write(13,*) '|  5. Rg2  - entire chain               = ', rg2tot
   write(13,*) '|  6. Rg2  - parallel to field          = ', rg2tot_para
   write(13,*) '|  7. Rg2  - perpendicular to field     = ', rg2tot_perp
   write(13,*) '|  8. Rend - entire chain               = ', rendtot
   write(13,*) '|  9. Rend - parallel to field          = ', rendtot_para
   write(13,*) '| 10. Rend - perpendicular to field     = ', rendtot_perp
   write(13,*) '| 11. Acceptance Percentage             = ', real(naccept)/real(nstep)*100
   write(13,*)'----- --------------------------------------------------------'
   
   close(77)
   close(13)   
   close(99)

END PROGRAM main
