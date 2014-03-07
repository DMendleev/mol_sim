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
!!!!!!!!!!!!!!!!!!!!!!   
   IMPLICIT NONE
   SAVE

   LOGICAL             :: linit, didimove, initrandom
   INTEGER             :: i, j, nstep, nbead,length, avlength
   REAL                :: t, energy, engmove, engmovetot, aveng, rm, rg2
   REAL                :: rend, rg, averg2, averend, endtoend, rg2tot, rendtot
   INTEGER             :: iounit, rx, ry, rz, check_condition, c, naccept, sample
   INTEGER             :: vec(3)

   open(unit=11,action='read',file='input.dat',status='old')
   open(unit=13,action='write',file='results.dat')
   open(unit=42,action='write',file='movie.xyz')
   open(unit=12,action='read',file='config_final.dat',status='old')

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
   close(iounit)
  
    
   L  = nbead
   sample = 200
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
         call slither(nbead, engmove, t, didimove)
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
         call properties(nbead, rg2, rend)      
         rg2tot = rg2tot + rg2
         rendtot = rendtot + rend
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
   write(77,*) engmove
   
   call PROPERTIES(nbead, rg2, rend)
   averg2  = rg2tot / (nstep)*sample
   averend = rendtot / (nstep)*sample
   aveng   = engmovetot/(nstep)*sample

   write(13,*) 'nseed = ', nseed
   write(13,*) 'Temperature [J/kb] = ', t
   write(13,*) 'nstep = ', nstep
   write(13,*) 'Energy = ', aveng
   write(13,*) 'RADIUSG= ', rg2
   write(13,*) 'endtoend= ', rend
   write(13,*) 'average rg2= ', averg2
   write(13,*) 'average rend ', averend
   write(13,*) 'Acceptance Percentage ', real(naccept)/real(nstep)*100
  
   close(77)
   close(13)   
   close(99)

END PROGRAM main
