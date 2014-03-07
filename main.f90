program main
use sumup
use variables
use moves
use util
use util_random
implicit none
save

   logical::linit,didimove,initrandom,writeintr,new_polymer,laminar
   integer::i,nstep,nbead,dimlatt,length,avlength,j,lcv
   real::t,energy,engmove,engmovetot,aveng,rm,rg2,rend,rg,averg2,averend,endtoend,rg2tot,rendtot
   integer::iounit, rx, ry, rz, check_condition, c, naccept, sample
   real :: c1, c2, c3
   integer, dimension(3) :: vec

   iounit=11   

   open(unit=iounit,action='read',file='input.dat',status='old')
   read(iounit,*) 
   read(iounit,*) nseed
   read(iounit,*)
   read(iounit,*) nstep  
   read(iounit,*)
   read(iounit,*) t
   read(iounit,*)
   read(iounit,*) nbead
   read(iounit,*)
   read(iounit,*) linit,new_polymer
   read(iounit,*) 
   read(iounit,*) writeintr 
   read(iounit,*) 
   read(iounit,*) initrandom
   read(iounit,*)
   read(iounit,*) pmreptation,pmpivot
   read(iounit,*) 
   read(iounit,*) mono, diblock, triblock
   close(iounit)
   
   call ranset(nseed)

   open(unit=13,action='write',file='results.dat')
   open(unit=42,action='write',file='movie.xyz')
   open(unit=99,action='write',file='intermediate.dat')
    
   !write(6,*) nstep,t,nbead,linit  
   if (linit) then
      !write(6,*) 'call initia'
      call initia(nbead,initrandom,new_polymer)
      length = nbead
      L = nbead
       !write(6,*) 'call totenergy'
      call totenergy(nbead,energy)
   else
      open(unit=12,action='read',file='config_final.dat',status='old')
      open(unit=21,action='read',file='polymer.dat',status='old')
      open(unit=78,action='read',file='interaction.dat',status='old')
      !allocate(moltyp(nbead),x(nbead),y(nbead),z(nbead),jij(nbead,nbead))
      allocate(x(nbead),y(nbead),z(nbead),jij(nbead,nbead))
       do i=1,nbead
         read(12,*) x(i),y(i),z(i)
       !  read(21,*) moltyp(i)
      enddo
!      do i=1,nbead
!         do j=i,nbead
!            read(78,*) jij(i,j)
!            jij(j,i) = jij(i,j)
!         enddo
!      enddo
      read(12,*) energy
      close(12)
      close(21)
      close(78)
   endif

         dir(:,1)  = (/1,0,0/) 
         dir(:,2)  = (/-1,0,0/)
         dir(:,3)  = (/0,1,0/)
         dir(:,4)  = (/0,-1,0/)
         dir(:,5)  = (/0,0,1/)
         dir(:,6)  = (/0,0,-1/)

   engmovetot = energy
   call writexyz(42,nbead,0)
    
   sample = 200
   naccept = 0
   do i=1,nstep
      rm=random()
      if (rm <= pmreptation) then
         call slither(nbead, engmove, t, didimove)
         if (didimove .eqv. .true.) then
             naccept = naccept + 1
!            call writexyz(42,nbead,i) 
!            write(99,*) i,energy/i,length/i
         endif
      else
         call pivot(nbead,engmovetot,t,didimove)
         if (didimove .eqv. .true.) then
             naccept = naccept + 1
!            call writexyz(42,nbead,i)
!            write(99,*) i,energy/i,length/i
         endif
      endif
      energy = energy + engmovetot

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

   enddo
   
   open(unit=77,action='write',file='config_final.dat') 
   do i=1,nbead
!      c1 = real(x(i))/real(L)                                                                                 
!      c2 = real(y(i))/real(L)                                                                                 
!      c3 = real(z(i))/real(L)                                                                                 
!      rx = x(i) - floor(c1)*L
!      ry = y(i) - floor(c2)*L
!      rz = z(i) - floor(c3)*L
!      write(77,*) rx,ry,rz
      write(77,*) x(i),y(i),z(i)
   enddo
   if (new_polymer) then
      open(unit=21,action='write',file='polymer.dat')
      do i=1,nbead
       !  write(21,*) moltyp(i)
      enddo
      close(21)
   endif
   write(77,*) engmove
   
   aveng = engmovetot/(nstep+1.0)
   call PROPERTIES(nbead, rg, endtoend)
   rg2 = rg
   rend =endtoend
   averg2 = rg2tot / (nstep)*sample
   averend = rendtot / (nstep)*sample

   write(13,*) 'nseed = ', nseed
   write(13,*) 'Temperature [J/kb] = ', t
   write(13,*) 'nstep = ', nstep
   write(13,*) 'Energy = ', aveng
   write(13,*) 'RADIUSG= ', rg2
   write(13,*) 'endtoend= ', rend
   write(13,*) 'average rg2= ', averg2
   write(13,*) 'average rend ', averend
   write(13,*) 'Acceptance Percentage ', real(naccept)/real(nstep)*100
   write(13,*) 'Acceptance Percentage ', naccept
  
   close(77)
   close(13)   
   close(99)

end program main
