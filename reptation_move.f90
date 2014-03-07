!*****************************************************************************************!
!                                MONTE CARLO MOVES                                        !
!*****************************************************************************************!

     MODULE MC_MOVES

         use shared_parameter
         use parameters
         use rng_module
         use latt_func

         IMPLICIT NONE

         CONTAINS

               !------------------------------------------------------------!
               !                SLITHERING SNAKE ALGORITHM                  !
               !------------------------------------------------------------!
!{{{
         SUBROUTINE slither(nbead,engmovetot,didimove, rg2, rend)
    
         IMPLICIT NONE

         INTEGER, ALLOCATABLE                         :: rx_trial, ry_trial, rz_trial
         INTEGER, DIMENSION(3)                        :: s 
         INTEGER                                      :: i, z, head, tail, c, v2, pnt, saw
         INTEGER                                      :: nbead
         LOGICAL                                      :: didimove
         REAL                                         :: engmovetot, engn

         allocate(rx_trial(nbead)) 
         allocate(ry_trial(nbead)) 
         allocate(rz_trial(nbead)) 

         rx_trial = 0                                  ! Dum vector is used to take the coordinate of one chain only              
         ry_trial = 0                                  ! and work on it for slithering alrogithm
         rz_trial = 0
         call totalenergy(nbead,engo)
        
         v2 = ceiling(random()*2)                                             ! Selecting one as head and other as a tail  
                                                       ! out of the twon ends od a chain
         if(v2==1)then
            head = 1                                                
            tail = nbead
         else
            head = nbead
            tail = 1
         endif

         pnt  = iuni(1,7)
         s(1) = dir(1,pnt)
         s(2) = dir(2,pnt)
         s(3) = dir(3,pnt)                                   ! Randomly selecting one of the 12 (for FCC) directions possible.
         rx_trial(head) = rx(head) + s(1)                  ! Storing the position of new head in the dum vector
         ry_trial(head) = ry(head) + s(2)
         rz_trial(head) = rz(head) + s(3)

         
           if(head==nbead)then
              rx_trial(1:nbead-1) = rx(2:nbead)              ! If head = N_dp, then 2:N_dp part of old chain will become 1:N_dp-1 part of new chain 
              ry_trial(1:nbead-1) = ry(2:nbead)
              rz_trial(1:nbead-1) = rz(2:nbead)

           else

              rx_trial(2:nbead) = rx(1:nbead-1)              ! If head = 1, then 1:N_dp-1 part of old chain will become 2:N_dp part of new chain 
              ry_trial(2:nbead) = ry(1:nbead-1)
              rz_trial(2:nbead) = rz(1:nbead-1)
            
           endif 

                                                              ! Checking Overlap. if saw  = 0 -  No overlap.
           saw = 0 
           do i=1,nbead
              do j=i+1,nbead    
                 if(rx_trial(i).eq.rx_trial(j).and.ry_trial(i).eq.ry_trial(j).and.rz_trial(i).eq.rz_trial(j))then
                    saw = 1
                 endif
              end do
           end do
    
           if(saw == 0)then
              call totalenergy(nbead,engn)
              P_acc  =exp(-(engn - engo))
              if(random().lt.P_acc)then
                 rx = rx_trial
                 ry = ry_trial
                 rz = rz_trial
                 didimove = .true.
                 engmovetot = engn               
              else
                 didimove = .false.
              endif
    
           else
              didimove = .false.
           endif
            

         END SUBROUTINE slither
!}}}
!*****************************************************************************************!

!*****************************************************************************************!
               !------------------------------------------------------------!
               !                  PIVOT ALGORITHM                           !
               !------------------------------------------------------------!


!subroutine pivot(nbead,engmovetot,lengthmove,t,didimove)
!
!   integer::i,j,k,nbead,beadrot,lengthmove, lcv
!   integer,allocatable::xo(:),yo(:),zo(:),xm(:),ym(:),zm(:),xmo(:),ymo(:),zmo(:)
!   real::t,diff,real,engo,engn,engmovetot,rm
!   logical::didimove,totheleft,clockw
!
!   allocate(xo(nbead)); allocate(yo(nbead)); allocate(zo(nbead))
!   allocate(xm(nbead)); allocate(ym(nbead)); allocate(zm(nbead))
!   allocate(xmo(nbead)); allocate(ymo(nbead)); allocate(zmo(nbead))
!   didimove = .true.
!!calculate energy in old config.
!   call totenergy(nbead,engo)
!
!   engmovetot = engo
!   lengthmove = maxlength(nbead)
!! choose bead at random
!   beadrot = int(random()*nbead + 1)
!   
!   if (beadrot == 1) then
!      totheleft = .false.
!   elseif (beadrot == nbead) then
!      totheleft = .true.
!   elseif (random() < 0.5) then
!      totheleft = .true.
!   else
!      totheleft = .false.
!   endif
!
!   do i=1,nbead
!      xo(i) = x(i)
!      yo(i) = y(i)
!      zo(i) = z(i)
!      xmo(i) = x(i) - x(beadrot)
!      ymo(i) = y(i) - y(beadrot)
!      zmo(i) = z(i) - z(beadrot)
!      xm(i) = x(i) - x(beadrot)
!      ym(i) = y(i) - y(beadrot)
!      zm(i) = z(i) - z(beadrot)
!   enddo
!
!!choose clock-wise or anti-clockwise
!   if (random() < 0.5) then
!      clockw = .true.
!   else
!      clockw = .false.
!   endif
!!rotate
!   rm = random()*3.0
!   if (rm < 1.0) then
!! in x,y-plane
!      if (totheleft .and. .not. clockw) then
!         do i=1,beadrot-1
!            xm(i) = -ymo(i)
!            ym(i) = xmo(i)
!            zm(i) = zmo(i)
!         enddo
!      elseif (totheleft .and. clockw) then
!         do i=1,beadrot-1
!            xm(i) = ymo(i)
!            ym(i) = -xmo(i)
!            zm(i) = zmo(i)
!         enddo
!      elseif (.not. totheleft .and. .not. clockw) then
!         do i=beadrot+1,nbead
!            xm(i) = -ymo(i)
!            ym(i) = xmo(i)
!            zm(i) = zmo(i)
!         enddo
!      else
!          do i=beadrot+1,nbead
!            xm(i) = ymo(i)
!            ym(i) = -xmo(i)
!            zm(i) = zmo(i)
!          enddo
!      endif
!   elseif (rm < 2.0) then
!! in x,z-plane
!      if (totheleft .and. .not. clockw) then
!         do i=1,beadrot-1
!            xm(i) = -zmo(i)
!            ym(i) = ymo(i)
!            zm(i) = xmo(i)
!         enddo
!      elseif (totheleft .and. clockw) then
!         do i=1,beadrot-1
!            xm(i) = zmo(i)
!            ym(i) = ymo(i)
!            zm(i) = -xmo(i)
!         enddo
!      elseif (.not. totheleft .and. .not. clockw) then
!         do i=beadrot+1,nbead
!            xm(i) = -zmo(i)
!            ym(i) = ymo(i)
!            zm(i) = xmo(i)
!         enddo
!      else
!          do i=beadrot+1,nbead
!            xm(i) = zmo(i)
!            ym(i) = ymo(i)
!            zm(i) = -xmo(i)
!          enddo
!      endif   
!   else
!! in y,z-plane
!      if (totheleft .and. .not. clockw) then
!         do i=1,beadrot-1
!            xm(i) = xmo(i)
!            ym(i) = -zmo(i)
!            zm(i) = ymo(i)
!         enddo
!      elseif (totheleft .and. clockw) then
!         do i=1,beadrot-1
!            xm(i) = xmo(i)
!            ym(i) = zmo(i)
!            zm(i) = -ymo(i)
!         enddo
!      elseif (.not. totheleft .and. .not. clockw) then
!         do i=beadrot+1,nbead
!            xm(i) = xmo(i)
!            ym(i) = -zmo(i)
!            zm(i) = ymo(i)
!         enddo
!      else
!          do i=beadrot+1,nbead
!            xm(i) = xmo(i)
!            ym(i) = zmo(i)
!            zm(i) = -ymo(i)
!          enddo
!      endif     
!   endif
!
!! new coordinates
!   do i=1,nbead
!      x(i) = xm(i) + x(beadrot)
!      y(i) = ym(i) + y(beadrot)
!      z(i) = zm(i) + z(beadrot)   
!   enddo
!   
!!check for overlap
!   do i = 1, nbead
!      do j = i+1, nbead
!      	 if (x(i) == x(j) .AND. y(i) == y(j) .AND. z(i) == z(j)) then
!            do k=1,nbead
!               x(k) = xo(k)
!               y(k) = yo(k)
!               z(k) = zo(k)
!            enddo
!            !write(6,*) 'rejected: overlap'
!            didimove = .false.
!            return
!         endif
!      enddo
!   enddo
!
!!calculate energy in new config.
!   call totenergy(nbead,engn)  
!   diff = engn - engo
!   if (diff .le. 0.0) then
!      !write(6,*) x(ibead),y(ibead),z(ibead)
!! accept
!   else if (exp(-diff/t) .gt. random()) then
!      !write(6,*) 'bead = ', ibead
!! accept
!   else
!      do i=1,nbead
!         x(i) = xo(i)
!         y(i) = yo(i)
!         z(i) = zo(i)
!      enddo
!      !write(6,*) 'rejected: boltz=',exp(-diff/t),'acc=',acc
!      didimove = .false.
!      return
!   endif
!  !write(6,*) 'accepted: overlap'  
!
!  engmovetot = diff + engo
!  lengthmove = maxlength(nbead)
!
!! recenter
!  call recenter(nbead)
!end subroutine pivot



     END MODULE MC_MOVES

!*****************************************************************************************!
