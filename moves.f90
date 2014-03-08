module moves
use sumup
use variables
use util_random,only:random
use util
implicit none
private
save
public:: initia,pivot,slither,random_init

contains
subroutine initia(nbead,initrandom)

   integer::nbead,i,j
   integer,dimension(1)::nseed
   logical::initrandom

   allocate(x(nbead))
   allocate(y(nbead))
   allocate(z(nbead))
   !allocate(moltyp(nbead))
   allocate(jij(nbead,nbead))

! Decide to initialize randomly or in line
 if(initrandom .eqv. .false.) then

! place beads in a line   
   !write(6,*) 'place beads in line'
   do i=1,nbead
      x(i) = nbead
      y(i) = i
      z(i) = nbead
   enddo
 else
  call random_init(nbead) 
   !write(6,*) 'placed beads randomly'
 endif
!  determine com
   call determineCOM(nbead) 
   
   !write(6,*) 'initialization complete'
   return
end subroutine initia

subroutine random_init(nbead)

   integer::nbead, i, j
   real::rwalk
   logical::overlap
   x(1) = nbead
   y(1) = 1
   z(1) = nbead
   overlap = .true.
   do i = 2, nbead
        !write(6,*) 'i=',i
        overlap = .true.
        do while (overlap .eqv. .true.)
         rwalk = random()*6.0
         !write(6,*) 'rand=',rwalk
         if (rwalk < 1) then
            x(i) = x(i-1) - 1
            y(i) = y(i-1)
            z(i) = z(i-1)
         else if (rwalk < 2)then
            x(i) = x(i-1) + 1
            y(i) = y(i-1)
            z(i) = z(i-1)
         else if (rwalk < 3) then
            y(i) = y(i-1) - 1
            x(i) = x(i-1)
            z(i) = z(i-1)
         else if (rwalk < 4) then
            y(i) = y(i-1) + 1
            x(i) = x(i-1)
            z(i) = z(i-1)
         else if (rwalk < 5) then
            z(i) = z(i-1) - 1
            y(i) = y(i-1)
            x(i) = x(i-1)
         else
            z(i) = z(i-1) + 1
            y(i) = y(i-1)
            x(i) = x(i-1)
         endif
         !write(6,*) 'going to ',x(i),y(i),z(i)
         do j = 1, i-1 
            !write(6,*) 'j=',j
         !Check for overlap
         !If overlap, redo point selection, overlap true 
         !If no overlap, overlap false, go to next point
            if (x(i) == x(j) .AND. y(i) == y(j) .AND. z(i) == z(j)) then
               !write(6,*) 'rejected: overlap'
               overlap = .true.
               exit 
            else
               !write(6,*) 'accepted'
               overlap = .false.
            endif
         enddo
         !write(6,*) 'overlap check complete'
      enddo
   enddo
end subroutine random_init

subroutine pivot(nbead,engmovetot,t,didimove)

   integer::i,j,k,nbead,beadrot, lcv, overlap
   integer,allocatable::xo(:),yo(:),zo(:),xm(:),ym(:),zm(:),xmo(:),ymo(:),zmo(:)
   real::t,diff,real,engo,engn,engmovetot,rm,rg,endtoend,rg2,rend
   logical::didimove,totheleft,clockw

   allocate(xo(nbead)); allocate(yo(nbead)); allocate(zo(nbead))
   allocate(xm(nbead)); allocate(ym(nbead)); allocate(zm(nbead))
   allocate(xmo(nbead)); allocate(ymo(nbead)); allocate(zmo(nbead))
   didimove = .true.
!calculate energy in old config.
   call totenergy(nbead,engo)

   engmovetot = engo
   
! choose bead at random
   beadrot = int(random()*nbead + 1)
   
   if (beadrot == 1) then
      totheleft = .false.
   elseif (beadrot == nbead) then
      totheleft = .true.
   elseif (random() < 0.5) then
      totheleft = .true.
   else
      totheleft = .false.
   endif

   do i=1,nbead
      xo(i) = x(i)
      yo(i) = y(i)
      zo(i) = z(i)
      xmo(i) = x(i) - x(beadrot)
      ymo(i) = y(i) - y(beadrot)
      zmo(i) = z(i) - z(beadrot)
      xm(i) = x(i) - x(beadrot)
      ym(i) = y(i) - y(beadrot)
      zm(i) = z(i) - z(beadrot)
   enddo

!choose clock-wise or anti-clockwise
   if (random() < 0.5) then
      clockw = .true.
   else
      clockw = .false.
   endif
!rotate
   rm = random()*3.0
   if (rm < 1.0) then
! in x,y-plane
      if (totheleft .and. .not. clockw) then
         do i=1,beadrot-1
            xm(i) = -ymo(i)
            ym(i) = xmo(i)
            zm(i) = zmo(i)
         enddo
      elseif (totheleft .and. clockw) then
         do i=1,beadrot-1
            xm(i) = ymo(i)
            ym(i) = -xmo(i)
            zm(i) = zmo(i)
         enddo
      elseif (.not. totheleft .and. .not. clockw) then
         do i=beadrot+1,nbead
            xm(i) = -ymo(i)
            ym(i) = xmo(i)
            zm(i) = zmo(i)
         enddo
      else
          do i=beadrot+1,nbead
            xm(i) = ymo(i)
            ym(i) = -xmo(i)
            zm(i) = zmo(i)
          enddo
      endif
   elseif (rm < 2.0) then
! in x,z-plane
      if (totheleft .and. .not. clockw) then
         do i=1,beadrot-1
            xm(i) = -zmo(i)
            ym(i) = ymo(i)
            zm(i) = xmo(i)
         enddo
      elseif (totheleft .and. clockw) then
         do i=1,beadrot-1
            xm(i) = zmo(i)
            ym(i) = ymo(i)
            zm(i) = -xmo(i)
         enddo
      elseif (.not. totheleft .and. .not. clockw) then
         do i=beadrot+1,nbead
            xm(i) = -zmo(i)
            ym(i) = ymo(i)
            zm(i) = xmo(i)
         enddo
      else
          do i=beadrot+1,nbead
            xm(i) = zmo(i)
            ym(i) = ymo(i)
            zm(i) = -xmo(i)
          enddo
      endif   
   else
! in y,z-plane
      if (totheleft .and. .not. clockw) then
         do i=1,beadrot-1
            xm(i) = xmo(i)
            ym(i) = -zmo(i)
            zm(i) = ymo(i)
         enddo
      elseif (totheleft .and. clockw) then
         do i=1,beadrot-1
            xm(i) = xmo(i)
            ym(i) = zmo(i)
            zm(i) = -ymo(i)
         enddo
      elseif (.not. totheleft .and. .not. clockw) then
         do i=beadrot+1,nbead
            xm(i) = xmo(i)
            ym(i) = -zmo(i)
            zm(i) = ymo(i)
         enddo
      else
          do i=beadrot+1,nbead
            xm(i) = xmo(i)
            ym(i) = zmo(i)
            zm(i) = -ymo(i)
          enddo
      endif     
   endif

! new coordinates
   do i=1,nbead
      x(i) = xm(i) + x(beadrot)
      y(i) = ym(i) + y(beadrot)
      z(i) = zm(i) + z(beadrot)   
   enddo
  
  overlap = 0 
!check for overlap
   do i = 1, nbead
      do j = i+1, nbead
      	 if (x(i) == x(j) .AND. y(i) == y(j) .AND. z(i) == z(j)) then
            do k=1,nbead
               x(k) = xo(k)
               y(k) = yo(k)
               z(k) = zo(k)
            enddo
            !write(6,*) 'rejected: overlap'
            didimove = .false.
            overlap = 1
            return
         endif
      enddo
   enddo

!calculate energy in new config.
if(overlap==0)then
   call totenergy(nbead,engn)  
   diff = engn - engo
   if (diff .le. 0.0) then
       engmovetot = engn
       didimove = .true.
   else if(exp(-diff/t) .gt. random()) then
       engmovetot = engn
       didimove = .true.
   else
      do i=1,nbead
         x(i) = xo(i)
         y(i) = yo(i)
         z(i) = zo(i)
      enddo
      !write(6,*) 'rejected: boltz=',exp(-diff/t),'acc=',acc
      didimove = .false.
       engmovetot = engo
      return
   endif
endif
  !write(6,*) 'accepted: overlap'  
  
! recenter
!  call recenter(nbead)
end subroutine pivot

 SUBROUTINE slither(nbead,engmovetot, t, didimove)

         IMPLICIT NONE

         INTEGER, ALLOCATABLE                         :: xo(:), yo(:), zo(:), rx_trial(:), ry_trial(:), rz_trial(:)
         INTEGER, DIMENSION(3)                        :: s
         INTEGER                                      :: i, j, head, tail, c, v2, pnt, saw
         INTEGER                                      :: nbead
         LOGICAL                                      :: didimove
         REAL                                         :: rn, t, engmovetot, engn, engo, P_acc

         allocate(xo(nbead))
         allocate(yo(nbead))
         allocate(zo(nbead))
         allocate(rx_trial(nbead))
         allocate(ry_trial(nbead))
         allocate(rz_trial(nbead))


         xo = x
         yo = y
         zo = z
         rx_trial = 0                                  ! Dum vector is used to take the coordinate of one chain only              
         ry_trial = 0                                  ! and work on it for slithering alrogithm
         rz_trial = 0
         call totenergy(nbead,engo)

         v2 = ceiling(random()*2)                                             ! Selecting one as head and other as a tail  
                                                       ! out of the twon ends od a chain
         if(v2==1)then
            head = 1
            tail = nbead
         else
            head = nbead
            tail = 1
         endif
       
         rn = random() 
       !  write(*,*)rn
         pnt  = ceiling(rn*6)
         s(1) = dir(1,pnt)
         s(2) = dir(2,pnt)
         s(3) = dir(3,pnt)                                   ! Randomly selecting one of the 12 (for FCC) directions possible.
         rx_trial(head) = x(head) + s(1)                  ! Storing the position of new head in the dum vector
         ry_trial(head) = y(head) + s(2)
         rz_trial(head) = z(head) + s(3)


           if(head==nbead)then
              rx_trial(1:nbead-1) = x(2:nbead)              ! If head = N_dp, then 2:N_dp part of old chain will become 1:N_dp-1 part of new chain 
              ry_trial(1:nbead-1) = y(2:nbead)
              rz_trial(1:nbead-1) = z(2:nbead)

           else

              rx_trial(2:nbead) = x(1:nbead-1)              ! If head = 1, then 1:N_dp-1 part of old chain will become 2:N_dp part of new chain 
              ry_trial(2:nbead) = y(1:nbead-1)
              rz_trial(2:nbead) = z(1:nbead-1)

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
      !        call totenergy(nbead,engn)
      !        P_acc  =exp(-(engn - engo)/t)
      !        if(random().lt.P_acc)then
      !           x = rx_trial
      !           y = ry_trial
      !           z = rz_trial
      !           didimove = .true.
      !           engmovetot = engn
      !        else
      !           didimove = .false.
      !        endif

                  x = rx_trial
                  y = ry_trial
                  z = rz_trial
                  call totenergy(nbead,engn)
                  didimove = .true.
                  engmovetot = engn
              if(exp(-(engn-engo)/t) .lt. random())then
                  x = xo
                  y = yo
                  z = zo
                  didimove = .false.
                  engmovetot = engo
              endif

           else
              didimove = .false.
           endif

!         rg2 = 0
!         rend = 0
!         call properties(nbead, rg2, rend)

         END SUBROUTINE slither
end module moves
