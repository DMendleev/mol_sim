MODULE sumup

   USE variables
   
   IMPLICIT NONE
   PRIVATE 
   SAVE
   
   public:: totenergy,PROPERTIES
   
   contains

!*****************************************************************************************!
!                                   CALCULATING ENERGY                                    !
!*****************************************************************************************!
   subroutine totenergy(totbead,toten)
      
      INTEGER      :: totbead, i, j, midbead, rx
      REAL         :: toten, period, c, Amp, midbox
     
      Amp   = -20.0 
      toten = 0.0
      
      if (mono) then

          toten = 0.0

      !For diblock energy divide polymer into two types of beads (A and B) and have them interact with a sine-wave potential
      else if (diblock) then
              midbead = totbead/2
              midbox  = real(L-1)/2
              period = 2*pi/(totbead -1)
  
              do i=1, totbead

                 if(i<=midbead)then
                    c = real(x(i))/real(L)
                    rx = x(i) - floor(c)*L
                    toten = toten + Amp*sin(rx*period)
                 else
                    c = real(x(i))/real(L)
                    rx = x(i) - floor(c)*L
                    toten = toten - Amp*sin(rx*period)
                endif

              enddo
              
      endif

      return

   end subroutine totenergy

!*****************************************************************************************!
!                                CALCULATING PROPERTIES                                   !
!*****************************************************************************************!

     SUBROUTINE PROPERTIES(nbead, Rg2, Rg2_para, Rg2_perp, Rend, Rend_para, Rend_perp)

         IMPLICIT NONE

         REAL, INTENT(OUT)                             :: Rg2, Rend
         REAL, INTENT(OUT)                             :: Rg2_para, Rg2_perp, Rend_para, Rend_perp
         INTEGER                                       :: i, j, nbead
         REAL,DIMENSION(3)                             :: R_cm, vec

         R_cm        = 0.0
         Rg2         = 0.0
         Rg2_para    = 0.0
         Rg2_perp    = 0.0
         Rend        = 0.0
         Rend_para   = 0.0
         Rend_perp   = 0.0


         R_cm(1) = sum(x(:))/real(nbead)
         R_cm(2) = sum(y(:))/real(nbead)
         R_cm(3) = sum(z(:))/real(nbead)
 
          do j=1,nbead
             Rg2 = Rg2 + (x(j)-R_cm(1))**2 &
                       + (y(j)-R_cm(2))**2 &
                       + (z(j)-R_cm(3))**2
             
             Rg2_para = Rg2_para + (x(j)-R_cm(1))**2
 
             Rg2_perp = Rg2_perp + ((y(j)-R_cm(2))**2 &
                                 + (z(j)-R_cm(3))**2)/2
          end do

         Rg2        = Rg2/real(nbead)
         Rg2_para   = Rg2_para/real(nbead)
         Rg2_perp   = Rg2_perp/real(nbead)

         vec(1)= x(nbead)-x(1)
         vec(2)= y(nbead)-y(1)
         vec(3)= z(nbead)-z(1)

         Rend      = dot_product(vec,vec)
         Rend_para = vec(1)*vec(1)
         Rend_perp = (vec(2)*vec(2) + vec(3)*vec(3))/2

         return

     END SUBROUTINE PROPERTIES


end module sumup
