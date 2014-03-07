module sumup
use variables
implicit none
private 
save
public:: totenergy,PROPERTIES

contains
   subroutine totenergy(totbead,toten)
      
      integer::totbead,i,j,midbead, rx
      real::toten, period, c, Amp
     
      Amp = 1.0 
      toten = 0.0
!Energy calculations not even needed for mono
      if (mono) then

          toten = 0.0

!For diblock energy divide polymer into two types of beads (A and B) and have them interact with a sine-wave potential
      else if (diblock) then
              midbead = totbead/2
              period = 2*pi/(totbead -1)
  
              do i=1, totbead

                 if(i<=midbead)then
                    c = real(x(i))/real(L)
                    rx = x(i) - floor(c)*L
                    if (rx<L/2)then
                        toten = toten + Amp*sin(x(i)*period)
                    else
                        toten = toten - Amp*sin(x(i)*period)
                    endif
                 
                 else

                    c = real(x(i))/real(L)
                    rx = x(i) - floor(c)*L
                    if (rx<L/2)then
                        toten = toten - Amp*sin(x(i)*period)
                    else
                        toten = toten + Amp*sin(x(i)*period)
                    endif
                endif
!                write (*,*) sin(x(i)*period)  
              enddo
              
       endif 
      return
   end subroutine totenergy
!*****************************************************************************************!
!                       CALCULATING PROPERTIES                                            !
!*****************************************************************************************!

     SUBROUTINE PROPERTIES(nbead,Rg2,Rend)


         IMPLICIT NONE

         REAL, INTENT(OUT)                                    :: Rg2, Rend
         REAL, DIMENSION(3)                                   :: R_cm, R_cm_A, R_cm_B
         INTEGER                                              :: i, j, nbead
         REAL,DIMENSION(3)                                    :: vec

         R_cm   = 0.0
         Rg2    = 0.0
         Rend   = 0.0


            R_cm(1) = sum(x(:))/real(nbead)                                      ! Center of Mass of chain i  
            R_cm(2) = sum(y(:))/real(nbead)
            R_cm(3) = sum(z(:))/real(nbead)

             do j=1,nbead
                Rg2 = Rg2 + (x(j)-R_cm(1))**2 &
                          + (y(j)-R_cm(2))**2 &
                          + (z(j)-R_cm(3))**2
             end do

         Rg2   = Rg2/real(nbead)

         vec(1)= x(nbead)-x(1)
         vec(2)= y(nbead)-y(1)
         vec(3)= z(nbead)-z(1)



         Rend  = dot_product(vec,vec)

         return
     END SUBROUTINE PROPERTIES

!*****************************************************************************************!

end module sumup
