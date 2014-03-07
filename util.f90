module util
use variables
implicit none
private 
save
public:: writexyz, determineCOM, recenter 

contains
   subroutine writexyz(iounit,nbead,msg)
      integer::i,iounit,nbead
      integer::msg

      write(iounit,'(I4)') nbead
      write(iounit,*) 'step ',msg
      do i=1,nbead
         !write(6,*) moltyp(i),aname(moltyp(i)
         if(i<=nbead/2)then
            write(iounit,'(A10,3I4)') 0,x(i),y(i),z(i)
         else
            write(iounit,'(A10,3I4)') 1,x(i),y(i),z(i)
         endif 
     enddo      
      return
   end subroutine writexyz

   subroutine determineCOM(nbead)
      integer::nbead
      icom = (nbead+1)/2
      xcom0 = x(icom)
      ycom0 = y(icom)
      zcom0 = z(icom)
   end subroutine determineCOM
   
   subroutine recenter(nbead)
      integer::i,nbead
      integer::xcom1,ycom1,zcom1
      xcom1=x(icom)
      ycom1=y(icom)
      zcom1=z(icom)
      do i=1,nbead
         x(i) = x(i) - (xcom1-xcom0)
         y(i) = y(i) - (ycom1-ycom0)
         z(i) = z(i) - (zcom1-zcom0) 
      enddo

   end subroutine recenter
end module util
