module variables
real(8),parameter :: pi = 3.141592653589793
integer,allocatable::x(:),y(:),z(:)
real,allocatable::jij(:,:)
character*50,dimension(20)::aname
integer::xcom0,ycom0,zcom0,icom,nseed, L
real::pmpivot,pmreptation
logical::mono,diblock,triblock
integer, dimension(3,6) :: dir
save
end module variables
