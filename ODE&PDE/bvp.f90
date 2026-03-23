program diff
implicit none

real*8 , parameter :: f_x=0.0d0 , end_x=1.0d0 
real*8 ,parameter :: f_y=0.0d0 , end_y=2.0d0
real*8 , parameter :: dx=0.01d0
integer ,parameter :: nop = int((end_x-f_x)/dx)

integer :: i,j,k,ll,cond
real*8 :: x(nop) , y(nop) , y_old(nop) , limit
real*8 , parameter :: dum1=1.0d0/(2.0d0-10.0d0*dx*dx) , dum2=1.0d0-2.50d0*dx , dum3=1.0d0+2.50d0*dx ,dum4 =-10.0d0*dx*dx

write(*,*) 'no of grid points' , nop
!initialization

limit=0.00010d0 !tolerance parameter
ll=0 !no. of iter to reach convergence
cond=0

open(80,file='diff.dat')

!giving the boundary conditions:
x(1)=0.0d0 ; x(nop)=1.0d0 ; y(1)=0.0d0; y(nop)=2.0d0

do i=2,nop-1
    x(i)=x(i-1)+dx
    y(i)=((end_y-f_y)/(end_x-f_x))*x(i)

end do 

do 
  ll=ll+1
  if (cond==1) exit
     y_old=y
     do i=2,nop-1
        y(i)=dum1*(dum2*y_old(i+1)+dum3*y(i-1)+dum4*x(i))

     end do 

     cond=1
     do i=2,nop-1
        if (abs(y_old(i)-y(i)).ge.limit) cond=0
     end do 
     
end do 

write(*,*) 'No. of iterations to achieve convergence =' ,ll

do i=1,nop
    write(080,*) x(i) , y(i)
    if ((nint(x(i)*10000.0d0)/10000.0d0).eq.0.8d0) then
        write(*,*) 'Value of y at x=0.80 is :' ,y(i)
    end if    

end do 

close(80)

end program diff


