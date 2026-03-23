program euler
implicit none

real(8):: y, k1,k2,k3,k4 ,x,dx=0.01
integer :: i ,niter

!initialisation::

x=0.0d0 ; y=0.0d0 ; niter= 1.550d0/dx
open (71,file="rk4.dat")

do i=1,niter+1
    k1=1.0d0+y*y
    k2=1.0d0+(y+dx/2.0d0*k1)**2.0d0
    k3=1.0d0+(y+dx/2.0d0*k2)**2.0d0
    k4=1.0d0+(y+dx*k3)**2.0d0

    y=y+dx/6.0d0*(k1+2.0d0*k2+2.0d0*k3+k4)


    write (71,*) i , dfloat(i)*dx , y
    
end do 

close(71)


end program euler    