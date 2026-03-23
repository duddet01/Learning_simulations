program couple
implicit none

real(8)::  k1x,k2x,k3x,k4x ,x,dx=0.01
real(8):: xtemp1 , xtemp2 , xtemp3 ,vtemp1 , vtemp2,vtemp3
real(8) :: k1v,k2v,k3v,k4v ,v, Time
integer :: i ,niter 

!initialisation::

x=0.0d0 ; v=2.40d0 ; Time=50.0d0; niter= int(Time/dx)
open (71,file="couple2.dat")

do i=1,niter+1
    k1x=v
    xtemp1=x+dx*0.50d0*k1x
    k1v=-sin(x)
    vtemp1=v+dx*0.50d0*k1v

    k2x=vtemp1
    xtemp2=x+dx*0.50d0*k2x
    k2v=-sin(xtemp1)
    vtemp2=v+dx*0.50d0*k2v

    k3x=vtemp2
    xtemp3=x+dx*k3x
    k3v=-sin(xtemp2)
    vtemp3=v+dx*k3v

    k4x=vtemp3
    k4v=-sin(xtemp3)
    

    x=x+dx/6.0d0*(k1x+2.0d0*k2x+2.0d0*k3x+k4x)
    v=v+dx/6.0d0*(k1v+2.0d0*k2v+2.0d0*k3v+k4v)


    write (71,*) i , dfloat(i)*dx , x , v    
end do 

close(71)


end program couple    