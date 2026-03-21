program ising
implicit none

  integer :: i,j,L,p,a,b,c,d,niter,time,mm,nn,N, earray(3)=(/26,30,34/) ,z
  real :: r,q,E,M,Mag,Ei,Ef,dE,u,h

  real:: T=3.3 , J_ising=1.00  !bascially Kb=1
  integer , dimension(:,:) , allocatable :: spin
  
  
   
  print *, 'Enter num of iterations'
  read * , niter
  
do z=1,size(earray)
    L=earray(z)

    if (allocated(spin)) then
        deallocate(spin)
    end if

    allocate(spin(L,L))
    E=0.0 !inistantaneoous energy of lattice
    M=0.0 !inst magnetization of lattice
    N=L*L  !total no f spins in lattice

    call random_seed

    !initialise lattice


    
    do i=1,L 
        do j=1,L 
            call random_number(r)
            spin(i,j)= +1
            
            ! if (r<0.5) then 
            !    spin(i,j)=-1
            ! else 
            !    spin(i,j)=+1
            ! end if
            
            

        end do 
    end do     




    !Calculate initial magnetization and energy

    do i=1,L 
        do j=1,L 
        a = mod(i, L) + 1
        b = mod(i-2+L, L) + 1
        c = mod(j, L) + 1
        d = mod(j-2+L, L) + 1  !getting the nearest neigbours;
            

            M=M+spin(i,j)
            E=E-J_ising*(spin(i,j)*(spin(i,d)+spin(i,c)+spin(a,j)+spin(a,d)+spin(b,j)+spin(b,c)))
            
        end do
    end do

    mag=M/(float(N))  !N=L*L: magnetization (instantaeous) per spin
    E=E*0.50d0   !to counteract double counting

    print*,'Initial energy E , Energy per spin : ' , E ,E/float(N)
    print*,'Initial magnetization M , M per spin : ' , M, M/float(N)

    !Done with init
    !--------------------------------------------------------

    !Evolve to reach equilibrium

    open(10,file='ising_a.dat')
    do time = 1,niter  !loop over MCS

        do mm=1,L 
            do nn= 1,L 
                call random_number(r) ;  i=int(r*float(L))+1 !choosing lattice site from 1 to L
                call random_number(r) ;  j=int(r*float(L))+1

            
                a = mod(i, L) + 1
                b = mod(i-2+L, L) + 1
                c = mod(j, L) + 1
                d = mod(j-2+L, L) + 1  !getting the nearest neigbours;
                
            

                Ei=-J_ising*(spin(i,j)*(spin(i,d)+spin(i,c)+spin(a,j)+spin(a,d)+spin(b,j)+spin(b,c)))

                spin(i,j)=-spin(i,j) !trial flip

                Ef=-J_ising*(spin(i,j)*(spin(i,d)+spin(i,c)+spin(a,j)+spin(a,d)+spin(b,j)+spin(b,c)))
                dE=Ef-Ei

                if (dE<0) then
                    E=E+dE  !reverseing the energy totally for the system 
                    M=M+(2.0*float(spin(i,j)))

                else 
                    u=exp(-dE/T)
                    call random_number(h)

                    if (h<u)then
                        E=E+dE  !reversing the energy totally for the system 
                        M=M+(2.0*float(spin(i,j)))
                    else
                        spin(i,j)=-spin(i,j) !not accepted-not updated

                    end if 
                end if  
            end do  
        end do  
        
        
        write(10,*) time , M/float(N) , E/float(N)  !writing E and M with no of iterations

    end do

end do  

close(10)   
            
            

end program ising