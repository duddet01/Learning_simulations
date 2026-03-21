program ising
implicit none

  integer :: i,j,L,p,a,b,c,d,niter,time,mm,nn,N,T_temp 
  real(8) :: r,q,E,M,mag,Ei,Ef,dE,u,h,mag1
  real(8) :: T, J_ising=1.00  !bascially Kb=1
  real(8) :: av_M , av_E , cv , av_E2 ,chi , av_M2 , av_M_N , av_E_N ,nmeans ,av_M1  
  real(8) :: av_M4 , binder
  integer , dimension(:,:) , allocatable :: spin
  integer :: seed, n_equil , n_stat
  !integer :: random_seedcharacter(len=30) :: charac_a , charac_b  ! char_b stores the name dump_pos

  seed=44859
  !charac_b= 'store_config'

  print *, 'Enter the no of lattice points in 1D'
  read *, L  
  print *, 'Enter no of iterations at each temp: '
  read * , niter

  allocate(spin(L,L))
  E=0.0 !inistantaneoous energy of lattice
  M=0.0 !inst magnetization of lattice
  N=L*L  !total no f spins in lattice
  N_equil=10000   !collect data after some equilibration time
  N_stat=10 !collect data after 10 steps
  call random_seed

!---------------------------------------------------------

  !initialise lattice

  !open(77,file='initialise_ising.dat')
  p=0
  do i=1,L 
     do j=1,L 
        call random_number(r)
        !spin(i,j)= -1
        
        if (r<0.5) then 
           spin(i,j)=-1
        else 
           spin(i,j)=+1
        end if
        
        write(77,*)  float(i),float(j),float(p),float(spin(i,j))

     end do 
 end do     

 !close(77)


!-----------------------------------------

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
E=E*0.5d0   !to counteract double counting

print*,'Initial energy E , Energy per spin : ' , E ,E/float(N)
print*,'Initial magnetization M , M per spin : ' , M, M/float(N)

 !Done with init
 !--------------------------------------------------------

 !Evolve to reach equilibrium

 !open(11,file='ising_T2_N40.dat')


open(17,file='ising_with_T.dat')

do T_temp=400,300,-1
    T=dfloat(T_temp)/100.0d0

    av_M=0.0d0 ; av_E=0.0d0  !for abs value
    av_M1=0.0d0  !for avg value fluctuation
    av_M_N=0.0d0 ; av_E_N=0.0d0
    av_E2=0.0d0 ; av_M2=0.0d0
    av_M4 = 0.0d0


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
                    E=E+dE  !reverseing the energy totally for the system 
                    M=M+(2.0*float(spin(i,j)))
                else
                    spin(i,j)=-spin(i,j) !not accepted-not updated

                end if 
            end if  
        end do  
    end do  
    
    
    !write(11,*) time , M/float(N) , E/float(N)  !writing E and M with no of iterations


!close(11)   

!After reaching equilibrium , collect data:

if (time>n_equil) then
    if (mod(time,n_stat).eq.0) then

        mag  = abs(M)/dfloat(N)
        mag1 = M/dfloat(N)

        av_M  = av_M + mag
        av_E  = av_E + E/dfloat(N)
        av_M1 = av_M1 + mag1

        av_M_N = av_M_N + M   !for entire lattice
        av_E_N = av_E_N + E

        av_M2 = av_M2 + M*M
        av_M4 = av_M4 + M*M*M*M
        av_E2 = av_E2 + E*E

    end if
end if


end do  ! time - 1 to niter
nmeans=dfloat((niter-n_equil)/n_stat)
av_M=av_M/(nmeans) ;av_E=av_E/(nmeans)
av_M1=av_M1/(nmeans)  !per spin E and M averaged
av_M2=av_M2/(nmeans) ;av_E2=av_E2/(nmeans) 
av_M_N=av_M_N/(nmeans) ;av_E_N=av_E_N/(nmeans)
av_M4 = av_M4/(nmeans)

cv= (av_E2-av_E_N**2)/T**2
chi= (av_M2-av_M_N**2)/T
binder = 1.0d0 - av_M4/(3.0d0*(av_M2**2))


write(17,*) T,av_M,av_E,cv,chi,av_M1,binder !Here av_M is abs M and av_m1 is average M


end do
        

close(17)


            

end program ising
