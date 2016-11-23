!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!These subroutines allow for an initialization Monte Carlo simulation of 
!a polymer chain to be performed
!Cases for parallel tempering (replica exchange) and without are included
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine wlcsim_init(wlc_p,wlc_d)
  
  use params
  use mersenne_twister

  implicit none

  !Structures for simulation
  type(wlcsim_params) :: wlc_p
  type(wlcsim_data) :: wlc_d
  
  !MPI status variables
  integer (kind = 4)  status(MPI_STATUS_SIZE)
  integer (kind = 4) error
  !Integers for commuinication among replicas
  integer rep
  integer dest
  integer source
  integer LK
  real(dp) Wr
  real(dp) eelas(4)
  integer repSTART,repEND
  real(dp) deEXCHANGE
  !Counting variables
  integer i,repIND

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Run  initialization Monte Carlo simulation
  !Cases:
  !parallel tempering on
  !parallel tempering off
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !If parallel tempering is on
  if (wlc_p%ptON) then

     !Cases:
     !Head node (id = 0)
     !Worker node (id /= 0)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Head node
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (wlc_d%id.eq.0) then
       
        !During replica exchange, replicas alternate between switching with
        !replica above and replica below. Initially, set the replica looking above
        !to the be first one
        repSTART = 1
        repEND = wlc_p%nLKs - 1

        !Begin replica exchange loop
        do repIND = 1,wlc_p%nREPLICAexchange
           
           !Send out Lk to the worker nodes to begin their simulations
           do rep = 1,wlc_p%nLKs
              LK = wlc_p%LKs(rep)
              dest = wlc_d%nodeNUMBER(rep)
              call MPI_Send (LK,1, MPI_INTEGER, dest,   0, &
                   MPI_COMM_WORLD,error )
           enddo

           !LK, Wr, and EELAS from the worker nodes to perform
           !replica exchange

           do rep = 1,wlc_p%nLKs
              LK = wlc_p%LKs(rep)
              source = wlc_d%nodeNUMBER(rep)
              call MPI_Recv(Wr,1, MPI_DOUBLE_PRECISION, source, 0, &
                   MPI_COMM_WORLD,status,error )
              call MPI_Recv(eelas,4, MPI_DOUBLE_PRECISION, source, 0, &
                   MPI_COMM_WORLD,status,error )
              wlc_d%Wrs(rep) = Wr
              wlc_d%eelasREPLICAS(rep,:) = eelas
           enddo

           !peform replica exchange
           call replicaEXCHANGE(wlc_p,wlc_d)
           

        enddo

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Worker nodes
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else
        !Begin replica exchange loop
        do repIND = 1,wlc_p%nREPLICAexchange
           !Receive LK from head node to begin simulation
           source = 0
           call MPI_Recv(LK,1, MPI_INTEGER, source,   0, &
                MPI_COMM_WORLD,status,error )
           wlc_p%LK = LK

           !Run a monte carlo simulation for NSTEPS
           call MCsim(wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%nINIT,wlc_p%intON,wlc_p%eb,&
                wlc_p%epar,wlc_p%eperp,wlc_p%gam,wlc_p%eta, wlc_p%lhc,wlc_p%vhc,wlc_p%lbox, &
                wlc_d%MCAMP,wlc_d%SUCCESS,wlc_d%successTOTAL,wlc_p%MOVEON,&
                wlc_d%WINDOW,wlc_p%RING,wlc_p%TWIST,wlc_p%Lk,wlc_p%LT,wlc_p%LP,wlc_p%L,wlc_d%rand_stat)
           
         
           !Recalculate structural quantities and energies
           call writhe(wlc_d%R,wlc_p%nB, wlc_d%Wr)
           call energy_elas(wlc_d%eelas,wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%eb,wlc_p%epar, &
                wlc_p%eperp,wlc_p%gam,wlc_p%eta,wlc_p%RING,wlc_p%twist,wlc_p%LK,wlc_p%lt,wlc_p%lp,wlc_p%l)
                    

           !Communicate with the head node for replica exchange
         
           !Send back writhe and elastic energy back to the head node
           call MPI_SEND(wlc_d%Wr,1, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,error )
           call MPI_SEND(wlc_d%eelas,4, MPI_DOUBLE_PRECISION, source, 0, &
                MPI_COMM_WORLD,error )
           
        enddo
     endif

  else


     !Run a MC simulation for NSTEPS
     call MCsim(wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%nINIT,wlc_p%intON,wlc_p%eb,&
          wlc_p%epar,wlc_p%eperp,wlc_p%gam,wlc_p%eta, wlc_p%lhc,wlc_p%vhc,wlc_p%lbox, &
          wlc_d%MCAMP,wlc_d%SUCCESS,wlc_d%successTOTAL,wlc_p%MOVEON,&
          wlc_d%WINDOW,wlc_p%RING,wlc_p%TWIST,wlc_p%Lk,wlc_p%LT,wlc_p%LP,wlc_p%L,wlc_d%rand_stat)
  end if
end subroutine wlcsim_init
        


  
  


