!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!wlc_ring simulation main
!This program performs a Monte Carlo simulation of a polymer modeled
!as a wormlike chain with topological constraints. The polymer can be
!a ring or linear, with the possibility for a linking constraint and
!an associated twist energy.

!The simulation can be performed with a single processor and a single
!linking number, or parallel tempering with respect to linking number
!can be performed using openMPI

!Brad Krajina and Andrew Spakowitz
!Last edited : 2016/11/22

program main

  use params
  use mpi
  implicit none

  !Simulation state
  type(wlcsim_params) :: wlc_p           !Structure containing static simulation parameters
  type(wlcsim_data) :: wlc_d             !Structure containing dynamic simulation data

  !I/O
  character(1024) infile !Input file
  !Counting variables
  integer ind,i

  !   variable for random number generator 
  integer Irand     ! Seed
  character*8 datedum  ! trash
  character*10 timedum ! trash
  character*5 zonedum  ! trash
  integer seedvalues(8) ! clock readings
 


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Begin reading inputs and initialization
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(infile,*)   'input/input'
  
  ! !Reaad and set simulation parameters
   call get_input_from_file(infile, wlc_p)

  !Initialize values
  call initialize_wlcsim_data(wlc_d, wlc_p)

  !Allocate arrays
  call MCvar_allocate(wlc_p,wlc_d)


  !Set initial chain configuration

  call initcond(wlc_d%R,wlc_d%U,wlc_p%nT,wlc_p%nB,wlc_p%nP,wlc_p%FRMfile,&
       wlc_p%gam,wlc_p%lbox,wlc_p%ring,wlc_d%rand_stat)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Perform an initialization Monte Carlo simulation
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (wlc_d%id.eq.0) then
     print *, 'run initialization monte carlo'
  endif

  call wlcsim_init(wlc_p,wlc_d)

  if (wlc_d%id.eq.0) then
     print *, 'finished initialization monte carlo'
     print *, '************************************'
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Loop over Monte Carlo simulations
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ind = 1,wlc_p%indMAX
     if (wlc_d%id.eq.0) then
        print *, '************************'
        print *, 'ind ', ind
        print *, '************************'
     endif
     
     !Perform the Monte Carlo simulation
     call wlcsim(wlc_p,wlc_d)
     
     !Save the polymer configuration
     call save_configuration(wlc_p,wlc_d,ind)
     
  enddo

  !Save statistics for structural quantities
  
  call save_stats(wlc_p,wlc_d)

  !Save the swapping probabilities for replicas
  if (wlc_d%id.eq.0)then
     open(unit =1, file = 'data/pSWAPup')
     do i = 1,wlc_p%nLKs - 1
        write(1,*) dble(wlc_d%nSWAPup(i))/dble(wlc_d%nTRIALup(i))
     enddo
     close(1)
  endif

  ! If parallel tempering is on, finalize MPI to exit properly

  if (wlc_p%ptON) then
     call MPI_Finalize(wlc_d%error)
     if (wlc_d%error.ne.0) then
        print*, "MPI_Finalize", wlc_d%error
     endif
              
  endif
  
end program main
  

  
