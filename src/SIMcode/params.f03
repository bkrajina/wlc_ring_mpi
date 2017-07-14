
module params
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: IEEE_ARITHMETIC
  use mersenne_twister
  use mpi
  IMPLICIT NONE

!  private
  public :: pi, wlcsim_params, wlcsim_data, get_input_from_file

  ! hardcoded params. will need to change if certain parts of code change
  integer, parameter :: nMoveTypes = 6 ! number of MC move types

  ! precision of simulations
  integer, parameter :: dp = real64 ! preferred over SELECTED_real_Kind(15,307)

  ! universal constants
  real(dp) :: pi = 4 * atan(1.0_dp) ! fully accurate, adaptive precision

  ! for all parameters that cannot change during individual simulations
  ! these are documented more thoroughly where they are read in (see the
  ! subroutine get_input_from_file), in the docs (TODO), and often the default
  ! values will help with understanding how the variable is used.

  type wlcsim_params

     !   Simulation parameters
     integer nT                ! Total number of beads  NT=NP*N*G
     integer nB                ! Number of beads in a polymer NB=N*G
     integer nP                ! Number of polymers
     real(dp) maxe2e   ! max distance between end to end points
     real(dp) mine2e   ! min distance between end to end points
     real(dp) lp       ! persistence length
     real(dp) lt       ! twist persistence length
     real(dp) b1       ! binding site 1
     real(dp) b2       ! binding site 2
     real(dp) l0       ! Equilibrium segment length (same as gam)
     real(dp) l        ! actual length of polymer
     real(dp) eb     ! Energy of bending
     real(dp) epar   ! energy "parallel" i.e. psring energy
     real(dp) eperp  ! energy "perpendicular" i.e. shear energy
     real(dp) gam    ! average equilibrium interbead spacing
     real(dp) eta    ! bend-shear coupling parameter
     real(dp) eps      ! number of kuhn lengths between beads
     real(dp) del      ! number of persistence lengths between beads
     real(dp) lhc      !hard-core diameter of lennard-jones repulsion
     real(dp) vhc      !hard-core repulsion strength
     integer LK                ! Linking number
     integer nLKs      !Number of linking number replicas
     integer, allocatable, dimension(:) ::  LKs    !Vector of linking numbers for replicas

     !   Monte Carlo Variables (for adaptation)
     integer movetypes
     real(dp) PDesire(nMoveTypes) ! desired hit rate
     real(dp) MAXWindoW(nMoveTypes)         ! Max Size of window for bead selection
     real(dp) MINWindoW(nMoveTypes)         ! Min Size of window for bead selection
     real(dp) MinAMP(nMoveTypes) ! minium amplitude
     real(dp) MaxAMP(nMoveTypes) ! maximum amplitude
     integer MOVEON(nMoveTypes)         ! Is the move active
     integer NADAPT(nMoveTypes) ! Nunber of steps between adapt

     !   Timing variables
     integer nINIT              ! number of initialization Monte Carlo steps to perform
     integer NPT                ! number of steps between parallel tempering
     integer indMAX             ! total number of save points
     integer NSTEP              ! steps per save point (no replica exchange) or steps between replica exchange
     integer nREPLICAexchange    !Number of replica exchanges between save points
     
    !   Switches
     integer ring              ! whether the polymer is a ring
     integer twist             ! whether to include twist (MC only for now)
     integer restrectedr       ! whether to include a restriction between site distances
     logical ptON              ! is parallel tempering on?
     integer intON             ! self interactions on ?
     logical FRMfile           ! read initial condition R from file
     logical FRMField          ! read initial field from file
     logical restart     ! whether we are restarting from a previous sim or not

    !   box size things
     real(dp) lbox(3)  ! Box length (approximate)

     
  end type wlcsim_params

     ! for variables that can change during the simulation
  type wlcsim_data

     real(dp), allocatable, dimension(:,:):: R   ! Conformation of polymer chains
     real(dp), allocatable, dimension(:,:):: U   ! Conformation of polymer chains
     real(dp), allocatable, dimension(:,:):: RP !Test Bead positions - only valid from IT1 to IT2
     real(dp), allocatable, dimension(:,:):: UP !Test target vectors - only valid from IT1 to IT2

     !Structural quantities
     real(dp) Wr                                 !Writhe of current structure
     
     !   Monte Carlo Variables (for adaptation)
     real(dp) MCAMP(nMoveTypes) ! Amplitude of random change
     integer  WindoW(nMoveTypes)         ! Size of window for bead selection
     integer SUCCESS(nMoveTypes)        ! Number of successes
     integer successTOTAL(nMoveTypes)               !Total number of successes
     real(dp) PHit(nMoveTypes) ! hit rate

     !   Energys
     real(dp) Eint     ! running interaction energy
     real(dp) EELAS(4) ! Elastic energies

     !Storged simulation trajectories
     real(dp), allocatable, dimension(:) :: wrTRAJ !vector of all writhe for the simulation
     real(dp), allocatable, dimension(:) :: rgsqTRAJ !vector of squared radii of gyration for the simulation

     !   Parallel tempering variables
     integer rep  ! which replica am I
     integer id   ! which thread am I
     integer error  ! MPI error
     integer p      !number of MPI processes
     real(dp), allocatable, dimension(:) :: Wrs !Vector of writhe for each replica
     real(dp), allocatable, dimension(:,:) :: eelasREPLICAS !elastic energies of replicas
     integer replicaSTART !index for replica to start with for exchange loop
     integer replicaEND   !index for replica to end at for exchange loop
     integer, allocatable, dimension(:) :: nTRIALup !number of times this replica has attempted to swap with replica above
     integer, allocatable, dimension(:) :: nTRIALdown !number of times this replica has attempted to swap with replica below
     integer, allocatable, dimension(:) :: nSWAPup !number of times this replica has swapped with replica above
     integer, allocatable, dimension(:) :: nSWAPdown !number of times this replica has swapped with replica below

     !Random number generator
     type(random_stat) rand_stat  !status of random number generator 
     integer, allocatable, dimension(:) :: nodeNUMBER !vector of replicas indices for nodes
     
  end type wlcsim_data
contains 
  subroutine set_param_defaults(wlc_p)
     implicit none
     type(wlcsim_params), intent(inout) :: wlc_p
     ! file IO
     wlc_p%FRMfile=.FALSE.      ! don't load initial bead positions from file
     wlc_p%restart=.FALSE.      ! don't restart from previously saved simulation
     
     ! geometry options
     wlc_p%NP  =1               ! one polymer
     wlc_p%nB  =200             ! 200 beads per polymer
     wlc_p%lp = 1                ! units of lp by default
     wlc_p%lt = 1                ! twist persistence length equals persistence length by default
     wlc_p%lk=0                  ! no linking number (lays flat) by default
     wlc_p%nLKs = 1              ! number of linking number replicas equal to 1
     wlc_p%lbox(1)=25.0_dp      ! arbitrary box size
     wlc_p%lbox(2)=25.0_dp
     wlc_p%lbox(3)=25.0_dp

     ! options
     wlc_p%movetypes=nMoveTypes
     wlc_p%ring=0    ! not a ring by default
     wlc_p%twist=0    ! don't include twist by default
     wlc_p%intON = 0   !don't include interactions

     ! polymer steric repulsion interaciton params
     wlc_p%lhc = 2.0_dp     !diameter of DNA by default
     wlc_p%vhc = 1.0_dp     !steric repulsion strength equal to thermal energy

     ! timing options
     wlc_p%nINIT=10000  ! number of initialization Monte Carlo steps to perform
     wlc_p%NStep=10000  ! number of steps to take per save point or replica exchnage
     wlc_p%indMAX=100    ! number of save points
     wlc_p%nREPLICAexchange = 1 !number of replica exchanges between save point
     ! replica options
     wlc_p%PTON=.FALSE.  !don't parallel temper by default
     wlc_p%NPT=100      ! 100 steps between parallel tempering is pretty frequent
     
     !move options
     wlc_p%moveON = 1  !Set moves on
     wlc_p%moveON(5) = 0
     wlc_p%moveON(6) = 0    !By default, set total chain moves off
     

   end subroutine set_param_defaults

subroutine read_from_file(infile, wlc_p)
    use INPUTparaMS, only : readLINE, readA, readF, readI, reado
    implicit none
    type(wlcsim_params), intent(inout) :: wlc_p
    character(1024), intent(in) :: infile
    LOGICAL :: fileend = .FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! parameter name currently being read in
    integer :: NITEMS ! number of items on the line in the parameter file
    integer pf
    pf = 1
    open(unit=PF,file=trim(adjustL(infile)),status='OLD')

    ! read in the keywords one line at a time
    do
       CALL readLINE(PF,fileend,NITEMS)
       if (fileend.and.nitems.eq.0) EXIT

       ! skip empty lines
       if (NITEMS.EQ.0) cycle

       ! read in the keyword for this line
       CALL readA(WORD,CASESET=1)

       ! Skip any empty lines or any comment lines
       if (WORD(1:1).EQ.'#') cycle

       SELECT CASE(WORD) ! pick which keyword
       CASE('FRMfile')
           call reado(wlc_p%FRMfile) ! read configuration from file
       CASE('TWIST')
           CALL readI(wlc_p%twist) ! whether to include twist energies in MC
       CASE('RING')
           CALL readI(wlc_p%ring) ! whether polymer is a ring or not
       CASE('LK')
           CALL readI(wlc_p%lk) ! linking number
       CASE('PTON')
           CALL reado(wlc_p%PTON) ! parallel Tempering on
       CASE('NB')
           Call readI(wlc_p%nb)  ! actual number of beads we want to simulate
       CASE('L')
           Call readF(wlc_p%l)  ! actual length in AU of polymer we want to simulate
       CASE ('RESTRICTEDR')
           Call readI(wlc_p%restrictedr) ! end to end distance restriction on or off
       CASE ('MAXEND2END')
           Call readF(wlc_p%maxe2e)  ! maximum distance between sites
       CASE ('MINEND2END')
           Call readF(wlc_p%mine2e)  ! minimum distance between sites
       CASE('LT')
           Call readF(wlc_p%lt)  ! twist persistence length
       CASE('LP')
           Call readF(wlc_p%lp)  ! persistence length
       CASE('LHC')
           Call readF(wlc_p%lp)  ! hard-core diameter for lennard jones 
       CASE('LBOX')
           Call readF(wlc_p%lbox(1)) ! side length of box
           wlc_p%lbox(2)=wlc_p%lbox(1)
           wlc_p%lbox(3)=wlc_p%lbox(1)
       CASE('LBOXX')
           Call readF(wlc_p%lbox(1)) ! side length of box in x direction
       CASE('LBOXY')
           Call readF(wlc_p%lbox(2)) ! side length of box in y direction
       CASE('LBOXZ')
           Call readF(wlc_p%lbox(3)) ! side length of box in z direction
       CASE('NP')
           CALL readI(wlc_p%NP)  ! Number of polymers
       CASE('INDMAX')
           Call readI(wlc_p%indMAX) ! total number of save points
       CASE('NSTEP')
           Call readI(wlc_p%NStep) ! steps per save point
       CASE('NINIT')
           Call readI(wlc_p%nINIT) ! number of initialization monte carlo steps
       CASE('NREPLICAEXCHANGE')
           Call readI(wlc_p%nREPLICAexchange) ! number of replica exchanges between save point
       CASE('NPT')
           Call readI(wlc_p%NPT) ! number of steps between parallel tempering
       CASE('CRANK_SHAFT_ON')
           Call readI(wlc_p%MOVEON(1)) ! is Crank shaft move on 1/0
       CASE('SLIDE_ON')
           Call readI(wlc_p%MOVEON(2)) ! is Slide move on 1/0
       CASE('PIVOT_ON')
           Call readI(wlc_p%MOVEON(3)) ! is Pivot move on 1/0
       CASE('ROTATE_ON')
           Call readI(wlc_p%MOVEON(4)) ! is single bead rotate on 1/0
       CASE('FULL_CHAIN_ROTATION_ON')
           Call readI(wlc_p%MOVEON(5)) ! is full chain rotate on 1/0
       CASE('FULL_CHAIN_SLIDE_ON')
           Call readI(wlc_p%MOVEON(6)) ! is full chain slide on 1/0
       CASE('MIN_CRANK_SHAFT_WIN')
           Call readF(wlc_p%MINWindoW(1)) ! min mean window size
       CASE('MIN_SLIDE_WIN')
           Call readF(wlc_p%MINWindoW(2))
       CASE('MIN_PIVOT_WIN')
           Call readF(wlc_p%MINWindoW(3))
       CASE('RESTART')
           call reado(wlc_p%restart) ! Restart from parallel tempering
       CASE('INTON')
           call readI(wlc_p%intON) ! self interactions on
       CASE('VHC')
           call readF(wlc_p%vhc) ! hard-core repulsion strength

       CASE DEFAULT
           print*, "Error in MCvar_setparams.  Unidentified keyword:", &
                   TRIM(WORD)
           stop 1
       endSELECT
    enddo
    close(pf)
end subroutine read_from_file

!Get Lks for parallel tempering from file
subroutine get_LKs_from_file(wlc_p)
  type(wlcsim_params), intent(inout) :: wlc_p
  integer nLKs !number of linking numbers
  integer IOstatus
  integer TempLk
  integer i
  nLKs = 0
  open (unit = 1, file = 'input/LKs')
  do
     read(unit = 1, fmt = *,iostat=IOstatus) TempLk
     if (IOstatus /= 0) exit
     nLKs = nLKs + 1
  end do
  close(unit = 1)

  wlc_p%nLKs = nLKs
  allocate(wlc_p%LKs(nLks))

  open(unit = 1, file = 'input/LKs')
  do i = 1, nLks
     read(unit = 1,fmt = *) wlc_p%Lks(i)
  enddo
  close(unit = 1)
end subroutine get_LKs_from_file
 
subroutine get_input_from_file(infile, wlc_p)
! Based on Elena's readkeys subroutine
    IMPLICIT NONE
    type(wlcsim_params), intent(inout) :: wlc_p
    character(1024), intent(in) :: infile  ! file with parameters
    integer :: PF   ! input file unit
    LOGICAL :: fileend = .FALSE. ! done reading file?
    CHARACTER*100 :: WORD ! parameter name currently being read in
    integer :: NITEMS ! number of items on the line in the parameter file

    call set_param_defaults(wlc_p)

    call read_from_file(infile, wlc_p)

    !If parallel tempering is on, read the Lks
    if (wlc_p%ptON) then
       call get_LKs_from_file(wlc_p)
    endif

    ! get derived parameters that aren't directly input from file
    call getpara(wlc_p)

end subroutine

!Perform checks to exit program if inputs not set up properly
subroutine idiot_checks(wlc_d,wlc_p)
  type(wlcsim_params), intent(in) :: wlc_p
  type(wlcsim_data), intent(in) :: wlc_d

  !Check that number of linking numbers is equal to number
  !of processes if parallel tempering is on
  if (wlc_p%ptON) then
     if (wlc_p%nLKs+1.ne.wlc_d%p) then
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print *, 'number of threads must be equal to number of replicas plus one'
        print *, 'exiting...'
        print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        stop
     endif
  endif
endsubroutine

!Initialize the random number generator
subroutine init_rng(wlc_d,wlc_p)
    implicit none

    type(wlcsim_params), intent(inout)  :: wlc_p
    type(wlcsim_data), intent(inout)   :: wlc_d

    !   variables for random number generator 
    integer Irand     ! Seed
    character*8 datedum  ! trash
    character*10 timedum ! trash
    character*5 zonedum  ! trash
    integer seedvalues(8) ! clock readings

    !Counting and MPI variables
    integer rep
    integer source,dest
    integer (kind = 4)  status(MPI_STATUS_SIZE)
    integer (kind = 4) error
    

    !Cases
    !parallel tempering on
    !parallel tempering off

    if (wlc_p%ptON) then
       !If parallel tempering is on, the head node generates a random seed and
       !sends it to the replicas for thread-safe random number generation

       !Cases
       !head node (id = 0)
       !worker node (id /= 0)

       if (wlc_d%id.eq.0) then
          
          !head node

          !Get the seed for the random number generator
          call date_and_time(datedum,timedum,zonedum,seedvalues)
          Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
               -seedvalues(7)*1E3-seedvalues(8))
          Irand=mod(Irand,10000)
          call random_setseed(Irand,wlc_d%rand_stat)

          !Loop over worker nodes and send them the seed
          do rep=1,wlc_p%nLKs
             dest=wlc_d%nodeNUMBER(rep)
             !Send out the random seed to the worker nodes
             call MPI_Send (Irand,1, MPI_INTEGER, dest,   0, &
                  MPI_COMM_WORLD,error )
          enddo
          
       else
          
          !worker nodes

          !receive the seed from the head node and modify it to start at 
          !a separate place in the RNG sequence
          
          source = 0

          call MPI_Recv(Irand,1, MPI_INTEGER, source,   0, &
               MPI_COMM_WORLD,status,error )

          call random_setseed(Irand*(wlc_d%id+1),wlc_d%rand_stat)

          

       endif

    else
       
       !Parallel tempering off

       !Get and set the seed for the random number generator
       
       call date_and_time(datedum,timedum,zonedum,seedvalues)
       Irand=int(-seedvalues(5)*1E7-seedvalues(6)*1E5 &
            -seedvalues(7)*1E3-seedvalues(8))
       Irand=mod(Irand,10000)
       call random_setseed(Irand,wlc_d%rand_stat)

    endif
end subroutine init_rng

subroutine initialize_wlcsim_data(wlc_d, wlc_p)
    implicit none
    type(wlcsim_params), intent(inout)  :: wlc_p
    type(wlcsim_data), intent(inout)   :: wlc_d
    integer i
    !MPI variables
    

    !Calculate total number of beads
    wlc_p%nT = wlc_p%nB*wlc_p%nP             

    ! initialize energies to zero
    wlc_d%EElas=0.0_dp
    wlc_d%Eint = 0.0_dp

    ! initialize windows to number of beads
    wlc_p%MAXWindoW = wlc_p%nB         ! Max Size of window for bead selection
    wlc_p% MINWindoW  = 1         ! Min Size of window for bead selection

    ! Window amplitudes
    wlc_p%MinAMP = 0.0_dp ! minium amplitude
    wlc_p%MinAMP(1) = 0.07_dp*pi
    wlc_p%MinAMP(2) = 0.01_dp*wlc_p%l/wlc_p%nB
    wlc_p%MaxAMP = 2.0_dp*pi
    wlc_p%MaxAMP(2) = wlc_p%lbox(1)
    wlc_p%MaxAMP(6) = wlc_p%lbox(1)

    ! Number of successes for MC
    wlc_d%SUCCESS = 0
    wlc_d%successTOTAL = 0
    wlc_d%pHIT = 0.0_dp

    !Initialize window size and move amplitudes
    wlc_d%Window = wlc_p%nB
    wlc_d%MCamp = 1.0_dp

    ! If ring is on, turn off the pivot move
    if (wlc_p%ring.eq.1) then
       wlc_p%moveON(3) = 0
    endif

    !If parallel tempering is on, initialize the nodeNumbers
    !and initialize MPI

    if (wlc_p%ptON) then
      
       !Allocate node numbers
       allocate(wlc_d%nodeNUMBER(wlc_p%nLKs))
       do i = 1,wlc_p%nLKs
          wlc_d%nodeNUMBER(i) = i
       enddo

       !Initially, replica start and replica end are the first and second to last replicas for even
       !nLKs and the first and second to last for odd nLKs
       if (mod(wlc_p%nLKs,2).eq.0) then
          wlc_d%replicaSTART = 1
          wlc_d%replicaEND = wlc_p%nLKs - 1
       else
          wlc_d%replicaSTART = 1
          wlc_d%replicaEND = wlc_p%nLKs - 2
       endif

       !Allocate the number of replica exchange trials and successes and initialize to zero
       allocate(wlc_d%nSWAPup(wlc_p%nLKs))
       allocate(wlc_d%nSWAPdown(wlc_p%nLKs))
       allocate(wlc_d%nTRIALup(wlc_p%nLKs))
       allocate(wlc_d%nTRIALdown(wlc_p%nLKs))

       wlc_d%nSWAPup = 0
       wlc_d%nSWAPdown = 0
       wlc_d%nTRIALup = 0
       wlc_d%nTRIALdown = 0

       !Initialize MPI
       call MPI_Init ( wlc_d%error )
       if (wlc_d%error.ne.0) then
          print*, "MPI_Init", wlc_d%error
       endif

      ! Get number of processes
       call MPI_Comm_size ( MPI_COMM_WORLD, wlc_d%p, wlc_d%error )
       if (wlc_d%error.ne.0) then
          print*, "MPI_Comm_size", wlc_d%error
       endif

      ! Get individual process id
       call MPI_Comm_rank ( MPI_COMM_WORLD, wlc_d%id, wlc_d%error )
       if (wlc_d%error.ne.0) then
          print*, "MPI_Comm_rank", wlc_d%error
       endif

       
    endif

    !Perform "idiot checks" for incompatible simulation parameters
    call idiot_checks(wlc_d,wlc_p)

    !Initialize random number generator (thread-safe)

    call init_rng(wlc_d,wlc_p)


  end subroutine initialize_wlcsim_data

subroutine MCvar_allocate(wlc_p,wlc_d)
  implicit none
  type(wlcsim_params), intent(inout) :: wlc_p
  type(wlcsim_data), intent(inout) :: wlc_d
  integer NT  ! total number of beads
  integer NP  ! total number of polymers

  NT = wlc_p%nT
  NP = wlc_p%nP

  !Allocate configuration vectors

  allocate(wlc_d%R(NT,3))
  allocate(wlc_d%U(NT,3))

  !Allocate vector of writhe and elastic energies for replicas
  if (wlc_p%ptON) then
     allocate(wlc_d%Wrs(wlc_p%nLKs))
     allocate(wlc_d%eelasREPLICAS(wlc_p%nLKs,4))
  endif

  !Allocate trajectory vectors
  allocate(wlc_d%wrTRAJ(wlc_p%indMAX))
  allocate(wlc_d%rgsqTRAJ(wlc_p%indMAX))

end subroutine

!subroutine for saving configurations
subroutine save_configuration(wlc_p,wlc_d,ind)
  implicit none
   type(wlcsim_params), intent(in) :: wlc_p
   type(wlcsim_data), intent(inout) :: wlc_d
   character(1024) savefile
   character(1024) lkSTRING
   character(1024) fileIND
   integer i,ind
   real(dp) rgsq
  !Cases
  !parallel tempering on
  !parallel tempering off

  !parallel tempering on
   if (wlc_p%ptON) then

      !!!!!!!!!!!!!!!!!!!!!!!!
      !Save r and u
      !!!!!!!!!!!!!!!!!!!!!!!!

      !Head node does not save configuration
      
      if (wlc_d%id.ne.0) then
         !Determine the folder to save in
         write(lkSTRING,*) wlc_p%lk
         write(savefile,*) 'data/LK_'//TRIM(ADJUSTL(lkSTRING))
         write(fileIND,*) ind
         !open file and write polymer configuration
         open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/r'//TRIM(ADJUSTL(fileIND)))
         open(unit = 2, file = TRIM(ADJUSTL(savefile))//'/u'//TRIM(ADJUSTL(fileIND)))

         do i = 1,wlc_p%nT
            write(1,*) wlc_d%r(i,:)
            write(2,*) wlc_d%u(i,:)
         enddo
         close(1)
         close(2)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Append writhe squared radius of gyration and energies to trajctory file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !calculate rgsq and polymer repulsive interaction energy (writhe and eelas already calculated
         !in wlcsim)

         call getRG(wlc_p,wlc_d,rgsq)

         if (wlc_p%inton.eq.1) then
            call ENERGY_SELF_CHAIN(wlc_d%Eint,wlc_d%R,wlc_p%nT,wlc_p%nB,wlc_p%LHC,wlc_p%VHC,wlc_p%RING)
         endif

         !save variables
         open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/wr', status = 'unknown', position = 'append')
         write(1,*) wlc_d%Wr
         close(1)

         open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/rgsq', status = 'unknown', position = 'append')
         write(1,*) rgsq
         close(1)

         open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/EELAS', status = 'unknown', position = 'append')
         write(1,*) wlc_d%EELAS
         close(1)


         if (wlc_p%inton.eq.1) then
            open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/Eint', status = 'unknown', position = 'append')
            write(1,*) wlc_d%Eint
            close(1)
         endif

      endif

  else
     write(savefile,*) 'data/'
     write(fileIND,*) ind
     !open file and write polymer configuration
     open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/r'//TRIM(ADJUSTL(fileIND)))
     open(unit = 2, file = TRIM(ADJUSTL(savefile))//'/u'//TRIM(ADJUSTL(fileIND)))

     do i = 1,wlc_p%nT
        write(1,*) wlc_d%r(i,:)
        write(2,*) wlc_d%u(i,:)
     enddo
     close(1)
     close(2)
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Append writhe squared radius of gyration and energies to trajctory file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !calculate rgsq and polymer repulsive interaction energy (writhe and eelas already calculated
     !in wlcsim)
     
     call getRG(wlc_p,wlc_d,rgsq)

     if (wlc_p%inton.eq.1) then
        call ENERGY_SELF_CHAIN(wlc_d%Eint,wlc_d%R,wlc_p%nT,wlc_p%nB,wlc_p%LHC,wlc_p%VHC,wlc_p%RING)
     endif

     !save variables
     open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/wr', status = 'unknown', position = 'append')
     write(1,*) wlc_d%Wr
     close(1)

     open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/rgsq', status = 'unknown', position = 'append')
     write(1,*) rgsq
     close(1)

     open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/EELAS', status = 'unknown', position = 'append')
     write(1,*) wlc_d%EELAS
     close(1)


     if (wlc_p%inton.eq.1) then
        open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/Eint', status = 'unknown', position = 'append')
        write(1,*) wlc_d%Eint
        close(1)
     endif


  endif

end subroutine save_configuration


!Save statistics for structural and energetic quantities
subroutine save_stats(wlc_p,wlc_d)
  implicit none
  type(wlcsim_params), intent(inout) :: wlc_p
  type(wlcsim_data), intent(inout) :: wlc_d
  !statistics
  real(dp) wrAVG
  real(dp) rgsqAVG
  real(dp) wrSTDEV
  real(dp) rgsqSTDEV
  real(dp) wrSTDER
  real(dp) rgsqSTDER
  !Variables for saving
  character(1024) savefile
  character(1024) lkSTRING
  !counting variables
  integer SaveInd

  !Load Wr and rgsq from their trajectory files

  !Cases
  !ptON = .TRUE. - head node does not load and replicas load from their directories
  !ptON = .FALSE. node loads directly from data directory

  if (wlc_p%ptON) then

     !only worker nodes get their trajectories
     if (wlc_d%id.ne.0 ) then
        !Determine the folder to load from
        
        write(lkSTRING,*) wlc_p%lk
        write(savefile,*) 'data/LK_'//TRIM(ADJUSTL(lkSTRING))

        !Read into writhe trajectory

        open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/wr', status = 'old')

        do SaveInd = 1,wlc_p%indMAX
           read(1,*) wlc_d%wrTRAJ(SaveInd)
        enddo

        close(1)

        !Read into rgsq trajectory

        open(unit = 1, file = TRIM(ADJUSTL(savefile))//'/rgsq', status = 'old')

        do SaveInd = 1,wlc_p%indMAX
           read(1,*) wlc_d%rgsqTRAJ(SaveInd)
        enddo

        close(1)
        
     endif

     
  else
     !Parallel tempering not on

        !Read into writhe trajectory

        open(unit = 1, file = 'data/wr', status = 'old')

        do SaveInd = 1,wlc_p%indMAX
           read(1,*) wlc_d%wrTRAJ(SaveInd)
        enddo

        close(1)

        !Read into rgsq trajectory

        open(unit = 1, file = 'data/rgsq', status = 'old')

        do SaveInd = 1,wlc_p%indMAX
           read(1,*) wlc_d%rgsqTRAJ(SaveInd)
        enddo

        close(1)
     

  endif

  !calculate average quantities
  wrAVG = sum(wlc_d%wrTRAJ)/dble(wlc_p%indMAX)
  rgsqAVG = sum(wlc_d%rgsqTRAJ)/dble(wlc_p%indMAX)

  !calculate standard deviations of quantities
  wrSTDEV = sqrt(sum((wlc_d%wrTRAJ - wrAVG)**2.0_dp)/dble(wlc_p%indMAX))
  rgsqSTDEV = sqrt(sum((wlc_d%rgsqTRAJ - rgsqAVG)**2.0_dp)/dble(wlc_p%indMAX))

  !calculate standard errors
  wrSTDER = wrSTDEV/sqrt(dble(wlc_p%indMAX))
  rgsqSTDER =rgsqSTDEV/sqrt(dble(wlc_p%indMAX))

  !Save statistics

  !Cases
  !parallel tempering on
  !parallel tempering off

  if (wlc_p%ptON) then

     !Head node does not write data
     if (wlc_d%id.ne.0) then
        !Get correct save folder for this lk
        write(lkSTRING,*) wlc_p%lk
        write(savefile,*) 'data/LK_'//trim(adjustL(lkSTRING))

     endif
     
  else
     write(savefile,*) 'data'

  endif
     
  !open files and save


  if ((wlc_p%ptON.and.wlc_d%id.ne.0).or.(.not.wlc_p%ptON)) then
     open(unit =1 , file = trim(adjustL(savefile))//'/wrAVG', status = 'replace')
     write(1,*) wrAVG
     close(1)

     open(unit =1 , file = trim(adjustL(savefile))//'/rgsqAVG', status = 'replace')
     write(1,*) rgsqAVG
     close(1)

     open(unit =1 , file = trim(adjustL(savefile))//'/wrSTDEV', status = 'replace')
     write(1,*) wrSTDEV
     close(1)

     open(unit =1 , file = trim(adjustL(savefile))//'/rgsqSTDEV', status = 'replace')
     write(1,*) rgsqSTDEV
     close(1)

     open(unit =1 , file = trim(adjustL(savefile))//'/wrSTDER', status = 'replace')
     write(1,*) wrSTDER
     close(1)

     open(unit =1 , file = trim(adjustL(savefile))//'/rgsqSTDER', status = 'replace')
     write(1,*) rgsqSTDER
     close(1)

  endif

endsubroutine

end module params
