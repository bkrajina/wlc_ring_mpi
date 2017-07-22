!! ---------------------------------------------------------------*
      
!     
!     This subroutine sets the initial condition for a chain
!     
!     Andrew Spakowitz and Brad Krajina
!     Updated 11/3/2016

      
SUBROUTINE initcond(R,U,NT,N,NP,FRMFILE,GAM,LBOX,RING,rand_stat,L,RESTRICTEDR,MAXEND2END,MINEND2END)
      
  use mersenne_twister

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: PI=3.141592654d0 ! Value of pi

  DOUBLE PRECISION R(NT,3)  ! Bead positions
  DOUBLE PRECISION U(NT,3)  ! Tangent vectors
  DOUBLE PRECISION L        ! Length of the chain
  INTEGER N,NP,NT           ! Number of beads
  DOUBLE PRECISION GAM      ! Equil bead separation
  DOUBLE PRECISION LBOX(3)  ! Box edge length
  INTEGER I,J,IB            ! Index Holders
  DOUBLE PRECISION S        ! Initial end to end distance
  LOGICAL FRMfile           ! Is conformation in file?
  INTEGER INPUT             ! Is input file set?
  INTEGER RING              ! Is polymer a ring?
  DOUBLE PRECISION RMIN
  DOUBLE PRECISION R0(3)
  type(random_stat) rand_stat  ! state of random number chain
  real urand(1)        !random number
  INTEGER RESTRICTEDR         ! Is rectriction between site distance on?
  DOUBLE PRECISION MAXEND2END ! Maximum allowed distance between sites   
  DOUBLE PRECISION MINEND2END ! Minimum allowed distance between sites    


  !     Input the conformation if FRMFILE=1

  if(FRMFILE)then
     OPEN (UNIT = 5, FILE = 'initial/rinit', STATUS = 'OLD')
     DO  I=1,N
        READ(5,*) R(I,:)
     ENDDO
     CLOSE(5)
     
     OPEN (UNIT = 5, FILE = 'initial/uinit', STATUS = 'OLD')
     DO  I=1,N
        READ(5,*) U(I,:)
     ENDDO
     CLOSE(5)
     

  !     Set the initial conformation to a straight chain if FRMFILE=0 and RING=0
  !     Set initial conformation to a circle if RING=1

  else
     IB=1
     DO  I=1,NP

        !Cases
        !NP > 1
        !NP = 1
        !If there are multiple polymers, initialize first bead randomly in the box
        !If there is one, put it in the middle of the box

        if (NP.GT.1) then
           call random_number(urand,rand_stat)
           R0(1)=urand(1)*lbox(1)
           call random_number(urand,rand_stat)
           R0(2)=urand(1)*lbox(2)
           call random_number(urand,rand_stat)
           R0(3)=urand(1)*lbox(3)
        else

           R0(1) = lbox(1)/(0.5d0)
           R0(2) = lbox(2)/(0.5d0)
           R0(3) = lbox(3)/(0.5d0)
        endif

        !Generate the initial chain configuration
        !Cases
        !RING = 0
        !RING = 1
        !RESTRICTEDR = 0 (Umbrella sampling off)
        !RESTRICTEDR = 1 (Umbrella sampling on)
        !If not a ring, initialize as a straight line
        !If a ring, initialize as a circle
        !If u.s. is on, initialize as a isosceles triangle
        
        S=(MAXEND2END+MINEND2END)/2  ! Initial distance between sites

        DO  J=1,N
           IF (RING.EQ.0) THEN
              IF (RESTRICTEDR.EQ.1) THEN
!               IF (GAM*(N-1).LT.L) THEN
!                R(IB,1)=R0(1)
!                R(IB,2)=R0(2)+GAM*(J-N/2.-0.5)
!                R(IB,3)=R0(3)
!                U(IB,1)=0.
!                U(IB,2)=1.
!                U(IB,3)=0.
!               ELSE
                IF (J.LT.N/2) THEN                          ! First half of triangle positive slope
                 R(IB,1)=R0(1)+J*S/(N-1)                    ! From similarity of triangle
                 R(IB,2)=R0(2)+J*sqrt((GAM*(N-1))**2-S**2)/(N-1) ! From similarity of triangle and Pythagorean theorem
                 R(IB,3)=R0(3)
                 U(IB,1)=S/(GAM*(N-1))
                 U(IB,2)=1-S**2/(GAM*(N-1))**2
                 U(IB,3)=0.
                ELSE                                                ! Second half of the triangle negative slope
                 R(IB,1)=R(IB-1,1)+S/(N-1)                        ! From similarity of triangle
                 R(IB,2)=R0(2)+(N-J)*sqrt((GAM*(N-1))**2-S**2)/(N-1)  ! From similarity of triangle and Pythagorean theorem 
                 R(IB,3)=R0(3)
                 U(IB,1)=S/(GAM*(N-1))
                 U(IB,2)=S**2/(GAM*(N-1))**2-1
                 U(IB,3)=0.
                ENDIF
!               ENDIF
              ELSE
              R(IB,1)=R0(1)
              R(IB,2)=R0(2)+GAM*(J-N/2.-0.5)
              R(IB,3)=R0(3)
              U(IB,1)=0.
              U(IB,2)=1.
              U(IB,3)=0.
             ENDIF
           ELSE
              R(IB,1)=R0(1)+((GAM*N)/(2*PI))*Cos(J*2*PI/N)
              R(IB,2)=R0(2)+((GAM*N)/(2*PI))*Sin(J*2*PI/N)
              R(IB,3)=0
              U(IB,1)=-Sin(J*2*PI/N)
              U(IB,2)=Cos(J*2*PI/N)
              U(IB,3)=0;
           ENDIF
           IB=IB+1
        ENDDO
     ENDDO

  endif

  RETURN     
END SUBROUTINE initcond

!---------------------------------------------------------------*
