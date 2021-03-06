!---------------------------------------------------------------*

! Find change in bead position for a crank-shaft type move
     
SUBROUTINE MC_move(R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,MCTYPE,MCAMP,WINDOW,RING,DIB,stat)
  use mersenne_twister
  IMPLICIT NONE
  
  
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0 ! Value of pi

  DOUBLE PRECISION R(NT,3)  ! Bead positions
  DOUBLE PRECISION U(NT,3)  ! Tangent vectors
  DOUBLE PRECISION RP(NT,3)  ! Bead positions
  DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
  INTEGER N,NP,NT           ! Number of beads

  INTEGER IP                ! Test polymer 
  INTEGER IB1               ! Test bead position 1
  INTEGER IT1               ! Index of test bead 1
  INTEGER IB2               ! Test bead position 2
  INTEGER IT2               ! Index of test bead 2
  INTEGER DIB               ! Number of beads between test bead 1 and test bead 2
  INTEGER I,J			! Test indices
  INTEGER RING              ! Is polymer a ring?

  ! Variables for the crank-shaft move

  DOUBLE PRECISION TA(3)    ! Axis of rotation
  DOUBLE PRECISION P1(3)    ! Point on rotation line
  DOUBLE PRECISION MAG      ! Magnitude of vector
  DOUBLE PRECISION ROT(4,4) ! Rotation matrix
  DOUBLE PRECISION ALPHA    ! Angle of move
  DOUBLE PRECISION BETA     ! Angle of move

  !     MC adaptation variables

  DOUBLE PRECISION MCAMP(6) ! Amplitude of random change 
  INTEGER WINDOW(6)         ! Window of MC move (max # of bead segments)
  INTEGER MCTYPE            ! Type of MC move
  DOUBLE PRECISION DR(3)    ! Displacement for slide move
  INTEGER TEMP

  !     Random number generator
  type(random_stat) stat  ! state of random number chain
  real urand(1)       ! random number
  !     Initialize new positions and orientations

  RP=R
  UP=U
  P1=0.

  !     Perform crank-shaft move (MCTYPE 1)
  
  if (MCTYPE.EQ.1) then
     call random_number(urand,stat)
     IP=nint(0.5+urand(1)*NP)
     call random_number(urand,stat)
     IB1=nint(0.5+urand(1)*N)
     call random_number(urand,stat)
     IB2=nint(0.5+urand(1)*N)

     IT1=N*(IP-1)+IB1
     IT2=N*(IP-1)+IB2

     if (RING.EQ.1) then                    !Polymer is a ring
        call random_number(urand,stat)
        DIB=nint(urand(1)*WINDOW(1))
        IF (IB1+DIB.LE.N) THEN
           IB2=IB1+DIB
        ELSE 
           IB2=DIB-(N-IB1)
        ENDIF
        IT2=N*(IP-1)+IB2

        if (IB1.EQ.IB2.AND.IB1.EQ.1) then
           TA(1)=R(IT1+1,1)-R(N*IP,1)
           TA(2)=R(IT1+1,2)-R(N*IP,2)
           TA(3)=R(IT1+1,3)-R(N*IP,3)
        elseif (IB1.EQ.IB2.AND.IB1.EQ.N) then
           TA(1)=R(N*(IP-1)+1,1)-R(IT1-1,1)
           TA(2)=R(N*(IP-1)+1,2)-R(IT1-1,2)
           TA(3)=R(N*(IP-1)+1,3)-R(IT1-1,3)
        elseif (IB1.EQ.IB2.AND.IB1.NE.1.AND.IB2.NE.N) then
           TA(1)=R(IT1+1,1)-R(IT1-1,1)
           TA(2)=R(IT1+1,2)-R(IT1-1,2)
           TA(3)=R(IT1+1,3)-R(IT1-1,3)
        else
           TA(1)=R(IT2,1)-R(IT1,1)
           TA(2)=R(IT2,2)-R(IT1,2)
           TA(3)=R(IT2,3)-R(IT1,3)
        endif
     else                                 !Polymer is not a ring
        if (IT1.GT.IT2) then
           TEMP=IT1
           IT1=IT2
           IT2=TEMP
           TEMP=IB1
           IB1=IB2
           IB2=TEMP
        endif
        DIB=IT2-IT1
        if (IB1.EQ.IB2.AND.IB1.EQ.1) then
           TA(1)=R(IT1+1,1)-R(IT1,1)
           TA(2)=R(IT1+1,2)-R(IT1,2)
           TA(3)=R(IT1+1,3)-R(IT1,3)
        elseif (IB1.EQ.IB2.AND.IB1.EQ.N) then
           TA(1)=R(N*IP,1)-R(N*IP-1,1)
           TA(2)=R(N*IP,2)-R(N*IP-1,2)
           TA(3)=R(N*IP,3)-R(N*IP-1,3)
        elseif (IB1.EQ.IB2.AND.IB1.NE.1.AND.IB2.NE.N) then
           TA(1)=R(IT1+1,1)-R(IT1-1,1)
           TA(2)=R(IT1+1,2)-R(IT1-1,2)
           TA(3)=R(IT1+1,3)-R(IT1-1,3)
        else
           TA(1)=R(IT2,1)-R(IT1,1)
           TA(2)=R(IT2,2)-R(IT1,2)
           TA(3)=R(IT2,3)-R(IT1,3)
        endif
     endif


     MAG=sqrt(TA(1)**2.+TA(2)**2.+TA(3)**2.)
     TA(1)=TA(1)/MAG
     TA(2)=TA(2)/MAG
     TA(3)=TA(3)/MAG
     P1(1)=R(IT1,1)
     P1(2)=R(IT1,2)
     P1(3)=R(IT1,3)	  
     call random_number(urand,stat)
     ALPHA=MCAMP(1)*(urand(1)-0.5)

     ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
     ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
     ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
          -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

     ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
     ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
     ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
          -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

     ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
     ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
     ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
     ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
          -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

     I=IT1

     DO  J=0,DIB
        if (I.EQ.(N*IP+1).AND.RING.EQ.1) then
           I=N*(IP-1)+1
        endif
        RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
        RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
        RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
        UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
        UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
        UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)

        I=I+1

     ENDDO

     !     Perform slide move (MCTYPE 2)

  elseif (MCTYPE.EQ.2) then
     call random_number(urand,stat)
     IP=nint(0.5+urand(1)*NP)
     call random_number(urand,stat)
     IB1=nint(0.5+urand(1)*N)
     call random_number(urand,stat)
     IB2=nint(0.5+urand(1)*N)

     IT1=N*(IP-1)+IB1
     IT2=N*(IP-1)+IB2

     if (RING.EQ.1) then
        call random_number(urand,stat)         
        DIB=nint(urand(1)*WINDOW(2))
        IF (IB1+DIB.LE.N) THEN
           IB2=IB1+DIB
        ELSE 
           IB2=DIB-(N-IB1)
        ENDIF
        IT2=N*(IP-1)+IB2

     else 
        if (IT1.GT.IT2) then
           TEMP=IT1
           IT1=IT2
           IT2=TEMP
           TEMP=IB1
           IB1=IB2
           IB2=TEMP
        endif
        DIB=IT2-IT1
     endif
     call random_number(urand,stat)
     DR(1)=MCAMP(2)*(urand(1)-0.5)
     call random_number(urand,stat)
     DR(2)=MCAMP(2)*(urand(1)-0.5)
     call random_number(urand,stat)
     DR(3)=MCAMP(2)*(urand(1)-0.5)


     I=IT1
     DO  J=0,DIB

        if (I.EQ.(N*IP+1).AND.RING.EQ.1) then
           I=N*(IP-1)+1
        endif

        RP(I,1)=R(I,1)+DR(1)
        RP(I,2)=R(I,2)+DR(2)
        RP(I,3)=R(I,3)+DR(3)
        UP(I,1)=U(I,1)
        UP(I,2)=U(I,2)
        UP(I,3)=U(I,3)
        I=I+1

     ENDDO


     !     Perform pivot move (MCTYPE 3)

  elseif (MCTYPE.EQ.3) then
     call random_number(urand,stat)
     IP=nint(0.5+urand(1)*NP)
     call random_number(urand,stat)
     IB1=nint(0.5+urand(1)*N)
     if (IB1.LT.(N/2.)) then
        IB2=IB1			
        IB1=1
        IT1=N*(IP-1)+IB1
        IT2=N*(IP-1)+IB2
        P1(1)=R(IT2,1)
        P1(2)=R(IT2,2)
        P1(3)=R(IT2,3)	  
     else
        IB2=N
        IT1=N*(IP-1)+IB1
        IT2=N*(IP-1)+IB2
        P1(1)=R(IT1,1)
        P1(2)=R(IT1,2)
        P1(3)=R(IT1,3)	  
     endif
     call random_number(urand,stat)             
     ALPHA=2.*PI*urand(1)
     call random_number(urand,stat)
     BETA=acos(2.*urand(1)-1.)
     TA(1)=sin(BETA)*cos(ALPHA)
     TA(2)=sin(BETA)*sin(ALPHA)
     TA(3)=cos(BETA)
     call random_number(urand,stat)
     ALPHA=MCAMP(3)*(urand(1)-0.5)

     ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
     ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
     ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
          -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

     ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
     ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
     ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
          -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

     ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
     ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
     ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
     ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
          -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

     DO  I=IT1,IT2
        RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
        RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
        RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
        UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
        UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
        UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)		 
     ENDDO

     !     Perform rotate move (MCTYPE 4)

  elseif (MCTYPE.EQ.4) then
     call random_number(urand,stat)
     IP=nint(0.5+urand(1)*NP)
     call random_number(urand,stat)
     IB1=nint(0.5+urand(1)*N)
     call random_number(urand,stat)
     IB2=IB1
     IT1=N*(IP-1)+IB1
     IT2=N*(IP-1)+IB2
     call random_number(urand,stat)
     ALPHA=2.*PI*urand(1)
     call random_number(urand,stat)
     BETA=acos(2.*urand(1)-1.)
     TA(1)=sin(BETA)*cos(ALPHA)
     TA(2)=sin(BETA)*sin(ALPHA)
     TA(3)=cos(BETA)
     call random_number(urand,stat)
     ALPHA=MCAMP(4)*(urand(1)-0.5)

     ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
     ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
     ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
          -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

     ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
     ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
     ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
          -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

     ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
     ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
     ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
     ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
          -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

     I=IT1
     UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
     UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
     UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)

     !     Perform a full chain rotation

  elseif (MCTYPE.EQ.5) then
     call random_number(urand,stat)
     IP=nint(0.5+urand(1)*NP)
     IB1=1
     IB2=N
     IT1=N*(IP-1)+IB1
     IT2=N*(IP-1)+IB2
     call random_number(urand,stat)
     ALPHA=2.*PI*urand(1)
     call random_number(urand,stat)
     BETA=acos(2.*urand(1)-1.)
     TA(1)=sin(BETA)*cos(ALPHA)
     TA(2)=sin(BETA)*sin(ALPHA)
     TA(3)=cos(BETA)
     call random_number(urand,stat)
     ALPHA=MCAMP(5)*(urand(1)-0.5)

     ROT(1,1)=TA(1)**2.+(TA(2)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(1,2)=TA(1)*TA(2)*(1.-cos(ALPHA))-TA(3)*sin(ALPHA)
     ROT(1,3)=TA(1)*TA(3)*(1.-cos(ALPHA))+TA(2)*sin(ALPHA)
     ROT(1,4)=(P1(1)*(1.-TA(1)**2.) &
          -TA(1)*(P1(2)*TA(2)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(2)*TA(3)-P1(3)*TA(2))*sin(ALPHA)

     ROT(2,1)=TA(1)*TA(2)*(1.-cos(ALPHA))+TA(3)*sin(ALPHA)
     ROT(2,2)=TA(2)**2.+(TA(1)**2.+TA(3)**2.)*cos(ALPHA)
     ROT(2,3)=TA(2)*TA(3)*(1.-cos(ALPHA))-TA(1)*sin(ALPHA)
     ROT(2,4)=(P1(2)*(1.-TA(2)**2.) &
          -TA(2)*(P1(1)*TA(1)+P1(3)*TA(3)))*(1.-cos(ALPHA))+(P1(3)*TA(1)-P1(1)*TA(3))*sin(ALPHA)

     ROT(3,1)=TA(1)*TA(3)*(1.-cos(ALPHA))-TA(2)*sin(ALPHA)
     ROT(3,2)=TA(2)*TA(3)*(1.-cos(ALPHA))+TA(1)*sin(ALPHA)
     ROT(3,3)=TA(3)**2.+(TA(1)**2.+TA(2)**2.)*cos(ALPHA)
     ROT(3,4)=(P1(3)*(1.-TA(3)**2.) &
          -TA(3)*(P1(1)*TA(1)+P1(2)*TA(2)))*(1.-cos(ALPHA))+(P1(1)*TA(2)-P1(2)*TA(1))*sin(ALPHA)

     DO  I=IT1,IT2
        RP(I,1)=ROT(1,4)+ROT(1,1)*R(I,1)+ROT(1,2)*R(I,2)+ROT(1,3)*R(I,3)
        RP(I,2)=ROT(2,4)+ROT(2,1)*R(I,1)+ROT(2,2)*R(I,2)+ROT(2,3)*R(I,3)
        RP(I,3)=ROT(3,4)+ROT(3,1)*R(I,1)+ROT(3,2)*R(I,2)+ROT(3,3)*R(I,3)
        UP(I,1)=ROT(1,1)*U(I,1)+ROT(1,2)*U(I,2)+ROT(1,3)*U(I,3)
        UP(I,2)=ROT(2,1)*U(I,1)+ROT(2,2)*U(I,2)+ROT(2,3)*U(I,3)
        UP(I,3)=ROT(3,1)*U(I,1)+ROT(3,2)*U(I,2)+ROT(3,3)*U(I,3)		 
     ENDDO

     !     Perform full chain slide move (MCTYPE 2)

  elseif (MCTYPE.EQ.6) then
     call random_number(urand,stat)
     IP=nint(0.5+urand(1)*NP)
     IB1=1
     IB2=N
     IT1=N*(IP-1)+IB1
     IT2=N*(IP-1)+IB2
     call random_number(urand,stat)
     DR(1)=MCAMP(6)*(urand(1)-0.5)
     call random_number(urand,stat)
     DR(2)=MCAMP(6)*(urand(1)-0.5)
     call random_number(urand,stat)
     DR(3)=MCAMP(6)*(urand(1)-0.5)

     DO  I=IT1,IT2
        RP(I,1)=R(I,1)+DR(1)
        RP(I,2)=R(I,2)+DR(2)
        RP(I,3)=R(I,3)+DR(3)
        UP(I,1)=U(I,1)
        UP(I,2)=U(I,2)
        UP(I,3)=U(I,3)
     ENDDO

  endif

  RETURN      
END SUBROUTINE MC_MOVE
      
!---------------------------------------------------------------!      
