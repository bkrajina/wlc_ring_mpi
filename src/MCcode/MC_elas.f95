!---------------------------------------------------------------*

! Calculate the change in the elastic energy after a monte carlo move.
! Elastic energy is calculated as that of a discrete stretchable-shearable worm-like chain
     
SUBROUTINE MC_eelas(DEELAS,R,U,RP,UP,NT,N,NP,IP,IB1,IB2,IT1,IT2,EB,EPAR,EPERP,GAM,ETA,RING,TWIST,Lk,lt,LP,L,MCTYPE,WR,WRP)
  IMPLICIT NONE

  REAL, PARAMETER :: PI = 3.141592654d0

  DOUBLE PRECISION R(NT,3)  ! Bead positions
  DOUBLE PRECISION U(NT,3)  ! Tangent vectors
  DOUBLE PRECISION RP(NT,3)  ! Bead positions (test polymer)
  DOUBLE PRECISION UP(NT,3)  ! Tangent vectors
  INTEGER  Lk               ! Linking number
  DOUBLE PRECISION L        ! Contour length 
  DOUBLE PRECISION lt       ! Twist persistence length
  DOUBLE PRECISION LP       ! Bend persistence length
  DOUBLE PRECISION Tw       ! Twist 
  DOUBLE PRECISION TWP      ! Twist of test structure
  DOUBLE PRECISION WR       ! Writhe
  DOUBLE PRECISION WRP      ! Writhe of test structure
  DOUBLE PRECISION DWR      ! Change in Writhe
  DOUBLE PRECISION WRM,WRMP ! Component of writhe affected by move
  INTEGER N,NP,NT           ! Number of beads
  DOUBLE PRECISION ETWIST

  INTEGER IP                ! Test polymer 
  INTEGER IB1               ! Test bead position 1
  INTEGER IT1               ! Index of test bead 1
  INTEGER IB2               ! Test bead position 2
  INTEGER IT2               ! Index of test bead 2
  INTEGER RING              ! Is polymer a ring?
  INTEGER TWIST             ! Include twist?
  INTEGER IT1M1
  INTEGER IT2P1
  INTEGER MCTYPE            ! MC move type

  DOUBLE PRECISION DEELAS   ! Change in ECOM      

  !     Polymer properties

  DOUBLE PRECISION EB,EPAR,EPERP
  DOUBLE PRECISION GAM,ETA

  !     Variables for force and torque calculations

  DOUBLE PRECISION DR(3),DRPAR,DRPERP(3),DRPERPM(3)
  DOUBLE PRECISION FI(3),TI(3)
  DOUBLE PRECISION U1U2,GI(3),DOTGU,HI(3)

! Setup parameters
      
  DEELAS=0.d0
      
!     Calculate the change in the energy
  
  
  if (RING.EQ.1) then

     if (IB1.EQ.1) then
        IT1M1=N*IP
     else
        IT1M1=IT1-1
     endif
  else 
     if (IB1.EQ.1) then
        GO TO 10 ! skip calculation for change in first bead energy
     else 
        IT1M1=IT1-1
     endif
  endif

  DR(1)=R(IT1,1)-R(IT1M1,1)
  DR(2)=R(IT1,2)-R(IT1M1,2)
  DR(3)=R(IT1,3)-R(IT1M1,3)
  DRPAR=DR(1)*U(IT1M1,1)+DR(2)*U(IT1M1,2)+DR(3)*U(IT1M1,3)

  DRPERP(1)=DR(1)-DRPAR*U(IT1M1,1)
  DRPERP(2)=DR(2)-DRPAR*U(IT1M1,2)
  DRPERP(3)=DR(3)-DRPAR*U(IT1M1,3)
  U1U2=U(IT1M1,1)*U(IT1,1)+U(IT1M1,2)*U(IT1,2)+U(IT1M1,3)*U(IT1,3)

  GI(1)=(U(IT1,1)-U(IT1M1,1)-ETA*DRPERP(1))
  GI(2)=(U(IT1,2)-U(IT1M1,2)-ETA*DRPERP(2))
  GI(3)=(U(IT1,3)-U(IT1M1,3)-ETA*DRPERP(3))

  DEELAS=DEELAS-0.5*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.) &
       -0.5*EPAR*(DRPAR-GAM)**2.-0.5*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

  DR(1)=RP(IT1,1)-R(IT1M1,1)
  DR(2)=RP(IT1,2)-R(IT1M1,2)
  DR(3)=RP(IT1,3)-R(IT1M1,3)
  DRPAR=DR(1)*U(IT1M1,1)+DR(2)*U(IT1M1,2)+DR(3)*U(IT1M1,3)

  DRPERP(1)=DR(1)-DRPAR*U(IT1M1,1)
  DRPERP(2)=DR(2)-DRPAR*U(IT1M1,2)
  DRPERP(3)=DR(3)-DRPAR*U(IT1M1,3)
  U1U2=U(IT1M1,1)*UP(IT1,1)+U(IT1M1,2)*UP(IT1,2)+U(IT1M1,3)*UP(IT1,3)

  GI(1)=(UP(IT1,1)-U(IT1M1,1)-ETA*DRPERP(1))
  GI(2)=(UP(IT1,2)-U(IT1M1,2)-ETA*DRPERP(2))
  GI(3)=(UP(IT1,3)-U(IT1M1,3)-ETA*DRPERP(3))

  DEELAS=DEELAS+0.5*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.) &
       +0.5*EPAR*(DRPAR-GAM)**2.+0.5*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)	 
10 if (RING.EQ.1) then
     if (IB2.EQ.N) then
        IT2P1=N*(IP-1)+1
     else
        IT2P1=IT2+1
     endif
  else 
     if (IB2.EQ.N) then
        GO TO 20 
     else
        IT2P1=IT2+1
     endif
  endif





  DR(1)=R(IT2P1,1)-R(IT2,1)
  DR(2)=R(IT2P1,2)-R(IT2,2)
  DR(3)=R(IT2P1,3)-R(IT2,3)
  DRPAR=DR(1)*U(IT2,1)+DR(2)*U(IT2,2)+DR(3)*U(IT2,3)

  DRPERP(1)=DR(1)-DRPAR*U(IT2,1)
  DRPERP(2)=DR(2)-DRPAR*U(IT2,2)
  DRPERP(3)=DR(3)-DRPAR*U(IT2,3)
  U1U2=U(IT2,1)*U(IT2P1,1)+U(IT2,2)*U(IT2P1,2)+U(IT2,3)*U(IT2P1,3)

  GI(1)=(U(IT2P1,1)-U(IT2,1)-ETA*DRPERP(1))
  GI(2)=(U(IT2P1,2)-U(IT2,2)-ETA*DRPERP(2))
  GI(3)=(U(IT2P1,3)-U(IT2,3)-ETA*DRPERP(3))

  DEELAS=DEELAS-0.5*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.) &
       -0.5*EPAR*(DRPAR-GAM)**2.-0.5*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

  DR(1)=R(IT2P1,1)-RP(IT2,1)
  DR(2)=R(IT2P1,2)-RP(IT2,2)
  DR(3)=R(IT2P1,3)-RP(IT2,3)
  DRPAR=DR(1)*UP(IT2,1)+DR(2)*UP(IT2,2)+DR(3)*UP(IT2,3)

  DRPERP(1)=DR(1)-DRPAR*UP(IT2,1)
  DRPERP(2)=DR(2)-DRPAR*UP(IT2,2)
  DRPERP(3)=DR(3)-DRPAR*UP(IT2,3)
  U1U2=UP(IT2,1)*U(IT2P1,1)+UP(IT2,2)*U(IT2P1,2)+UP(IT2,3)*U(IT2P1,3)

  GI(1)=(U(IT2P1,1)-UP(IT2,1)-ETA*DRPERP(1))
  GI(2)=(U(IT2P1,2)-UP(IT2,2)-ETA*DRPERP(2))
  GI(3)=(U(IT2P1,3)-UP(IT2,3)-ETA*DRPERP(3))

  DEELAS=DEELAS+0.5*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.) &
       +0.5*EPAR*(DRPAR-GAM)**2.+0.5*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)


  IF (RING.EQ.1.AND.TWIST.EQ.1) THEN
     IF (MCTYPE.EQ.1) THEN
        CALL WRITHECRANK(R,IT1,IT2,N,WRM)
        CALL WRITHECRANK(RP,IT1,IT2,N,WRMP)
        DWR=WRMP-WRM
     ELSEIF (MCTYPE.EQ.2) THEN
        CALL WRITHESLIDE(R,IT1,IT2,N,WRM)
        CALL WRITHESLIDE(RP,IT1,IT2,N,WRMP)
        DWR=WRMP-WRM
     ELSE
        DWR=0.
     ENDIF
     WRP=WR+DWR
     TW=REAL(LK)-WR
     TWP=REAL(LK)-WRP
     DEELAS=DEELAS+(((2.*pi*TWP)**2.)*LT/(2.*L))-(((2.*pi*TW)**2.)*LT/(2.*L))


  ENDIF



20 RETURN      
ENDSUBROUTINE MC_eelas
      
!---------------------------------------------------------------*
