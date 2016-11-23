!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This subroutine calculates the elastic energy
!for a discrete-stretchable shearable wormlike
!chain with twist. 
!Andrew J. Spakowitz and Brad A. Krajina
!Last edited: 2016/11/3
      
SUBROUTINE energy_elas(EELAS,R,U,NT,N,NP,EB,EPAR,EPERP,GAM,ETA,RING,TWIST,Lk,lt,LP,L)
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0

  DOUBLE PRECISION EELAS(4) ! Elastic force
  DOUBLE PRECISION R(NT,3)  ! Bead positions
  DOUBLE PRECISION U(NT,3)  ! Tangent vectors
  DOUBLE PRECISION B(NT)  ! Tangent vectors
  DOUBLE PRECISION UR(NT,3)  ! Tangent vectors
  DOUBLE PRECISION KAP,EPS  ! Elastic props
  DOUBLE PRECISION L0       ! Bead separation
  DOUBLE PRECISION FCOM(3)  ! Compress force
  DOUBLE PRECISION FBEND(3) ! Bend force
  INTEGER I,J,IB,IBP1       ! Index holders
  INTEGER N,NT,NP           ! Number of bead

  !     Polymer properties
  DOUBLE PRECISION EB,EPAR,EPERP
  DOUBLE PRECISION GAM,ETA
  DOUBLE PRECISION L       ! Contour length
  DOUBLE PRECISION LT      ! Twist persistence length
  DOUBLE PRECISION LP      ! Bend persistence length
  INTEGER LK               ! Linking number
  INTEGER RING             ! Is polymer a ring?
  INTEGER TWIST            ! Include twist?
  DOUBLE PRECISION Tw      ! Twist
  DOUBLE PRECISION Wr      ! Writhe

  !     Variables for force and torque calculations

  DOUBLE PRECISION DR(3),DRPAR,DRPERP(3)
  DOUBLE PRECISION FI(3),TI(3)
  DOUBLE PRECISION U1U2,GI(3),DOTGU,HI(3)

  !     Calculate the forces and torques


  EELAS(1)=0.
  EELAS(2)=0.
  EELAS(3)=0.
  EELAS(4)=0.
  IB=1
  DO  I=1,NP
     DO  J=1,N
        IF (RING.EQ.1) THEN
           IF (J.EQ.N) THEN
              IBP1=1+(I-1)*N
           ELSE
              IBP1=IB+1
           ENDIF
        ELSEIF (RING.EQ.0.AND.J.EQ.N) THEN
           GOTO 40
        ELSE
           IBP1=IB+1
        ENDIF



        DR(1)=R(IBP1,1)-R(IB,1)
        DR(2)=R(IBP1,2)-R(IB,2)
        DR(3)=R(IBP1,3)-R(IB,3)
        DRPAR=DR(1)*U(IB,1)+DR(2)*U(IB,2)+DR(3)*U(IB,3)

        DRPERP(1)=DR(1)-DRPAR*U(IB,1)
        DRPERP(2)=DR(2)-DRPAR*U(IB,2)
        DRPERP(3)=DR(3)-DRPAR*U(IB,3)
        U1U2=U(IB,1)*U(IBP1,1)+U(IB,2)*U(IBP1,2)+U(IB,3)*U(IBP1,3)

        GI(1)=(U(IBP1,1)-U(IB,1)-ETA*DRPERP(1))
        GI(2)=(U(IBP1,2)-U(IB,2)-ETA*DRPERP(2))
        GI(3)=(U(IBP1,3)-U(IB,3)-ETA*DRPERP(3))

        EELAS(1)=EELAS(1)+0.5*EB*(GI(1)**2.+GI(2)**2.+GI(3)**2.)
        EELAS(2)=EELAS(2)+0.5*EPAR*(DRPAR-GAM)**2.
        EELAS(3)=EELAS(3)+0.5*EPERP*(DRPERP(1)**2.+DRPERP(2)**2.+DRPERP(3)**2.)

        IB=IB+1		 
40      CONTINUE
     ENDDO
     IB=IB+1
30   CONTINUE
  ENDDO

  ! Get Twist Energy
  IF (TWIST.EQ.1) THEN
     call WRITHE(R,N,Wr)
     Tw=Lk-Wr
     EELAS(4)=((2*PI*Tw)**2)*LT/(2*L)
  ENDIF

  RETURN

END SUBROUTINE energy_elas

!---------------------------------------------------------------*
