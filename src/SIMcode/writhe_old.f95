  !Calculate the Writhe of a curve using the method of Klenin (2000)

  SUBROUTINE WRITHEOLD(R,N,Wr)
  IMPLICIT NONE

  REAL, PARAMETER :: PI=3.141592654 ! Value of pi

  
  INTEGER, INTENT(IN) :: N                 ! Number of beads
  DOUBLE PRECISION, INTENT(IN) :: R(N,3)  ! Positions
  !Geometric variables
  DOUBLE PRECISION  r1(3)                  ! Position bead 1
  DOUBLE PRECISION  r2(3)                  ! Position bead 2
  DOUBLE PRECISION  r12(3)                 ! Relative position vector
  DOUBLE PRECISION  s1                     ! Length of segment 1
  DOUBLE PRECISION  s2                     ! Length of segment 2
  DOUBLE PRECISION  beta                   ! Angle between tangents
  DOUBLE PRECISION  e1(3)                  ! Tangent of first segment
  DOUBLE PRECISION  e2(3)                  ! Tangent of second segment
  DOUBLE PRECISION  cosB
  DOUBLE PRECISION  sin2B
  DOUBLE PRECISION  e3(3)
  DOUBLE PRECISION  OMEGA(N,N)             ! Matrix of solid angles
  DOUBLE PRECISION  a0
  DOUBLE PRECISION  a1
  DOUBLE PRECISION  a2
  DOUBLE PRECISION  a3
 
  !Variables for writhe integral
  DOUBLE PRECISION  t1
  DOUBLE PRECISION  t2
  DOUBLE PRECISION  F1
  DOUBLE PRECISION  F2
  DOUBLE PRECISION  F3
  DOUBLE PRECISION  F4
  
  !Counter variables
  INTEGER  I,IP1
  INTEGER  J,JP1

  DOUBLE PRECISION, INTENT(OUT) :: Wr      ! Writhe
  

  OMEGA=0
  Wr=0

  DO  I=2,N
     DO J=1,I-1
        
        IF (I.EQ.N) THEN
           IP1=1
        ELSE 
           IP1=I+1
        ENDIF

        JP1=J+1
        
        r1=R(i,:)
        r2=R(j,:)
        r12=r2-r1
        s1=SQRT(SUM((R(IP1,:)-R(I,:))**2))
        s2=SQRT(SUM((R(JP1,:)-R(J,:))**2))
        e1=(R(IP1,:)-R(I,:))/s1
        e2=(R(JP1,:)-R(J,:))/s2
        cosB=DOT_PRODUCT(e1,e2)
         
        sin2B=1-(cosB**2)

        e3(1)=e1(2)*e2(3)-e1(3)*e2(2)
        e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
        e3(3)=e1(1)*e2(2)-e1(2)*e2(1)

        a1=DOT_PRODUCT(r12,e2*cosB-e1)/(sin2B)
        a2=DOT_PRODUCT(r12,e2-e1*cosB)/sin2B
        a0=DOT_PRODUCT(r12,e3)/sin2B
        IF(A0.NE.0..AND.SIN2B.NE.0.) THEN
           t1=a1+s1
           t2=a2+s2
           F1=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4*PI)
           t1=a1+s2
           t2=a2
           F2=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4*PI)
           t1=a1
           t2=a2+s2
           F3=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4*PI)
           t1=a1
           t2=a2
           F4=-ATAN((t1*t2+(a0**2)*cosB)/(a0*SQRT((t1**2+t2**2-2*t1*t2*cosB+(a0**2)*sin2B))))/(4*PI)

           Wr=Wr+F1-F2-F3+F4
        ENDIF 
     ENDDO
  ENDDO

  Wr=Wr*2

  RETURN
  END SUBROUTINE





        
        
  

  
