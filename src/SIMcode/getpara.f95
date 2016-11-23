! *---------------------------------------------------------------*
!     This subroutine determines the effective elastic parameters for 
!     a discrete-stretchable shearable worm-like chain, given the discretization
!     length, DEL

      
subroutine getpara(wlc_p)

  use params
  implicit none

  !Simulation variables
  type(wlcsim_params) :: wlc_p  !Structure containing static simulation parameters
  real(dp) del !number of persistence lengths per chain inter-bead segment

  !Counting variables
  integer i,ind,crs

  !Variables for interpolating tabulated values
  real(dp) PVEC(679,8)
  real(dp) M



  !Determine the number of persistence lengths per  segment

  wlc_p%L = wlc_p%L/wlc_p%lp

  IF (wlc_p%ring.EQ.0) THEN
     del=wlc_p%L/(dble(wlc_p%nB)-1.0_dp)
  ELSE
     del=wlc_p%L/(dble(wlc_p%nB))
  ENDIF

  !     Load the tabulated parameters

  OPEN (UNIT=5,FILE='input/dssWLCparams',STATUS='OLD')
  DO  I=1,679
     READ(5,*) PVEC(I,1),PVEC(I,2),PVEC(I,3),PVEC(I,4),PVEC(I,5),PVEC(I,6),PVEC(I,7),PVEC(I,8)
  ENDDO
  CLOSE(5)

  if (DEL.LT.PVEC(1,1)) then
     DEL=PVEC(1,1)
  endif
  if (DEL.GT.PVEC(679,1)) then
     DEL=PVEC(679,1)
  endif

  CRS=0
  IND=1
  do while (CRS.EQ.0)
     if (DEL.LE.PVEC(IND,1)) then
        CRS=1
     else
        IND=IND+1
     endif
  enddo

  !     Perform linear interpolations 

  I=2 
  M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
  wlc_p%eb=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

  I=3 
  M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
  wlc_p%gam=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

  I=4
  M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
  wlc_p%epar=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

  I=5
  M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
  wlc_p%EPERP=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

  I=6
  M=(PVEC(IND,I)-PVEC(IND-1,I))/(PVEC(IND,1)-PVEC(IND-1,1))
  wlc_p%ETA=M*(DEL-PVEC(IND,1))+PVEC(IND,I)

  wlc_p%EB=wlc_p%LP*wlc_p%EB/(DEL*wlc_p%LP)
  wlc_p%EPAR=wlc_p%EPAR/(DEL*wlc_p%LP*wlc_p%LP)
  wlc_p%EPERP=wlc_p%EPERP/(DEL*wlc_p%LP*wlc_p%LP)
  wlc_p%GAM=DEL*wlc_p%LP*wlc_p%GAM
  wlc_p%ETA=wlc_p%ETA/wlc_p%LP
  
  wlc_p%L=wlc_p%L*wlc_p%LP

  RETURN     
END SUBROUTINE getpara

!---------------------------------------------------------------*
