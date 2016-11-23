!This subroutine calculates the center of mass and radius of gyrations of 
!a single polymer

!Brad Krajina
!updated 2016/11/22

SUBROUTINE getRG(wlc_p,wlc_d,RGYRSQ)
  use params

  implicit none

  type(wlcsim_params), intent(in) :: wlc_p
  type(wlcsim_data), intent(in) :: wlc_d

  DOUBLE PRECISION  RCOM(3) ! Center of mass matrix
  DOUBLE PRECISION  RGYRSQ ! Radius of gyration squared
  DOUBLE PRECISION XGYR(wlc_p%nB,3) ! Bead positions relative to center of mass
  
  RCOM=0.0_dp
  RGYRSQ=0._dp
 
  !Calculate center of mass and radius of gyration squared
  RCOM(:)=SUM(wlc_d%r,DIM = 1)/DBLE(wlc_p%nB)
  XGYR(:,1)=wlc_d%r(:,1)-RCOM(1)
  XGYR(:,2)=wlc_d%r(:,2)-RCOM(2)
  XGYR(:,3)=wlc_d%r(:,3)-RCOM(3)
  RGYRSQ=SUM(XGYR**2.0_dp)/DBLE(wlc_p%nB)

  RETURN

END SUBROUTINE getRG








