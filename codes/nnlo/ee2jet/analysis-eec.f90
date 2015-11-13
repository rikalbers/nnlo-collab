! This file holds the user analysis routine and related functions
! and subroutines:
module analysis_supp
implicit none
!
  character (len=1) , dimension(9) :: cn = &
    (/'1','2','3','4','5','6','7','8','9'/)
!
end module analysis_supp
!
subroutine init_analysis
use histo
use analysis_supp
implicit none
!
!
!
!
! We have to initialize the histograms, we allow for 100 bins in each:
  call init_hist(100)
!
! EEC:
  call bookup_hist("EEC",0.02d0,-1d0,1d0)
!
end subroutine init_analysis
!
subroutine analysis(p,wgt)
use histo
use particles
use observables
use analysis_supp
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:) , intent(in) :: wgt
!
  integer :: ipart,jpart
  integer :: npart
  real(kind(1d0)) :: coschi,EEC
!
  npart = size(p)
!
  do ipart=3,npart-1
    if (abs(p(ipart)%flv).gt.5) cycle
    do jpart=ipart+1,npart
      if (abs(p(jpart)%flv).gt.5) cycle
      EEC = CalcEEC(ipart,jpart,5,p)
      coschi = get_costheta(p(ipart),p(jpart))
      call fill_hist("EEC",coschi,EEC*wgt)
    end do
  end do
!
!
end subroutine analysis
