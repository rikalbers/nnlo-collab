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
! We have to initialize the histograms:
  call init_hist(100)
!
! sigma:
  call bookup_hist("sigma",1d0,0d0,1d0)
!
! C parameter:
  call bookup_hist("Cpar",0.02d0,0d0,1d0)
! C moments:
  call bookup_hist("Cmoments",1d0,0.5d0,10.5d0)
! sigma with C:
  call bookup_hist("sigmaC",1d0,0d0,1d0)
!
!  call print_hist
!
!  stop
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
  integer :: ipart,imom
  integer :: npart
!
  real(kind(1d0)) :: Cpar,Dpar,S
!
  call fill_hist("sigma",0.5d0,wgt)
!
! Calculate the parameters:
  call CDpars(p,5,Cpar,Dpar)
  call Sphericity(p,5,S)
!
  call fill_hist("Cpar",Cpar,wgt)
  call fill_hist("sigmaC",0.5d0,Cpar*wgt)
!
  do imom=1,10
    call fill_hist("Cmoments",dble(imom),Cpar**imom*wgt)
  end do
!
end subroutine analysis
