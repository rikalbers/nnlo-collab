! This file holds the user analysis routine and related functions
! and subroutines:
module analysis_supp
implicit none
!
  integer , parameter :: nycut = 5
  real(kind(1d0)) , dimension(nycut) :: ycuts = &
    (/0.001d0,0.005d0,0.01d0,0.05d0,0.1d0/)
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
  integer :: ijet
!
!
! We have to initialize the histograms:
  call init_hist(100)
!
! sigma:
  call bookup_hist("sigma",1d0,0d0,1d0)
!
! Cross sections with different y cuts using JADE:
  do ijet=3,5
    call bookup_hist("jade ycut njet>"//cn(ijet),1d0,0.5d0,5.5d0)
  end do
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
  integer :: ipart,ijet,itrack,iycut
  integer :: npart
  integer :: ntrack,njet
  integer , parameter :: maxjet = 20
  integer , dimension(maxjet) :: jetvec
  real(kind(1d0)) , dimension(4,maxjet) :: pjet
  real(kind(1d0)) , dimension(4,maxjet) :: ptrack
!
  logical , dimension(nycut) :: passed
!
  call fill_hist("sigma",0.5d0,wgt)
!
  ntrack = 0
  njet = 0
  jetvec = 0
!
! Number of final state particles:
  npart = size(p)
! We have to construct the tracks:
  do ipart=3,npart
    if (abs(p(ipart)%flv).gt.5) cycle
    ntrack = ntrack + 1
    ptrack(1,ntrack) = p(ipart)%p%px
    ptrack(2,ntrack) = p(ipart)%p%py
    ptrack(3,ntrack) = p(ipart)%p%pz
    ptrack(4,ntrack) = p(ipart)%p%E
  end do
!
  do iycut=1,nycut
    call JadeMine(ptrack,ntrack,ycuts(iycut),njet)
    do ijet=3,5
    if (njet.ge.ijet) call fill_hist("jade ycut njet>"//cn(ijet), &
                                     dble(iycut),wgt)
    end do
  end do
!
end subroutine analysis
