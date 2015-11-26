! This file holds the user analysis routine and related functions
! and subroutines:
module analysis_supp
implicit none
!
  integer , parameter :: nycut = 48
  character (len=4) , dimension(nycut) :: cycut = &
    (/'.001','.002','.003','.004','.005','.006','.007','.008','.009', &
      '.01 ','.02 ','.03 ','.04 ','.05 ','.06 ','.07 ','.08 ','.09 ', &
      '.10 ','.11 ','.12 ','.13 ','.14 ','.15 ','.16 ','.17 ','.18 ','.19 ', &
      '.20 ','.21 ','.22 ','.23 ','.24 ','.25 ','.26 ','.27 ','.28 ','.29 ', &
      '.30 ','.31 ','.32 ','.33 ','.34 ','.35 ','.36 ','.37 ','.38 ','.39 '/)
  real(kind(1d0)) , dimension(nycut) :: ycut = &
    (/0.001d0,0.002d0,0.003d0,0.004d0,0.005d0,0.006d0,0.007d0,0.008d0,0.009d0, &
      0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.06d0,0.07d0,0.08d0,0.09d0, &
      0.10d0,0.11d0,0.12d0,0.13d0,0.14d0,0.15d0,0.16d0,0.17d0,0.18d0,0.19d0, &
      0.20d0,0.21d0,0.22d0,0.23d0,0.24d0,0.25d0,0.26d0,0.27d0,0.28d0,0.29d0, &
      0.30d0,0.31d0,0.32d0,0.33d0,0.34d0,0.35d0,0.36d0,0.37d0,0.38d0,0.39d0/)
  real(kind(1d0)) , dimension(nycut+1) :: ycutbins = &
    (/0.0005d0,0.0015d0,0.0025d0,0.0035d0,0.0045d0,0.0055d0,0.0065d0,0.0075d0,0.0085d0, &
      0.0095d0,0.015d0,0.025d0,0.035d0,0.045d0,0.055d0,0.065d0,0.075d0,0.085d0,0.095d0, &
      0.105d0,0.115d0,0.125d0,0.135d0,0.145d0,0.155d0,0.165d0,0.175d0,0.185d0,0.195d0, &
      0.205d0,0.215d0,0.225d0,0.235d0,0.245d0,0.255d0,0.265d0,0.275d0,0.285d0,0.295d0, &
      0.305d0,0.315d0,0.325d0,0.335d0,0.345d0,0.355d0,0.365d0,0.375d0,0.385d0,0.395d0/)
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
  do ijet=2,3
    call bookup_hist("jade ycut njet>"//cn(ijet),nycut,ycutbins)
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
    call JadeMine(ptrack,ntrack,ycut(iycut),njet)
    do ijet=2,3
    if (njet.eq.ijet) call fill_hist("jade ycut njet>"//cn(ijet), &
                                     ycut(iycut),wgt)
    end do
  end do
!
end subroutine analysis
