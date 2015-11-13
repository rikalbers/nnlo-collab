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
!
!
! We have to initialize the histograms, we allow for 100 bins in each:
  call init_hist(100)
!
! 1 - Thrust:
  call bookup_hist("1-T",0.005d0,0d0,0.4d0)
! Heavy jet mass:
  call bookup_hist("rho",0.005d0,0d0,0.4d0)
! Total jet broadening:
  call bookup_hist("BT",0.005d0,0d0,0.4d0)
! Wide jet broadening:
  call bookup_hist("BW",0.005d0,0d0,0.4d0)
! C parameter:
  call bookup_hist("Cpar",0.01d0,0d0,1d0)
! jet transition variable: Y3:
  call bookup_hist("Y3",0.1d0,0d0,10d0)
! The histos corresponding to moments are commented out,
! as far as we can make distributions moments are not
! interesting at all:
! <(1-T)^n>:
!  call bookup_hist("1-Tn",1d0,0.5d0,10.5d0)
! <rho^n>:
!  call bookup_hist("rhon",1d0,0.5d0,10.5d0)
! <BT^n>:
!  call bookup_hist("BTn",1d0,0.5d0,10.5d0)
! <BW^n>:
!  call bookup_hist("BWn",1d0,0.5d0,10.5d0)
! <C^n>:
!  call bookup_hist("Cn",1d0,0.5d0,10.5d0)
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
  integer :: ipart,imom,iycut
  integer :: npart,ntrack,njet
  integer , parameter :: maxjet = 5
  real(kind(1d0)) , dimension(4,maxjet) :: ptrack
  real(kind(1d0)) , dimension(20) :: EvntShapes
  real(kind(1d0)) , dimension(4,3) :: ThrustDirs
  integer :: err
!
  real(kind(1d0)) :: T,Tmaj,Tmin,Oblat,Cpar,Dpar,MH2,md2,BW,BT, &
                     Y3,rho
!
  ntrack = 0
  njet   = 0
!
  npart = size(p)
! Constructing the tracks:
  do ipart=3,npart
    if (abs(p(ipart)%flv).gt.5) cycle
    ntrack = ntrack + 1
    ptrack(1,ntrack) = p(ipart)%p%px
    ptrack(2,ntrack) = p(ipart)%p%py
    ptrack(3,ntrack) = p(ipart)%p%pz
    ptrack(4,ntrack) = p(ipart)%p%E
  end do
!
  call EvShapes(ptrack,ntrack,.false.,EvntShapes,ThrustDirs,err)
!
  T     = EvntShapes(1)
  Tmaj  = EvntShapes(2)
  Tmin  = EvntShapes(3)
  Oblat = EvntShapes(4)
  Cpar  = EvntShapes(5)
  MH2   = EvntShapes(6)
  md2   = EvntShapes(7)
  BT    = EvntShapes(8)
  BW    = EvntShapes(9)
  Y3    = EvntShapes(10)
  Dpar  = EvntShapes(11)
!
  rho   = MH2
!
  call fill_hist("1-T",1d0 - T,(1d0-T)*wgt)
  call fill_hist("rho",rho,rho*wgt)
  call fill_hist("BT",BT,BT*wgt)
  call fill_hist("BW",BW,BW*wgt)
  call fill_hist("Cpar",Cpar,Cpar*wgt)
  call fill_hist("Y3",-log(Y3),Y3*wgt)
!
!  do imom=1,10
!    call fill_hist("1-Tn",dble(imom),(1d0-T)**imom*wgt)
!    call fill_hist("rhon",dble(imom),rho**imom*wgt)
!    call fill_hist("BTn",dble(imom),BT**imom*wgt)
!    call fill_hist("BWn",dble(imom),BW**imom*wgt)
!    call fill_hist("Cn",dble(imom),Cpar**imom*wgt)
!  end do
!
end subroutine analysis
