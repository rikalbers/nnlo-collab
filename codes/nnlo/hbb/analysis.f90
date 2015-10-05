! This file holds the user analysis routine and related functions
! and subroutines:
subroutine init_analysis
use histo
implicit none
!
!
  integer , parameter :: nptbins = 19
  real(kind(1d0)) , dimension(nptbins+1) :: ptbins =   &
    (/0d0,10d0,20d0,30d0,40d0,50d0,60d0,70d0,80d0,90d0, &
      100d0,120d0,140d0,160d0,180d0,200d0,250d0,        &
      300d0,400d0,500d0/)
  integer :: ipart,jpart
  integer , parameter :: npart = 5
  character (len=1) , dimension(9) :: cn = &
    (/'1','2','3','4','5','6','7','8','9'/)
!
! We have to initialize the histograms:
  call init_hist(100)
!
! sigma:
  call bookup_hist("sigma",1d0,0d0,1d0)
!
! Jet pts: 
  do ipart=1,npart
    call bookup_hist("pt j"//cn(ipart),1d0,0d0,100d0)
  end do
!
! Jet rapidities:
  do ipart=1,npart
    call bookup_hist("y j"//cn(ipart),0.1d0,-4d0,4d0)
  end do
!
! invariant masses:
  do ipart=1,npart-1
    do jpart=ipart+1,npart
      call bookup_hist("m j"//cn(ipart)//"j"//cn(jpart),1d0,0d0,100d0)
    end do
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
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:) , intent(in) :: wgt
!
  integer :: ijet,jjet,ipart
  integer :: npart
  integer :: ntrack,njet
  integer , parameter :: maxjet = 20
  integer , dimension(maxjet) :: jetvec
  real(kind(1d0)) , dimension(4,maxjet) :: pjet
  real(kind(1d0)) , dimension(4,maxjet) :: ptrack
  real(kind(1d0)) , dimension(maxjet) :: ptj
  real(kind(1d0)) , dimension(maxjet) :: yj
  real(kind(1d0)) , dimension(maxjet,maxjet) :: mjj
  type(mom) , dimension(maxjet) :: pjets
!
  real(kind(1d0)) , parameter :: PtMin = 5d0
  real(kind(1d0)) , parameter :: R = 0.4d0
  real(kind(1d0)) , parameter :: palg = 0d0
!
  character (len=1) , dimension(9) :: cn = &
    (/'1','2','3','4','5','6','7','8','9'/)
!
  call fill_hist("sigma",0.5d0,wgt)
!
  ntrack = 0
  njet = 0
  jetvec = 0
!
! Number of final state particles:
  npart = size(p) - 2
! We have to construct the tracks:
  do ipart=3,npart+2
    if (abs(p(ipart)%flv).gt.5) cycle
    ntrack = ntrack + 1
    ptrack(1,ntrack) = p(ipart)%p%px
    ptrack(2,ntrack) = p(ipart)%p%py
    ptrack(3,ntrack) = p(ipart)%p%pz
    ptrack(4,ntrack) = p(ipart)%p%E
  end do
!
!  call fastjetppgenkt(ptrack,ntrack,R,palg,PtMin,pjet,njet,jetvec)
  call fastjeteegenjade(ptrack,ntrack,0.01d0,2d0,pjet,njet,jetvec)
!
  if (njet.lt.3) return
!
  do ijet=1,njet
    pjets(ijet)%E  = pjet(4,ijet)
    pjets(ijet)%px = pjet(1,ijet)
    pjets(ijet)%py = pjet(2,ijet)
    pjets(ijet)%pz = pjet(3,ijet)
  end do
!
! We calculate the pts and ys for all the jets:
  do ijet=1,njet
    ptj(ijet) = get_pt(pjets(ijet))
    yj(ijet)  = get_rapidity(pjets(ijet))
  end do
!
  do ijet=1,njet
    ptj(ijet) = get_pt(pjets(ijet))
    yj(ijet)  = get_rapidity(pjets(ijet))
    do jjet=ijet+1,njet
      mjj(ijet,jjet) = get_invm(pjets(ijet) + pjets(jjet))
    end do
  end do
!
  do ijet=1,njet
    call fill_hist("pt j"//cn(ijet),ptj(ijet),wgt)
    call fill_hist("y j"//cn(ijet),yj(ijet),wgt)
  end do
!
  do ijet=1,njet-1
    do jjet=ijet+1,njet
      call fill_hist("m j"//cn(ijet)//"j"//cn(jjet), &
                     mjj(ijet,jjet),wgt)
    end do
  end do
!
end subroutine analysis
