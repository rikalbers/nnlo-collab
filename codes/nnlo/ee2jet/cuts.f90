! This source contains routines related to physical cuts:
subroutine apply_cuts(p,icut)
use observables
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  integer , intent(out) :: icut
!
  integer :: npart
  integer :: ipart,jpart
  integer , parameter :: maxjet = 20
  integer , dimension(maxjet) :: jetvec
  integer :: ntrack,njet
  real(kind(1d0)) , dimension(4,maxjet) :: pjet
  real(kind(1d0)) , dimension(4,maxjet) :: ptrack
!
  real(kind(1d0)) , parameter :: PtMin = 5d0
  real(kind(1d0)) , parameter :: R = 0.4d0
  real(kind(1d0)) , parameter :: palg = 0d0
!
  logical , save :: init = .true.
! 
  if (init) then
    print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    print *,"@@@@@@@@@@@@@@@@@@@@@@@        $@@@@@@@@@@@@@@@@@@@@@@@@@@"
    print *,"@@@@@@@@@@@@@@@ Cuts are being used... @@@@@@@@@@@@@@@@@@@"
    print *,"@@@@@@@@@@@@@@@@@@@@@@@        @@@@@@@@@@@@@@@@@@@@@@@@@@@"
    print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
    init = .false.
  end if
! 
  icut = 0
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
  call JadeMine(ptrack,ntrack,0.01d0,njet)
!  call fastjeteegenjade(ptrack,ntrack,0.2d0,2d0,pjet,njet,jetvec)
!
!  do ipart=3,npart+2
!    print *,"ipart: ",ipart," pt: ",get_pt(p(ipart))
!    if (get_pt(p(ipart)).lt.PtMin) return
!    do jpart=ipart+1,npart+2
!      print *,"jpart: ",jpart, " dRij: ",get_dR(p(ipart),p(jpart))
!      if (get_dR(p(ipart),p(jpart)).lt.R) return
!    end do
!  end do
!
! ######################## Place cuts here: ###########################
  if (njet.lt.2) return
! #####################################################################
!
! If this line is reached the event passed all the cuts:
  icut = 1
!
end subroutine apply_cuts
!
subroutine cutfunc(p,cfunc)
use particles
use observables
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(out) :: cfunc
!
  real(kind(1d0)) :: Cpar,Dpar,EEC,T
!
!
  cfunc = 0d0
!
! ###################### Place function here: #########################
!  call CDpars(p,5,Cpar,Dpar)
!  cfunc = Cpar
!
!  EEC = CalcEEC(0,0,5,p)
!  cfunc = EEC
!
!  T = CalcThrust(p,5)
!  cfunc = 1 - T
  cfunc = 1d0
!
! #####################################################################
!
end subroutine cutfunc
!
subroutine JadeMine(ptrack,ntrack,ycut,njet)
implicit none
!
  integer , intent(in) :: ntrack
  real(kind(1d0)) , dimension(4,ntrack) , intent(in) :: ptrack
  real(kind(1d0)) , intent(in) :: ycut
  integer , intent(out) :: njet
!
  integer :: itrack,i,j,itrk,jtrk
  real(kind(1d0)) :: Q2,ymin,yij,piabs,pjabs
  real(kind(1d0)) , dimension(4) :: Q,p_tmp
  real(kind(1d0)) , dimension(4,ntrack) :: pjet
!
!
  njet = 0
!
  Q = 0
  do itrack=1,ntrack
    Q = Q + ptrack(:,itrack)
  end do
  Q2 = Q(4)**2 - Q(1)**2 - Q(2)**2 - Q(3)**2
!  print *,"Q2: ",Q2
!  print *,"sqrt: ",sqrt(Q2)
!
  pjet = ptrack
!
  ymin = 1d99
  do while (.true.)
    do i=1,ntrack-1
      if (sum(pjet(:,i)).eq.0) cycle
      piabs = sqrt(pjet(1,i)**2 + pjet(2,i)**2 + pjet(3,i)**2)
      do j=i+1,ntrack
        if (sum(pjet(:,j)).eq.0) cycle
        pjabs = sqrt(pjet(1,j)**2 + pjet(2,j)**2 + pjet(3,j)**2)
        yij = pjet(4,i)*pjet(4,j) &
            * (1 - (pjet(1,i)*pjet(1,j) + pjet(2,i)*pjet(2,j) &
            + pjet(3,i)*pjet(3,j))/(piabs*pjabs))/Q2
        if (yij.lt.ymin) then
          itrk = i
          jtrk = j
          ymin = yij
!          print *,"itrk,jtrk: ",itrk,jtrk
!          print *,"ymin changed: ",ymin
        end if
      end do
    end do
    if (ymin.gt.ycut) exit
    p_tmp = pjet(:,itrk) + pjet(:,jtrk)
    pjet(:,itrk) = p_tmp
    pjet(:,jtrk) = 0
    ymin = 1d99
  end do
!
  njet = 0
  do i=1,ntrack
    if (sum(pjet(:,i)).ne.0) njet = njet + 1
  end do
!
end subroutine JadeMine
