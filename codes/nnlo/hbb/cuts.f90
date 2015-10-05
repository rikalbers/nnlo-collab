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
  call fastjeteegenjade(ptrack,ntrack,0.01d0,2d0,pjet,njet,jetvec)
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
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(out) :: cfunc
!
!
  cfunc = 0d0
!
! ###################### Place function here: #########################
! #####################################################################
!
end subroutine cutfunc
