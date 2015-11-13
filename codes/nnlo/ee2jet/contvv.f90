! This source contains routines related to the Doubly-Virtual contribution...
module VVirt_data
use particles
implicit none
!
  logical :: ini_VV = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: pvv
!
  integer :: iptrn_saved
  integer :: jleg_saved
  real(kind(1d0)) :: muref
  real(kind(1d0)) , dimension(-4:2) :: VVLaurentAtMuRef
  real(kind(1d0)) , dimension(-4:2) :: VLaurentAtMuRef
  real(kind(1d0)) , dimension(-4:2) :: BLaurentAtMuRef
!
contains
!
subroutine init_VVirt()
use process
use reshuffle_data
implicit none
!
!
  integer :: istat
!
  if (.not.ini_VV) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
! For e+ e- -> q q~ g we have only two subprocesses:
  numptrns = 2
!
  allocate(ptrns(nleg_born,numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with patterns allocation in double-virtual..."
    stop
  end if
!
  allocate(prefacts(numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with prefacts allocation in double-virtual..."
    stop
  end if
!
! Filling up...
! Note that the labeling for massless partons is simplified.
! The first massless quark flavor starts with 1 and if the 
! next is different than it is 2 otherwise 1:
!
! e+ e- -> d d~ g
  ptrns(1,1) = -11 ; ptrns(2,1) = 11 ; ptrns(3,1) =  1
  ptrns(4,1) =  -1 ; ptrns(5,1) =  0
! e+ e- -> u u~ g
  ptrns(1,2) = -11 ; ptrns(2,2) = 11 ; ptrns(3,2) =  2
  ptrns(4,2) =  -2 ; ptrns(5,2) =  0
!
  prefacts = 1d0
!
  print *,"The following patterns were created for the Double-Virtual: "
  call PrintPatterns(nleg_born,numptrns,ptrns)
!
! We also allocate an array which will hold Born momenta and 
! flavor information for the Virtual:
  allocate(pvv(nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"pvv cannot be allocated..."
    stop
  end if
!
  ini_VV = .false.
!
end subroutine init_VVirt
!
end module VVirt_data
!
subroutine CalcVV(parts,VVirt,VVirtLaurent)
use process
use particles
use VVirt_data
use math
use flags
use scales
use QCDparams
!
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: VVirt
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VVirtLaurent
!
  integer :: iptrn
  real(kind(1d0)) :: prefact
  real(kind(1d0)) :: L
  real(kind(1d0)) , dimension(-4:2) :: Laurent_tmp
!
  interface
    subroutine reshufflemomud(n,p,pout,numptrn,ptrn,returnptrn, &
                              prefactor)
    use particles
    implicit none
!
      integer , intent(in) :: n,numptrn
      type(particle) , intent(in) , dimension(:) :: p
      type(particle) , intent(out) , dimension(:) :: pout
      integer , intent(in) , dimension(:,:) :: ptrn
      integer , intent(out) :: returnptrn
      real(kind(1d0)) , intent(out) :: prefactor
!
    end subroutine reshufflemomud
!
    subroutine VVirtSME(iptrn,pvv,mur,VVLaurent,VLaurent,BLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvv
      real(kind(1d0)) , intent(in) :: mur
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: &
        VVLaurent,VLaurent,BLaurent
!
    end subroutine VVirtSME
  end interface
!
!
  VVirt = 0d0
!
  if (flg_scaledep) then
! xi: mur2 / muref2:
    L = log(mur**2/muref**2)
!
    Laurent_tmp = 0
    Laurent_tmp(-4) = 0
    Laurent_tmp(-3) = 2*VVLaurentAtMuRef(-4)*L
    Laurent_tmp(-2) = 2*VVLaurentAtMuRef(-4)*L**2       &
                    + 2*VVLaurentAtMuRef(-3)*L          &
                    + 1*qcd_beta0*VLaurentAtMuRef(-2)*L
    Laurent_tmp(-1) = 4d0/3d0*VVLaurentAtMuRef(-4)*L**3          &
                    + 2*VVLaurentAtMuRef(-3)*L**2                &
                    + 3d0/2d0*qcd_beta0*VLaurentAtMuRef(-2)*L**2 &
                    + 2*VVLaurentAtMuRef(-2)*L                   &
                    + 1*qcd_beta0*VLaurentAtMuRef(-1)*L
    Laurent_tmp( 0) = 2d0/3d0*VVLaurentAtMuRef(-4)*L**4            &
                    + 4d0/3d0*VVLaurentAtMuRef(-3)*L**3            & 
                    + 7d0/6d0*qcd_beta0*VLaurentAtMuRef(-2)*L**3   &
                    + 2*VVLaurentAtMuRef(-2)*L**2                  &
                    + 3d0/2d0*qcd_beta0*VLaurentAtMuRef(-1)*L**2   &
                    + 0.25d0*qcd_beta0**2*BLaurentAtMuRef( 0)*L**2 &
                    + 2*VVLaurentAtMuRef(-1)*L                     &
                    + 1*qcd_beta0*VLaurentAtMuRef( 0)*L            &
                    + qcd_beta1/2d0*BLaurentAtMuRef( 0)*L 
!
    Laurent_tmp = Laurent_tmp + VVLaurentAtMuRef
!
    VVirt = Laurent_tmp(0) * cGamma * (4d0*pi)**5
!
    if (present(VVirtLaurent)) then
      VVirtLaurent = Laurent_tmp * cGamma * (4d0*pi)**5
    end if
!
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvv,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvv)
!
! Calling the double-virtual routine to calculate the full Laurent series of it at
! muref and the Taylor expansion of the Born, too:
  iptrn_saved = iptrn
  muref = mur
  call VVirtSME(iptrn,pvv,muref, &
                VVLaurentAtMuRef,VLaurentAtMuRef,BLaurentAtMuRef)
!
  VVirt = VVLaurentAtMuRef(0)
  VVirt = VVirt * cGamma
  VVirt = VVirt * (4d0*pi)**5
!
  if (present(VVirtLaurent)) then
    VVirtLaurent = VVLaurentAtMuRef * cGamma * (4d0*pi)**5
  end if
!
end subroutine CalcVV
