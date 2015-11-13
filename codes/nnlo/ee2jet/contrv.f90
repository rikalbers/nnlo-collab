! This source contains routines related to the Virtual contribution...
module RVirt_data
use particles
implicit none
!
  logical :: ini_RV = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: prv
!
  integer :: iptrn_saved
  real(kind(1d0)) :: smeRV_muindep_saved
!
contains
!
subroutine init_RVirt()
use process
use reshuffle_data
implicit none
!
!
  integer :: istat
!
  if (.not.ini_RV) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
  numptrns = 7
!
  allocate(ptrns(nleg_born+1,numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with patterns allocation in real-virtual..."
    stop
  end if
!
  allocate(prefacts(numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with prefacts allocation in real-virtual..."
    stop
  end if
!
! Filling up...
! Note that the labeling for massless partons is simplified.
! The first massless quark flavor starts with 1 and if the 
! next is different than it is 2 otherwise 1:
!
! e+ e- -> d d~ g g
  ptrns(1,1) = -11 ; ptrns(2,1) = 11 ; ptrns(3,1) =  1
  ptrns(4,1) =  -1 ; ptrns(5,1) =  0 ; ptrns(6,1) =  0
! e+ e- -> u u~ g g
  ptrns(1,2) = -11 ; ptrns(2,2) = 11 ; ptrns(3,2) =  2
  ptrns(4,2) =  -2 ; ptrns(5,2) =  0 ; ptrns(6,2) =  0
! e+ e- -> u u~ d d~
  ptrns(1,3) = -11 ; ptrns(2,3) = 11 ; ptrns(3,3) =  2
  ptrns(4,3) =  -2 ; ptrns(5,3) =  1 ; ptrns(6,3) = -1
! e+ e- -> u u~ u u~
  ptrns(1,4) = -11 ; ptrns(2,4) = 11 ; ptrns(3,4) =  2
  ptrns(4,4) =  -2 ; ptrns(5,4) =  2 ; ptrns(6,4) = -2
! e+ e- -> d d~ d d~
  ptrns(1,5) = -11 ; ptrns(2,5) = 11 ; ptrns(3,5) =  1
  ptrns(4,5) =  -1 ; ptrns(5,5) =  1 ; ptrns(6,5) = -1
! e+ e- -> u u~ c c~
  ptrns(1,6) = -11 ; ptrns(2,6) = 11 ; ptrns(3,6) =  2
  ptrns(4,6) =  -2 ; ptrns(5,6) =  4 ; ptrns(6,6) = -4
! e+ e- -> d d~ s s~
  ptrns(1,7) = -11 ; ptrns(2,7) = 11 ; ptrns(3,7) =  1
  ptrns(4,7) =  -1 ; ptrns(5,7) =  3 ; ptrns(6,7) = -3
!
  prefacts = 1d0
!
  print *,"The following patterns were created for the Real-Virtual: "
  call PrintPatterns(nleg_born+1,numptrns,ptrns)
!
! We also allocate an array which will hold Born momenta and 
! flavor information for the Virtual:
  allocate(prv(nleg_born+1),stat=istat)
  if (istat.ne.0) then
    print *,"prv cannot be allocated..."
    stop
  end if
!
  ini_RV = .false.
!
end subroutine init_RVirt
!
end module RVirt_data
!
subroutine CalcRV(parts,RVirt,RVLaurent)
use process
use particles
use RVirt_data
use math
use flags
use scales
!
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: RVirt
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: RVLaurent
!
  integer :: iptrn
  real(kind(1d0)) :: prefact
  real(kind(1d0)) :: RVirt_dep,RVirt_indep
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
    subroutine RVirtSME(iptrn,prv,mur,mode,RVirt,RVirt_dep,RVirt_indep, &
                        RVLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: prv
      real(kind(1d0)) , intent(in) :: mur
      character , intent(in) :: mode
      real(kind(1d0)) , intent(out) :: RVirt
      real(kind(1d0)) , intent(out) :: RVirt_dep,RVirt_indep
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: RVLaurent
!
    end subroutine RVirtSME
  end interface
!
!
  RVirt = 0d0
!
! If flg_scaledep is true we are having the same PS point with the
! same flavor configuration hence we only have to reevaluate the mu
! dependent term only:
  if (flg_scaledep) then
! pvirt is declared inside the module hence still usable, the 
! pattern number is also saved:
    iptrn = iptrn_saved
    call RVirtSME(iptrn,prv,mur,'d',RVirt,RVirt_dep,RVirt_indep)
! The total real-virtual is built up from the dependent and independent
! terms:
    RVirt = smeRV_muindep_saved + RVirt_dep
! Normalization factors:
! c_Gamma:
    RVirt = RVirt * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    RVirt = RVirt * (4d0*pi)**5
! For debugging purposes it is possible that the full Laurent series
! is needed for the virtual if this is the case we simply calculate
! everything:
    if (present(RVLaurent)) then
      call RVirtSME(iptrn,prv,mur,'f',RVirt,RVirt_dep,RVirt_indep,RVLaurent)
      RVirt = RVirt * cGamma
      RVLaurent = RVLaurent * cGamma
      RVirt = RVirt * (4d0*pi)**5
      RVLaurent = RVLaurent * (4d0*pi)**5
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born+1,parts,prv,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(prv)
!
! For the real-virtual part we also have to supply the renormalization 
! scale to the routine responsible for calculating the SME:
  if (present(RVLaurent)) then
    call RVirtSME(iptrn,prv,mur,'f',RVirt,RVirt_dep,RVirt_indep,RVLaurent)
  else
! If we have a new PS point we reevaluate the whole real-virtual part,
! that is the mu dependent and independent part as well:
    if (.not.flg_scaledep) then
      call RVirtSME(iptrn,prv,mur,'b',RVirt,RVirt_dep,RVirt_indep)
! We stash away the independent part:
      smeRV_muindep_saved = RVirt_indep
! We also stash away the pattern used:
      iptrn_saved = iptrn
    end if
  end if
!
! Normalization coming from c_Gamma:
  RVirt = RVirt * cGamma
  if (present(RVLaurent)) RVLaurent = RVLaurent * cGamma
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  RVirt = RVirt * (4d0*pi)**5
  if (present(RVLaurent)) RVLaurent = RVLaurent * (4d0*pi)**5
!
end subroutine CalcRV
