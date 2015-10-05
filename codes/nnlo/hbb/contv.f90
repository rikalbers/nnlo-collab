! This source contains routines related to the Virtual contribution...
module Virt_data
use particles
implicit none
!
  logical :: ini = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: pvirt
!
  integer :: iptrn_saved
  real(kind(1d0)) :: smeV_muindep_saved
!
contains
!
subroutine init_Virt()
use process
use reshuffle_data
implicit none
!
!
  integer :: istat
!
  if (.not.ini) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
! For e+ e- -> b b~ we have only one subprocess:
  numptrns = 1
!
  allocate(ptrns(nleg_born,numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with patterns allocation in born..."
    stop
  end if
!
  allocate(prefacts(numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with prefacts allocation in born..."
    stop
  end if
!
! Filling up...
! Note that the labeling for massless partons is simplified.
! The first massless quark flavor starts with 1 and if the 
! next is different than it is 2 otherwise 1:
!
! e+ e- -> b b~
  ptrns(1,1) = -11 ; ptrns(2,1) = 11 ; ptrns(3,1) =  5
  ptrns(4,1) =  -5
!
  prefacts = 1d0
!
  print *,"The following patterns were created for the Virtual: "
  call PrintPatterns(nleg_born,numptrns,ptrns)
!
! We also allocate an array which will hold Born momenta and 
! flavor information for the Virtual:
  allocate(pvirt(nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"pborn cannot be allocated..."
    stop
  end if
!
  ini = .false.
!
end subroutine init_Virt
!
end module Virt_data
!
subroutine CalcV(parts,Virt,VirtLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: Virt
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
!
  integer :: iptrn
  real(kind(1d0)) :: prefact
  real(kind(1d0)) :: Virt_dep,Virt_indep
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
    subroutine VirtSME(iptrn,pvirt,mur,mode,Virt,Virt_dep,Virt_indep, &
                       VirtLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      character , intent(in) :: mode
      real(kind(1d0)) , intent(out) :: Virt
      real(kind(1d0)) , intent(out) :: Virt_dep,Virt_indep
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
!
    end subroutine VirtSME
  end interface
!
!
  Virt = 0d0
  if (present(VirtLaurent)) VirtLaurent = 0d0
!
! If flg_scaledep is true we are having the same PS point with the
! same flavor configuration hence we only have to reevaluate the mu
! dependent term only:
  if (flg_scaledep) then
! pvirt is declared inside the module hence still usable, the 
! pattern number is also saved:
    iptrn = iptrn_saved
    call VirtSME(iptrn,pvirt,mur,'d',Virt,Virt_dep,Virt_indep)
! The total virtual is built up from the dependent and independent
! terms:
    Virt = smeV_muindep_saved + Virt_dep
! For debugging purposes it is possible that the full Laurent series
! is needed for the virtual if this is the case we simply calculate
! everything:
    if (present(VirtLaurent)) then
      call VirtSME(iptrn,pvirt,mur,'f',Virt,Virt_dep,Virt_indep,VirtLaurent)
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvirt,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvirt)
!
! For the virtual part we also have to supply the renormalization scale to the
! routine responsible for calculating the SME:
  if (present(VirtLaurent)) then
    call VirtSME(iptrn,pvirt,mur,'f',Virt,Virt_dep,Virt_indep,VirtLaurent)
  else
! If we have a new PS point we reevaluate the whole virtual part,
! that is the mu dependent and independent part as well:
    if (.not.flg_scaledep) then
      call VirtSME(iptrn,pvirt,mur,'b',Virt,Virt_dep,Virt_indep)
! We stash away the independent part:
      smeV_muindep_saved = Virt_indep
! We also stash away the pattern used:
      iptrn_saved = iptrn
    end if
  end if
!
! All normalizations with an alpha_S factored out are in the Virtual part
!
  Virt = -Virt
!
end subroutine CalcV
