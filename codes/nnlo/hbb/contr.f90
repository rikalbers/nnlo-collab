! This source contains routines related to the Real contribution...
module Real_data
use particles
implicit none
!
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: preal
!
  real(kind(1d0)) :: smeR_saved
!
contains
!
subroutine init_Real()
use process
use reshuffle_data
implicit none
!
!
  integer :: istat
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
  numptrns = 1
!
  allocate(ptrns(nleg_born+1,numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with patterns allocation in real..."
    stop
  end if
!
  allocate(prefacts(numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with prefacts allocation in real..."
    stop
  end if
!
! Filling up...
! Note that the labeling for massless partons is simplified.
! The first massless quark flavor starts with 1 and if the 
! next is different than it is 2 otherwise 1:
!
! e+ e- -> b b~ g
  ptrns(1,1) = -11 ; ptrns(2,1) = 11 ; ptrns(3,1) =  5
  ptrns(4,1) =  -5 ; ptrns(5,1) =  0
!
  prefacts = 1d0
!
  print *,"The following patterns were created for the Real: "
  call PrintPatterns(nleg_born+1,numptrns,ptrns)
!
! We also allocate an array which will hold Real momenta and 
! flavor information:
  allocate(preal(nleg_born+1),stat=istat)
  if (istat.ne.0) then
    print *,"preal cannot be allocated..."
    stop
  end if
!
end subroutine init_Real
!
end module Real_data
!
subroutine CalcR(parts,smeR)
use process
use particles
use Real_data
use math
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: smeR
!
  integer :: iptrn
  real(kind(1d0)) :: prefact
!
  interface
    subroutine reshufflemom(n,p,pout,numptrn,ptrn,returnptrn,prefactor)
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
    end subroutine reshufflemom
!
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
    subroutine RealSME(iptrn,preal,smeR)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: preal
      real(kind(1d0)) , intent(out) :: smeR
!
    end subroutine RealSME
  end interface
!
!
  smeR = 0d0
!
  if (flg_scaledep) then
    smeR = smeR_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  print *,"Before: "
!  call PrintParts(parts)
  call reshufflemomud(nleg_born+1,parts,preal,numptrns,ptrns,iptrn,prefact)
!  print *,"After: "
!  call PrintParts(preal)
!
! We call the SME routine:
  call RealSME(iptrn,preal,smeR)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
! The real SME is such that only \alpha_S is factored out:
  smeR = smeR
!
  smeR_saved = smeR
!
!  print *,"smeR: ",smeR
!
end subroutine CalcR
