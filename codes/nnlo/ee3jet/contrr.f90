! This source contains routines related to the Doubly-Real contribution...
module RReal_data
use particles
implicit none
!
  logical :: ini_RR = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: prr
!
  real(kind(1d0)) :: smeRR_saved
!
contains
!
subroutine init_RReal()
use process
use reshuffle_data
implicit none
!
!
  integer :: istat
!
  if (.not.ini_RR) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
  numptrns = 7
!
  allocate(ptrns(nleg_born+2,numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with patterns allocation in RR..."
    stop
  end if
!
  allocate(prefacts(numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with prefacts allocation in RR..."
    stop
  end if
!
! Filling up...
! Note that the labeling for massless partons is simplified.
! The first massless quark flavor starts with 1 and if the 
! next is different than it is 2 otherwise 1:
!
! e+ e- -> d d~ g g g
  ptrns(1,1) = -11 ; ptrns(2,1) = 11 ; ptrns(3,1) =  1
  ptrns(4,1) =  -1 ; ptrns(5,1) =  0 ; ptrns(6,1) =  0 
  ptrns(7,1) =   0
! e+ e- -> u u~ g g g
  ptrns(1,2) = -11 ; ptrns(2,2) = 11 ; ptrns(3,2) =  2
  ptrns(4,2) =  -2 ; ptrns(5,2) =  0 ; ptrns(6,2) =  0
  ptrns(7,2) =   0
! e+ e- -> u u~ d d~ g
  ptrns(1,3) = -11 ; ptrns(2,3) = 11 ; ptrns(3,3) =  2
  ptrns(4,3) =  -2 ; ptrns(5,3) =  1 ; ptrns(6,3) = -1
  ptrns(7,3) =   0
! e+ e- -> u u~ u u~ g
  ptrns(1,4) = -11 ; ptrns(2,4) = 11 ; ptrns(3,4) =  2
  ptrns(4,4) =  -2 ; ptrns(5,4) =  2 ; ptrns(6,4) = -2
  ptrns(7,4) =   0
! e+ e- -> d d~ d d~ g
  ptrns(1,5) = -11 ; ptrns(2,5) = 11 ; ptrns(3,5) =  1
  ptrns(4,5) =  -1 ; ptrns(5,5) =  1 ; ptrns(6,5) = -1
  ptrns(7,5) =   0
! e+ e- -> u u~ c c~ g
  ptrns(1,6) = -11 ; ptrns(2,6) = 11 ; ptrns(3,6) =  2
  ptrns(4,6) =  -2 ; ptrns(5,6) =  4 ; ptrns(6,6) = -4
  ptrns(7,6) =   0
! e+ e- -> d d~ s s~ g
  ptrns(1,7) = -11 ; ptrns(2,7) = 11 ; ptrns(3,7) =  1
  ptrns(4,7) =  -1 ; ptrns(5,7) =  3 ; ptrns(6,7) = -3
  ptrns(7,7) =   0
!
  prefacts = 1d0
!
  print *,"The following patterns were created for the Real: "
  call PrintPatterns(nleg_born+2,numptrns,ptrns)
!
! We also allocate an array which will hold RR momenta and 
! flavor information:
  allocate(prr(nleg_born+2),stat=istat)
  if (istat.ne.0) then
    print *,"prr cannot be allocated..."
    stop
  end if
!
  ini_RR = .false.
!
end subroutine init_RReal
!
end module RReal_data
!
subroutine CalcRR(parts,smeRR)
use process
use particles
use RReal_data
use math
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: smeRR
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
    subroutine RRealSME(iptrn,prr,smeRR)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: prr
      real(kind(1d0)) , intent(out) :: smeRR
!
    end subroutine RRealSME
  end interface
!
!
  smeRR = 0d0
!
  if (flg_scaledep) then
    smeRR = smeRR_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  print *,"Before: "
!  call PrintParts(parts)
  call reshufflemomud(nleg_born+2,parts,prr,numptrns,ptrns,iptrn,prefact)
!  print *,"After: "
!  call PrintParts(prr)
!  print *,"pattern: ",iptrn
!  read(*,*)
!
! We call the SME routine:
  call RRealSME(iptrn,prr,smeRR)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  smeRR = smeRR * (4d0*pi)**5
!
  smeRR_saved = smeRR
!
!  print *,"smeRR: ",smeRR
!  read(*,*)
!
end subroutine CalcRR
