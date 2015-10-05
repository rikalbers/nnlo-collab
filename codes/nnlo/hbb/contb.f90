! This source contains routines related to the Born contribution...
module Born_data
use particles
implicit none
!
  logical :: ini = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: pborn
!
  real(kind(1d0)) :: smeB_saved
  real(kind(1d0)) , dimension(0:3,0:3) :: smeBmunu_saved
  real(kind(1d0)) , dimension(:,:) , allocatable :: smeBij_saved
!
contains
!
subroutine init_Born()
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
! For e+ e- -> b b~ we have only two subprocesses:
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
  print *,"The following patterns were created for the Born: "
  call PrintPatterns(nleg_born,numptrns,ptrns)
!
! We also allocate an array which will hold Born momenta and 
! flavor information:
  allocate(pborn(nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"pborn cannot be allocated..."
    stop
  end if
!
  allocate(smeBij_saved(nleg_born,nleg_born),stat=istat)
  if (istat.ne.0) then
    print *,"smeBij_saved cannot be allocated..."
    stop
  end if
!
  ini = .false.
!
end subroutine init_Born
!
end module Born_data
!
subroutine CalcB(parts,Born)
use process
use particles
use Born_data
use math
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: Born
!
  integer :: iptrn
  real(kind(1d0)) :: prefact
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
    subroutine BornSME(iptrn,pborn,Born)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pborn
      real(kind(1d0)) , intent(out) :: Born
!
    end subroutine BornSME
  end interface
!
!
  Born = 0d0
!
! If flg_scaledep is true we do not have to reevaluate the Born provided
! by sitting at the same PS point and there is no \mu dependent terms
! inside:
  if (flg_scaledep) then
    Born = smeB_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  call BornSME(iptrn,pborn,Born)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
! but this time the Born is independent of \alpha_S and \alpha_{EM}:
  Born = Born
!
! We also store the Born which can be useful for scale dependency studies:
  smeB_saved = Born
!
!  print *,"Born: ",Born
!  read(*,*)
!
end subroutine CalcB
!
subroutine CalcBmunu(ileg,parts,Bmunu)
use process
use particles
use Born_data
use math
use flags
implicit none
!
  integer , intent(in) :: ileg
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
  integer :: iptrn,jleg
  real(kind(1d0)) :: prefact
!
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
    subroutine BmunuSME(iptrn,ileg,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn,ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine BmunuSME
  end interface
!
  Bmunu = 0d0
!
  if (flg_scaledep) then
    Bmunu = smeBmunu_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  print *,"For this process there is no Bmunu..."
  stop
!
end subroutine CalcBmunu
!
subroutine CalcBij(parts,Bij)
use process
use particles
use Born_data
use math
use misc
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
  integer :: iptrn,jleg
  real(kind(1d0)) :: prefact
!
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
    subroutine BijSME(iptrn,p,Bij)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
!
    end subroutine BijSME
  end interface
!
  Bij = 0d0
!
  if (flg_scaledep) then
    Bij = smeBij_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  call BijSME(iptrn,pborn,Bij)
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Bij:
  call SwapColor(nleg_born,parts,pborn,Bij)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
! In the SMEs alpha_S is properly factored out:
  Bij = Bij
!
  smeBij_saved = Bij
!
end subroutine CalcBij
