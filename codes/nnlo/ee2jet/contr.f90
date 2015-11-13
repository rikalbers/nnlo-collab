! This source contains routines related to the Real contribution...
module Real_data
use particles
implicit none
!
  logical :: ini_R = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: preal
!
  real(kind(1d0)) :: smeR_saved
  real(kind(1d0)) , dimension(0:3,0:3) :: smeRmunu_saved
  real(kind(1d0)) , dimension(:,:) , allocatable :: smeRij_saved
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
  if (.not.ini_R) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
  numptrns = 7
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
  allocate(smeRij_saved(nleg_born+1,nleg_born+1),stat=istat)
  if (istat.ne.0) then
    print *,"smeRij_saved cannot be allocated..."
    stop
  end if
!
  ini_R = .false.
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
  smeR = smeR * (4d0*pi)**4
!
  smeR_saved = smeR
!
!  print *,"smeR: ",smeR
!
end subroutine CalcR
!
subroutine CalcRmunu(ileg,parts,Rmunu)
use process
use particles
use Real_data
use math
use flags
implicit none
!
  integer , intent(in) :: ileg
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
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
    subroutine RmunuSME(iptrn,ileg,p,Rmunu)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn,ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
    end subroutine RmunuSME
  end interface
!
  Rmunu = 0d0
!
  if (flg_scaledep) then
    Rmunu = smeRmunu_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born+1,parts,preal,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(preal)
!
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do jleg=1,size(parts)
! We only care about gluons:
    if (preal(jleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.preal(jleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.preal(jleg)%p%px).and. &
        (parts(ileg)%p%py.eq.preal(jleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.preal(jleg)%p%pz)) then
      exit
    end if
  end do
! We include a check to be sure that we really found the leg:
  if (jleg.eq.size(parts)+1) then
    print *,"Error in CalcRmunu..."
    print *,"ileg,jleg: ",ileg,jleg
    call PrintParts(parts)
    call PrintParts(preal)
    stop
  end if
!
  call RmunuSME(iptrn,jleg,preal,Rmunu)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Rmunu = Rmunu * (4d0*pi)**4
!
! We store Bmunu:
  smeRmunu_saved = Rmunu
!
end subroutine CalcRmunu
!
subroutine CalcRij(parts,Rij)
use process
use particles
use Real_data
use math
use misc
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
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
    subroutine RijSME(iptrn,p,Rij)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
    end subroutine RijSME
  end interface
!
  Rij = 0d0
!
  if (flg_scaledep) then
    Rij = smeRij_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born+1,parts,preal,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(preal)
!
  call RijSME(iptrn,preal,Rij)
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Rij:
  call SwapColor(nleg_born+1,parts,preal,Rij)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Rij = Rij * (4d0*pi)**4
!
  smeRij_saved = Rij
!
end subroutine CalcRij
