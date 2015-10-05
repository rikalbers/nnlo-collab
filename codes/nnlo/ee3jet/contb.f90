! This source contains routines related to the Born contribution...
module Born_data
use particles
implicit none
!
  logical :: ini_B = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: pborn
!
  real(kind(1d0)) :: smeB_saved
  real(kind(1d0)) , dimension(-4:2) :: smeBLaurent_saved
  real(kind(1d0)) , dimension(0:3,0:3) :: smeBmunu_saved
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) :: smeBalbemunu_saved
  real(kind(1d0)) , dimension(:,:) , allocatable :: smeBij_saved
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: smeBijLaurent_saved
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: smeBijk_saved
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: smeBijkl_saved
  real(kind(1d0)) , dimension(:,:,:,:,:) , allocatable :: smeBijklLaurent_saved
  real(kind(1d0)) , dimension(:,:,:,:) , allocatable :: smeBmunuij_saved
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
  if (.not.ini_B) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
! For e+ e- -> q q~ g we have only two subprocesses:
  numptrns = 2
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
! e+ e- -> d d~ g
  ptrns(1,1) = -11 ; ptrns(2,1) = 11 ; ptrns(3,1) =  1
  ptrns(4,1) =  -1 ; ptrns(5,1) =  0
! e+ e- -> u u~ g
  ptrns(1,2) = -11 ; ptrns(2,2) = 11 ; ptrns(3,2) =  2
  ptrns(4,2) =  -2 ; ptrns(5,2) =  0
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
  allocate(smeBijLaurent_saved(nleg_born,nleg_born,-4:2),stat=istat)
  if (istat.ne.0) then
    print *,"smeBijLaurent_saved cannot be allocated..."
    stop
  end if
!
! For this process the next color-correlated Born is not needed:
!  allocate(smeBijk_saved(nleg_born,nleg_born,nleg_born),stat=istat)
!  if (istat.ne.0) then
!    print *,"smeBijk_saved cannot be allocated..."
!    stop
!  end if
!
  allocate(smeBijkl_saved(nleg_born,nleg_born,nleg_born,nleg_born), &
           stat=istat)
  if (istat.ne.0) then
    print *,"smeBijkl_saved cannot be allocated..."
    stop
  end if
!
  allocate(smeBijklLaurent_saved(nleg_born, &
                                 nleg_born, &
                                 nleg_born, &
                                 nleg_born,-4:2), &
                                 stat=istat)
  if (istat.ne.0) then
    print *,"smeBijklLaurent_saved cannot be allocated..."
    stop
  end if
!
  allocate(smeBmunuij_saved(0:3,0:3,nleg_born,nleg_born), &
           stat=istat)
  if (istat.ne.0) then
    print *,"smeBijkl_saved cannot be allocated..."
    stop
  end if
!
  ini_B = .false.
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
  Born = Born * (4d0*pi)**3
!
! We also store the Born which can be useful for scale dependency studies:
  smeB_saved = Born
!
!  print *,"Born: ",Born
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
! Debug!
  integer :: mu,nu
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu_tmp
!
!
  interface
    subroutine reshufflemomud(n,p,pout,numptrn,ptrn,returnptrn, &
                              prefactor)
    use particles
    implicit none
!
      integer , intent(in) :: n,numptrn
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: pout
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
!
    subroutine BmunuCS(p,Bmunu)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine BmunuCS
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
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do jleg=1,size(parts)
! We only care about gluons:
    if (pborn(jleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pborn(jleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pborn(jleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pborn(jleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pborn(jleg)%p%pz)) then
      exit
    end if
  end do
! We include a check to be sure that we really found the leg:
  if (jleg.eq.size(parts)+1) then
    print *,"Error in CalcBmunu..."
    print *,"ileg,jleg: ",ileg,jleg
    call PrintParts(parts)
    call PrintParts(pborn)
    stop
  end if
!
  call BmunuSME(iptrn,jleg,pborn,Bmunu)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bmunu = Bmunu * (4d0*pi)**3
!
! We store Bmunu:
  smeBmunu_saved = Bmunu
!
end subroutine CalcBmunu
!
subroutine CalcBmunu_q1q2(ileg,q1,q2,parts,Bmunu)
use process
use particles
use Born_data
use math
use flags
implicit none
!
  integer , intent(in) :: ileg
  type(mom) , intent(in) :: q1,q2
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
    subroutine BmunuSME_q1q2(iptrn,ileg,q1,q2,p,Bmunu)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn,ileg
      type(mom) , intent(in) :: q1,q2
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
    end subroutine BmunuSME_q1q2
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
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do jleg=1,size(parts)
! We only care about gluons:
    if (pborn(jleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pborn(jleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pborn(jleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pborn(jleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pborn(jleg)%p%pz)) then
      exit
    end if
  end do
! We include a check to be sure that we really found the leg:
  if (jleg.eq.size(parts)+1) then
    print *,"Error in CalcBmunu_q1q2..."
    print *,"ileg,jleg: ",ileg,jleg
    call PrintParts(parts)
    call PrintParts(pborn)
    stop
  end if
!
  call BmunuSME_q1q2(iptrn,jleg,q1,q2,pborn,Bmunu)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bmunu = Bmunu * (4d0*pi)**3
!
! We store Bmunu:
  smeBmunu_saved = Bmunu
!
end subroutine CalcBmunu_q1q2
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
  integer :: iptrn
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
  Bij = Bij * (4d0*pi)**3
!
  smeBij_saved = Bij
!
end subroutine CalcBij
!
subroutine CalcBijk(parts,Bijk)
use process
use particles
use Born_data
use math
use misc
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:,:) , intent(out) :: Bijk
!
  integer :: iptrn
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
!    subroutine BijkSME(iptrn,p,Bijk)
!    use particles
!    implicit none
!!
!      integer , intent(in) :: iptrn
!      type(particle) , dimension(:) , intent(in) :: p
!      real(kind(1d0)) , dimension(:,:,:) , intent(out) :: Bijk
!!
!    end subroutine BijkSME
  end interface
!
  Bijk = 0d0
!
  print *,"We entered CalcBijk..."
  print *,"For this process this cannot happen..."
  stop
!
  if (flg_scaledep) then
    Bijk = smeBijk_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
!  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
!  call BijkSME(iptrn,pborn,Bijk)
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Bijk:
! TODO extend the functionality of SwapColor to be usable even with Bijk:
!  call SwapColor(nleg_born,parts,pborn,Bijk)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bijk = Bijk * (4d0*pi)**3
!
  smeBijk_saved = Bijk
!
end subroutine CalcBijk
!
subroutine CalcBijkl(parts,Bijkl)
use process
use particles
use Born_data
use math
use misc
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
  integer :: iptrn
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
    subroutine BijklSME(iptrn,p,Bijkl)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
    end subroutine BijklSME
  end interface
!
  Bijkl = 0d0
!
  if (flg_scaledep) then
    Bijkl = smeBijkl_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  call BijklSME(iptrn,pborn,Bijkl)
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Bijk:
  call SwapColor(nleg_born,parts,pborn,Bijkl)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bijkl = Bijkl * (4d0*pi)**3
!
  smeBijkl_saved = Bijkl
!
end subroutine CalcBijkl
!
subroutine CalcBmunuij(ileg,parts,Bmunuij)
use process
use particles
use Born_data
use math
use misc
use flags
use QCDparams
implicit none
!
  integer , intent(in) :: ileg
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
  integer :: iptrn
  integer :: jleg
  real(kind(1d0)) :: prefact
!
  integer :: mu,nu,i,j
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu,Bmunu_tmp
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
    subroutine BmunuijSME(iptrn,ileg,p,Bmunuij)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      integer , intent(in) :: ileg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
    end subroutine BmunuijSME
  end interface
!
  if (flg_scaledep) then
    Bmunuij = smeBmunuij_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do jleg=1,size(parts)
! We only care about gluons:
    if (pborn(jleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pborn(jleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pborn(jleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pborn(jleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pborn(jleg)%p%pz)) then
      exit
    end if
  end do
! We include a check to be sure that we really found the leg:
  if (jleg.eq.size(parts)+1) then
    print *,"Error in CalcBmunuij..."
    print *,"ileg,jleg: ",ileg,jleg
    call PrintParts(parts)
    call PrintParts(pborn)
    stop
  end if
!
  call BmunuijSME(iptrn,jleg,pborn,Bmunuij)
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Bijk:
  call SwapColormunuij_part_part(nleg_born,parts,pborn,Bmunuij)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bmunuij = Bmunuij * (4d0*pi)**3
!
  smeBmunuij_saved = Bmunuij
!
end subroutine CalcBmunuij
!
subroutine CalcBalbemunu(ileg,jleg,parts,Balbemunu)
use process
use particles
use Born_data
use math
use flags
implicit none
!
  integer , intent(in) :: ileg,jleg
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(out) :: Balbemunu
!
  integer :: iptrn,kleg,lleg
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
    subroutine BalbemunuSME(iptrn,ileg,jleg,p,Balbemunu)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn,ileg,jleg
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(0:3,0:3,0:3,0:3) , intent(out) :: Balbemunu
!
    end subroutine BalbemunuSME
  end interface
!
  Balbemunu = 0d0
!
  print *,"This routine is not implemented yet..."
  stop "CalcBalbemunu"
!
  if (flg_scaledep) then
    Balbemunu = smeBalbemunu_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do kleg=1,size(parts)
! We only care about gluons:
    if (pborn(kleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pborn(kleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pborn(kleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pborn(kleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pborn(kleg)%p%pz)) then
      exit
    end if
  end do
  do lleg=1,size(parts)
! We only care about gluons:
    if (pborn(lleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pborn(lleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pborn(lleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pborn(lleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pborn(lleg)%p%pz)) then
      exit
    end if
  end do
! We include a check to be sure that we really found the leg:
  if (kleg.eq.size(parts)+1) then
    print *,"Error in CalcBalbemunu..."
    print *,"ileg,kleg: ",ileg,kleg
    call PrintParts(parts)
    call PrintParts(pborn)
    stop
  end if
  if (lleg.eq.size(parts)+1) then
    print *,"Error in CalcBalbemunu..."
    print *,"ileg,lleg: ",ileg,lleg
    call PrintParts(parts)
    call PrintParts(pborn)
    stop
  end if
!
!  call BalbemunuSME(iptrn,kleg,lleg,pborn,Balbemunu)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Balbemunu = Balbemunu * (4d0*pi)**3
!
! We store Bmunu:
  smeBalbemunu_saved = Balbemunu
!
end subroutine CalcBalbemunu
!
subroutine CalcBmunu_hel(ileg,parts,Balbe)
use process
use particles
use Born_data
use math
use flags
implicit none
!
  integer , intent(in) :: ileg
  type(particle) , dimension(:) , intent(in) :: parts
  complex(kind(1d0)) , dimension(2,2) , intent(out) :: Balbe
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
      type(particle) , dimension(:) , intent(in) :: p
      type(particle) , dimension(:) , intent(out) :: pout
      integer , intent(in) , dimension(:,:) :: ptrn
      integer , intent(out) :: returnptrn
      real(kind(1d0)) , intent(out) :: prefactor
!
    end subroutine reshufflemomud
!
    subroutine BalbeSME(iptrn,ileg,p,Balbe)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn,ileg
      type(particle) , dimension(:) , intent(in) :: p
      complex(kind(1d0)) , dimension(2,2) , intent(out) :: Balbe
!
    end subroutine BalbeSME
  end interface
!
  Balbe = 0d0
!
!  if (flg_scaledep) then
!    Balbe = smeBalbe_saved
!    return
!  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do jleg=1,size(parts)
! We only care about gluons:
    if (pborn(jleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pborn(jleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pborn(jleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pborn(jleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pborn(jleg)%p%pz)) then
      exit
    end if
  end do
! We include a check to be sure that we really found the leg:
  if (jleg.eq.size(parts)+1) then
    print *,"Error in CalcBalbe..."
    print *,"ileg,jleg: ",ileg,jleg
    call PrintParts(parts)
    call PrintParts(pborn)
    stop
  end if
!
  call BalbeSME(iptrn,jleg,pborn,Balbe)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Balbe = Balbe * (4d0*pi)**3
!
! We store Balbe:
!  smeBalbe_saved = Balbe
!
end subroutine CalcBmunu_hel
!
! This is a carbon-copy of the CalcB routine with the difference
! that this one calls a d-dimensional SME:
subroutine CalcBddim(parts,Born,BornLaurent)
use process
use particles
use Born_data
use math
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: Born
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
    BornLaurent
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
    subroutine BornSMEddim(iptrn,pborn,Born,BornLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pborn
      real(kind(1d0)) , intent(out) :: Born
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        BornLaurent
!
    end subroutine BornSMEddim
  end interface
!
!
  Born = 0
  if (present(BornLaurent)) BornLaurent = 0
!
! If flg_scaledep is true we do not have to reevaluate the Born provided
! by sitting at the same PS point and there is no \mu dependent terms
! inside:
  if (flg_scaledep) then
    Born = smeB_saved
    if (present(BornLaurent)) BornLaurent = smeBLaurent_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  if (.not.present(BornLaurent)) then
    call BornSMEddim(iptrn,pborn,Born)
  else
    call BornSMEddim(iptrn,pborn,Born,BornLaurent)
    BornLaurent = BornLaurent * (4d0*pi)**3
    smeBLaurent_saved = BornLaurent
  end if
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Born = Born * (4d0*pi)**3
!
! We also store the Born which can be useful for scale dependency studies:
  smeB_saved = Born
!
!  print *,"Born: ",Born
!
end subroutine CalcBddim
!
subroutine CalcBijddim(parts,Bij,BijLaurent)
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
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: &
    BijLaurent
!
  integer :: iptrn
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
    subroutine BijSMEddim(iptrn,p,Bij,BijLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: &
        BijLaurent
!
    end subroutine BijSMEddim
  end interface
!
  Bij = 0
  if (present(BijLaurent)) BijLaurent = 0
!
  if (flg_scaledep) then
    Bij = smeBij_saved
    if (present(BijLaurent)) BijLaurent = smeBijLaurent_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  if (.not.present(BijLaurent)) then
    call BijSMEddim(iptrn,pborn,Bij)
  else
    call BijSMEddim(iptrn,pborn,Bij,BijLaurent)
    BijLaurent = BijLaurent * (4d0*pi)**3
    call SwapColorLaurent(nleg_born,parts,pborn,BijLaurent)
    smeBijLaurent_saved = BijLaurent
  end if
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Bij:
  call SwapColor(nleg_born,parts,pborn,Bij)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bij = Bij * (4d0*pi)**3
!
  smeBij_saved = Bij
!
end subroutine CalcBijddim
!
subroutine CalcBijklddim(parts,Bijkl,BijklLaurent)
use process
use particles
use Born_data
use math
use misc
use flags
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
  real(kind(1d0)) , optional , dimension(:,:,:,:,-4:) , intent(out) :: &
    BijklLaurent
!
  integer :: iptrn
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
    subroutine BijklSMEddim(iptrn,p,Bijkl,BijklLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
      real(kind(1d0)) , optional , dimension(:,:,:,:,-4:) , intent(out) :: &
        BijklLaurent
!
    end subroutine BijklSMEddim
  end interface
!
  Bijkl = 0d0
  if (present(BijklLaurent)) BijklLaurent = 0
!
  if (flg_scaledep) then
    Bijkl = smeBijkl_saved
    if (present(BijklLaurent)) BijklLaurent = smeBijklLaurent_saved
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pborn,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pborn)
!
  if (.not.present(BijklLaurent)) then
    call BijklSMEddim(iptrn,pborn,Bijkl)
  else
    call BijklSMEddim(iptrn,pborn,Bijkl,BijklLaurent)
    BijklLaurent = BijklLaurent * (4d0*pi)**3
    call SwapColorLaurent(nleg_born,parts,pborn,BijklLaurent)
    smeBijklLaurent_saved = BijklLaurent
  end if
!
! If a few of the momenta were reshuffled we have to change the position of
! entries even in Bijk:
  call SwapColor(nleg_born,parts,pborn,Bijkl)
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Bijkl = Bijkl * (4d0*pi)**3
!
  smeBijkl_saved = Bijkl
!
end subroutine CalcBijklddim
