! This source contains routines related to the Virtual contribution...
module Virt_data
use particles
implicit none
!
  logical :: ini_V = .true.
  integer :: numptrns
  integer , dimension(:,:) , allocatable :: ptrns
  real(kind(1d0)) , dimension(:) , allocatable :: prefacts
  type(particle) , dimension(:) , allocatable :: pvirt
!
  integer :: iptrn_saved
  integer :: jleg_saved
  real(kind(1d0)) :: V_muindep_saved
  real(kind(1d0)) :: muref
  real(kind(1d0)) , dimension(-4:2) :: VirtLaurentAtMuRef
  real(kind(1d0)) , dimension(-4:2) :: BornLaurentAtMuRef
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: VijLaurentAtMuRef
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: BijLaurentAtMuRef
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: MijLaurent_tmp
  real(kind(1d0)) , dimension(0:3,0:3) :: Vmunu_muindep_saved
  real(kind(1d0)) , dimension(:,:) , allocatable :: Vij_muindep, &
                                                    Vij_mudep
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
  if (.not.ini_V) return
!
! We have to allocate the array which will hold the patterns
! and we have to fill these arrays up too:
! For e+ e- -> q q~ we have only two subprocesses:
  numptrns = 2
!
  allocate(ptrns(nleg_born,numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with patterns allocation in virtual..."
    stop
  end if
!
  allocate(prefacts(numptrns),stat=istat)
  if (istat.ne.0) then
    print *,"Problem with prefacts allocation in virtual..."
    stop
  end if
!
! Filling up...
! Note that the labeling for massless partons is simplified.
! The first massless quark flavor starts with 1 and if the 
! next is different than it is 2 otherwise 1:
!
! e+ e- -> d d~
  ptrns(1,1) = -11 ; ptrns(2,1) = 11
  ptrns(3,1) =  1 ; ptrns(4,1) =  -1
! e+ e- -> u u~
  ptrns(1,2) = -11 ; ptrns(2,2) = 11
  ptrns(3,2) =  2 ; ptrns(4,2) =  -2
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
    print *,"pvirt cannot be allocated..."
    stop
  end if
!
  allocate(VijLaurentAtMuRef(nleg_born,nleg_born,-4:2), &
           BijLaurentAtMuRef(nleg_born,nleg_born,-4:2), &
           MijLaurent_tmp(nleg_born,nleg_born,-4:2) ,   &
           stat=istat)
  if (istat.ne.0) then
    print *,"?ijLaurentAtMuRef cannot be allocated..."
    stop
  end if
!
  allocate(Vij_mudep(nleg_born,nleg_born), &
           Vij_muindep(nleg_born,nleg_born) , &
           stat=istat)
  if (istat.ne.0) then
    print *,"Vij_mu??dep cannot be allocated..."
    stop
  end if
!
  ini_V = .false.
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
!
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
    Virt = V_muindep_saved + Virt_dep
! Normalization factors:
! c_Gamma:
!    Virt = Virt * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    Virt = Virt * (4d0*pi)**2
! For debugging purposes it is possible that the full Laurent series
! is needed for the virtual if this is the case we simply calculate
! everything:
    if (present(VirtLaurent)) then
      call VirtSME(iptrn,pvirt,mur,'f',Virt,Virt_dep,Virt_indep,VirtLaurent)
!      Virt = Virt * cGamma
!      VirtLaurent = VirtLaurent * cGamma
      Virt = Virt * (4d0*pi)**2
      VirtLaurent = VirtLaurent * (4d0*pi)**2
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
      V_muindep_saved = Virt_indep
! We also stash away the pattern used:
      iptrn_saved = iptrn
    end if
  end if
!
! Normalization coming from c_Gamma:
!  Virt = Virt * cGamma
!  if (present(VirtLaurent)) VirtLaurent = VirtLaurent * cGamma
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Virt = Virt * (4d0*pi)**2
  if (present(VirtLaurent)) VirtLaurent = VirtLaurent * (4d0*pi)**2
!
!  print *,"Virt: ",Virt
end subroutine CalcV
!
subroutine CalcVmunu(ileg,parts,Vmunu,VmunuLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
!
implicit none
!
  integer , intent(in) :: ileg
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu
  real(kind(1d0)) , optional , dimension(0:3,0:3,-4:2) , intent(out) :: VmunuLaurent
!
  integer :: iptrn
  integer :: jleg
  real(kind(1d0)) :: prefact
  real(kind(1d0)) , dimension(0:3,0:3) :: Vmunu_dep,Vmunu_indep
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
    subroutine VmunuSME(iptrn,ileg,pvirt,mur,mode,Vmunu,Vmunu_dep,Vmunu_indep, &
                        VmunuLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn,ileg
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      character , intent(in) :: mode
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu
      real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu_dep,Vmunu_indep
      real(kind(1d0)) , optional , dimension(0:3,0:3,-4:2) , intent(out) :: VmunuLaurent
!
    end subroutine VmunuSME
  end interface
!
!
  Vmunu = 0d0
!
! If flg_scaledep is true we are having the same PS point with the
! same flavor configuration hence we only have to reevaluate the mu
! dependent term only:
  if (flg_scaledep) then
! pvirt is declared inside the module hence still usable, the 
! pattern number is also saved:
    iptrn = iptrn_saved
    jleg  = jleg_saved
    call VmunuSME(iptrn,jleg,pvirt,mur,'d',Vmunu,Vmunu_dep,Vmunu_indep)
! The total spin-correlated virtual is built up from the dependent and independent
! terms:
    Vmunu = Vmunu_muindep_saved + Vmunu_dep
! Normalization factors:
! c_Gamma:
    Vmunu = Vmunu * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    Vmunu = Vmunu * (4d0*pi)**2
! If the Laurent expansion is also needed we simply calculate the full 
! contribution along with the Laruent expansion, since it is only used
! in the testing phase not in production code:
    if (present(VmunuLaurent)) then
      call VmunuSME(iptrn,jleg,pvirt,mur,'f',Vmunu,Vmunu_dep,Vmunu_indep, &
                    VmunuLaurent)
! Vmunu has to be normalized again, since it is calculated again:
      Vmunu = Vmunu * cGamma
      VmunuLaurent = VmunuLaurent * cGamma
      Vmunu = Vmunu * (4d0*pi)**2
      VmunuLaurent = VmunuLaurent * (4d0*pi)**2
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvirt,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvirt)
!
! The reshuffling can change the location of ileg, hence we loop over the
! partons and try to find the leg corresponding to the momentum of p(ileg):
  do jleg=1,size(parts)
! We only care about gluons:
    if (pvirt(jleg)%flv.ne.0) cycle
    if ((parts(ileg)%p%E.eq.pvirt(jleg)%p%E).and.   &
        (parts(ileg)%p%px.eq.pvirt(jleg)%p%px).and. &
        (parts(ileg)%p%py.eq.pvirt(jleg)%p%py).and. &
        (parts(ileg)%p%pz.eq.pvirt(jleg)%p%pz)) then
      exit
    end if
  end do
!  print *,"ileg,jleg: ",ileg,jleg
! We include a check to be sure that we really found the leg:
  if (jleg.eq.size(parts)+1) then
    print *,"Error in CalcVmunu..."
    print *,"ileg,jleg: ",ileg,jleg
    call PrintParts(parts)
    call PrintParts(pvirt)
    stop
  end if
  jleg_saved = jleg
!
! If we have a new PS point we reevaluate the whole virtual part,
! that is the mu dependent and independent part as well:
  if (present(VmunuLaurent)) then
    call VmunuSME(iptrn,jleg,pvirt,mur,'f',Vmunu,Vmunu_dep,Vmunu_indep, &
                  VmunuLaurent)
  else
    if (.not.flg_scaledep) then
      call VmunuSME(iptrn,jleg,pvirt,mur,'b',Vmunu,Vmunu_dep,Vmunu_indep)
! We stash away the independent part:
      Vmunu_muindep_saved = Vmunu_indep
! We also stash away the pattern used:
      iptrn_saved = iptrn
    end if
  end if
!
! Normalization coming from c_Gamma:
  Vmunu = Vmunu * cGamma
  if (present(VmunuLaurent)) VmunuLaurent = VmunuLaurent * cGamma
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Vmunu = Vmunu * (4d0*pi)**2
  if (present(VmunuLaurent)) VmunuLaurent = VmunuLaurent * (4d0*pi)**2
!
end subroutine CalcVmunu
!
subroutine CalcVij(parts,Vij,VijLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
!
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
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
    subroutine VijSME(iptrn,pvirt,mur,mode,Vij,Vij_dep,Vij_indep,VijLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      character , intent(in) :: mode
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij_dep,Vij_indep
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine VijSME
  end interface
!
!
  Vij = 0d0
!
! If flg_scaledep is true we are having the same PS point with the
! same flavor configuration hence we only have to reevaluate the mu
! dependent term only:
  if (flg_scaledep) then
! pvirt is declared inside the module hence still usable, the 
! pattern number is also saved:
    iptrn = iptrn_saved
    call VijSME(iptrn,pvirt,mur,'d',Vij,Vij_mudep,Vij_muindep)
! The total spin-correlated virtual is built up from the dependent and independent
! terms:
    Vij = Vij_muindep + Vij_mudep
! Normalization factors:
! c_Gamma:
    Vij = Vij * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    Vij = Vij * (4d0*pi)**2
    if (present(VijLaurent)) then
      call VijSME(iptrn,pvirt,mur,'f',Vij,Vij_mudep,Vij_muindep,VijLaurent)
      Vij = Vij * cGamma
      VijLaurent = VijLaurent * cGamma
      Vij = Vij * (4d0*pi)**2
      VijLaurent = VijLaurent * (4d0*pi)**2
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvirt,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvirt)
!
! If we have a new PS point we reevaluate the whole virtual part,
! that is the mu dependent and independent part as well:
  if (present(VijLaurent)) then
    call VijSME(iptrn,pvirt,mur,'f',Vij,Vij_mudep,Vij_muindep,VijLaurent)
  else
    if (.not.flg_scaledep) then
      call VijSME(iptrn,pvirt,mur,'b',Vij,Vij_mudep,Vij_muindep)
! We stash away the pattern used:
      iptrn_saved = iptrn
    end if
  end if
!
! Normalization coming from c_Gamma:
  Vij = Vij * cGamma
  if (present(VijLaurent)) VijLaurent = VijLaurent * cGamma
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Vij = Vij * (4d0*pi)**2
  if (present(VijLaurent)) VijLaurent = VijLaurent * (4d0*pi)**2
!
end subroutine CalcVij
!
subroutine CalcVddim(parts,Virt,VirtLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
!
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
    subroutine VirtSMEddim(iptrn,pvirt,mur,mode,Virt,Virt_dep,Virt_indep, &
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
    end subroutine VirtSMEddim
  end interface
!
!
  Virt = 1d0
  if (present(VirtLaurent)) then
    VirtLaurent = 1d0
  end if
!
end subroutine CalcVddim
!
subroutine CalcVijddim(parts,Vij,VijLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
!
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
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
    subroutine VijSMEddim(iptrn,pvirt,mur,mode,Vij,Vij_dep,Vij_indep,VijLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      character , intent(in) :: mode
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
      real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij_dep,Vij_indep
      real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
    end subroutine VijSMEddim
  end interface
!
!
  Vij = 0d0
!
! If flg_scaledep is true we are having the same PS point with the
! same flavor configuration hence we only have to reevaluate the mu
! dependent term only:
  if (flg_scaledep) then
! pvirt is declared inside the module hence still usable, the 
! pattern number is also saved:
    iptrn = iptrn_saved
    call VijSMEddim(iptrn,pvirt,mur,'d',Vij,Vij_mudep,Vij_muindep)
! The total spin-correlated virtual is built up from the dependent and independent
! terms:
    Vij = Vij_muindep + Vij_mudep
! Normalization factors:
! c_Gamma:
    Vij = Vij * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    Vij = Vij * (4d0*pi)**2
    if (present(VijLaurent)) then
      call VijSMEddim(iptrn,pvirt,mur,'f',Vij,Vij_mudep,Vij_muindep,VijLaurent)
      Vij = Vij * cGamma
      VijLaurent = VijLaurent * cGamma
      Vij = Vij * (4d0*pi)**2
      VijLaurent = VijLaurent * (4d0*pi)**2
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvirt,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvirt)
!
! If we have a new PS point we reevaluate the whole virtual part,
! that is the mu dependent and independent part as well:
  if (present(VijLaurent)) then
    call VijSMEddim(iptrn,pvirt,mur,'f',Vij,Vij_mudep,Vij_muindep,VijLaurent)
  else
    if (.not.flg_scaledep) then
      call VijSMEddim(iptrn,pvirt,mur,'b',Vij,Vij_mudep,Vij_muindep)
! We stash away the pattern used:
      iptrn_saved = iptrn
    end if
  end if
!
! Normalization coming from c_Gamma:
  Vij = Vij * cGamma
  if (present(VijLaurent)) VijLaurent = VijLaurent * cGamma
!
! The matrix elements are defined such 4\pi\alpha_s = 4\pi\alpha_{EM} = 1:
  Vij = Vij * (4d0*pi)**2
  if (present(VijLaurent)) VijLaurent = VijLaurent * (4d0*pi)**2
!
end subroutine CalcVijddim
!
subroutine CalcVijddim_new(parts,Vij,VijLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
use QCDparams
!
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
  integer :: iptrn
  real(kind(1d0)) :: prefact
  real(kind(1d0)) :: L
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
    subroutine VijSMEddim_new(iptrn,pvirt,mur,VijLaurent,BijLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      real(kind(1d0)) , dimension(:,:,-4:) , intent(out) :: &
        VijLaurent,BijLaurent
!
    end subroutine VijSMEddim_new
  end interface
!
!
  Vij = 0
!
! If flg_scaledep is true we evolve the contributions from muref to mur:
  if (flg_scaledep) then
    MijLaurent_tmp = 0
! xi: mur2 / muref2:
    L = log(mur**2/muref**2)
! Evolving the virtual from muref to mur:
    MijLaurent_tmp(:,:,-2) = 0
    MijLaurent_tmp(:,:,-1) = VijLaurentAtMuRef(:,:,-2)*L
    MijLaurent_tmp(:,:, 0) = 0.5d0*VijLaurentAtMuRef(:,:,-2)*L**2  &
                           + (VijLaurentAtMuRef(:,:,-1)            &
                           + qcd_beta0*BijLaurentAtMuRef(:,:,0))*L
    MijLaurent_tmp(:,:, 1) = 1d0/6d0*VijLaurentAtMuRef(:,:,-2)*L**3   &
                           + 0.5d0*(VijLaurentAtMuRef(:,:,-1)         &
                           + qcd_beta0*BijLaurentAtMuRef(:,:,0))*L**2 &
                           + (VijLaurentAtMuRef(:,:, 0)               &
                           + qcd_beta0*BijLaurentAtMuRef(:,:,1))*L
    MijLaurent_tmp(:,:, 2) = 1d0/24d0*VijLaurentAtMuRef(:,:,-2)*L**4  &
                           + 1d0/6d0*(VijLaurentAtMuRef(:,:,-1)       &
                           + qcd_beta0*BijLaurentAtMuRef(:,:,0))*L**3 &
                           + 0.5d0*(VijLaurentAtMuRef(:,:, 0)         &
                           + qcd_beta0*BijLaurentAtMuRef(:,:,1))*L**2 &
                           + (VijLaurentAtMuRef(:,:, 1)               &
                           + qcd_beta0*BijLaurentAtMuRef(:,:,2))*L
!
    MijLaurent_tmp = MijLaurent_tmp + VijLaurentAtMuRef
!
    Vij = MijLaurent_tmp(:,:,0)
! Normalization factors:
! c_Gamma:
    Vij = Vij * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    Vij = Vij * (4d0*pi)**2
! For debugging purposes it is possible that the full Laurent series
! is needed for the virtual if this is the case we simply calculate
! everything:
    if (present(VijLaurent)) then
      VijLaurent = MijLaurent_tmp * cGamma
      VijLaurent = VijLaurent * (4d0*pi)**2
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvirt,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvirt)
!
! Calling the virtual routine to calculate the full Laurent series of it at
! muref and the Taylor expansion of the Born, too:
  iptrn_saved = iptrn
  muref = mur
  call VijSMEddim_new(iptrn,pvirt,muref,VijLaurentAtMuRef,BijLaurentAtMuRef)
!
  Vij = VijLaurentAtMuRef(:,:,0)
  Vij = Vij * cGamma
  Vij = Vij * (4d0*pi)**2
!
  if (present(VijLaurent)) then
    VijLaurent = VijLaurentAtMuRef
    VijLaurent = VijLaurent * cGamma
    VijLaurent = VijLaurent * (4d0*pi)**2
  end if
!
end subroutine CalcVijddim_new
!
subroutine CalcVddim_new(parts,Virt,VirtLaurent)
use process
use particles
use Virt_data
use math
use flags
use scales
use QCDparams
!
implicit none
!
  type(particle) , dimension(:) , intent(in) :: parts
  real(kind(1d0)) , intent(out) :: Virt
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
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
    subroutine VirtSMEddim_new(iptrn,pvirt,mur,VirtLaurent,BornLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: pvirt
      real(kind(1d0)) , intent(in) :: mur
      real(kind(1d0)) , dimension(-4:2) , intent(out) :: VirtLaurent,BornLaurent
!
    end subroutine VirtSMEddim_new
  end interface
!
!
  Virt = 0
!
! If flg_scaledep is true we evolve the contributions from muref to mur:
  if (flg_scaledep) then
    Laurent_tmp = 0
! xi: mur2 / muref2:
    L = log(mur**2/muref**2)
! Evolving the virtual from muref to mur:
    Laurent_tmp(-2) = 0
    Laurent_tmp(-1) = VirtLaurentAtMuRef(-2)*L
    Laurent_tmp( 0) = 0.5d0*VirtLaurentAtMuRef(-2)*L**2  &
                    + (VirtLaurentAtMuRef(-1)            &
                    + qcd_beta0*BornLaurentAtMuRef(0))*L
    Laurent_tmp( 1) = 1d0/6d0*VirtLaurentAtMuRef(-2)*L**3   &
                    + 0.5d0*(VirtLaurentAtMuRef(-1)         &
                    + qcd_beta0*BornLaurentAtMuRef(0))*L**2 &
                    + (VirtLaurentAtMuRef( 0)               &
                    + qcd_beta0*BornLaurentAtMuRef(1))*L
    Laurent_tmp( 2) = 1d0/24d0*VirtLaurentAtMuRef(-2)*L**4  &
                    + 1d0/6d0*(VirtLaurentAtMuRef(-1)       &
                    + qcd_beta0*BornLaurentAtMuRef(0))*L**3 &
                    + 0.5d0*(VirtLaurentAtMuRef( 0)         &
                    + qcd_beta0*BornLaurentAtMuRef(1))*L**2 &
                    + (VirtLaurentAtMuRef( 1)               &
                    + qcd_beta0*BornLaurentAtMuRef(2))*L
!
    Laurent_tmp = Laurent_tmp + VirtLaurentAtMuRef
!
    Virt = Laurent_tmp(0)
! Normalization factors:
! c_Gamma:
    Virt = Virt * cGamma
! We factor out a tower of \alpha_S and \alpha_{EM}:
    Virt = Virt * (4d0*pi)**2
! For debugging purposes it is possible that the full Laurent series
! is needed for the virtual if this is the case we simply calculate
! everything:
    if (present(VirtLaurent)) then
      VirtLaurent = Laurent_tmp * cGamma
      VirtLaurent = VirtLaurent * (4d0*pi)**2
    end if
    return
  end if
!
! We have to find the correct pattern corresponding to parts:
!  call PrintParts(parts)
  call reshufflemomud(nleg_born,parts,pvirt,numptrns,ptrns,iptrn,prefact)
!  call PrintParts(pvirt)
!
! Calling the virtual routine to calculate the full Laurent series of it at
! muref and the Taylor expansion of the Born, too:
  iptrn_saved = iptrn
  muref = mur
  call VirtSMEddim_new(iptrn,pvirt,muref,VirtLaurentAtMuRef,BornLaurentAtMuRef)
!
  Virt = VirtLaurentAtMuRef(0)
  Virt = Virt * cGamma
  Virt = Virt * (4d0*pi)**2
!
  if (present(VirtLaurent)) then
    VirtLaurent = VirtLaurentAtMuRef * cGamma * (4d0*pi)**2
  end if
!
end subroutine CalcVddim_new
