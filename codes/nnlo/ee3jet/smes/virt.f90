! mode is a character variable:
! mode = 'f' the full virtual part is needed,
! mode = 'i' the mu independent part is needed,
! mode = 'd' the mu dependent part is needed,
! mode = 'b' both the mu dependent and indepent part is calculated.
subroutine VirtSME(iptrn,pin,mur,mode,Virt,Virt_dep,Virt_indep, &
                   VirtLaurent)
use process
use particles
use QCDparams
use my_model
use math
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(in) :: mur
  character , intent(in) :: mode
  real(kind(1d0)) , intent(out) :: Virt
  real(kind(1d0)) , intent(out) :: Virt_dep,Virt_indep
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  complex(kind(1d0)) , dimension(2,2,2) :: C
  complex(kind(1d0)) , dimension(2) :: Cax
  common/COUPLINGS/ C &
        /AXIALCOUPLINGS/ Cax &
        /QCD/XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  real(kind(1d0)) :: COMU
  real(kind(1d0)) :: Q2
  common/dotproducts/S &
        /spinorproducts/A,B &
        /scales/ COMU &
        /qsquared/ Q2
!
  real(kind(1d0)) , external :: PSI2d1gBorn,PSI2u1gBorn, &
                                PSI2d1gVirtNLO,PSI2u1gVirtNLO, &
                                PSI2d1gVirtNLOem2,PSI2d1gVirtNLOem1, &
                                PSI2u1gVirtNLOem2,PSI2u1gVirtNLOem1, &
                                PSI2d1gVirtNLOmudep,PSI2d1gVirtNLOmuindep, &
                                PSI2u1gVirtNLOmudep,PSI2u1gVirtNLOmuindep
!
! We initialize the QCD and COUPLINGS blocks for each and every
! contribution, safety first...
  if (init) then
    init = .false.
!
    XNf = qcd_nf
    XNu = qcd_nu
    XNd = qcd_nd
    XNc = qcd_nc
    XNa = qcd_nc**2 - 1d0
!
    C(:,:,:) = couplings
    Cax(:)   = axialcouplings
!
  end if
!
! We have to setup the dot and spinor products:
! To calculate these products the momenta have to be given in
! an ordinary array:
  P = 0d0
  do ipart=1,nleg_born
! Note that we are working in an all-outgoing scheme, hence the
! incoming momenta should acquire an additional minus sign:
! Note, too, that the position for the incoming particles are at the
! 7th and 8th positions:
    if (ipart.le.2) then
      P(:,ipart+6) = -pin(ipart)%p
! Final state momenta should be accordingly shifted to the left:
    else
      P(:,ipart-2) =  pin(ipart)%p
    end if
  end do
!
!  do ipart=1,9
!    print *,P(:,ipart)
!  end do
!
  call getdotproducts(P,9,S)
  call getspinorproducts(P,9,A,B)
!
! Renormalization scale:
  comu = mur
  Q2   = 2d0*pin(1)%p*pin(2)%p
!
  Virt       = 0d0
  Virt_indep = 0d0
  Virt_dep   = 0d0
  if (present(VirtLaurent)) VirtLaurent = 0d0
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
!    print *,"e+ e- -> d d~ g"
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born = PSI2d1gBorn(1,3,2,8,7)
    end if
! The full virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
      Virt = PSI2d1gVirtNLO(1,3,2,8,7)
      Virt = Virt - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      Virt_indep = PSI2d1gVirtNLOmuindep(1,3,2,8,7)
! We sweep the renormalization into the independent part:
      Virt_indep = Virt_indep - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      Virt_dep = PSI2d1gVirtNLOmudep(1,3,2,8,7)
    end if
! We construct the total virtual part out of the two contributions:
    if (mode.eq.'b') Virt = Virt_dep + Virt_indep
! If the VirtLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(VirtLaurent)) then
      VirtLaurent = 0d0
      VirtLaurent(-2) = PSI2d1gVirtNLOem2(1,3,2,8,7)
      VirtLaurent(-1) = PSI2d1gVirtNLOem1(1,3,2,8,7) &
                      - qcd_beta0*Born
    end if
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g"
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born = PSI2u1gBorn(1,3,2,8,7)
    end if
! The full virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
      Virt = PSI2u1gVirtNLO(1,3,2,8,7)
      Virt = Virt - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      Virt_indep = PSI2u1gVirtNLOmuindep(1,3,2,8,7)
! We sweep the renormalization into the independent part:
      Virt_indep = Virt_indep - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      Virt_dep = PSI2u1gVirtNLOmudep(1,3,2,8,7)
    end if
! We construct the total virtual part out of the two contributions:
    if (mode.eq.'b') Virt = Virt_dep + Virt_indep
! If the VirtLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(VirtLaurent)) then
      VirtLaurent = 0d0
      VirtLaurent(-2) = PSI2u1gVirtNLOem2(1,3,2,8,7)
      VirtLaurent(-1) = PSI2u1gVirtNLOem1(1,3,2,8,7) &
                      - qcd_beta0*Born
    end if
!    do ipart=1,5
!      print *,pin(ipart)%p
!    end do
!    print *,"virt: ",virt
!    print *,"virtem2: ",PSI2u1gVirtNLOem2(1,3,2,8,7)
!    print *,"virtem1: ",PSI2u1gVirtNLOem1(1,3,2,8,7) &
!                      - qcd_beta0*Born
!    print *,"born: ",born
!    print *,"Q: ",sqrt(Q2)
!    print *,"mur: ",mur
!    stop "VirtSME..."
    return
  end if
!
end subroutine VirtSME
!
! Note that the quark sitting in position 3 is used as the reference
! vector for constructing the gluon polarizations:
subroutine VmunuSME(iptrn,ileg,pin,mur,mode,Vmunu,Vmunu_dep,Vmunu_indep, &
                    VmunuLaurent)
use process
use particles
use QCDparams
use my_model
use math
use misc
implicit none
!
  integer , intent(in) :: iptrn,ileg
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(in) :: mur
  character , intent(in) :: mode
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Vmunu_dep,Vmunu_indep
  real(kind(1d0)) , optional , dimension(0:3,0:3,-4:2) , intent(out) :: VmunuLaurent
!
  integer :: ipart,qtype
  logical , save :: init = .true.
  real(kind(1d0)) :: Born,Virt,Virt_dep,Virt_indep,Virtem2,Virtem1
  complex(kind(1d0)) , dimension(2,2) :: m2ijV,m2ijB
  complex(kind(1d0)) , dimension(2,2) :: m2ijV_dep,m2ijV_indep
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  complex(kind(1d0)) , dimension(2,2,2) :: C
  complex(kind(1d0)) , dimension(2) :: Cax
  common/COUPLINGS/ C &
        /AXIALCOUPLINGS/ Cax &
        /QCD/XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  real(kind(1d0)) :: COMU
  real(kind(1d0)) :: Q2
  common/dotproducts/S &
        /spinorproducts/A,B &
        /scales/ COMU &
        /qsquared/ Q2
!
  real(kind(1d0)) , external :: PSI2q1gVmunu, &
                                PSI2d1gBornV,PSI2u1gBornV, &
                                PSI2q1gVmunu_mudep, &
                                PSI2q1gVmunu_muindep, &
                                PSI2q1gVmunuem1, &
                                PSI2q1gVmunuem2
!
! We initialize the QCD and COUPLINGS blocks for each and every
! contribution, safety first...
  if (init) then
    init = .false.
!
    XNf = qcd_nf
    XNu = qcd_nu
    XNd = qcd_nd
    XNc = qcd_nc
    XNa = qcd_nc**2 - 1d0
!
    C(:,:,:) = couplings
    Cax(:)   = axialcouplings
!
  end if
!
! We have to setup the dot and spinor products:
! To calculate these products the momenta have to be given in
! an ordinary array:
  P = 0d0
  do ipart=1,nleg_born
! Note that we are working in an all-outgoing scheme, hence the
! incoming momenta should acquire an additional minus sign:
! Note, too, that the position for the incoming particles are at the
! 7th and 8th positions:
    if (ipart.le.2) then
      P(:,ipart+6) = -pin(ipart)%p
! Final state momenta should be accordingly shifted to the left:
    else
      P(:,ipart-2) =  pin(ipart)%p
    end if
  end do
!
!  do ipart=1,9
!    print *,P(:,ipart)
!  end do
!
  call getdotproducts(P,9,S)
  call getspinorproducts(P,9,A,B)
!
! Renormalization scale:
  comu = mur
  Q2   = 2d0*pin(1)%p*pin(2)%p
!
  Vmunu       = 0d0
  Vmunu_indep = 0d0
  Vmunu_dep   = 0d0
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
!    print *,"e+ e- -> d d~ g"
    qtype = 2
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born = PSI2d1gBornV(1,3,2,8,7,m2ijB)
    end if
! The full virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
      Virt = PSI2q1gVmunu(1,3,2,8,7,qtype,m2ijV)
      Virt = Virt - 2d0*qcd_cf*Born
      m2ijV = m2ijV - 2d0*qcd_cf*m2ijB
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      Virt_indep = PSI2q1gVmunu_muindep(1,3,2,8,7,qtype,m2ijV_indep)
! We sweep the renormalization into the independent part:
      Virt_indep = Virt_indep - 2d0*qcd_cf*Born
      m2ijV_indep = m2ijV_indep - 2d0*qcd_cf*m2ijB
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      Virt_dep = PSI2q1gVmunu_mudep(1,3,2,8,7,qtype,m2ijV_dep)
    end if
! We construct the total virtual part out of the two contributions:
    if (mode.eq.'b') then
      m2ijV = m2ijV_dep + m2ijV_indep
      Virt = Virt_dep + Virt_indep
    end if
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    if (mode.eq.'f') then
! When the full contribution is asked for only m2ijV is transformed
! to Lorentz basis:
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,Vmunu)
    elseif (mode.eq.'b') then
! When both contributions are requested for the sum and the
! separate contributions are transformed to Lorentz basis:
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,Vmunu)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_dep,Vmunu_dep)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_indep,Vmunu_indep)
! Only the independent part is treated:
    elseif (mode.eq.'i') then
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_indep,Vmunu_indep)
! Only the dependent part is treated:
    elseif (mode.eq.'d') then
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_dep,Vmunu_dep)
    end if
! If the poles are also needed additional calculations have to 
! be done:
    if (present(VmunuLaurent)) then
      VmunuLaurent = 0d0
! Calculation of the O(e^-2) term:
      Virtem2 = PSI2q1gVmunuem2(1,3,2,8,7,qtype,m2ijV)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,VmunuLaurent(:,:,-2))
! Calculation of the O(e^-1) term:
      Virtem1 = PSI2q1gVmunuem1(1,3,2,8,7,qtype,m2ijV)
! Renormalization:
      m2ijV   = m2ijV - qcd_beta0*m2ijB
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,VmunuLaurent(:,:,-1))
    end if
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g"
    qtype = 1
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born = PSI2u1gBornV(1,3,2,8,7,m2ijB)
    end if
! The full virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
      Virt = PSI2q1gVmunu(1,3,2,8,7,qtype,m2ijV)
      Virt = Virt - 2d0*qcd_cf*Born
      m2ijV = m2ijV - 2d0*qcd_cf*m2ijB
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      Virt_indep = PSI2q1gVmunu_muindep(1,3,2,8,7,qtype,m2ijV_indep)
! We sweep the renormalization into the independent part:
      Virt_indep = Virt_indep - 2d0*qcd_cf*Born
      m2ijV_indep = m2ijV_indep - 2d0*qcd_cf*m2ijB
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      Virt_dep = PSI2q1gVmunu_mudep(1,3,2,8,7,qtype,m2ijV_dep)
    end if
! We construct the total virtual part out of the two contributions:
    if (mode.eq.'b') then
      m2ijV = m2ijV_dep + m2ijV_indep
      Virt = Virt_dep + Virt_indep
    end if
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    if (mode.eq.'f') then
! When the full contribution is asked for only m2ijV is transformed
! to Lorentz basis:
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,Vmunu)
    elseif (mode.eq.'b') then
! When both contributions are requested for the sum and the
! separate contributions are transformed to Lorentz basis:
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,Vmunu)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_dep,Vmunu_dep)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_indep,Vmunu_indep)
! Only the independent part is treated:
    elseif (mode.eq.'i') then
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_indep,Vmunu_indep)
! Only the dependent part is treated:
    elseif (mode.eq.'d') then
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV_dep,Vmunu_dep)
    end if
! If the poles are also needed additional calculations have to 
! be done:
    if (present(VmunuLaurent)) then
      VmunuLaurent = 0d0
! Calculation of the O(e^-2) term:
      Virtem2 = PSI2q1gVmunuem2(1,3,2,8,7,qtype,m2ijV)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,VmunuLaurent(:,:,-2))
! Calculation of the O(e^-1) term:
      Virtem1 = PSI2q1gVmunuem1(1,3,2,8,7,qtype,m2ijV)
! Renormalization:
      m2ijV   = m2ijV - qcd_beta0*m2ijB
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ijV,VmunuLaurent(:,:,-1))
    end if
    return
  end if
!
end subroutine VmunuSME
!
subroutine VijSME(iptrn,pin,mur,mode,Vij,Vij_dep,Vij_indep,VijLaurent)
use particles
use QCDparams
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(in) :: mur
  character , intent(in) :: mode
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij_dep,Vij_indep
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
  logical , save :: init = .true.
  real(kind(1d0)) :: Virt,Virt_dep,Virt_indep
  real(kind(1d0)) , dimension(-4:2) :: VirtLaurent
  real(kind(1d0)) , dimension(5,5) , save :: Vij_tmp
!
!
  interface
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
  if (init) then
! A matrix template is constructed for the color-correlated
! virtual:
    Vij_tmp = 0d0
!
    Vij_tmp(3,3) = qcd_cf
    Vij_tmp(4,4) = qcd_cf
    Vij_tmp(5,5) = qcd_ca
    Vij_tmp(3,4) = 0.5d0*(Vij_tmp(5,5) - Vij_tmp(3,3) - Vij_tmp(4,4))
    Vij_tmp(3,5) = 0.5d0*(Vij_tmp(4,4) - Vij_tmp(3,3) - Vij_tmp(5,5))
    Vij_tmp(4,5) = 0.5d0*(Vij_tmp(3,3) - Vij_tmp(4,4) - Vij_tmp(5,5))
    Vij_tmp(4,3) = Vij_tmp(3,4)
    Vij_tmp(5,3) = Vij_tmp(3,5)
    Vij_tmp(5,4) = Vij_tmp(4,5)
!
    init = .false.
  end if
!
  Vij       = 0d0
  Vij_indep = 0d0
  Vij_dep   = 0d0
!
! The color-correlated virtual is trivial for 3 parton
! processes, we obtain the virtual then construct the 
! matrix.
! An extra conditional environment is introduced to tackle the
! case when the Laurent series is asked for:
  if (.not.present(VijLaurent)) then
    call VirtSME(iptrn,pin,mur,mode,Virt,Virt_dep,Virt_indep)
  else
    call VirtSME(iptrn,pin,mur,mode,Virt,Virt_dep,Virt_indep,VirtLaurent)
  end if
!
! The mode flag determines which contribution is to be calculated:
! color-correlated virtual is constructed according:
! Only the full contribution is calculated:
  if (mode.eq.'f') then
    Vij = Vij_tmp * Virt
! Only the \mu_R independent part is calculated:
  elseif (mode.eq.'i') then
    Vij_indep = Vij_tmp * Virt_indep
! Only the \mu_R dependent part is calculated:
  elseif (mode.eq.'d') then
    Vij_dep = Vij_tmp * Virt_dep
! Both the \mu_R dependent and independent part are calculated:
  elseif (mode.eq.'b') then
    Vij       = Vij_tmp * Virt
    Vij_dep   = Vij_tmp * Virt_dep
    Vij_indep = Vij_tmp * Virt_indep
  else
    print *,"Wrong mode selector is given to VijSME..."
    print *,"mode: ",mode
    stop
  end if
!
  if (present(VijLaurent)) then
    VijLaurent = 0d0
    VijLaurent(:,:,-2) = Vij_tmp*VirtLaurent(-2)
    VijLaurent(:,:,-1) = Vij_tmp*VirtLaurent(-1)
  end if
!
end subroutine VijSME
!
subroutine VirtSMEddim(iptrn,p,mur,mode,Virt,Virt_dep,Virt_indep, &
                       VirtLaurent)
use process
use particles
use QCDparams
use my_model
!use math
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: mur
  character , intent(in) :: mode
  real(kind(1d0)) , intent(out) :: Virt
  real(kind(1d0)) , intent(out) :: Virt_dep,Virt_indep
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: VirtLaurent
!
  integer :: i
  real(kind(1d0)) :: Q2,y12,y13,y23
  real(kind(1d0)) :: Born,xi
  real(kind(1d0)) ::  NC,nf,TR,nfg,pi,z3,eg
  common/piez/ pi,z3
  common/const/ NC,nf,TR,eg
  real(kind(1d0)) , dimension(-4:2) :: BLaurent,VLaurent,VLaurent1,Vlaurent2
  complex(kind(1d0)) , dimension(-4:2) :: f1yz,f1zy,f2yz,f2zy
  type(mom) :: Q
!
!
  interface 
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
  pi=3.14159265358979324d0
  z3=1.20205690315959429d0 !zeta(3)
  eg=0.57721566490153286d0 !eulergamma
!
  Virt       = 0
  Virt_indep = 0
  Virt_dep   = 0
  if (present(VirtLaurent)) VirtLaurent = 0
!
! Calculation of invariants:
!
  Q = p(1)%p + p(2)%p
  Q2 = Q**2
!
  y12 = 2*p(3)%p*p(4)%p/Q2
  y13 = 2*p(3)%p*p(5)%p/Q2
  y23 = 2*p(4)%p*p(5)%p/Q2
!
! To construct a Virtual with a d-dimensional Born we need
! the complete Laurent series of both:
  call BornSMEddim(iptrn,p,Born,BLaurent)
!
  VLaurent = 0
!
! Calculation of mu independent part:
  if ((mode.eq.'f').or.(mode.eq.'b').or.(mode.eq.'i')) then
    call f1(f1yz,y13,y23)
    call f1(f1zy,y23,y13)
    call f2(f2yz,y13,y23)
    call f2(f2zy,y23,y13)
    do i=-4,2
      VLaurent(i) = 8*qcd_nc*qcd_cf*real(qcd_nc*(f1yz(i)+f1zy(i)) &
                  + 1d0/qcd_nc*(f2yz(i)+f2zy(i)))
    end do
! Pattern 1 is for d d~ g and pattern 2 is for u u~ g:
    if (iptrn.eq.2) VLaurent = 4*VLaurent
    VLaurent = VLaurent/Q2/27d0
! Renormalization:
! Note the extra factor of 1/2 compensation for a factor of 2
! difference between the beta0 definition used in the 
! framework and in the SME:
    do i=-4,1
      VLaurent(i) = VLaurent(i) - qcd_beta0*BLaurent(i+1)
    end do
    Virt_indep = VLaurent(0)
    if (present(VirtLaurent)) VirtLaurent = VLaurent
  end if
!
! Calculation of the mu dependent part:
  if ((mode.eq.'f').or.(mode.eq.'b').or.(mode.eq.'d')) then
    Virt_dep = 0
!
! xi = mur**2/Q**2:
    xi = mur**2/Q2
!
! This corresponds to the original scale dependence of sqmbeta.f
    VLaurent1 = 0
!
    VLaurent1(-2) = 0
    VLaurent1(-1) = VLaurent(-2)*log(xi)
    VLaurent1( 0) = 0.5d0*VLaurent(-2)*log(xi)**2 &
                  + (VLaurent(-1) + qcd_beta0*BLaurent(0))*log(xi)
    VLaurent1( 1) = 1d0/6d0*VLaurent(-2)*log(xi)**3 &
                  + 0.5d0*(VLaurent(-1) + qcd_beta0*BLaurent(0))*log(xi)**2 &
                  + (VLaurent( 0) + qcd_beta0*BLaurent(1))*log(xi)
    VLaurent1( 2) = 1d0/24d0*VLaurent(-2)*log(xi)**4 &
                  + 1d0/6d0*(VLaurent(-1) + qcd_beta0*BLaurent(0))*log(xi)**3 &
                  + 0.5d0*(VLaurent( 0) + qcd_beta0*BLaurent(1))*log(xi)**2 &
                  + (VLaurent( 1) + qcd_beta0*BLaurent(2))*log(xi)
!
    VLaurent = VLaurent + VLaurent1
!
    Virt_dep = VLaurent1(0)
!
    if (present(VirtLaurent)) VirtLaurent = VirtLaurent + VLaurent1
!
  end if
!
  Virt = Virt_indep + Virt_dep
!
end subroutine VirtSMEddim
!
subroutine VijSMEddim(iptrn,p,mur,mode,Vij,Vij_dep,Vij_indep,VijLaurent)
use particles
use QCDparams
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: mur
  character , intent(in) :: mode
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Vij_dep,Vij_indep
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: VijLaurent
!
  logical , save :: init = .true.
  real(kind(1d0)) :: Virt,Virt_dep,Virt_indep
  real(kind(1d0)) , dimension(-4:2) :: VirtLaurent
  real(kind(1d0)) , dimension(5,5) , save :: Vij_tmp
!
!
  interface
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
  if (init) then
! A matrix template is constructed for the color-correlated
! virtual:
    Vij_tmp = 0d0
!
    Vij_tmp(3,3) = qcd_cf
    Vij_tmp(4,4) = qcd_cf
    Vij_tmp(5,5) = qcd_ca
    Vij_tmp(3,4) = 0.5d0*(Vij_tmp(5,5) - Vij_tmp(3,3) - Vij_tmp(4,4))
    Vij_tmp(3,5) = 0.5d0*(Vij_tmp(4,4) - Vij_tmp(3,3) - Vij_tmp(5,5))
    Vij_tmp(4,5) = 0.5d0*(Vij_tmp(3,3) - Vij_tmp(4,4) - Vij_tmp(5,5))
    Vij_tmp(4,3) = Vij_tmp(3,4)
    Vij_tmp(5,3) = Vij_tmp(3,5)
    Vij_tmp(5,4) = Vij_tmp(4,5)
!
    init = .false.
  end if
!
  Vij       = 0d0
  Vij_indep = 0d0
  Vij_dep   = 0d0
!
! The color-correlated virtual is trivial for 3 parton
! processes, we obtain the virtual then construct the 
! matrix.
! An extra conditional environment is introduced to tackle the
! case when the Laurent series is asked for:
  if (.not.present(VijLaurent)) then
    call VirtSMEddim(iptrn,p,mur,mode,Virt,Virt_dep,Virt_indep)
  else
    call VirtSMEddim(iptrn,p,mur,mode,Virt,Virt_dep,Virt_indep,VirtLaurent)
  end if
!
! The mode flag determines which contribution is to be calculated:
! color-correlated virtual is constructed according:
! Only the full contribution is calculated:
  if (mode.eq.'f') then
    Vij = Vij_tmp * Virt
! Only the \mu_R independent part is calculated:
  elseif (mode.eq.'i') then
    Vij_indep = Vij_tmp * Virt_indep
! Only the \mu_R dependent part is calculated:
  elseif (mode.eq.'d') then
    Vij_dep = Vij_tmp * Virt_dep
! Both the \mu_R dependent and independent part are calculated:
  elseif (mode.eq.'b') then
    Vij       = Vij_tmp * Virt
    Vij_dep   = Vij_tmp * Virt_dep
    Vij_indep = Vij_tmp * Virt_indep
  else
    print *,"Wrong mode selector is given to VijSME..."
    print *,"mode: ",mode
    stop
  end if
!
  if (present(VijLaurent)) then
    VijLaurent = 0d0
    VijLaurent(:,:,-2) = Vij_tmp*VirtLaurent(-2)
    VijLaurent(:,:,-1) = Vij_tmp*VirtLaurent(-1)
    VijLaurent(:,:, 0) = Vij_tmp*VirtLaurent( 0)
    VijLaurent(:,:, 1) = Vij_tmp*VirtLaurent( 1)
    VijLaurent(:,:, 2) = Vij_tmp*VirtLaurent( 2)
  end if
!
end subroutine VijSMEddim
!
subroutine VijSMEddim_new(iptrn,p,mur,VijLaurent,BijLaurent)
use particles
use QCDparams
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: mur
  real(kind(1d0)) , dimension(:,:,-4:) , intent(out) :: VijLaurent,BijLaurent
!
  logical , save :: init = .true.
  integer :: i
  real(kind(1d0)) , dimension(-4:2) :: VirtLaurent,BornLaurent
  real(kind(1d0)) , dimension(5,5) , save :: Vij_tmp
!
!
  interface
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
  if (init) then
! A matrix template is constructed for the color-correlated
! virtual:
    Vij_tmp = 0d0
!
    Vij_tmp(3,3) = qcd_cf
    Vij_tmp(4,4) = qcd_cf
    Vij_tmp(5,5) = qcd_ca
    Vij_tmp(3,4) = 0.5d0*(Vij_tmp(5,5) - Vij_tmp(3,3) - Vij_tmp(4,4))
    Vij_tmp(3,5) = 0.5d0*(Vij_tmp(4,4) - Vij_tmp(3,3) - Vij_tmp(5,5))
    Vij_tmp(4,5) = 0.5d0*(Vij_tmp(3,3) - Vij_tmp(4,4) - Vij_tmp(5,5))
    Vij_tmp(4,3) = Vij_tmp(3,4)
    Vij_tmp(5,3) = Vij_tmp(3,5)
    Vij_tmp(5,4) = Vij_tmp(4,5)
!
    init = .false.
  end if
!
  VijLaurent = 0
  BijLaurent = 0
!
  call VirtSMEddim_new(iptrn,p,mur,VirtLaurent,BornLaurent)
!
  do i=-2,2
    VijLaurent(:,:,i) = Vij_tmp(:,:)*VirtLaurent(i)
    BijLaurent(:,:,i) = Vij_tmp(:,:)*BornLaurent(i)
  end do
!
end subroutine VijSMEddim_new
!
subroutine VirtSMEddim_new(iptrn,p,mur,VirtLaurent,BornLaurent)
use process
use particles
use QCDparams
use my_model
!use math
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(in) :: mur
  real(kind(1d0)) , dimension(-4:2) , intent(out) :: VirtLaurent,BornLaurent
!
  integer :: i
  real(kind(1d0)) :: Q2,y12,y13,y23
  real(kind(1d0)) :: Born,xi
  real(kind(1d0)) ::  NC,nf,TR,nfg,pi,z3,eg
  real(kind(1d0)) :: L
  common/piez/ pi,z3
  common/const/ NC,nf,TR,eg
  real(kind(1d0)) , dimension(-4:2) :: VLaurent
  complex(kind(1d0)) , dimension(-4:2) :: f1yz,f1zy,f2yz,f2zy
  type(mom) :: Q
!
!
  interface 
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
  pi=3.14159265358979324d0
  z3=1.20205690315959429d0 !zeta(3)
  eg=0.57721566490153286d0 !eulergamma
!
  VirtLaurent = 0
!
! Calculation of invariants:
!
  Q = p(1)%p + p(2)%p
  Q2 = Q**2
!
  y12 = 2*p(3)%p*p(4)%p/Q2
  y13 = 2*p(3)%p*p(5)%p/Q2
  y23 = 2*p(4)%p*p(5)%p/Q2
!
! To construct a Virtual with a d-dimensional Born we need
! the complete Laurent series of both:
  call BornSMEddim(iptrn,p,Born,BornLaurent)
!
! Calculation of mu independent part:
  call f1(f1yz,y13,y23)
  call f1(f1zy,y23,y13)
  call f2(f2yz,y13,y23)
  call f2(f2zy,y23,y13)
  do i=-4,2
    VirtLaurent(i) = 8*qcd_nc*qcd_cf*real(qcd_nc*(f1yz(i)+f1zy(i)) &
                   + 1d0/qcd_nc*(f2yz(i)+f2zy(i)))
  end do
! Pattern 1 is for d d~ g and pattern 2 is for u u~ g:
  if (iptrn.eq.2) VirtLaurent = 4*VirtLaurent
  VirtLaurent = VirtLaurent/Q2/27d0
! Renormalization:
! Note the extra factor of 1/2 compensation for a factor of 2
! difference between the beta0 definition used in the 
! framework and in the SME:
  do i=-4,1
    VirtLaurent(i) = VirtLaurent(i) - qcd_beta0*BornLaurent(i+1)
  end do
!
! Calculation of the mu dependent part:
!
! xi = mur**2/Q**2:
  L = log(mur**2/Q2)
!
  VLaurent = 0
!
  VLaurent(-2) = 0
  VLaurent(-1) = VirtLaurent(-2)*L
  VLaurent( 0) = 0.5d0*VirtLaurent(-2)*L**2 &
               + (VirtLaurent(-1) &
               + qcd_beta0*BornLaurent(0))*L
  VLaurent( 1) = 1d0/6d0*VirtLaurent(-2)*L**3 &
               + 0.5d0*(VirtLaurent(-1) &
               + qcd_beta0*BornLaurent(0))*L**2 &
               + (VirtLaurent( 0) &
               + qcd_beta0*BornLaurent(1))*L
  VLaurent( 2) = 1d0/24d0*VirtLaurent(-2)*L**4 &
               + 1d0/6d0*(VirtLaurent(-1) &
               + qcd_beta0*BornLaurent(0))*L**3 &
               + 0.5d0*(VirtLaurent( 0) &
               + qcd_beta0*BornLaurent(1))*L**2 &
               + (VirtLaurent( 1) &
               + qcd_beta0*BornLaurent(2))*L
!
  VirtLaurent = VLaurent + VirtLaurent
!
end subroutine VirtSMEddim_new
