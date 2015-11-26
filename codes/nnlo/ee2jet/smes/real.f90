subroutine RealSME(iptrn,pin,smeR)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(out) :: smeR
!
  integer :: ipart
  logical , save :: init = .true.
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  complex(kind(1d0)) , dimension(2,2,2) :: C
  common/COUPLINGS/ C &
        /QCD/XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  common/dotproducts/S &
        /spinorproducts/A,B
!
  real(kind(1d0)) , external :: PSI2d1gBorn,PSI2u1gBorn
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
!
  end if
!
! We have to setup the dot and spinor products:
! To calculate these products the momenta have to be given in
! an ordinary array:
  P = 0d0
  do ipart=1,nleg_born+1
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
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
    smeR = PSI2d1gBorn(1,3,2,8,7)
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    smeR = PSI2u1gBorn(1,3,2,8,7)
    return
  end if
!
end subroutine RealSME
!
subroutine RmunuSME(iptrn,ileg,pin,Rmunu)
use process
use particles
use QCDparams
use my_model
use misc
implicit none
!
  integer , intent(in) :: iptrn,ileg
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Rmunu
!
  integer :: ipart,qtype
  logical , save :: init = .true.
  real(kind(1d0)) :: smeR
  complex(kind(1d0)) , dimension(2,2) :: m2ij
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  complex(kind(1d0)) , dimension(2,2,2) :: C
  common/COUPLINGS/ C &
        /QCD/XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  common/dotproducts/S &
        /spinorproducts/A,B
!
  real(kind(1d0)) , external :: PSI2q2gBornV
  real(kind(1d0)) , external :: PSI2q2gBAKV
!
  real(kind(1d0)) , external :: PSI2d2gBorn,PSI2u2gBorn, &
                                PSI2u2d,PSI2u2c,PSI2d2s, &
                                PSI4u,PSI4usl,PSI4d,PSI4dsl
!
  print *,"Using 3jet version of RmunuSME"
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
!
  end if
!
! We have to setup the dot and spinor products:
! To calculate these products the momenta have to be given in
! an ordinary array:
  P = 0d0
  do ipart=1,nleg_born+1
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
  Rmunu = 0d0
!
! e+ e- -> d d~ g g
  if (iptrn.eq.1) then
!    print *,"e+ e- -> d d~ g g"
    qtype = 2
    if (ileg.eq.5) then
      smeR = PSI2q2gBAKV(1,3,4,2,8,7,qtype,m2ij)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Rmunu)
      return
    elseif (ileg.eq.6) then
      smeR = PSI2q2gBAKV(1,4,3,2,8,7,qtype,m2ij)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Rmunu)
      return
    end if
! e+ e- -> u u~ g g
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g g"
    qtype = 1
!
    if (ileg.eq.5) then
      smeR = PSI2q2gBAKV(1,3,4,2,8,7,qtype,m2ij)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Rmunu)
      return
    elseif (ileg.eq.6) then
      smeR = PSI2q2gBAKV(1,4,3,2,8,7,qtype,m2ij)
      call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Rmunu)
      return
    end if
    print *,"Error in RmunuSME..."
    print *,"iptrn,ileg: ",iptrn,ileg
    stop
! e+ e- -> u u~ d d~
  elseif (iptrn.eq.3) then
!    print *,"e+ e- -> u u~ d d~"
    print *,"No gluon present, Rmunu..."
    print *,"iptrn: ",iptrn
    stop
    return
! e+ e- -> u u~ u u~
  elseif (iptrn.eq.4) then
!    print *,"e+ e- -> u u~ u u~"
    print *,"No gluon present, Rmunu..."
    print *,"iptrn: ",iptrn
    stop
    return
! e+ e- -> d d~ d d~
  elseif (iptrn.eq.5) then
!    print *,"e+ e- -> d d~ d d~"
    print *,"No gluon present, Rmunu..."
    print *,"iptrn: ",iptrn
    stop
    return
! e+ e- -> u u~ c c~
  elseif (iptrn.eq.6) then
!    print *,"e+ e- -> u u~ c c~"
    print *,"No gluon present, Rmunu..."
    print *,"iptrn: ",iptrn
    stop
    return
! e+ e- -> d d~ s s~
  elseif (iptrn.eq.7) then
!    print *,"e+ e- -> d d~ s s~"
    print *,"No gluon present, Rmunu..."
    print *,"iptrn: ",iptrn
    stop
    return
  else 
    print *,"unknown real SCSME is asked for..."
    print *,"iptrn: ",iptrn
    stop
  end if
!
end subroutine RmunuSME
!
subroutine RijSME(iptrn,pin,Rij)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Rij
!
  integer :: ipart
  integer :: qtype 
  logical , save :: init = .true.
  real(kind(1d0)) :: smeR
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  complex(kind(1d0)) , dimension(2,2,2) :: C
  common/COUPLINGS/ C &
        /QCD/XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  common/dotproducts/S &
        /spinorproducts/A,B
!
  real(kind(1d0)) , external :: CCM12INT,CCM13INT,CCM14INT, &
                                CCM23INT,CCM24INT,CCM34INT, &
                                PSI4u,PSI4d,PSI4usl,PSI4dsl, &
                                PSI2u2d,PSI2u2c,PSI2d2s, &
                                PSI4qBorn12,PSI4qBorn13,PSI4qBorn14, &
                                PSI2d2gBorn,PSI2u2gBorn, &
                                PSI2u2dBorn12,PSI2u2dBorn13,PSI2u2dBorn14, &
                                PSI2u2cBorn12,PSI2u2cBorn13,PSI2u2cBorn14, &
                                PSI2d2sBorn12,PSI2d2sBorn13,PSI2d2sBorn14
!
  print *,"Using 3jet version of RijSME"
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
!
  end if
!
! We have to setup the dot and spinor products:
! To calculate these products the momenta have to be given in
! an ordinary array:
  P = 0d0
  do ipart=1,nleg_born+1
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
  Rij = 0d0
!
! e+ e- -> d d~ g g
  if (iptrn.eq.1) then
!    print *,"e+ e- -> d d~ g g"
    qtype = 2
    smeR = PSI2d2gBorn(1,3,4,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! g(1)-g(1)
    Rij(5,5) = qcd_ca * smeR
! g(2)-g(2)
    Rij(6,6) = Rij(5,5)
! q-q~
    Rij(3,4) = CCM14INT(1,4,3,2,8,7,qtype)
    Rij(4,3) = Rij(3,4)
! q-g(1)
    Rij(3,5) = CCM12INT(1,4,3,2,8,7,qtype)
    Rij(5,3) = Rij(3,5)
! q-g(2)
    Rij(3,6) = CCM13INT(1,4,3,2,8,7,qtype)
    Rij(6,3) = Rij(3,6)
! q~-g(1)
    Rij(4,5) = CCM24INT(1,4,3,2,8,7,qtype)
    Rij(5,4) = Rij(4,5)
! q~-g(2)
    Rij(4,6) = CCM34INT(1,4,3,2,8,7,qtype)
    Rij(6,4) = Rij(4,6)
! g(1)-g(2)
    Rij(5,6) = - Rij(5,5) - Rij(3,5) - Rij(4,5)
    Rij(6,5) = Rij(5,6)
    return
! e+ e- -> u u~ g g
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g g"
    qtype = 1
    smeR = PSI2u2gBorn(1,3,4,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! g(1)-g(1)
    Rij(5,5) = qcd_ca * smeR
! g(2)-g(2)
    Rij(6,6) = Rij(5,5)
! q-q~
    Rij(3,4) = CCM14INT(1,4,3,2,8,7,qtype)
    Rij(4,3) = Rij(3,4)
! q-g(1)
    Rij(3,5) = CCM12INT(1,4,3,2,8,7,qtype)
    Rij(5,3) = Rij(3,5)
! q-g(2)
    Rij(3,6) = CCM13INT(1,4,3,2,8,7,qtype)
    Rij(6,3) = Rij(3,6)
! q~-g(1)
    Rij(4,5) = CCM24INT(1,4,3,2,8,7,qtype)
    Rij(5,4) = Rij(4,5)
! q~-g(2)
    Rij(4,6) = CCM34INT(1,4,3,2,8,7,qtype)
    Rij(6,4) = Rij(4,6)
! g(1)-g(2)
    Rij(5,6) = CCM23INT(1,4,3,2,8,7,qtype)
    Rij(6,5) = Rij(5,6)
    return
! e+ e- -> u u~ d d~
  elseif (iptrn.eq.3) then
!    print *,"e+ e- -> u u~ d d~"
    smeR = PSI2u2d(1,4,3,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! Q-Q
    Rij(5,5) = Rij(3,3)
! Q~-Q~
    Rij(6,6) = Rij(3,3)
! q-q~
    Rij(3,4) = PSI2u2dBorn12(1,4,3,2,8,7)
    Rij(4,3) = Rij(3,4)
! q-Q
    Rij(3,5) = Psi2u2dBorn13(1,4,3,2,8,7)
    Rij(5,3) = Rij(3,5)
! q-Q~
    Rij(3,6) = Psi2u2dBorn14(1,4,3,2,8,7)
    Rij(6,3) = Rij(3,6)
! q~-Q
    Rij(4,5) = 0.5d0*(Rij(6,6) - Rij(3,3) - Rij(4,4) - Rij(5,5)) &
             - Rij(3,4) - Rij(3,5)
    Rij(5,4) = Rij(4,5)
! q~-Q~
    Rij(4,6) = Rij(3,5)
    Rij(6,4) = Rij(4,6)
! Q-Q~
    Rij(5,6) = Rij(3,4)
    Rij(6,5) = Rij(5,6)
    return
! e+ e- -> u u~ u u~
  elseif (iptrn.eq.4) then
!    print *,"e+ e- -> u u~ u u~"
    qtype = 1
    smeR = PSI4u(1,4,3,2,8,7) + PSI4usl(1,4,3,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! q-q
    Rij(5,5) = Rij(3,3)
! q~-q~
    Rij(6,6) = Rij(3,3)
! q-q~
    Rij(3,4) = PSI4qBorn12(1,4,3,2,8,7,qtype)
    Rij(4,3) = Rij(3,4)
! q-Q
    Rij(3,5) = PSI4qBorn13(1,4,3,2,8,7,qtype)
    Rij(5,3) = Rij(3,5)
! q-Q~
    Rij(3,6) = PSI4qBorn14(1,4,3,2,8,7,qtype)
    Rij(6,3) = Rij(3,6)
! The remaining color connections are worked out following the
! appendix of arXiv:hep-ph/9605323, cf. (A.7):
! q~-Q
    Rij(4,5) = 0.5d0*(Rij(6,6) - Rij(3,3) - Rij(4,4) - Rij(5,5)) &
             - Rij(3,4) - Rij(3,5)
    Rij(5,4) = Rij(4,5)
! q~-Q~
    Rij(4,6) = Rij(3,5)
    Rij(6,4) = Rij(4,6)
! Q-Q~
    Rij(5,6) = Rij(3,4)
    Rij(6,5) = Rij(5,6)
    return
! e+ e- -> d d~ d d~
  elseif (iptrn.eq.5) then
!    print *,"e+ e- -> d d~ d d~"
    qtype = 2
    smeR = PSI4d(1,4,3,2,8,7) + PSI4dsl(1,4,3,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! q-q
    Rij(5,5) = Rij(3,3)
! q~-q~
    Rij(6,6) = Rij(3,3)
! q-q~
    Rij(3,4) = PSI4qBorn12(1,4,3,2,8,7,qtype)
    Rij(4,3) = Rij(3,4)
! q-Q
    Rij(3,5) = PSI4qBorn13(1,4,3,2,8,7,qtype)
    Rij(5,3) = Rij(3,5)
! q-Q~
    Rij(3,6) = PSI4qBorn14(1,4,3,2,8,7,qtype)
    Rij(6,3) = Rij(3,6)
! q~-Q
    Rij(4,5) = 0.5d0*(Rij(6,6) - Rij(3,3) - Rij(4,4) - Rij(5,5)) &
             - Rij(3,4) - Rij(3,5)
    Rij(5,4) = Rij(4,5)
! q~-Q~
    Rij(4,6) = Rij(3,5)
    Rij(6,4) = Rij(4,6)
! Q-Q~
    Rij(5,6) = Rij(3,4)
    Rij(6,5) = Rij(5,6)
    return
! e+ e- -> u u~ c c~
  elseif (iptrn.eq.6) then
!    print *,"e+ e- -> u u~ c c~"
    smeR = PSI2u2c(1,4,3,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! Q-Q
    Rij(5,5) = Rij(3,3)
! Q~-Q~
    Rij(6,6) = Rij(3,3)
! q-q~
    Rij(3,4) = PSI2u2cBorn12(1,4,3,2,8,7)
    Rij(4,3) = Rij(3,4)
! q-Q
    Rij(3,5) = Psi2u2cBorn13(1,4,3,2,8,7)
    Rij(5,3) = Rij(3,5)
! q-Q~
    Rij(3,6) = Psi2u2cBorn14(1,4,3,2,8,7)
    Rij(6,3) = Rij(3,6)
! q~-Q
    Rij(4,5) = 0.5d0*(Rij(6,6) - Rij(3,3) - Rij(4,4) - Rij(5,5)) &
             - Rij(3,4) - Rij(3,5)
    Rij(5,4) = Rij(4,5)
! q~-Q~
    Rij(4,6) = Rij(3,5)
    Rij(6,4) = Rij(4,6)
! Q-Q~
    Rij(5,6) = Rij(3,4)
    Rij(6,5) = Rij(5,6)
    return
! e+ e- -> d d~ s s~
  elseif (iptrn.eq.7) then
!    print *,"e+ e- -> d d~ s s~"
    smeR = PSI2d2s(1,4,3,2,8,7)
! q-q
    Rij(3,3) = qcd_cf * smeR
! q~-q~
    Rij(4,4) = Rij(3,3)
! Q-Q
    Rij(5,5) = Rij(3,3)
! Q~-Q~
    Rij(6,6) = Rij(3,3)
! q-q~
    Rij(3,4) = PSI2d2sBorn12(1,4,3,2,8,7)
    Rij(4,3) = Rij(3,4)
! q-Q
    Rij(3,5) = Psi2d2sBorn13(1,4,3,2,8,7)
    Rij(5,3) = Rij(3,5)
! q-Q~
    Rij(3,6) = Psi2d2sBorn14(1,4,3,2,8,7)
    Rij(6,3) = Rij(3,6)
! q~-Q
    Rij(4,5) = 0.5d0*(Rij(6,6) - Rij(3,3) - Rij(4,4) - Rij(5,5)) &
             - Rij(3,4) - Rij(3,5)
    Rij(5,4) = Rij(4,5)
! q~-Q~
    Rij(4,6) = Rij(3,5)
    Rij(6,4) = Rij(4,6)
! Q-Q~
    Rij(5,6) = Rij(3,4)
    Rij(6,5) = Rij(5,6)
    return
  else 
    print *,"unknown RijSME is asked for..."
    print *,"iptrn: ",iptrn
    stop
  end if
!
end subroutine RijSME
