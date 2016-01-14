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
  integer :: ipart
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
  real(kind(1d0)) , external :: PSI2d1gBornV,PSI2u1gBornV
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
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
    smeR = PSI2d1gBornV(1,3,2,8,7,m2ij)
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Rmunu)
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    smeR = PSI2u1gBornV(1,3,2,8,7,m2ij)
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Rmunu)
    return
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
  Rij = 0d0
!
  Rij(3,3) = qcd_cf
  Rij(4,4) = qcd_cf
  Rij(5,5) = qcd_ca
  Rij(3,4) = 0.5d0*(Rij(5,5) - Rij(3,3) - Rij(4,4))
  Rij(3,5) = 0.5d0*(Rij(4,4) - Rij(3,3) - Rij(5,5))
  Rij(4,5) = 0.5d0*(Rij(3,3) - Rij(4,4) - Rij(5,5))
  Rij(4,3) = Rij(3,4)
  Rij(5,3) = Rij(3,5)
  Rij(5,4) = Rij(4,5)
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
    Rij  = Rij * smeR
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    smeR = PSI2u1gBorn(1,3,2,8,7)
    Rij  = Rij * smeR
    return
  end if
!
end subroutine RijSME
