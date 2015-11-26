subroutine BornSME(iptrn,pin,Born)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(out) :: Born
!
  integer :: ipart
  logical , save :: init = .true.
!
! In what follows is ugly as hell but we still have some
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
  real(kind(1d0)) , external :: PSI2dBorn, PSI2uBorn
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
  do ipart=1, nleg_born
! Note that we are working in an all-outgoing scheme, hence the
! incoming momenta should acquire an additional minus sign:
! Note, too that the position for the incoming particles are at the 
! 7th and 8th positions:
! p1,p2 => -p7,-p8
  if (ipart.le.2) then
    P(:,ipart+6) = -pin(ipart)%p
! Final state momenta should be accordingly shifted to the left:
! p3,p4 => p1,p2
  else
    P(:,ipart-2) = pin(ipart)%p
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
! Ordering is important here! From the q2g1 case:
! Probably we have: q,qb,e+,e- => q,qb,e-,e+
!   Note that the position of the quark and the antiquark is interchanged.
!   To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
!   to 1,3,2,8,7
!   7e+ 8e- 1q 2qb 3g
!   The ordering along momenta: q,qb,g,e+,e-
!   The variable called iptrn determines which contribution we should
!   calculate:
! e+ e- -> d d~
  if (iptrn.eq.1) then
    Born = PSI2dBorn(1,2,8,7)
    return
! e+ e- -> u u~
  elseif (iptrn.eq.2) then
    Born = PSI2uBorn(1,2,8,7)
    return
  end if
!
end subroutine BornSME
!
subroutine BmunuSME(iptrn,ileg,pin,Bmunu)
use process
use particles
use QCDparams
use my_model
use misc
implicit none
!
  integer , intent(in) :: iptrn,ileg
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
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
  print *,"Using 3jet version of BmunuSME"
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
  Bmunu = 0d0
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
    Born = PSI2d1gBornV(1,3,2,8,7,m2ij)
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Bmunu)
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    Born = PSI2u1gBornV(1,3,2,8,7,m2ij)
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    call CastHelicityToLorentz(pin(ileg),pin(3),m2ij,Bmunu)
    return
  end if
!
end subroutine BmunuSME
!
subroutine BmunuSME_q1q2(iptrn,ileg,q1,q2,pin,Bmunu)
use process
use particles
use QCDparams
use my_model
use misc
implicit none
!
  integer , intent(in) :: iptrn,ileg
  type(mom) , intent(in) :: q1,q2
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
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
  print *,"Using 3jet version of BmunuSME_q1q2"
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
  Bmunu = 0d0
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
    Born = PSI2d1gBornV(1,3,2,8,7,m2ij)
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    call CastHelicityToLorentz(pin(ileg),q1,q2,m2ij,Bmunu)
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    Born = PSI2u1gBornV(1,3,2,8,7,m2ij)
! We give the quark as a reference particle which is needed to
! obtain the polarization vectors:
    call CastHelicityToLorentz(pin(ileg),q1,q2,m2ij,Bmunu)
    return
  end if
!
end subroutine BmunuSME_q1q2
!
subroutine BijSME(iptrn,pin,Bij)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
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
  common/COUPLINGS/ C &
        /QCD/XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  common/dotproducts/S &
        /spinorproducts/A,B
!
  real(kind(1d0)) , external :: PSI2dBorn,PSI2uBorn
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
  Bij = 0d0
!
! 3q 4qb 5g
  Bij(3,3) = qcd_cf
  Bij(4,4) = qcd_cf
  Bij(3,4) = -qcd_cf
  Bij(4,3) = Bij(3,4)
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,7,8
! to 1,2,8,7.
! The ordering among momenta: q,qb,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~
  if (iptrn.eq.1) then
    Born = PSI2dBorn(1,2,8,7)
    Bij  = Bij * Born
    return
! e+ e- -> u u~
  elseif (iptrn.eq.2) then
    Born = PSI2uBorn(1,2,8,7)
    Bij  = Bij * Born
    return
  end if
!
end subroutine BijSME
!
subroutine BijklSME(iptrn,pin,Bijkl)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
!
  integer :: i,ipart,jpart,kpart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
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
  print *,"Using 3jet version of BijklSME"
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
  Bijkl = 0d0
!
! no. 5 particle is the gluon:
! Filling up (TiTi)(TiTi):
  do ipart=3,5
    if (ipart.ne.5) Bijkl(ipart,ipart,ipart,ipart) = qcd_cf**2
    if (ipart.eq.5) Bijkl(ipart,ipart,ipart,ipart) = qcd_ca**2
!    print *,"(TiTi)(TiTi): ",ipart
  end do
! Filling up (TiTi)(TjTj) , i != j:
  do ipart=3,4
    do jpart=ipart+1,5
      if (jpart.ne.5) Bijkl(ipart,ipart,jpart,jpart) = qcd_cf**2
      if (jpart.eq.5) Bijkl(ipart,ipart,jpart,jpart) = qcd_cf*qcd_ca
      Bijkl(jpart,jpart,ipart,ipart) = Bijkl(ipart,ipart,jpart,jpart)
!    print *,"(TiTi)(TjTj): ",ipart,jpart
    end do
  end do
! Filling up (TiTi)(TjTk):
  do ipart=3,5
    do jpart=3,4
      do kpart=jpart+1,5
        if (ipart.ne.5) Bijkl(ipart,ipart,jpart,kpart) = qcd_cf
        if (ipart.eq.5) Bijkl(ipart,ipart,jpart,kpart) = qcd_ca
        if (kpart.ne.5) then
          Bijkl(ipart,ipart,jpart,kpart) = Bijkl(ipart,ipart,jpart,kpart) &
                                         * 0.5d0*(qcd_ca - 2*qcd_cf)
        else
          Bijkl(ipart,ipart,jpart,kpart) = Bijkl(ipart,ipart,jpart,kpart) &
                                         * (-0.5d0*qcd_ca)
        end if
        Bijkl(ipart,ipart,kpart,jpart) = Bijkl(ipart,ipart,jpart,kpart)
        Bijkl(jpart,kpart,ipart,ipart) = Bijkl(ipart,ipart,jpart,kpart)
        Bijkl(kpart,jpart,ipart,ipart) = Bijkl(ipart,ipart,jpart,kpart)
      end do
!    print *,"(TiTi)(TjTk): ",ipart,jpart,kpart
    end do
  end do
! Filling up (TiTj)(TiTj) , i != j:
  do ipart=3,4
    do jpart=ipart+1,5
      if (jpart.ne.5) Bijkl(ipart,jpart,ipart,jpart) = 0.25d0/qcd_nc**2
      if (jpart.eq.5) Bijkl(ipart,jpart,ipart,jpart) = 0.25d0*qcd_nc**2
      Bijkl(jpart,ipart,ipart,jpart) = Bijkl(ipart,jpart,ipart,jpart)
      Bijkl(ipart,jpart,jpart,ipart) = Bijkl(ipart,jpart,ipart,jpart)
      Bijkl(jpart,ipart,jpart,ipart) = Bijkl(ipart,jpart,ipart,jpart)
!    print *,"(TiTj)(TiTj): ",ipart,jpart
    end do
  end do
! Filling up (TiTj)(TiTk) , i != j , i != k:
  do i=1,3
    ipart=i+2 ; jpart = mod(i,3)+3 ; kpart = mod(i+1,3)+3
!    print *,"(TiTj)(TiTk): ",ipart,jpart,kpart
! (TiTj)(TiTk)
    if (((ipart+jpart).eq.7).or.((ipart+kpart).eq.7)) then
! The configuration is (T1T2)(T1T3) or (T2T3)(T2T1)
      Bijkl(ipart,jpart,ipart,kpart) = -0.25d0
    else 
! The configuration is (T3T1)(T3T2)
      Bijkl(ipart,jpart,ipart,kpart) = 0.25d0*qcd_nc**2
    end if
! (TjTi)(TiTk)
    Bijkl(jpart,ipart,ipart,kpart) = Bijkl(ipart,jpart,ipart,kpart)
! (TiTj)(TkTi)
    Bijkl(ipart,jpart,kpart,ipart) = Bijkl(ipart,jpart,ipart,kpart)
! (TjTi)(TkTi)
    Bijkl(jpart,ipart,kpart,ipart) = Bijkl(ipart,jpart,ipart,kpart)
! (TiTk)(TiTj)
    Bijkl(ipart,kpart,ipart,jpart) = Bijkl(ipart,jpart,ipart,kpart)
! (TkTi)(TiTj)
    Bijkl(kpart,ipart,ipart,jpart) = Bijkl(ipart,jpart,ipart,kpart)
! (TiTk)(TjTi)
    Bijkl(ipart,kpart,jpart,ipart) = Bijkl(ipart,jpart,ipart,kpart)
! (TkTi)(TjTi)
    Bijkl(kpart,ipart,jpart,ipart) = Bijkl(ipart,jpart,ipart,kpart)
  end do
!
  Bijkl = 2d0*Bijkl
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
    Born = PSI2d1gBorn(1,3,2,8,7)
    Bijkl  = Bijkl * Born
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    Born = PSI2u1gBorn(1,3,2,8,7)
    Bijkl  = Bijkl * Born
    return
  end if
!
end subroutine BijklSME
!
subroutine BmunuijSME(iptrn,ileg,pin,Bmunuij)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  integer , intent(in) :: ileg
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , dimension(0:,0:,:,:) , intent(out) :: Bmunuij
!
  integer :: mu,nu
  real(kind(1d0)) , dimension(0:3,0:3) :: Bmunu
!
!
  interface
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
  print *,"Using 3jet version of BmunuijSME"
!
  Bmunuij = 0d0
!
! The color-correlated SME is trivial and the color of the Born factorizes
! Bmunuij can be directly obtained from Bmunu:
  call BmunuSME(iptrn,ileg,pin,Bmunu)
!
  do mu=0,3
    do nu=0,3
      Bmunuij(mu,nu,3,3) = qcd_cf*Bmunu(mu,nu)
      Bmunuij(mu,nu,4,4) = qcd_cf*Bmunu(mu,nu)
      Bmunuij(mu,nu,5,5) = qcd_ca*Bmunu(mu,nu)
      Bmunuij(mu,nu,3,4) = 0.5d0*(qcd_ca - 2*qcd_cf)*Bmunu(mu,nu)
      Bmunuij(mu,nu,4,3) = 0.5d0*(qcd_ca - 2*qcd_cf)*Bmunu(mu,nu)
      Bmunuij(mu,nu,3,5) = -0.5d0*qcd_ca*Bmunu(mu,nu)
      Bmunuij(mu,nu,4,5) = -0.5d0*qcd_ca*Bmunu(mu,nu)
      Bmunuij(mu,nu,5,3) = -0.5d0*qcd_ca*Bmunu(mu,nu)
      Bmunuij(mu,nu,5,4) = -0.5d0*qcd_ca*Bmunu(mu,nu)
    end do
  end do
!
end subroutine BmunuijSME
!
! This routine implements the Bmunu as it can be found in 
! arXiv:hep-ph/9605323:
subroutine BmunuCS(p,Bmunu)
use particles
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
  integer :: ipart,mu,nu
  real(kind(1d0)) :: smeB
  real(kind(1d0)) :: x1,x2,Q2
  real(kind(1d0)) , dimension(0:3) :: p1a,p2a,p3a
  real(kind(1d0)) , dimension(0:3,0:3) :: gmunu
  type(mom) :: Q,p1,p2,p3
!
!
  interface
    subroutine CalcB(p,smeB)
    use particles
    implicit none
!
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: smeB
!
    end subroutine CalcB
  end interface
!
  print *,"Using 3jet version of BmunuCS"
!
  Bmunu = 0
!
  gmunu(0,:) = (/1,0,0,0/)
  gmunu(1,:) = (/0,-1,0,0/)
  gmunu(2,:) = (/0,0,-1,0/)
  gmunu(3,:) = (/0,0,0,-1/)
! We have to obtain the SME first:
  call CalcB(p,smeB)
  print *,"smeB: ",smeB
!
  print *,"gmunu: "
  do mu=0,3
    print *,(gmunu(mu,nu),nu=0,3)
  end do 
  call PrintParts(p)
!
! we have to identify the various partons in the final state
! and assign value to the appropriate arrays:
  do ipart=3,5
    if (p(ipart)%flv.gt.0) p1 = p(ipart)%p
    if (p(ipart)%flv.lt.0) p2 = p(ipart)%p
    if (p(ipart)%flv.eq.0) p3 = p(ipart)%p
  end do
  print *,"quark: "
  call PrintMom(p1)
  print *,"antiquark: "
  call PrintMom(p2)
  print *,"gluon: "
  call PrintMom(p3)
!
  p1a(0) = p1%E
  p1a(1) = p1%px
  p1a(2) = p1%py
  p1a(3) = p1%pz
!
  p2a(0) = p2%E
  p2a(1) = p2%px
  p2a(2) = p2%py
  p2a(3) = p2%pz
!
  p3a(0) = p3%E
  p3a(1) = p3%px
  p3a(2) = p3%py
  p3a(3) = p3%pz
!
  Q  = p(1)%p + p(2)%p
  Q2 = Q**2
!
  x1 = 2*p1*Q/Q2
  x2 = 2*p2*Q/Q2
!
  do mu=0,3
    do nu=0,3
      Bmunu(mu,nu) = 2*p1a(mu)*p2a(nu)/Q2 + 2*p2a(mu)*p1a(nu)/Q2 &
                   - 2*(1 - x1)*p1a(mu)*p1a(nu)/((1 - x2)*Q2) &
                   - 2*(1 - x2)*p2a(mu)*p2a(nu)/((1 - x1)*Q2) &
                   - (1 - x1 - x2 + x2**2)/((1 - x2)*Q2) &
                   * (p1a(mu)*p3a(nu) + p3a(mu)*p1a(nu)) &
                   - (1 - x2 - x1 + x1**2)/((1 - x1)*Q2) &
                   * (p2a(mu)*p3a(nu) + p3a(mu)*p2a(nu)) &
                   + (1 + 0.5d0*x1**2 + 0.5d0*x2**2 - x1 - x2) &
                   * gmunu(mu,nu)
    end do
  end do
  Bmunu = -Bmunu * smeB / (x1**2 + x2**2)
!
end subroutine BmunuCS
!
subroutine BalbeSME(iptrn,ileg,pin,Balbe)
use process
use particles
use QCDparams
use my_model
use misc
implicit none
!
  integer , intent(in) :: iptrn,ileg
  type(particle) , dimension(:) , intent(in) :: pin
  complex(kind(1d0)) , dimension(2,2) , intent(out) :: Balbe
!
  integer :: ipart
  logical , save :: init = .true.
  real(kind(1d0)) :: Born
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
  print *,"Using 3jet version of BalbeSME"
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
  Balbe = 0d0
!
! Note that the position of the quark and the antiquark is interchanged.
! To get agreement with HELAC we had to change the ordering from 1,2,3,7,8
! to 1,3,2,8,7.
! The ordering among momenta: q,qb,g,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g
  if (iptrn.eq.1) then
    Born = PSI2d1gBornV(1,3,2,8,7,m2ij)
    Balbe = m2ij
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
    Born = PSI2u1gBornV(1,3,2,8,7,m2ij)
    Balbe = m2ij
    return
  end if
!
end subroutine BalbeSME
!
subroutine BornSME_AK(p,smeB)
use particles
use QCDparams
use math
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(out) :: smeB
!
  integer :: ipart
  real(kind(1d0)) :: qq
  real(kind(1d0)) :: Q2,s,t,u,x1,x2
  type(mom) :: p1,p2,p3,p4,p5,Q
!
  print *,"Using 3jet version of BornSME_AK"
!
  p1 = p(1)%p ; p2 = p(2)%p
!
  call PrintParts(p)
!
! we have to identify the various partons in the final state
! and assign value to the appropriate arrays:
  do ipart=3,5
    if (p(ipart)%flv.gt.0) p3 = p(ipart)%p
    if (p(ipart)%flv.lt.0) p4 = p(ipart)%p
    if (p(ipart)%flv.eq.0) p5 = p(ipart)%p
  end do
  print *,"quark: "
  call PrintMom(p3)
  print *,"antiquark: "
  call PrintMom(p4)
  print *,"gluon: "
  call PrintMom(p5)
! 
! Determining the change of the quark:
  if (mod(max(p(3)%flv,p(4)%flv,p(5)%flv),2).eq.0) then
    qq = 2d0/3d0
  else
    qq = 1d0/3d0
  end if
!
  Q = p1 + p2
  s = Q**2
  Q2 = s
! The minus signs appear due to crossing from the all outgoing channel:
  t = -2*p1*p3
  u = -2*p1*p4
!
  x1 = 2*p3*Q/Q2
  x2 = 2*p4*Q/Q2
!
  smeB = 2*t**2 + 2*u**2 + 2*s*t*x1 + 2*s*u*x2 + s**2*(x1**2 + x2**2)
  smeB = smeB/(s**3*(1 - x1)*(1 - x2))
  smeB = (4*pi)**3*2*2*qcd_nc*qcd_cf*qq**2*smeB
!
end subroutine BornSME_AK
!
subroutine Bmunu_AK(p,Bmunu)
use particles
use math
use QCDparams
implicit none
!
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(0:3,0:3) , intent(out) :: Bmunu
!
  integer :: ipart,mu,nu
  real(kind(1d0)) , dimension(0:3) :: p1a,p3a,p4a,p5a
  real(kind(1d0)) , dimension(0:3,0:3) :: gmunu
  real(kind(1d0)) :: qq
  real(kind(1d0)) :: Q2,s,t,u,x1,x2
  type(mom) :: p1,p2,p3,p4,p5,Q
!
  print *,"Using 3jet version of Bmunu_AK"
!
  p1 = p(1)%p ; p2 = p(2)%p
!
  Bmunu = 0
!
  gmunu(0,:) = (/1,0,0,0/)
  gmunu(1,:) = (/0,-1,0,0/)
  gmunu(2,:) = (/0,0,-1,0/)
  gmunu(3,:) = (/0,0,0,-1/)
!
! we have to identify the various partons in the final state
! and assign value to the appropriate arrays:
  do ipart=3,5
    if (p(ipart)%flv.gt.0) p3 = p(ipart)%p
    if (p(ipart)%flv.lt.0) p4 = p(ipart)%p
    if (p(ipart)%flv.eq.0) p5 = p(ipart)%p
  end do
!
! p1 changes sign due to crossing:
  p1a(0) = -p1%E
  p1a(1) = -p1%px
  p1a(2) = -p1%py
  p1a(3) = -p1%pz
!
  p3a(0) = p3%E
  p3a(1) = p3%px
  p3a(2) = p3%py
  p3a(3) = p3%pz
!
  p4a(0) = p4%E
  p4a(1) = p4%px
  p4a(2) = p4%py
  p4a(3) = p4%pz
!
  p5a(0) = p5%E
  p5a(1) = p5%px
  p5a(2) = p5%py
  p5a(3) = p5%pz
!
! Determining the change of the quark:
  if (mod(max(p(3)%flv,p(4)%flv,p(5)%flv),2).eq.0) then
    qq = 2d0/3d0
  else
    qq = 1d0/3d0
  end if
!
  Q = p1 + p2
  s = Q**2
  Q2 = s
! The minus signs appear due to crossing from the all outgoing channel:
  t = -2*p1*p3
  u = -2*p1*p4
!
  x1 = 2*p3*Q/Q2
  x2 = 2*p4*Q/Q2
!
  do mu=0,3
    do nu=0,3
      Bmunu(mu,nu) = &
      (qcd_cf*qcd_nc*((16*(s**2*(1-x1+x1**2-x2)-2*s*t*(-1+x2)+  &
      2*t*(t+u-u*x1-t*x2))*p4a(nu)*p5a(mu))/(s**4*(-1+  &
      x1)**2*(-1+x2))+p1a(nu)*((32*(s+2*u)*p3a(mu))/(s**3*(-1+x2))+  &
      (32*(s+2*t)*p4a(mu))/(s**3*(-1+x1))+(32*(u*x1+t*x2+s*(-1+x1+  &
      x2))*p5a(mu))/(s**3*(-1+x1)*(-1+x2)))+  &
      p3a(nu)*((-32*(2*t*u+s*(t+u))*p4a(mu))/(s**4*(-1+x1)*(-1+x2))-  &
      (16*(2*s*u*(-1+x1)+2*u*(u*(-1+x1)+t*(-1+x2))+s**2*(-1+x1+x2-  &
      x2**2))*p5a(mu))/(s**4*(-1+x1)*(-1+x2)**2))+  &
      p1a(mu)*((-64*p1a(nu))/s**2+(32*(s+2*u)*p3a(nu))/(s**3*(-1+x2))  &
      +(32*(s+2*t)*p4a(nu))/(s**3*(-1+x1))+(32*(u*x1+t*x2+s*(-1+x1+  &
      x2))*p5a(nu))/(s**3*(-1+x1)*(-1+x2)))+  &
      p4a(mu)*((-32*(s**2+2*s*t+2*t**2)*p4a(nu))/(s**4*(-1+x1)**2)+  &
      (16*(s**2*(1-x1+x1**2-x2)-2*s*t*(-1+x2)+2*t*(t+u-u*x1-  &
      t*x2))*p5a(nu))/(s**4*(-1+x1)**2*(-1+x2)))+  &
      p3a(mu)*((-32*(s**2+2*s*u+2*u**2)*p3a(nu))/(s**4*(-1+x2)**2)-  &
      (32*(2*t*u+s*(t+u))*p4a(nu))/(s**4*(-1+x1)*(-1+x2))-  &
      (16*(2*s*u*(-1+x1)+2*u*(u*(-1+x1)+t*(-1+x2))+s**2*(-1+x1+x2-  &
      x2**2))*p5a(nu))/(s**4*(-1+x1)*(-1+x2)**2))+  &
      (8*(-8*t*u*p5a(mu)*p5a(nu)+s*(2*t**2+2*u**2+  &
      2*s*u*(x1+x2)+s**2*(x1**2+x2**2)+2*t*(2*u+s*(x1+  &
      x2)))*gmunu(mu,nu)))/(s**4*(-1+x1)*(-1+x2))))/4.d0
    end do
  end do
!
  Bmunu = -Bmunu * (4*pi)**3*qq**2
!
!  print *,"Bmunu: "
!  do mu=0,3
!    print *,(Bmunu(mu,nu),nu=0,3)
!  end do
!
  print *,"smeB from Bmunu: ",Bmunu(1,1)+Bmunu(2,2)+Bmunu(3,3)-Bmunu(0,0)
!
end subroutine Bmunu_AK
!
subroutine BornSMEddim(iptrn,p,Born,BornLaurent)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , intent(out) :: Born
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
    BornLaurent
!
  real(kind(1d0)) :: y12,y13,y23,Q2
  real(kind(1d0)) , dimension(-4:2) :: Laurent
  type(mom) :: Q
!
  print *,"Using 3jet version of BornSMEddim"
!
  Born = 1d0
!
end subroutine BornSMEddim
!
subroutine BijSMEddim(iptrn,p,Bij,BijLaurent)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:) , intent(out) :: Bij
  real(kind(1d0)) , optional , dimension(:,:,-4:) , intent(out) :: &
    BijLaurent
!
  integer :: i
  real(kind(1d0)) :: Born
  real(kind(1d0)) , dimension(-4:2) :: BornLaurent
!
  interface
    subroutine BornSMEddim(iptrn,p,Born,BornLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        BornLaurent
!
    end subroutine BornSMEddim
  end interface
!
  print *,"Using 3jet version of BijSMEddim"
!
  Bij = 1d0
!
  if (present(BijLaurent)) then
    BijLaurent = 1d0
  end if
!
end subroutine BijSMEddim
!
subroutine BijklSMEddim(iptrn,p,Bijkl,BijklLaurent)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: p
  real(kind(1d0)) , dimension(:,:,:,:) , intent(out) :: Bijkl
  real(kind(1d0)) , optional , dimension(:,:,:,:,-4:) , intent(out) :: &
    BijklLaurent
!
  integer :: i,ipart,jpart,kpart
  real(kind(1d0)) :: Born
  real(kind(1d0)) , dimension(-4:2) :: BornLaurent
!
  interface
    subroutine BornSMEddim(iptrn,p,Born,BornLaurent)
    use particles
    implicit none
!
      integer , intent(in) :: iptrn
      type(particle) , dimension(:) , intent(in) :: p
      real(kind(1d0)) , intent(out) :: Born
      real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: &
        BornLaurent
!
    end subroutine BornSMEddim
  end interface
!
  print *,"Using 3jet version of BijklSMEddim"
!
  Bijkl = 1d0
!
  if (present(BijklLaurent)) then
    BijklLaurent(:,:,:,:,i) = 1d0
  end if
!
end subroutine BijklSMEddim
