subroutine RRealSME(iptrn,pin,smeRR)
use process
use particles
use QCDparams
use my_model
implicit none
!
  integer , intent(in) :: iptrn
  type(particle) , dimension(:) , intent(in) :: pin
  real(kind(1d0)) , intent(out) :: smeRR
!
  integer :: ipart
  logical , save :: init = .true.
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  real(kind(1d0)) Xnx,Xny,Xnz,XC3
  complex(kind(1d0)) , dimension(2,2,2) :: C
  common/COUPLINGS/ C &
        /QCD/XNf,XNu,XNd,XNc,XNa
!        /QCDTWO/Xnx,Xny,Xnz,XC3
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  common/dotproducts/S &
        /spinorproducts/A,B
!
  real(kind(1d0)) , external :: PSI2d2gBorn,PSI2u2gBorn, &
                                PSI2u2d,PSI2u2c,PSI2d2s, &
                                PSI4u,PSI4usl,PSI4d,PSI4dsl

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
  do ipart=1,nleg_born+2
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
! To get agreement with HELAC we had to change the ordering from 1,2,3,4,7,8
! to 1,2,3,4,8,7.
! The ordering among momenta: q,g,g,qb,e-,e+
! The ordering among momenta: q,Qb,Q,qb,e-,e+
! The variable called iptrn determines which contribution we should 
! calculate:
! e+ e- -> d d~ g g
  if (iptrn.eq.1) then
!    print *,"e+ e- -> d d~ g g"
    smeRR = PSI2d2gBorn(1,3,4,2,8,7)
    return
! e+ e- -> u u~ g g
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g g"
    smeRR = PSI2u2gBorn(1,3,4,2,8,7)
    return
! e+ e- -> u u~ d d~
  elseif (iptrn.eq.3) then
!    print *,"e+ e- -> u u~ d d~"
    smeRR = PSI2u2d(1,4,3,2,8,7)
    return
! e+ e- -> u u~ u u~
  elseif (iptrn.eq.4) then
!    print *,"e+ e- -> u u~ u u~"
    smeRR = PSI4u(1,4,3,2,8,7) + PSI4usl(1,4,3,2,8,7)
    return
! e+ e- -> d d~ d d~
  elseif (iptrn.eq.5) then
!    print *,"e+ e- -> d d~ d d~"
    smeRR = PSI4d(1,4,3,2,8,7) + PSI4dsl(1,4,3,2,8,7)
    return
! e+ e- -> u u~ c c~
  elseif (iptrn.eq.6) then
!    print *,"e+ e- -> u u~ c c~"
    smeRR = PSI2u2c(1,4,3,2,8,7)
    return
! e+ e- -> d d~ s s~
  elseif (iptrn.eq.7) then
!    print *,"e+ e- -> d d~ s s~"
    smeRR = PSI2d2s(1,4,3,2,8,7)
    return
  else 
    print *,"unknown real SME is asked for..."
    print *,"iptrn: ",iptrn
    stop
  end if
!
end subroutine RRealSME
