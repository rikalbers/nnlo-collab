! mode is a character variable:
! mode = 'f' the full virtual part is needed,
! mode = 'i' the mu independent part is needed,
! mode = 'd' the mu dependent part is needed,
! mode = 'b' both the mu dependent and indepent part is calculated.
subroutine RVirtSME(iptrn,pin,mur,mode,RVirt,RVirt_dep,RVirt_indep, &
                    RVLaurent)
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
  real(kind(1d0)) , intent(out) :: RVirt
  real(kind(1d0)) , intent(out) :: RVirt_dep,RVirt_indep
  real(kind(1d0)) , optional , dimension(-4:2) , intent(out) :: RVLaurent
!
  integer :: ipart,qtype,rtype
  logical , save :: init = .true.
  real(kind(1d0)) :: Born, smeR
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
! Renormalization scale:
  comu = mur
  Q2   = 2d0*pin(1)%p*pin(2)%p
!
  RVirt       = 0d0
  RVirt_indep = 0d0
  RVirt_dep   = 0d0
  if (present(RVLaurent)) RVLaurent = 0d0
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
! The Real part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      smeR = PSI2d1gBorn(1,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
      RVirt = PSI2d1gVirtNLO(1,3,2,8,7)
      RVirt = RVirt - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*smeR
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2d1gVirtNLOmuindep(1,3,2,8,7)
! We sweep the renormalization into the independent part:
      RVirt_indep = RVirt_indep - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*smeR
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2d1gVirtNLOmudep(1,3,2,8,7)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the real-virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      RVLaurent(-2) = PSI2d1gVirtNLOem2(1,3,2,8,7)
      RVLaurent(-1) = PSI2d1gVirtNLOem1(1,3,2,8,7) &
                      - qcd_beta0*smeR
    end if
    return
! e+ e- -> u u~ g
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g"
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      smeR = PSI2u1gBorn(1,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
      RVirt = PSI2u1gVirtNLO(1,3,2,8,7)
      RVirt = RVirt - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*smeR
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2u1gVirtNLOmuindep(1,3,2,8,7)
! We sweep the renormalization into the independent part:
      RVirt_indep = RVirt_indep - (2d0*(qcd_cf)-0d0*4d0*qcd_ca/6d0)*smeR
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2u1gVirtNLOmudep(1,3,2,8,7)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the real-virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      RVLaurent(-2) = PSI2u1gVirtNLOem2(1,3,2,8,7)
      RVLaurent(-1) = PSI2u1gVirtNLOem1(1,3,2,8,7) &
                      - qcd_beta0*smeR
    end if
!    do ipart=1,5
!      print *,pin(ipart)%p
!    end do
!    print *,"rvirt: ",rvirt
!    print *,"rvirtem2: ",PSI2u1gVirtNLOem2(1,3,2,8,7)
!    print *,"rvirtem1: ",PSI2u1gVirtNLOem1(1,3,2,8,7) &
!                      - qcd_beta0*Born
!    print *,"real: ",smeR
!    print *,"Q: ",sqrt(Q2)
!    print *,"mur: ",mur
!    stop "RVirtSME..."
    return
  end if
!
end subroutine RVirtSME
