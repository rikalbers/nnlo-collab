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
  real(kind(1d0)) :: Born
!
! In what follows is ugly as hell, but we still have some
! f77 legacy code around...
  real(kind(1d0)) , dimension(4,9) :: P
  real(kind(1d0)) XNf,XNu,XNd,XNc,XNa
  complex(kind(1d0)) , dimension(2,2,2) :: C
  complex(kind(1d0)) , dimension(2) :: Cax
  complex(kind(1d0)) , dimension(2) :: Cfsum
  common/COUPLINGS/ C &
        /AXIALCOUPLINGS/ Cax &
        /QCD/XNf,XNu,XNd,XNc,XNa &
        /FSUMCOUPLING/Cfsum
  real(kind(1d0)) , dimension(11,11) :: S
  complex(kind(1d0)) , dimension(11,11) :: A,B
  real(kind(1d0)) :: COMU
  real(kind(1d0)) :: Q2
  real(kind(1d0)) :: xmt,mtsq
  common/dotproducts/S &
        /spinorproducts/A,B &
        /scales/ COMU &
        /qsquared/ Q2 &
        /topmass/xmt,mtsq
!
  real(kind(1d0)) :: beta0
!
  real(kind(1d0)) , external :: PSI2q2gvirtNLOe0, &
                                PSI2q2gvirtNLOem2, &
                                PSI2q2gvirtNLOem1, &
                                PSI2d2gBorn,PSI2u2gBorn, &
                                PSI2u2d,PSI2u2c,PSI2d2s, &
                                PSI4uBorn,PSI4dBorn, &
                                PSI2q2gvirtNLOmudep, &
                                PSI2q2gvirtNLOmuindep, &
                                PSI2q2QvirtNLOem2, &
                                PSI2q2QvirtNLOem1, &
                                PSI2q2QvirtNLO, &
                                PSI2q2QvirtNLOmudep, &
                                PSI2q2QvirtNLOmuindep, &
                                PSI4qvirtNLOem2, &
                                PSI4qvirtNLOem1, &
                                PSI4qvirtNLO, &
                                PSI4qvirtNLOmudep, &
                                PSI4qvirtNLOmuindep
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
    Cfsum(1:2) = 0.5d0*(XNu*(C(1:2,1,1)+C(1:2,2,1)) &
               +        XNd*(C(1:2,1,2)+C(1:2,2,2)))
!
    xmt = phys_Tmass
    mtsq = xmt**2
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
! To get agreement with HELAC we had to change the ordering from 1,2,3,4,7,8
! to 1,3,4,2,8,7.
! The ordering among momenta: q,g,g,qb,e+,e-
! The variable called iptrn determines which contribution we should 
! calculate:
!********************
! e+ e- -> d d~ g g *
!********************
  if (iptrn.eq.1) then
!    print *,"e+ e- -> d d~ g g"
    qtype = 2
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born  = PSI2d2gBorn(1,3,4,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI2q2gVirtNLOe0(1,3,4,2,8,7,qtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born - (qcd_nc**2 - 1d0)/qcd_nc*born
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt - (qcd_nc**2 - 1d0)/qcd_nc*born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2q2gVirtNLOmuindep(1,3,4,2,8,7,qtype)
! We sweep the renormalization into the independent part:
      RVirt_indep = RVirt_indep - (qcd_nc**2 - 1d0)/qcd_nc*born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2q2gVirtNLOmudep(1,3,4,2,8,7,qtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI2q2gVirtNLOem2(1,3,4,2,8,7,qtype)
      RVLaurent(-1) = PSI2q2gVirtNLOem1(1,3,4,2,8,7,qtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*born
    end if
    return
!********************
! e+ e- -> u u~ g g *
!********************
  elseif (iptrn.eq.2) then
!    print *,"e+ e- -> u u~ g g"
    qtype = 1
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born  = PSI2u2gBorn(1,3,4,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI2q2gVirtNLOe0(1,3,4,2,8,7,qtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born - (qcd_nc**2 - 1d0)/qcd_nc*born
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt - (qcd_nc**2 - 1d0)/qcd_nc*born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2q2gVirtNLOmuindep(1,3,4,2,8,7,qtype)
! We sweep the renormalization into the independent part:
      RVirt_indep = RVirt_indep - (qcd_nc**2 - 1d0)/qcd_nc*born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2q2gVirtNLOmudep(1,3,4,2,8,7,qtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI2q2gVirtNLOem2(1,3,4,2,8,7,qtype)
      RVLaurent(-1) = PSI2q2gVirtNLOem1(1,3,4,2,8,7,qtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*born
    end if
    return
!*********************
! e+ e- -> u u~ d d~ *
!*********************
  elseif (iptrn.eq.3) then
!    print *,"e+ e- -> u u~ d d~"
    qtype = 1
    rtype = 2
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born  = PSI2u2d(1,4,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI2q2QvirtNLO(1,4,3,2,8,7,qtype,rtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 4d0*qcd_nc/3d0*Born - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born 
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt + 2d0*qcd_nc/3d0*Born &
            - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2q2QVirtNLOmuindep(1,4,3,2,8,7,qtype,rtype)
! We sweep the renormalization and scheme convertion into the 
! independent part:
      RVirt_indep = RVirt_indep + 2d0*qcd_nc/3d0*Born &
                  - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2q2QVirtNLOmudep(1,4,3,2,8,7,qtype,rtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI2q2QVirtNLOem2(1,4,3,2,8,7,qtype,rtype)
      RVLaurent(-1) = PSI2q2QVirtNLOem1(1,4,3,2,8,7,qtype,rtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*Born
    end if
    return
!*********************
! e+ e- -> u u~ u u~ *
!*********************
  elseif (iptrn.eq.4) then
!    print *,"e+ e- -> u u~ u u~"
    qtype = 1
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born = PSI4uBorn(1,4,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI4qVirtNLO(1,4,3,2,8,7,qtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 4d0*qcd_nc/3d0*Born &
!             - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born 
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt + 2d0*qcd_nc/3d0*Born &
            - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI4qVirtNLOmuindep(1,4,3,2,8,7,qtype)
! We sweep the renormalization and scheme convertion into the 
! independent part:
      RVirt_indep = RVirt_indep + 2d0*qcd_nc/3d0*Born &
                  - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI4qVirtNLOmudep(1,4,3,2,8,7,qtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI4qVirtNLOem2(1,4,3,2,8,7,qtype)
      RVLaurent(-1) = PSI4qVirtNLOem1(1,4,3,2,8,7,qtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*Born
    end if
    return
!*********************
! e+ e- -> d d~ d d~ *
!*********************
  elseif (iptrn.eq.5) then
!    print *,"e+ e- -> d d~ d d~"
    qtype = 2
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born = PSI4dBorn(1,4,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI4qVirtNLO(1,4,3,2,8,7,qtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 4d0*qcd_nc/3d0*Born &
!             - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born 
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt + 2d0*qcd_nc/3d0*Born &
            - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI4qVirtNLOmuindep(1,4,3,2,8,7,qtype)
! We sweep the renormalization and scheme convertion into the 
! independent part:
      RVirt_indep = RVirt_indep + 2d0*qcd_nc/3d0*Born &
                  - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI4qVirtNLOmudep(1,4,3,2,8,7,qtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI4qVirtNLOem2(1,4,3,2,8,7,qtype)
      RVLaurent(-1) = PSI4qVirtNLOem1(1,4,3,2,8,7,qtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*Born
    end if
    return
!*********************
! e+ e- -> u u~ c c~ *
!*********************
  elseif (iptrn.eq.6) then
!    print *,"e+ e- -> u u~ c c~"
    qtype = 1
    rtype = 1
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born  = PSI2u2c(1,4,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI2q2QvirtNLO(1,4,3,2,8,7,qtype,rtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 4d0*qcd_nc/3d0*Born - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born 
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt + 2d0*qcd_nc/3d0*Born &
            - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2q2QVirtNLOmuindep(1,4,3,2,8,7,qtype,rtype)
! We sweep the renormalization and scheme convertion into the 
! independent part:
      RVirt_indep = RVirt_indep + 2d0*qcd_nc/3d0*Born &
                  - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2q2QVirtNLOmudep(1,4,3,2,8,7,qtype,rtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI2q2QVirtNLOem2(1,4,3,2,8,7,qtype,rtype)
      RVLaurent(-1) = PSI2q2QVirtNLOem1(1,4,3,2,8,7,qtype,rtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*Born
    end if
    return
!*********************
! e+ e- -> d d~ s s~ *
!*********************
  elseif (iptrn.eq.7) then
!    print *,"e+ e- -> d d~ s s~"
    qtype = 2
    rtype = 2
! The Born part is needed because of renormalization, it is always 
! calculated but when the dependent part is needed:
    if ((mode.eq.'f').or.(mode.eq.'i').or.(mode.eq.'b')) then
      Born  = PSI2d2s(1,4,3,2,8,7)
    end if
! The full real-virtual part is only calculated if it is explicitly needed
    if (mode.eq.'f') then
! First we obtain the unrenormalized result in FDH:
      RVirt = PSI2q2QvirtNLO(1,4,3,2,8,7,qtype,rtype)
! The result in FDH, but renormalized:
!      RVirt = RVirt + 2d0*qcd_nc/3d0*born
! This is the result in tHV but with an extra renormalization term, only
! needed if a comparison is meant with HELAC1LOOP, in an actual calculation
! this _SHOULD_ not be used!!!!!!!!!!!!!!!
!      RVirt = RVirt + 4d0*qcd_nc/3d0*Born - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born 
! This is the correct formula to convert the FDH, unrenormalized result
! to the renormalized one in CDR:
      RVirt = RVirt + 2d0*qcd_nc/3d0*Born &
            - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
! The mu independent part is calculated if it is asked for or both
! contributions are needed:
    if ((mode.eq.'i').or.(mode.eq.'b')) then
      RVirt_indep = PSI2q2QVirtNLOmuindep(1,4,3,2,8,7,qtype,rtype)
! We sweep the renormalization and scheme convertion into the 
! independent part:
      RVirt_indep = RVirt_indep + 2d0*qcd_nc/3d0*Born &
                  - 2d0*(qcd_nc**2 - 1d0)/qcd_nc*Born
    end if
    if ((mode.eq.'d').or.(mode.eq.'b')) then
      RVirt_dep = PSI2q2QVirtNLOmudep(1,4,3,2,8,7,qtype,rtype)
    end if
! We construct the total real-virtual part out of the two contributions:
    if (mode.eq.'b') RVirt = RVirt_dep + RVirt_indep
! If the RVLaurent variable is present in the argument we give back
! the Laurent series of the virtual part not just the finite piece:
    if (present(RVLaurent)) then
      RVLaurent = 0d0
      beta0 = (11d0*qcd_nc - 2d0*qcd_nf)/3d0
      RVLaurent(-2) = PSI2q2QVirtNLOem2(1,4,3,2,8,7,qtype,rtype)
      RVLaurent(-1) = PSI2q2QVirtNLOem1(1,4,3,2,8,7,qtype,rtype)
      RVLaurent(-1) = RVLaurent(-1) - 2d0*beta0*Born
    end if
    return
  else
    print *,"Wrong pattern is given to RVirtSME..."
    print *,"iptrn: ",iptrn
    stop
  end if
!
end subroutine RVirtSME
