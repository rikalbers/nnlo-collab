! This source contains routines used by the Kaleu PS generator, but
! process dependent:
subroutine Kaleu_initmyprocess
use KaleuPS
use kaleu_particles
use my_model
implicit none
!
!
  integer :: istat
!
! We initialize particles:
  call init_particles
! We allocate the arrays related:
  allocate(labels(nfirst:nlast),stat=istat)
  if (istat.ne.0) then
    print *,"Problem ocurred during allocation of labels..."
    stop
  end if
  allocate(masses(nfirst:nlast),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocation of masses..."
    stop
  end if
  allocate(widths(nfirst:nlast),stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocation of widths..."
    stop
  end if
! The minimal value for all the invariants is set to one:
  smin = 1d-4
! No PDF for this process:
  PDF_option = 0
!
  onlyqcd   = .false.
  withqcd   = .true.
  withhiggs = .false.
!
  masses = 0d0
  widths = 0d0
!
! We initialize all the particles in the model:
  labels( gluon)='g '; masses( gluon)=0d0        ; widths( gluon)=0d0
  labels(photon)='A '; masses(photon)=0d0        ; widths(photon)=0d0
  labels(wboson)='W '; masses(wboson)=phys_Wmass ; widths(wboson)=phys_Wwidth
  labels(zboson)='Z '; masses(zboson)=phys_Zmass ; widths(zboson)=phys_Zwidth
  labels( higgs)='H '; masses( higgs)=phys_Hmass ; widths( higgs)=phys_Hwidth
  labels(uquark)='u '; masses(uquark)=0d0        ; widths(uquark)=0d0      
  labels(dquark)='d '; masses(dquark)=0d0        ; widths(dquark)=0d0      
  labels(cquark)='c '; masses(cquark)=0d0        ; widths(cquark)=0d0      
  labels(squark)='s '; masses(squark)=0d0        ; widths(squark)=0d0      
  labels(tquark)='t '; masses(tquark)=phys_Tmass ; widths(tquark)=phys_Twidth
  labels(bquark)='b'; masses(bquark)=0d0         ; widths(bquark)=0d0      
  labels(eleneu)='ne'; masses(eleneu)=0d0        ; widths(eleneu)=0d0      
  labels(electr)='el'; masses(electr)=0d0        ; widths(electr)=0d0      
  labels(muneu )='nm'; masses(muneu )=0d0        ; widths(muneu )=0d0      
  labels(muon  )='mu'; masses(muon  )=0d0        ; widths(muon  )=0d0      
  labels(tauneu)='nt'; masses(tauneu)=0d0        ; widths(tauneu)=0d0      
  labels(tauon )='ta'; masses(tauon )=0d0        ; widths(tauon )=0d0      
!
end subroutine Kaleu_initmyprocess
