! This source contains routines used by the Rambo PS generator, but
! process dependent:
subroutine Rambo_initmyprocess
use RamboPS
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
! The minimal value for the invariants is set here:
  ymin = 1d-8
! No PDF for this process:
  PDF_option = 0
!
  masses = 0d0
!
! We initialize all the particles in the model:
  labels( gluon)='g '; masses( gluon)=0d0
  labels(photon)='A '; masses(photon)=0d0
  labels(wboson)='W '; masses(wboson)=phys_Wmass
  labels(zboson)='Z '; masses(zboson)=phys_Zmass
  labels( higgs)='H '; masses( higgs)=phys_Hmass
  labels(uquark)='u '; masses(uquark)=0d0
  labels(dquark)='d '; masses(dquark)=0d0
  labels(cquark)='c '; masses(cquark)=0d0
  labels(squark)='s '; masses(squark)=0d0
  labels(tquark)='t '; masses(tquark)=phys_Tmass
  labels(bquark)='b'; masses(bquark)=0d0
  labels(eleneu)='ne'; masses(eleneu)=0d0
  labels(electr)='el'; masses(electr)=0d0
  labels(muneu )='nm'; masses(muneu )=0d0
  labels(muon  )='mu'; masses(muon  )=0d0
  labels(tauneu)='nt'; masses(tauneu)=0d0
  labels(tauon )='ta'; masses(tauon )=0d0
!
end subroutine Rambo_initmyprocess
