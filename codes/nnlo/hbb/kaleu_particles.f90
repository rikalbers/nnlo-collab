module kaleu_particles
!
  integer ,parameter :: nfirst=1,nlast=17 &
   ,eleneu=1 ,electr=2 ,uquark=3 ,dquark=4 &
   , muneu=5 ,  muon=6 ,cquark=7 ,squark=8 &
   ,tauneu=9 , tauon=10,tquark=11,bquark=12 &
   ,wboson=13,photon=14,zboson=15, higgs=16 ,gluon=17
!
! logical ,protected :: is_vector(nfirst:nlast)= .false. ! with protected
! logical ,protected :: is_spinor(nfirst:nlast)= .false. ! with protected
! logical ,protected :: is_higgs(nfirst:nlast)= .false.  ! with protected
! logical ,protected :: is_qcd(nfirst:nlast)= .false.    ! with protected
! logical ,protected :: is_quark(nfirst:nlast)= .false.  ! with protected
! logical ,protected :: is_gluon(nfirst:nlast)= .false.  ! with protected
! logical ,protected :: is_chlept(nfirst:nlast)= .false. ! with protected
! logical ,protected :: is_neutri(nfirst:nlast)= .false. ! with protected
! logical ,protected :: is_photon(nfirst:nlast)= .false. ! with protected
!
  logical :: is_vector(nfirst:nlast)= .false. ! without protected
  logical :: is_spinor(nfirst:nlast)= .false. ! without protected
  logical :: is_higgs(nfirst:nlast)= .false.  ! without protected
  logical :: is_qcd(nfirst:nlast)= .false.    ! without protected
  logical :: is_quark(nfirst:nlast)= .false.  ! without protected
  logical :: is_gluon(nfirst:nlast)= .false.  ! without protected
  logical :: is_chlept(nfirst:nlast)= .false. ! without protected
  logical :: is_neutri(nfirst:nlast)= .false. ! without protected
  logical :: is_photon(nfirst:nlast)= .false. ! without protected
!
contains
!
  subroutine init_particles
  implicit none
  is_neutri(eleneu)=.true.; is_spinor(eleneu)=.true. 
  is_neutri( muneu)=.true.; is_spinor( muneu)=.true.
  is_neutri(tauneu)=.true.; is_spinor(tauneu)=.true.
  is_chlept(electr)=.true.; is_spinor(electr)=.true.
  is_chlept(  muon)=.true.; is_spinor(  muon)=.true.
  is_chlept( tauon)=.true.; is_spinor( tauon)=.true.
   is_quark(uquark)=.true.; is_spinor(uquark)=.true.; is_qcd(uquark)=.true.
   is_quark(dquark)=.true.; is_spinor(dquark)=.true.; is_qcd(dquark)=.true.
   is_quark(cquark)=.true.; is_spinor(cquark)=.true.; is_qcd(cquark)=.true.
   is_quark(squark)=.true.; is_spinor(squark)=.true.; is_qcd(squark)=.true.
   is_quark(tquark)=.true.; is_spinor(tquark)=.true.; is_qcd(tquark)=.true.
   is_quark(bquark)=.true.; is_spinor(bquark)=.true.; is_qcd(bquark)=.true.
   is_gluon( gluon)=.true.; is_vector( gluon)=.true.; is_qcd( gluon)=.true.
  is_photon(photon)=.true.; is_vector(photon)=.true.;
                            is_vector(wboson)=.true.;
                            is_vector(zboson)=.true.;
   is_higgs( higgs)=.true.
  end subroutine
!
end module
