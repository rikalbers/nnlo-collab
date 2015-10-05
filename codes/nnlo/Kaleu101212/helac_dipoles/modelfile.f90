  subroutine helac_kaleu_model &
             ( model,vertx ,hel2kal ,parmas,parwid ,ihiggs,oqcd,wqcd )
!********************************************************************
!* Prepared for ILC-physics with Helac.
!* All leptons massless, so need only one pair.
!* All quarks except the top-quark and possibly the botton-quark
!* massless, so need two pairs.
!********************************************************************
  use avh_kaleu_model
  implicit none
  type(model_type) ,intent(out) :: model
  type(vertx_type) ,intent(out) :: vertx
  integer :: gluon,photon,wboson,zboson,higgs
  integer :: uquark,dquark,cquark,squark,tquark,bquark
  integer :: electr,eleneu,muon,muneu,tauon,tauneu
  integer :: npartic
  logical :: withhiggs,onlyqcd,withqcd
!
  data   onlyqcd/.false./
  data   withqcd/.true./
  data withhiggs/.false./
!
  parameter(  gluon=1  )
  parameter( photon=2  )
  parameter( wboson=3  )
  parameter( zboson=4  )
  parameter(  higgs=5  )
  parameter( uquark=6  )
  parameter( dquark=7  )
  parameter( cquark=8  )
  parameter( squark=9  )
  parameter( tquark=10 )
  parameter( bquark=11 )
  parameter( eleneu=12 )
  parameter( electr=13 )
  parameter(  muneu=14 )
  parameter(   muon=15 )
  parameter( tauneu=16 )
  parameter(  tauon=17 )
!
  parameter( npartic=17 )
!
  character(2) symbol(npartic)
  data symbol/'g ','A ','W ','Z ','H ','u ','d ','c ','s ','t ','b ' ,&
              'En','E ','Mn','M ','Tn','T '/
!
! Start Helac-specific stuff ******************************************
!
  integer         ,intent(out) :: hel2kal(-12:41)
  integer         ,intent(in)  :: ihiggs
  logical         ,intent(in)  :: oqcd,wqcd
  real(kind(1d0)) ,intent(in)  :: parmas(-12:41),parwid(-12:41)
  integer :: k2h(npartic),ii
!
    onlyqcd = oqcd
    withqcd = wqcd
  withhiggs = (ihiggs.eq.1)
!
  hel2kal( 1) = eleneu
  hel2kal( 2) = electr
  hel2kal( 3) = uquark
  hel2kal( 4) = dquark
  hel2kal( 5) = muneu
  hel2kal( 6) = muon
  hel2kal( 7) = cquark
  hel2kal( 8) = squark
  hel2kal( 9) = tauneu
  hel2kal(10) = tauon
  hel2kal(11) = tquark
  hel2kal(12) = bquark
  hel2kal(31) = photon
  hel2kal(32) = zboson
  hel2kal(33) = wboson
  hel2kal(34) = wboson
  hel2kal(35) = gluon
  hel2kal(41) = higgs
  do ii=1,12
    hel2kal(-ii) = hel2kal(ii)
  enddo
  k2h(  gluon ) = 35
  k2h( photon ) = 31
  k2h( wboson ) = 33
  k2h( zboson ) = 32
  k2h(  higgs ) = 41
  k2h( uquark ) = 3
  k2h( dquark ) = 4
  k2h( cquark ) = 7
  k2h( squark ) = 8
  k2h( tquark ) = 11
  k2h( bquark ) = 12
  k2h( eleneu ) = 1
  k2h( electr ) = 2
  k2h( muneu  ) = 5
  k2h( muon   ) = 6
  k2h( tauneu ) = 9
  k2h( tauon  ) = 10
!
! Put masses using Helac-specific arrays parmas and parwid
  do ii=1,npartic
    call addparticle &
            ( model ,ii ,symbol(ii) ,parmas(k2h(ii)) ,parwid(k2h(ii)) )
  enddo
!
! End Helac-specific stuff ********************************************
!
  if (onlyqcd.or.withqcd) then
    call addvertex( vertx , gluon,gluon,gluon )
    call addvertex( vertx ,uquark,uquark,gluon )
    call addvertex( vertx ,dquark,dquark,gluon )
    call addvertex( vertx ,cquark,cquark,gluon )
    call addvertex( vertx ,squark,squark,gluon )
    call addvertex( vertx ,tquark,tquark,gluon )
    call addvertex( vertx ,bquark,bquark,gluon )
  endif
!
  if (.not.onlyqcd) then
    call addvertex( vertx ,uquark,uquark,photon )
    call addvertex( vertx ,dquark,dquark,photon )
    call addvertex( vertx ,cquark,cquark,photon )
    call addvertex( vertx ,squark,squark,photon )
    call addvertex( vertx ,tquark,tquark,photon )
    call addvertex( vertx ,bquark,bquark,photon )
!  
    call addvertex( vertx ,uquark,uquark,zboson )
    call addvertex( vertx ,dquark,dquark,zboson )
    call addvertex( vertx ,cquark,cquark,zboson )
    call addvertex( vertx ,squark,squark,zboson )
    call addvertex( vertx ,tquark,tquark,zboson )
    call addvertex( vertx ,bquark,bquark,zboson )
!  
    call addvertex( vertx ,uquark,dquark,wboson )
    call addvertex( vertx ,cquark,squark,wboson )
    call addvertex( vertx ,tquark,bquark,wboson )
!  
    call addvertex( vertx ,wboson,wboson,zboson )
    call addvertex( vertx ,wboson,wboson,photon )
!  
    call addvertex( vertx ,electr,electr,photon )
    call addvertex( vertx ,electr,electr,zboson )
    call addvertex( vertx ,eleneu,eleneu,zboson )
    call addvertex( vertx ,electr,eleneu,wboson )
!  
    call addvertex( vertx ,  muon,  muon,photon )
    call addvertex( vertx ,  muon,  muon,zboson )
    call addvertex( vertx , muneu, muneu,zboson )
    call addvertex( vertx ,  muon, muneu,wboson )
!  
    call addvertex( vertx , tauon, tauon,photon )
    call addvertex( vertx , tauon, tauon,zboson )
    call addvertex( vertx ,tauneu,tauneu,zboson )
    call addvertex( vertx , tauon,tauneu,wboson )
!  
    if (withhiggs) then
      call addvertex( vertx , higgs, higgs,higgs )
      call addvertex( vertx ,zboson,zboson,higgs )
      call addvertex( vertx ,wboson,wboson,higgs )
      call addvertex( vertx ,tquark,tquark,higgs )
      call addvertex( vertx ,bquark,bquark,higgs )
      call addvertex( vertx ,cquark,cquark,higgs )
      call addvertex( vertx ,squark,squark,higgs )
      call addvertex( vertx , tauon, tauon,higgs )
      call addvertex( vertx ,  muon,  muon,higgs )
!      call addvertex( vertx ,dquark,dquark,higgs )
!      call addvertex( vertx ,uquark,uquark,higgs )
    endif
  endif
!
  end
