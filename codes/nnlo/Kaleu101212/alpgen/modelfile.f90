  subroutine user_model( model,vertx ,apar ,oqcd,wqcd,whiggs )
!********************************************************************
!* The Standard Model.
!* Realize that you might need to include less fermions if you make
!* many of them massless. See the lhc example.
!********************************************************************
  use avh_kaleu_model
  implicit none
  type(model_type) ,intent(out) :: model
  type(vertx_type) ,intent(out) :: vertx
  real(kind(1d0))  ,intent(in)  :: apar(100)
  logical          ,intent(in)  :: oqcd,wqcd,whiggs
  integer :: gluon,photon,wboson,zboson,higgs
  integer :: uquark,dquark,cquark,squark,tquark,bquark
  integer :: electr,eleneu,muon,muneu,tauon,tauneu
!
  logical :: onlyqcd = .false.
  logical :: withqcd = .true.
  logical :: withhiggs = .false.
!
  parameter(  gluon=1  )
  parameter( photon=2  )
  parameter( wboson=3  )
  parameter( zboson=4  )
  parameter(  higgs=5  )
  parameter( uquark=6  )
  parameter( dquark=7  )
  parameter( cquark=8 )
  parameter( squark=9 )
  parameter( tquark=10 )
  parameter( bquark=11 )
  parameter( eleneu=12 )
  parameter( electr=13 )
  parameter( muneu =14 )
  parameter( muon  =15 )
  parameter( tauneu=16 )
  parameter( tauon =17 )
!
  onlyqcd   = oqcd
  withqcd   = wqcd
  withhiggs = whiggs
!
  call addparticle( model , gluon ,'g ' ,0d0     ,0d0 )
  call addparticle( model ,photon ,'A ' ,0d0     ,0d0 )
  call addparticle( model ,wboson ,'W ' ,apar( 3),apar(4) )
  call addparticle( model ,zboson ,'Z ' ,apar( 1),apar(2) )
  call addparticle( model , higgs ,'H ' ,apar( 5),apar(6) )
  call addparticle( model ,uquark ,'u ' ,apar(12),0d0 )
  call addparticle( model ,dquark ,'d ' ,apar(11),0d0 )
  call addparticle( model ,cquark ,'c ' ,apar(14),0d0 )
  call addparticle( model ,squark ,'s ' ,apar(13),0d0 )
  call addparticle( model ,tquark ,'t ' ,apar(16),0d0 )
  call addparticle( model ,bquark ,'b ' ,apar(15),0d0 )
  call addparticle( model ,eleneu ,'En' ,apar(22),0d0 )
  call addparticle( model ,electr ,'E ' ,apar(21),0d0 )
  call addparticle( model ,muneu  ,'Mn' ,apar(24),0d0 )
  call addparticle( model ,muon   ,'M ' ,apar(23),0d0 )
  call addparticle( model ,tauneu ,'Tn' ,apar(26),0d0 )
  call addparticle( model ,tauon  ,'T ' ,apar(25),0d0 )
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
      call addvertex( vertx , tauon, tauon,higgs )
!      call addvertex( vertx ,squark,squark,higgs )
!      call addvertex( vertx ,  muon,  muon,higgs )
!      call addvertex( vertx ,dquark,dquark,higgs )
!      call addvertex( vertx ,uquark,uquark,higgs )
!      call addvertex( vertx ,electr,electr,higgs )
    endif
  endif
!
  end
