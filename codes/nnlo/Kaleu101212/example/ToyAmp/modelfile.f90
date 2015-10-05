  subroutine simple_model( model,vertx ,masses,widths,labels    &
                                       ,onlyqcd,withqcd,withhiggs )
!********************************************************************
!* The Standard Model.
!********************************************************************
  use particles ! This is not part of ToyAmp, should come from outside
  use avh_toyamp_model
  implicit none
  type(model_type) ,intent(out) :: model
  type(vertx_type) ,intent(out) :: vertx
!
  real(kind(1d0)) ,intent(in) :: masses(nfirst:nlast),widths(nfirst:nlast)
  character(20)   ,intent(in) :: labels(nfirst:nlast)
  logical         ,intent(in) :: onlyqcd,withqcd,withhiggs
!
  integer :: ii
!
  do ii=nfirst,nlast
    call addparticle( model ,ii ,labels(ii) ,masses(ii) ,widths(ii) )
  enddo
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
!      call addvertex( vertx ,electr,electr,higgs )
    endif
  endif
!
  end
