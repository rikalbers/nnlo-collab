module QQ_kaleu
  use avh_kaleu
  use avh_kaleu_model
  use avh_kaleu_strf
  logical          ,save :: yes = .false.
  integer          ,save :: nfin
  type(model_type) ,save :: mdl     ! Masses and widths
  type(kaleu_type) ,save :: obj(16) ! Instances of the phase space generator
  type(strf_type)  ,save :: str(16) ! Concerns the generation of x1,x2
end module
  

  subroutine alpgen_kaleu_init( ikaleu )
  use QQ_kaleu
  yes = (ikaleu.eq.1)
  end subroutine

  function alpgen_kaleu_yes() result(value)
  use QQ_kaleu
  logical value
  value = yes
  end function
    

  subroutine QQ_kaleu_init( ihvy, njets ,roots ,apar &
                         ,ptjmin,ptQmin,etajmax,etaQmax,drjmin,drQmin &
                         ,nopt,niter  ) 
!*********************************************************************
!*********************************************************************
  use QQ_kaleu
  implicit none
  integer         ,intent(in) :: ihvy,njets,niter
  real(kind(1d0)) ,intent(in) :: roots,apar(100),nopt &
                         ,ptjmin,ptQmin,etajmax,etaQmax,drjmin,drQmin
  type(vertx_type) :: vtx
  real(kind(1d0))  :: smin(-2:17,-2:17)
  integer :: process(-2:17),ii,nproc,ihq
!
  call kaleu_hello
!
!                                ,onlyqcd,withqcd,withhiggs
  call user_model( mdl,vtx ,apar ,.false.,.true. ,.true.   )
!
  if     (ihvy.eq.4) then
    ihq = 8
  elseif (ihvy.eq.5) then
    ihq = 11
  elseif (ihvy.eq.6) then
    ihq = 10
  else
    write(*,*) 'ERROR in QQ_kaleu_init: heavy quark',ihvy,' not defined'
    stop
  endif
  nfin = njets+2
!
  smin(-1,1     ) = -2d0*ptQmin*ptQmin*dexp(-etaQmax)
  smin(-1,2     ) = smin(-1,1)
  smin(-1,3:nfin) = -2d0*ptjmin*ptjmin*dexp(-etajmax)
  smin(-2,1:nfin) = smin(-1,1:nfin)
  smin(1:nfin,-2:-1) = smin(-2:-1,1:nfin)
!
  smin(1     ,2     ) = 2d0*ptQmin*ptQmin*(1d0-dcos(drQmin))
  smin(1     ,3:nfin) = 2d0*ptQmin*ptjmin*(1d0-dcos(drjmin)) ! not sure
  smin(2     ,3:nfin) = smin(1,3:nfin)
  smin(3:nfin,1:2   ) = smin(1:2,3:nfin)
  smin(3:nfin,3:nfin) = 2d0*ptjmin*ptjmin*(1d0-dcos(drjmin))
!
  process(1:2   ) = ihq
  process(3:nfin) = 1
!---Q Qbar final state (+ up to 4 jets)
  nproc = 2
!
!  1  g g    -> Q Qbar (+ gluons)
  process(-2) = 1
  process(-1) = 1
  call kaleu_put_process( mdl,vtx,obj(1) ,process ,nfin ,roots )
!  2  q qbar -> Q Qbar (+ gluons)   q=u,d,c,s
!     qbar q -> Q Qbar (+ gluons)   q=u,d,c,s
  process(-2) = 6
  process(-1) = 6
  call kaleu_put_process( mdl,vtx,obj(2) ,process ,nfin ,roots )
!
!---Q Qbar + 1 q jet final state
  if (njets.ge.1) then
  nproc = 4
!
!  3  g q    -> Q Qbar q    (+ gluons)
!     g qbar -> Q Qbar qbar (+ gluons)
  process(-2) = 6
  process(-1) = 1
  process( 3) = 6
  call kaleu_put_process( mdl,vtx,obj(3) ,process ,nfin ,roots )
!  4  q g    -> Q Qbar q    (+ gluons)
!     qbar g -> Q Qbar qbar (+ gluons)
  process(-2) = 1
  process(-1) = 6
  process( 3) = 6
  call kaleu_put_process( mdl,vtx,obj(4) ,process ,nfin ,roots )
!
!---Q Qbar + 2 q jets final state
  if (njets.ge.2) then
  nproc = 10
!
!  5  g g    -> Q Qbar q qbar   (+ gluons)
  process(-2) = 1
  process(-1) = 1
  process( 3) = 6
  process( 4) = 6
  call kaleu_put_process( mdl,vtx,obj(5) ,process ,nfin ,roots )
!  6  q qbar -> Q Qbar q qbar   (+ gluons)   q=u,d,c,s
!     qbar q -> Q Qbar q qbar   (+ gluons)   q=u,d,c,s
  process(-2) = 6
  process(-1) = 6
  process( 3) = 6
  process( 4) = 6
  call kaleu_put_process( mdl,vtx,obj(6) ,process ,nfin ,roots )
!  7  q qbar -> Q Qbar q' qbar' (+ gluons)   q=u,d,c,s
!     qbar q -> Q Qbar q' qbar' (+ gluons)   q=u,d,c,s
  process(-2) = 6
  process(-1) = 6
  process( 3) = 7
  process( 4) = 7
  call kaleu_put_process( mdl,vtx,obj(7) ,process ,nfin ,roots )
!  8  q q'       -> Q Qbar q q'           (+ gluons)    q=u,d,c,s
!     qbar q'bar -> Q Qbar qbar q'bar     (+ gluons)    q=u,d,c,s
  process(-2) = 7
  process(-1) = 6
  process( 3) = 6
  process( 4) = 7
  call kaleu_put_process( mdl,vtx,obj(8) ,process ,nfin ,roots )
!  9  q qbar'    -> Q Qbar q qbar'        (+ gluons)    q=u,d,c,s
!     qbar q'    -> Q Qbar qbar q'        (+ gluons)    q=u,d,c,s
  process(-2) = 7
  process(-1) = 6
  process( 3) = 6
  process( 4) = 7
  call kaleu_put_process( mdl,vtx,obj(9) ,process ,nfin ,roots )
! 10  q q        -> Q Qbar q q            (+ gluons)    q=u,d,c,s
!     qbar qbar  -> Q Qbar qbar qbar      (+ gluons)    q=u,d,c,s
  process(-2) = 6
  process(-1) = 6
  process( 3) = 6
  process( 4) = 6
  call kaleu_put_process( mdl,vtx,obj(10) ,process ,nfin ,roots )
!
!---Q Qbar + 3 q jets final state
  if (njets.ge.3) then
  nproc = 14
!
! 11  g q    -> Q Qbar q q' qbar'    (+ gluon)
!     g qbar -> Q Qbar qbar q' qbar' (+ gluon)
  process(-2) = 6
  process(-1) = 1
  process( 3) = 6
  process( 4) = 7
  process( 5) = 7
  call kaleu_put_process( mdl,vtx,obj(11) ,process ,nfin ,roots )
! 12  q g    -> Q Qbar q q' qbar'    (+ gluon)
!     qbar g -> Q Qbar qbar q' qbar' (+ gluon)
  process(-2) = 1
  process(-1) = 6
  process( 3) = 6
  process( 4) = 7
  process( 5) = 7
  call kaleu_put_process( mdl,vtx,obj(12) ,process ,nfin ,roots )
! 13  g q    -> Q Qbar q q qbar      (+ gluon)
!     g qbar -> Q Qbar qbar q qbar   (+ gluon)
  process(-2) = 6
  process(-1) = 1
  process( 3) = 6
  process( 4) = 6
  process( 5) = 6
  call kaleu_put_process( mdl,vtx,obj(13) ,process ,nfin ,roots )
! 14  q g    -> Q Qbar q q qbar      (+ gluon)
!     qbar g -> Q Qbar qbar q qbar   (+ gluon)
  process(-2) = 1
  process(-1) = 6
  process( 3) = 6
  process( 4) = 6
  process( 5) = 6
  call kaleu_put_process( mdl,vtx,obj(14) ,process ,nfin ,roots )
!
!---Q Qbar + 4 q jets final state
  if (njets.ge.4) then
  nproc = 16
!
! 15  g g    -> Q Qbar q qbar q qbar   (+ gluons)  q=u,d,c,s
  process(-2) = 1
  process(-1) = 1
  process( 3) = 6
  process( 4) = 6
  process( 5) = 6
  process( 6) = 6
  call kaleu_put_process( mdl,vtx,obj(15) ,process ,nfin ,roots )
! 16  g g    -> Q Qbar q qbar q' qbar' (+ gluons)
  process(-2) = 1
  process(-1) = 1
  process( 3) = 6
  process( 4) = 6
  process( 5) = 7
  process( 6) = 7
  call kaleu_put_process( mdl,vtx,obj(16) ,process ,nfin ,roots )
!
  endif
  endif
  endif
  endif
!
  do ii=1,nproc
    call kaleu_init_strf( str(ii),obj(ii) ,0d0 )
    call kaleu_init_adapt( mdl,obj(ii) ,int(nopt/nproc),niter,1d-3 )
    call kaleu_updt_smin( obj(ii) ,smin )
    call kaleu_updt_strf( str(ii),obj(ii) ,smin )
  enddo
  end


  subroutine QQ_kaleu_gnrt( jproc ,x1,x2 ,pcm ,lnot )
!*********************************************************************
!*********************************************************************
  use QQ_kaleu
  implicit none
  integer          ,intent(in)  :: jproc
  real(kind(1d0))  ,intent(out) :: x1,x2
  real(kind(1d0))  ,intent(out) :: pcm(0:3,10)
  integer          ,intent(out) :: lnot
  real(kind(1d0)) :: pkaleu(0:3,-2:17)
  logical         :: discard
  lnot = 0
  call kaleu_gnrt_strf_lab( str(jproc),obj(jproc) ,discard ,x1,x2 )
  if (discard) then
    lnot = 1
    return
  endif
  call kaleu_gnrt( mdl,obj(jproc) ,discard ,pkaleu )
  if (discard) then
    lnot = 1
  else
    pcm(0:3,2) = -pkaleu(0:3,-2)
    pcm(0:3,1) = -pkaleu(0:3,-1)
    pcm(0:3,3:2+nfin) = pkaleu(0:3,1:nfin)
  endif
  end


  subroutine alpgen_kaleu_wght( jproc ,weight )
!*********************************************************************
!*********************************************************************
  use QQ_kaleu
  implicit none
  integer         ,intent(in)  :: jproc
  real(kind(1d0)) ,intent(out) :: weight
  real(kind(1d0)) :: ww
  call kaleu_wght( mdl,obj(jproc) ,weight )
  call strf_wght( str(jproc) ,ww )
  weight = weight*ww
  end


  subroutine alpgen_kaleu_collect( jproc ,weight )
!*********************************************************************
!*********************************************************************
  use QQ_kaleu
  implicit none
  integer         ,intent(in) :: jproc
  real(kind(1d0)) ,intent(in) :: weight
  if (weight.gt.0d0) call strf_collect( str(jproc) ,weight )
  call kaleu_collect( obj(jproc) ,weight )
  end
