module QQQQ_kaleu
  use avh_kaleu
  use avh_kaleu_model
  use avh_kaleu_strf
  logical          ,save :: yes = .false.
  integer          ,save :: nfin
  type(model_type) ,save :: mdl     ! Masses and widths
  type(kaleu_type) ,save :: obj(6) ! Instances of the phase space generator
  type(strf_type)  ,save :: str(6) ! Concerns the generation of x1,x2
end module
  

  subroutine alpgen_kaleu_init( ikaleu )
  use QQQQ_kaleu
  yes = (ikaleu.eq.1)
  end subroutine

  function alpgen_kaleu_yes() result(value)
  use QQQQ_kaleu
  logical value
  value = yes
  end function
    

  subroutine QQQQ_kaleu_init( ihvy,ihvy2, njets ,roots ,apar &
                             ,ptjmin ,ptQmin ,ptQ2min   &
                             ,etajmax,etaQmax,etaQ2max  &
                             ,drjmin ,drQmin ,drQ2min   &
                             ,nopt,niter  ) 
!*********************************************************************
!*********************************************************************
  use QQQQ_kaleu
  implicit none
  integer         ,intent(in) :: ihvy,ihvy2,njets,niter
  real(kind(1d0)) ,intent(in) :: roots,apar(100),nopt &
  ,ptjmin,ptQmin,ptQ2min,etajmax,etaQmax,etaQ2max,drjmin,drQmin,drQ2min
  type(vertx_type) :: vtx
  real(kind(1d0))  :: smin(-2:17,-2:17)
  integer :: process(-2:17),ii,nproc,ihq,ihq2
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
    write(*,*) 'ERROR in QQQQ_kaleu_init: heavy quark',ihvy,' not defined'
    stop
  endif
  if     (ihvy2.eq.4) then
    ihq2 = 8
  elseif (ihvy2.eq.5) then
    ihq2 = 11
  elseif (ihvy2.eq.6) then
    ihq2 = 10
  else
    write(*,*) 'ERROR in QQQQ_kaleu_init: heavy quark',ihvy2,' not defined'
    stop
  endif
  nfin = njets+4
!
  smin(-2:nfin,-2:nfin) = 0d0
  smin(-1,1     ) = -2d0*ptQmin*ptQmin*dexp(-etaQmax)
  smin(-1,2     ) = smin(-1,1)
  smin(-1,3     ) = -2d0*ptQ2min*ptQ2min*dexp(-etaQ2max)
  smin(-1,4     ) = smin(-1,3)
  smin(-1,5:nfin) = -2d0*ptjmin*ptjmin*dexp(-etajmax)
  smin(-2,1:nfin) = smin(-1,1:nfin)
  smin(1:nfin,-2:-1) = smin(-2:-1,1:nfin)
!
  smin(1,2) = 2d0*ptQmin*ptQmin*(1d0-dcos(drQmin))
  smin(3,4) = 2d0*ptQ2min*ptQ2min*(1d0-dcos(drQ2min))
  smin(1,3) = 2d0*ptQmin*ptQ2min*(1d0-dcos(min(drQmin,drQ2min))) 
  smin(1,4) = smin(1,3)
  smin(2,3) = smin(1,3)
  smin(2,4) = smin(1,3)
!
  smin(2,1) = smin(1,2)
  smin(4,3) = smin(3,4)
  smin(3,1) = smin(1,3)
  smin(4,1) = smin(1,4)
  smin(3,2) = smin(2,3)
  smin(4,2) = smin(2,4)
!
  smin(1,5:nfin) = 2d0*ptQmin*ptjmin*(1d0-dcos(drjmin)) ! not sure
  smin(2,5:nfin) = smin(1,5:nfin)
  smin(3,5:nfin) = 2d0*ptQ2min*ptjmin*(1d0-dcos(drjmin)) ! not sure
  smin(4,5:nfin) = smin(3,5:nfin)
  smin(5:nfin,1:4   ) = smin(1:4,5:nfin)
  smin(5:nfin,5:nfin) = 2d0*ptjmin*ptjmin*(1d0-dcos(drjmin))
!
  process(1:2   ) = ihq
  process(3:4   ) = ihq2
  process(5:nfin) = 1
!---Q1 Q1bar Q2 Q2bar final state
  nproc = 2
!
!  1  g g    -> Q1 Q1bar Q2 Q2bar (+ gluons)
  process(-2) = 1
  process(-1) = 1
  call kaleu_put_process( mdl,vtx,obj(1) ,process ,nfin ,roots )
!  2  q qbar -> Q1 Q1bar Q2 Q2bar (+ gluons)   
!     qbar q -> Q1 Q1bar Q2 Q2bar (+ gluons)   
  process(-2) = 6
  process(-1) = 6
  call kaleu_put_process( mdl,vtx,obj(2) ,process ,nfin ,roots )
!
!---Q1 Q1bar Q2 Q2bar + 1 q jet final state
  if (njets.ge.1) then
  nproc = 4
!
!  3  g q    -> Q1 Q1bar Q2 Q2bar q    (+ gluon)
!     g qbar -> Q1 Q1bar Q2 Q2bar qbar (+ gluon)
  process(-2) = 6
  process(-1) = 1
  process( 5) = 6
  call kaleu_put_process( mdl,vtx,obj(3) ,process ,nfin ,roots )
!  4  q g    -> Q1 Q1bar Q2 Q2bar q    (+ gluon)
!     qbar g -> Q1 Q1bar Q2 Q2bar qbar (+ gluon)
  process(-2) = 1
  process(-1) = 6
  process( 5) = 6
  call kaleu_put_process( mdl,vtx,obj(4) ,process ,nfin ,roots )
!
!---Q1 Q1bar Q2 Q2bar + 2 q jets final state
  if (njets.ge.2) then
  nproc = 6
!
!  5  g g    -> Q1 Q1bar Q2 Q2bar q qbar
  process(-2) = 1
  process(-1) = 1
  process( 5) = 6
  process( 6) = 6
  call kaleu_put_process( mdl,vtx,obj(5) ,process ,nfin ,roots )
!  6  g g    -> Q1 Q1bar Q2 Q2bar Q1 Q1bar
  process(-2) = 1
  process(-1) = 1
  process( 5) = ihq
  process( 6) = ihq
  call kaleu_put_process( mdl,vtx,obj(6) ,process ,nfin ,roots )
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


  subroutine QQQQ_kaleu_gnrt( jproc ,x1,x2 ,pcm ,lnot )
!*********************************************************************
!*********************************************************************
  use QQQQ_kaleu
  implicit none
  integer          ,intent(in)  :: jproc
  real(kind(1d0))  ,intent(out) :: x1,x2
  real(kind(1d0))  ,intent(out) :: pcm(0:3,10)
  integer          ,intent(out) :: lnot
  real(kind(1d0)) :: pkaleu(0:3,-2:17)
  logical         :: discard
  integer         :: ii
  ii = jproc
  lnot = 0
  call kaleu_gnrt_strf_lab( str(ii),obj(ii) ,discard ,x1,x2 )
  if (discard) then
    lnot = 1
    return
  endif
  call kaleu_gnrt( mdl,obj(ii) ,discard ,pkaleu )
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
  use QQQQ_kaleu
  implicit none
  integer         ,intent(in)  :: jproc
  real(kind(1d0)) ,intent(out) :: weight
  real(kind(1d0)) :: ww
  integer         :: ii
  ii = jproc
  call kaleu_wght( mdl,obj(ii) ,weight )
  call strf_wght( str(ii) ,ww )
  weight = weight*ww
  end


  subroutine alpgen_kaleu_collect( jproc ,weight )
!*********************************************************************
!*********************************************************************
  use QQQQ_kaleu
  implicit none
  integer         ,intent(in) :: jproc
  real(kind(1d0)) ,intent(in) :: weight
  integer         :: ii
  ii = jproc
  if (weight.gt.0d0) call strf_collect( str(ii) ,weight )
  call kaleu_collect( obj(ii) ,weight )
  end
