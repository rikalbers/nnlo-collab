module multi_kaleu
  use avh_kaleu
  use avh_kaleu_model
  use avh_kaleu_kinem
  use avh_kaleu_strf
  implicit none
!
  integer , parameter :: nmax = 20
!
  type(model_type) , save :: mdl
  type(kinem_type) , save :: kin(nmax)
  type(kaleu_type) , save :: obj(nmax)
  type(strf_type) , save :: str(nmax)
!
  integer , save :: iinst
end module multi_kaleu

subroutine multi_kaleu_init(iproc,process,nfinst,ecm ,option, &
                      masses,widths,labels, &
                      onlyqcd,withqcd,withhiggs, &
                      nbatch,nstep,thrs,smin)
  use multi_kaleu
  use particles
  implicit none
!
  integer , intent(in) :: iproc,process(-2:17),nfinst
  real(kind(1d0)) , intent(in) :: ecm
  integer , intent(in) :: option
  real(kind(1d0)) , intent(in) :: masses(nfirst:nlast),widths(nfirst:nlast)
  character(20) , intent(in) :: labels(nfirst:nlast)
  logical , intent(in) :: onlyqcd,withqcd,withhiggs
  integer , intent(in) :: nbatch,nstep
  real(kind(1d0)) , intent(in) :: thrs,smin(-2:17,-2:17)
  integer :: prcss(-2:17)
  logical :: cancel
!
  type(vertx_type) :: vtx
!
  call user_model(mdl,vtx,masses,widths,labels, &
                  onlyqcd,withqcd,withhiggs)
!
  prcss(-2:nfinst) = abs(process(-2:nfinst))
  call kaleu_put_process(mdl,vtx,obj(iproc),prcss,nfinst,ecm,cancel)
  if (cancel) stop
!
  iinst = option
  if (iinst.eq.1) call kaleu_init_strf(str(iproc),obj(iproc),0d0)
!
  call kaleu_init_adapt(mdl,obj(iproc),nbatch,nstep,thrs)
  call kaleu_updt_smin(obj(iproc),smin)
  call kaleu_print_smin(mdl,obj(iproc),prcss,nfinst,ecm,smin,6)
  if (iinst.eq.1) call kaleu_updt_strf(str(iproc),obj(iproc),smin)
!
end subroutine multi_kaleu_init
!
subroutine multi_kaleu_gnrt(iproc,discard,x1kaleu,x2kaleu,pkaleu)
  use multi_kaleu
  implicit none
  integer , intent(in) :: iproc
  logical , intent(out)   :: discard
  real(kind(1d0)) , intent(out)   :: x1kaleu,x2kaleu
  real(kind(1d0)) , intent(inout) :: pkaleu(0:3,-2:17)
  x1kaleu = 1d0
  x2kaleu = 1d0
  if (iinst.eq.1) then
    call kaleu_gnrt_strf(str(iproc),obj(iproc),discard,x1kaleu,x2kaleu)
  elseif (iinst.eq.2) then
    call kaleu_inst(obj(iproc),discard,pkaleu)
  endif
  call kaleu_gnrt(mdl,obj(iproc),discard,pkaleu)
end subroutine multi_kaleu_gnrt
!
subroutine multi_kaleu_wght(iproc,weight)
  use multi_kaleu
  implicit none
!
  integer , intent(in) :: iproc
  real(kind(1d0)) ,intent(out) :: weight
  real(kind(1d0)) :: ww
!
  call kaleu_wght(mdl,obj(iproc),weight)
  if (iinst.eq.1) then
    call strf_wght(str(iproc),ww)
    weight = weight*ww
  endif
end subroutine multi_kaleu_wght
!
subroutine multi_kaleu_collect(iproc,weight)
  use multi_kaleu
  implicit none
!
  integer , intent(in) :: iproc
  real(kind(1d0)) ,intent(in) :: weight
!
  if (iinst.eq.1.and.weight.gt.0d0) call strf_collect(str(iproc),weight)
  call kaleu_collect(obj(iproc),weight)
end subroutine multi_kaleu_collect
!
subroutine multi_kaleu_close(iproc)
  use multi_kaleu
  implicit none
!
  integer , intent(in) :: iproc
!
  call kaleu_close(obj(iproc))
  if (iinst.eq.1) call strf_close(str(iproc))
end subroutine multi_kaleu_close
