! This version of the Kaleu interface allows to use dipole
! kinematics:
module multi_kaleu
use avh_kaleu
use avh_kaleu_model
use avh_kaleu_dipol
use avh_kaleu_strf
implicit none
!
  type(model_type) :: mdl
  type(dipol_type) , dimension(:) , allocatable :: obj
  type(strf_type)  , dimension(:) , allocatable :: str
!
  integer :: iinst
!
end module multi_kaleu
!
subroutine InitKaleuArrays(numproc)
use multi_kaleu
implicit none
!
  integer , intent(in) :: numproc
!
  integer :: istat
!
  allocate(obj(numproc),stat=istat)
  if (istat.ne.0) then
    print *,"Kaleu array obj is not allocated properly..."
    stop
  end if
  allocate(str(numproc),stat=istat)
  if (istat.ne.0) then
    print *,"Kaleu array str is not allocated properly..."
    stop
  end if
!
end subroutine InitKaleuArrays
!
subroutine multi_kaleu_init(iproc,process,nfinst,ecm ,option, &
                            masses,widths,labels,             &
                            onlyqcd,withqcd,withhiggs,        &
                            nbatch,nstep,thrs,smin,           &
                            aff,afi,aif,aii,                  &
                            ndip,dip)
use multi_kaleu
use kaleu_particles
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
  integer , intent(in) :: ndip
  integer , dimension(3,ndip) , intent(in) :: dip
  real(kind(1d0)) , intent(in) :: aff,afi,aif,aii
!
  integer :: prcss(-2:17)
  integer :: noffset,nbatch_grid,nbatch_cut
  logical :: cancel
!
  type(vertx_type) :: vtx
!
! Loading the model with its vertices:
  call user_model(mdl,vtx,masses,widths,labels, &
                  onlyqcd,withqcd,withhiggs)
!
! We are not interested in the sign of our particles:
  prcss(-2:nfinst) = abs(process(-2:nfinst))
! More optimization parameters:
  noffset = nbatch*nstep
  nbatch_grid = 1000
  nbatch_cut  = 1000
!
  call dipol_init( mdl,vtx,obj(iproc),prcss,nfinst,ecm, &
                   smin,                                &
                   nbatch,nstep,thrs,                   &
                   noffset,nbatch_grid,nbatch_cut,      &
                   dip,ndip,gluon,                      &
                   aff,afi,aif,aii)
  iinst = option
  if (iinst.eq.1) call dipol_init_strf(str(iproc),obj(iproc))
!
end subroutine multi_kaleu_init
!
subroutine multi_kaleu_gnrt(iproc,discard,x1kaleu,x2kaleu,pkaleu)
use multi_kaleu
implicit none
!
  integer , intent(in) :: iproc
  logical , intent(out)   :: discard
  real(kind(1d0)) , intent(out)   :: x1kaleu,x2kaleu
  real(kind(1d0)) , intent(inout) :: pkaleu(0:3,-2:17)
!
  x1kaleu = 1d0
  x2kaleu = 1d0
  if (iinst.eq.1) then
    call dipol_gnrt_strf(str(iproc),obj(iproc),discard,x1kaleu,x2kaleu)
  elseif (iinst.eq.2) then
    call dipol_inst(obj(iproc),pkaleu)
  endif
  call dipol_gnrt(mdl,obj(iproc),discard,pkaleu)
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
  call dipol_wght(mdl,obj(iproc),weight)
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
  call dipol_collect(obj(iproc),weight)
end subroutine multi_kaleu_collect
!
subroutine multi_kaleu_close(iproc)
use multi_kaleu
implicit none
!
  integer , intent(in) :: iproc
!
  call dipol_close(obj(iproc))
  if (iinst.eq.1) call strf_close(str(iproc))
end subroutine multi_kaleu_close
!
subroutine multi_kaleu_plotgrids(iproc,iunit)
use multi_kaleu
implicit none
!
  integer , intent(in) :: iproc,iunit
!
!
!
  call dipol_plotgrids(obj(iproc),iunit)
  if (iinst.eq.1) call strf_plot(str(iproc),iunit)
!
end subroutine multi_kaleu_plotgrids
