module helac_kaleu
  use avh_kaleu
  use avh_kaleu_model
  use avh_kaleu_strf
  type(model_type) ,save :: mdl ! Masses and widths
  type(kaleu_type) ,save :: obj ! Instance of the phase space generator
  type(strf_type)  ,save :: str ! Concerns the generation of x1,x2
  logical          ,save :: strf_yes ! Should Kaleu generate x1,x2 ?
end module
      

  subroutine helac_kaleu_init( parmas,parwid,ihiggs,onlyqcd,withqcd ,&
                               ifl,nn ,&
                               ecm,istruc ,&
                               nbatch,nstep ,&
                               smin ) 
!*********************************************************************
!* Everything is input.
!*   parmas: particle masses
!*   parwid: particle widths
!*   ihiggs: should the higgs be included?
!*  onlyqcd: speaks for itself
!*  withqcd: speaks for itself
!*      ifl: the particles involved in the process
!*       nn: the number of particles involved in the process
!*      ecm: the centre-of-mass energy
!*   istruc: should Kaleu generate x1,x2 and the inst-momenta?
!*   nbatch: number of accepted ps-points before optimization step
!*    nstep: number of optimization steps
!*     smin: matrix of lower limits on positive 2-particle invariants,
!*           and upper limits on negative invariants
!*********************************************************************
  use helac_kaleu
  implicit none
  integer         ,intent(in) :: ihiggs,nn,istruc,nbatch,nstep
  integer         ,intent(in) :: ifl(nn)
  logical         ,intent(in) :: onlyqcd,withqcd
  real(kind(1d0)) ,intent(in) :: parmas(-12:41),parwid(-12:41),ecm
  real(kind(1d0)) ,intent(in) :: smin(-2:17,-2:17)
  integer :: ii,hel2kal(-12:41),process(-2:17)
  logical :: cancel
! The list vtx of possible interaction vertices is temporary. It is
! only used to create the tree stored inside obj
  type(vertx_type) :: vtx
!
  real(kind(1d0)) ,parameter :: small = 1d-6 ! for channel removal
!
  call kaleu_hello
!
  call helac_kaleu_model &
      ( mdl,vtx ,hel2kal ,parmas,parwid ,ihiggs,onlyqcd,withqcd )
  process(-2) = hel2kal(ifl(2))
  process(-1) = hel2kal(ifl(1))
  do ii=1,nn-2
    process(ii) = hel2kal(ifl(ii+2))
  enddo
  call kaleu_put_process( mdl,vtx,obj ,process ,nn-2 ,ecm ,cancel )
  if (cancel) stop
!
  strf_yes = (istruc.eq.1)
  if (strf_yes) call kaleu_init_strf( str,obj ,0d0 )
!
  call kaleu_init_adapt( mdl,obj ,nbatch,nstep,small )
!
  call kaleu_updt_smin( obj ,smin )
  if (strf_yes) call kaleu_updt_strf( str,obj ,smin )
!
  end


  subroutine helac_kaleu_gnrt( discard ,x1kaleu,x2kaleu ,pkaleu )
!*********************************************************************
!*********************************************************************
  use helac_kaleu
  implicit none
  logical          ,intent(out) :: discard
  real(kind(1d0))  ,intent(out) :: x1kaleu,x2kaleu,pkaleu(0:3,-2:17)
  if (strf_yes) then
    call kaleu_gnrt_strf( str,obj ,discard ,x1kaleu,x2kaleu )
  else
    x1kaleu = 1d0
    x2kaleu = 1d0
  endif
  call kaleu_gnrt( mdl,obj ,discard ,pkaleu )
  end


  subroutine helac_kaleu_wght( weight )
!*********************************************************************
!*********************************************************************
  use helac_kaleu
  implicit none
  real(kind(1d0)) ,intent(out) :: weight
  real(kind(1d0)) :: ww
  call kaleu_wght( mdl,obj ,weight )
  if (strf_yes) then
    call strf_wght( str ,ww )
    weight = weight*ww
  endif
  end


  subroutine helac_kaleu_collect( weight )
!*********************************************************************
!*********************************************************************
  use helac_kaleu
  implicit none
  real(kind(1d0)) ,intent(in) :: weight
  if (strf_yes.and.weight.gt.0d0) call strf_collect( str ,weight )
  call kaleu_collect( obj ,weight )
  end


  subroutine helac_kaleu_close
!*********************************************************************
!*********************************************************************
  use helac_kaleu
  implicit none
  call kaleu_close( obj )
  if (strf_yes) call strf_close( str )
  end
