  subroutine kaleu_init_adapt( mdl ,obj ,nbatch,nstep,small )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(kaleu_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: nbatch,nstep
  real(kind(1d0))  ,intent(in)    :: small
!
  call mulcha_init_adapt( obj%mch ,nbatch,nstep,small )
!
  end subroutine
