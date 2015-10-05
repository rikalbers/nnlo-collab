  subroutine kaleu_init_adapt( mdl,obj ,nbatch,nstep,small &
                                       ,nbatch_grid,nbatch_cut )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(kaleu_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: nbatch,nstep
  integer ,optional,intent(in)    :: nbatch_grid,nbatch_cut
  real(kind(1d0))  ,intent(in)    :: small
  integer :: ii
!
  call mulcha_init_adapt( obj%mch ,nbatch,nstep,small )
!
  if (present(nbatch_cut)) call ranvar_global( nbatch_cut )
!
  end subroutine
