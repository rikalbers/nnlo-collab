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
  integer         ,parameter :: nchmax_grid = 200
  real(kind(1d0)) ,parameter :: fusefrac_grid = 2d0/(2*nchmax_grid-1)
!
  call mulcha_init_adapt( obj%mch ,nbatch,nstep,small )
!
  if     (present(nbatch_grid).and.present(nbatch_cut)) then
    call ranvar_global( nbatch_grid ,nchmax_grid ,fusefrac_grid ,nbatch_cut )
  elseif (present(nbatch_grid)                        ) then
    call ranvar_global( nbatch_grid ,nchmax_grid ,fusefrac_grid )
  elseif (                         present(nbatch_cut)) then
    call ranvar_global( nchmax_g=nchmax_grid ,fusefrac_g=fusefrac_grid &
                       ,nbatch_c=nbatch_cut )
  else
    call ranvar_global( nchmax_g=nchmax_grid ,fusefrac_g=fusefrac_grid )
  endif
!
  do ii=1,obj%mch%nv
    call vertex_init_adapt( obj%vtx(ii) ,ii )
  enddo
!
  end subroutine
