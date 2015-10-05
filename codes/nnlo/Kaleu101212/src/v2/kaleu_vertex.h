  subroutine vertex_init_adapt( obj ,iv )
!********************************************************************
!********************************************************************
  implicit none
  type(vertex_type) ,intent(inout) :: obj
  integer           ,intent(in)    :: iv
  character(8) :: label
  label(1:8) = 'S000V000'
  write(label(6:8),'(i3.3)') iv
  if (obj%sch1) then
    write(label(2:4),'(i3.3)') obj%j1
    call ranvar_init_grid( obj%rv1 ,label )
  endif
  if (obj%sch2) then
    write(label(2:4),'(i3.3)') obj%j2
    call ranvar_init_grid( obj%rv2  ,label )
  endif
  if (obj%iz.gt.0) then
    write(label(2:4),'(i3.3)') obj%j1+obj%iz
    call ranvar_init_grid( obj%rvt ,label )
  endif
  end subroutine
!
!
  subroutine gnrt_s( mdl ,rvr ,kin,j0,jj,jprev ,discard )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: mdl
  type(ranvar_type) ,intent(inout) :: rvr
  type(kinem_type)  ,intent(inout) :: kin
  integer           ,intent(in)    :: j0,jj,jprev
  logical           ,intent(out)   :: discard
  real(kind(1d0)) :: smin,smax
  if (jprev.gt.0) then
    smax = ( wrkm(j0) - wrkm(jprev) )**2
  else
    smax = wrks(j0)
  endif
  smin = kin%smin(jj)
  call ranvar_gnrt( mdl ,rvr ,discard,wrks(jj) ,kin%techcut,smin,smax )
  if (discard) return
  if (wrks(jj).le.0d0) then
    if (unitg.gt.0) write(unitg,*) &
      'WARNING from vertex_gnrt_s: s =',wrks(jj),', discard event.' 
    discard = .true.
    return
  endif
  wrkm(jj) = dsqrt(wrks(jj))
  end subroutine
!
  subroutine wght_s( mdl ,rvr ,kin,j0,jj,jprev ,weight )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: mdl
  type(ranvar_type) ,intent(inout) :: rvr
  type(kinem_type)  ,intent(in)    :: kin
  integer           ,intent(in)    :: j0,jj,jprev
  real(kind(1d0))   ,intent(out)   :: weight
  real(kind(1d0)) :: smin,smax
  weight = 0d0
  if (jprev.gt.0) then
    smax = ( wrkm(j0) - wrkm(jprev) )**2
  else
    smax = wrks(j0)
  endif
  smin = kin%smin(jj)
  call ranvar_wght( mdl ,rvr ,weight,wrks(jj) ,kin%techcut,smin,smax )
  end subroutine
!
!
  subroutine gnrt_t( mdl ,rvr ,kin,jj ,tmin,tmax ,discard )
!********************************************************************
!********************************************************************
  type(model_type)  ,intent(in)    :: mdl
  type(ranvar_type) ,intent(inout) :: rvr
  type(kinem_type)  ,intent(in)    :: kin
  integer           ,intent(in)    :: jj
  real(kind(1d0))   ,intent(in)    :: tmin,tmax
  logical           ,intent(out)   :: discard
  real(kind(1d0)) :: t0,t1
  if     (tmax.gt.dabs(tmin)) then
!    t0 = max( tmin ,kin%smin(jj) )
    t0 = kin%smin(jj)
    t1 = tmax
    call ranvar_gnrt( mdl ,rvr ,discard,wrks(jj) ,kin%techcut,t0,t1 )
  elseif (tmin.lt.-dabs(tmax)) then
    t0 = tmin
!    t1 = min( tmax ,kin%smin(jj) )
    t1 = kin%smin(jj)
    call ranvar_gnrt( mdl ,rvr ,discard,wrks(jj) ,t0,t1,-kin%techcut )
  endif
  end subroutine
!
!
  subroutine wght_t( mdl ,rvr ,kin,jj ,tmin,tmax ,weight )
!********************************************************************
!********************************************************************
  type(model_type)  ,intent(in)    :: mdl
  type(ranvar_type) ,intent(inout) :: rvr
  type(kinem_type)  ,intent(in)    :: kin
  integer           ,intent(in)    :: jj
  real(kind(1d0))   ,intent(in)    :: tmin,tmax
  real(kind(1d0))   ,intent(out)   :: weight
  real(kind(1d0)) :: t0,t1
  if     (tmax.gt.dabs(tmin)) then
!    t0 = max( tmin ,kin%smin(jj) )
    t0 = kin%smin(jj)
    t1 = tmax
    call ranvar_wght( mdl ,rvr ,weight,wrks(jj) ,kin%techcut,t0,t1 )
  elseif (tmin.lt.-dabs(tmax)) then
    t0 = tmin
!    t1 = min( tmax ,kin%smin(jj) )
    t1 = kin%smin(jj)
    call ranvar_wght( mdl ,rvr ,weight,wrks(jj) ,t0,t1,-kin%techcut )
  endif
  end subroutine
