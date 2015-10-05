  subroutine model_gnrt_above( model,ii ,discard,ss ,xx,slow,supp )
!********************************************************************
!* Construct invariant ss from input xx \in (0,1)
!********************************************************************
  implicit none
  type(model_type) ,intent(in)   :: model
  integer          ,intent(in)   :: ii
  logical         ,intent(inout) :: discard
  real(kind(1d0)) ,intent(out)   :: ss
  real(kind(1d0)) ,intent(in)    :: xx,slow,supp
  real(kind(1d0)) :: mgam,mas2,xupp,xlow,ex
!
  discard = .false.
  if (slow.ge.supp) then
    discard = .true.
    return
  endif
!
  if (model%pa(ii)%haswidth) then
    mgam = model%pa(ii)%mass * model%pa(ii)%width
    mas2 = model%pa(ii)%mass**2
    xupp = datan( (supp-mas2)/mgam )
    xlow = datan( (slow-mas2)/mgam )
    ss = mas2 + mgam*dtan( (xupp-xlow)*xx + xlow )
!
  elseif (model%pa(ii)%nowidth) then
!  elseif (.false.) then
    ex = model%pa(ii)%expon
    xupp = dabs(supp)**ex
    xlow = dabs(slow)**ex
    if     (slow.ge.0d0) then
      ss = ( (xupp-xlow)*xx + xlow )**(1d0/ex)
    elseif (supp.le.0d0) then
      ss =-( (xlow-xupp)*xx + xupp )**(1d0/ex)
    else!if(slow.lt.0d0.and.supp.gt.0d0)
      ss = (xupp+xlow)*xx - xlow
      if (ss.ge.0d0) then
        ss = ss**(1d0/ex)
      else
        ss =-(-ss)**(1d0/ex)
      endif
    endif
!
  else
    ss = ( supp-slow )*xx + slow
!
  endif
  end subroutine
!
  subroutine model_wght_above( model,ii ,ww,xx ,ss,slow,supp )
!********************************************************************
!* Re-construct xx\in(0,1) from invariant ss and return jacobian ww
!********************************************************************
  implicit none
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  real(kind(1d0)) ,intent(out) :: ww,xx
  real(kind(1d0)) ,intent(in)  :: ss,slow,supp
  real(kind(1d0)) :: mgam,mas2,xupp,xlow,ex
!
  if (ss.le.slow.or.supp.le.ss) then
    ww = 0d0
    xx = 0d0
    return
  endif
!
  if (model%pa(ii)%haswidth) then
    mgam = model%pa(ii)%mass * model%pa(ii)%width
    mas2 = model%pa(ii)%mass**2
    xupp = datan( (supp-mas2)/mgam )
    xlow = datan( (slow-mas2)/mgam )
    ww = xupp-xlow
    xx = ( datan((ss-mas2)/mgam) - xlow )/ww
    ww = ww*( ((ss-mas2)**2)/mgam + mgam )
!
  elseif (model%pa(ii)%nowidth) then
    ex = model%pa(ii)%expon
    xupp = dabs(supp)**ex
    xlow = dabs(slow)**ex
    if     (slow.ge.0d0) then
      ww = xupp-xlow
      xx = (ss**ex - xlow)/ww
      ww = ww*ss**(1d0-ex)/ex
    elseif (supp.le.0d0) then
      ww = xlow-xupp
      xx = ((-ss)**ex - xupp)/ww
      ww = ww*(-ss)**(1d0-ex)/ex
    else!if(slow.lt.0d0.and.supp.gt.0d0)
      ww = xlow+xupp
      if (ss.ge.0d0) then
        xx = (ss**ex + xlow)/ww
      else!if(ss.lt.0d0) then
        xx = (xlow - (-ss)**ex)/ww
      endif
      ww = ww*dabs(ss)**(1d0-ex)/ex
    endif
!
  else
    ww = supp-slow
    xx = (ss-slow)/ww
!
  endif
  end subroutine
!
  subroutine model_gnrt_below( model,ii ,discard,ss ,xx,slow,supp )
!********************************************************************
!* generate p(s) = 1/( s + frac*supp ) for  slow < s < supp
!********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: model
  integer          ,intent(in)    :: ii
  logical          ,intent(inout) :: discard
  real(kind(1d0))  ,intent(out)   :: ss
  real(kind(1d0))  ,intent(in)    :: xx,slow,supp
  real(kind(1d0)) :: hh
  discard = .false.
  if (slow.ge.supp) then
    discard = .true.
    return
  endif
!
  if (model%pa(ii)%nowidth) then
    hh = frac*(supp+slow)
    ss = (slow+hh)*( (supp+hh)/(slow+hh) )**xx - hh
!
  else
    ss = (supp-slow)*xx + slow
!
  endif
  end subroutine
!
  subroutine model_wght_below( model,ii ,ww,xx ,ss,slow,supp )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type) ,intent(in)  :: model
  integer          ,intent(in)  :: ii
  real(kind(1d0))  ,intent(out) :: ww,xx
  real(kind(1d0))  ,intent(in)  :: ss,slow,supp
  real(kind(1d0)) :: hh
  if (ss.le.slow.or.supp.le.ss) then
    ww = 0d0
    xx = 0d0
    return
  endif
!
  if (model%pa(ii)%nowidth) then
    hh = frac*(supp+slow)
    ww = dlog( (supp+hh)/(slow+hh) )
    xx = dlog( (ss  +hh)/(slow+hh) )/ww
    ww = ww*dabs(ss+hh)
!
  else
    ww = supp-slow
    xx = (ss-slow)/ww
!
  endif
  end subroutine
