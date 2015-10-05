module avh_kaleu_vertex
  use avh_kaleu_model
  use avh_kaleu_tree ,only : tree_type
  use avh_kaleu_ranvar
  use avh_kaleu_kinem
!
  private
  public :: vertex_type,vertex_init,vertex_init_adapt,vertex_close &
           ,vertex_gnrt,vertex_wght,vertex_collect,vertex_plotgrids
!
! The number 2*pi
  real(kind(1d0)) ,parameter :: twopi = 6.2831853071795864769252867665590d0
! Unit messages are send to
  integer ,parameter :: unitg=0
  integer ,parameter :: unitw=6
!
!
  type :: vertex_type
    private
    integer :: typ ,i0,i1,i2 ,iz=0 ,j0,j1,j2
    logical :: sch1=.false.,sch2=.false.
    type(ranvar_type) :: rv1,rv2,rvt
    real(kind(1d0)) :: prob=0.5d0,av1=0d0,av2=0d0,dn1=0d0,dn2=0d0
    integer         :: idat=0,istp=0
  end type
!
!
contains
!
!
!********************************************************************
  include 'kaleu_vertex.h'
!********************************************************************
!
!
  subroutine vertex_init( tree,iv ,obj )
!********************************************************************
!********************************************************************
  use avh_bint
  implicit none
  type(tree_type)   ,intent(in)    :: tree
  type(vertex_type) ,intent(inout) :: obj
  integer           ,intent(in)    :: iv
  integer :: f1,f2,f3
  obj%typ = tree%typ(iv)
  obj%iz  = tree%zax(iv)
  obj%i0 = tree%po( tree%ov(0,iv) )
  obj%i1 = tree%po( tree%ov(1,iv) )
  obj%i2 = tree%po( tree%ov(2,iv) )
  f1 = tree%fo( tree%ov(1,iv) )
  f2 = tree%fo( tree%ov(2,iv) )
  f3 = tree%fo( tree%ov(3,iv) )
  if     (obj%typ.eq.0) then
    obj%j0 = obj%i0
    obj%j1 = obj%i1
    obj%j2 = obj%i2
    obj%sch1 = (avh_bint_l(obj%j1).gt.1)
    obj%sch2 = (avh_bint_l(obj%j2).gt.1)
  elseif (obj%typ.eq.2) then ! s-channel vertex
    obj%j0 = obj%i0
    obj%j1 = obj%i1
    obj%j2 = obj%i2
    obj%sch1 = (avh_bint_l(obj%j1).gt.1)
    obj%sch2 = (avh_bint_l(obj%j2).gt.1)
    if (obj%sch1) call ranvar_init( obj%rv1 ,f1 )
    if (obj%sch2) call ranvar_init( obj%rv2 ,f2 )
  elseif (obj%typ.eq.1) then ! t-variable is s(j1+iz)=s(j1+1)=s(i1)
    obj%j0 = obj%i0-1
    obj%j1 = obj%i1-1
    obj%j2 = obj%i2
    obj%sch1 = (avh_bint_l(obj%j1).gt.1)
    obj%sch2 = (avh_bint_l(obj%j2).gt.1)
    if (obj%sch1) call ranvar_init( obj%rv1 ,f3 )
    if (obj%sch2) call ranvar_init( obj%rv2 ,f2 )
                  call ranvar_init( obj%rvt ,f1 )
  elseif (obj%typ.eq.3) then ! t-variable is s(j1+iz)=s(i1-1+iz)
    obj%j0 = obj%i0-1
    obj%j1 = obj%i1-1
    obj%j2 = obj%i2
    obj%sch1 = (avh_bint_l(obj%j1).gt.1)
    obj%sch2 = (avh_bint_l(obj%j2).gt.1)
    if (obj%sch1) call ranvar_init( obj%rv1 ,f3 )
    if (obj%sch2) call ranvar_init( obj%rv2 ,f2 )
                  call ranvar_init( obj%rvt , 1 ) ! overrule it to a gluon
  elseif (obj%typ.eq.5) then ! t-variable is s(j1+iz)=s(j2+iz)
    obj%j0 = obj%i0-1
    obj%j1 = obj%i2
    obj%j2 = obj%i1-1
    obj%sch1 = (avh_bint_l(obj%j1).gt.1)
    obj%sch2 = (avh_bint_l(obj%j2).gt.1)
    if (obj%sch1) call ranvar_init( obj%rv1 ,f2 )
    if (obj%sch2) call ranvar_init( obj%rv2 ,f3 )
                  call ranvar_init( obj%rvt , 1 ) ! overrule it to a gluon
  elseif (obj%typ.eq.4) then ! t-variable is s(j1+iz)=s(i1+iz)
    obj%j0 = obj%i0
    obj%j1 = obj%i1
    obj%j2 = obj%i2
    obj%sch1 = (avh_bint_l(obj%j1).gt.1)
    obj%sch2 = (avh_bint_l(obj%j2).gt.1)
    if (obj%sch1) call ranvar_init( obj%rv1 ,f1 )
    if (obj%sch2) call ranvar_init( obj%rv2 ,f2 )
                  call ranvar_init( obj%rvt , 1 ) ! overrule it to a gluon
  endif
  end subroutine
!
!
  subroutine vertex_close( obj )
!********************************************************************
!********************************************************************
  implicit none
  type(vertex_type) ,intent(inout) :: obj
  if (obj%sch1)    call ranvar_close( obj%rv1 )
  if (obj%sch2)    call ranvar_close( obj%rv2 )
  if (obj%iz.gt.0) call ranvar_close( obj%rvt )
  obj%sch1 = .false.
  obj%sch2 = .false.
  obj%iz = 0
  end subroutine
!
!
  subroutine vertex_gnrt( mdl ,kin ,obj ,discard )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: mdl
  type(vertex_type) ,intent(inout) :: obj
  type(kinem_type)  ,intent(inout) :: kin
  logical           ,intent(out)   :: discard
  integer :: j0,j1,j2,iz,ityp,jj
  real(kind(1d0)) :: xx,lami,lamo,phi,ct,st,tmin,tmax,Ez,E1,hz,h1
  real(kind(1d0)) :: pz0,pz1,pz2,pz3,p0,p1,p2,p3,absz,abs1
  logical :: above
!
  discard = .false.
!
  j0 = obj%j0
  j1 = obj%j1
  j2 = obj%j2
  ityp = obj%typ
  iz = obj%iz                              
!
  if (ityp.eq.0) then
!No generation, only happens as first vertex for pure S-channel
    wrkp(0,obj%i2) = wrkp(0,obj%i0) - wrkp(0,obj%i1)
    wrkp(1,obj%i2) = wrkp(1,obj%i0) - wrkp(1,obj%i1)
    wrkp(2,obj%i2) = wrkp(2,obj%i0) - wrkp(2,obj%i1)
    wrkp(3,obj%i2) = wrkp(3,obj%i0) - wrkp(3,obj%i1)
!
  elseif (iz.gt.0) then
    call gnrt_s1s2( mdl ,kin ,obj ,discard )
    if (discard) return
    if (wrkm(j0).lt.wrkm(j1)+wrkm(j2)) then
      if (unitg.gt.0) write(unitg,*) &
        'WARNING from vertex_gnrt T-chan: m0-m1-m2 < 0, discard event.' 
      discard = .true.
      return
    endif
    Ez  = wrkp(0,iz)
    pz1 = wrkp(1,iz)
    pz2 = wrkp(2,iz)
    pz3 = wrkp(3,iz)
    if (Ez.lt.0d0) then ! pz may have negative energy,
      Ez  = -Ez         ! put it positive for this construction.
      pz1 = -pz1
      pz2 = -pz2
      pz3 = -pz3
    endif
    call tsoob( Ez,pz1,pz2,pz3 ,&
                wrkp(0,j0),wrkp(1,j0),wrkp(2,j0),wrkp(3,j0),wrkm(j0) )
    absz = dsqrt( pz1**2 + pz2**2 + pz3**2 )
    h1 = (wrkm(j0)-wrkm(j1)-wrkm(j2))*(wrkm(j0)+wrkm(j1)+wrkm(j2)) &
        *(wrkm(j0)+wrkm(j1)-wrkm(j2))*(wrkm(j0)-wrkm(j1)+wrkm(j2))
    if (h1.lt.0d0) then
      if (unitg.gt.0) write(unitg,*) &
        'WARNING from vertex_gnrt T-chan: lam =',h1,', discard event.' 
      discard = .true.
      return
    endif
    abs1 = dsqrt(h1)/( 2*wrkm(j0) )
    E1 = abs1*dsqrt( 1d0 + 4*wrks(j0)*wrks(j1)/h1 )
    if (iz.eq.1.or.iz.eq.kin%b(kin%n1)-1) then
      ct = wrks(j1) + wrks(iz) - 2*E1*Ez
    else
      ct = wrks(j1) + wrks(iz) + 2*E1*Ez ! never tested for |wrks(j1/iz)|>0
    endif
    st = 2*abs1*absz
    tmin = ct-st
    tmax = ct+st
    jj = j1+iz
    if (jj.gt.kin%b(kin%n1)) jj = iz-j1
    call gnrt_t( mdl ,obj%rvt ,kin,jj ,tmin,tmax ,discard )
    if (discard) return
    call avh_random( phi )
    phi = twopi*phi
!Given the variables, construct the momenta 
    if (iz.eq.1.or.iz.eq.kin%b(kin%n1)-1) then
      ct = ( wrks(jj) - wrks(j1) - wrks(iz) + 2*Ez*E1 )/(2*absz)
    else
      ct = (-wrks(jj) + wrks(j1) + wrks(iz) + 2*Ez*E1 )/(2*absz)
    endif
    st = (abs1+ct)*(abs1-ct)
    if (st.lt.0d0) then
      if (unitg.gt.0) write(unitg,*) &
        'WARNING from vertex_gnrt T-chan: 1-ct^2 =',st/abs1**2,', discard event.' 
      discard = .true.
      return
    endif
    st = dsqrt(st)
    p1 = st*dsin(phi)
    p2 = st*dcos(phi)
    p3 = ct
    call rot3( p1,p2,p3 ,pz1,pz2,pz3,absz )
    call boost( E1,p1,p2,p3 ,&
                wrkp(0,j0),wrkp(1,j0),wrkp(2,j0),wrkp(3,j0),wrkm(j0) )
    wrkp(0,j1) = E1
    wrkp(1,j1) = p1
    wrkp(2,j1) = p2
    wrkp(3,j1) = p3
    wrkp(0,j2) = wrkp(0,j0) - wrkp(0,j1)
    wrkp(1,j2) = wrkp(1,j0) - wrkp(1,j1)
    wrkp(2,j2) = wrkp(2,j0) - wrkp(2,j1)
    wrkp(3,j2) = wrkp(3,j0) - wrkp(3,j1)
    if (obj%i0.ne.j0) then
      if (j1.eq.obj%i2) j1 = j2
      wrkp(0,obj%i1) = wrkp(0,j1) + wrkp(0,1)
      wrkp(1,obj%i1) = wrkp(1,j1) + wrkp(1,1)
      wrkp(2,obj%i1) = wrkp(2,j1) + wrkp(2,1)
      wrkp(3,obj%i1) = wrkp(3,j1) + wrkp(3,1)
      if (ityp.ne.1) then
        wrks(obj%i1) = wrkp(0,obj%i1)**2 &
                     - wrkp(1,obj%i1)**2 &
                     - wrkp(2,obj%i1)**2 &
                     - wrkp(3,obj%i1)**2
      endif
    endif
!
  else
!S-channel generation (without correlated cos(theta))
    call gnrt_s1s2( mdl ,kin ,obj ,discard )
    if (discard) return
    if (wrkm(j0).lt.wrkm(j1)+wrkm(j2)) then
      if (unitg.gt.0) write(unitg,*) &
        'WARNING from vertex_gnrt S-chan: m0-m1-m2 < 0, discard event.' 
      discard = .true.
      return
    endif
    lamo = (wrkm(j0)-wrkm(j1)-wrkm(j2))*(wrkm(j0)+wrkm(j1)+wrkm(j2))  &
          *(wrkm(j0)+wrkm(j1)-wrkm(j2))*(wrkm(j0)-wrkm(j1)+wrkm(j2))
    if (lamo.lt.0d0) then
      if (unitg.gt.0) write(unitg,*) &
        'WARNING from vertex_gnrt S-chan: lamo =',lamo,', discard event.'
      discard = .true.
      return
    endif
    lamo = dsqrt(lamo)
    call avh_random( xx )
    ct = 2d0*xx - 1d0
    call avh_random( phi )
    phi = twopi*phi
!Given the variables, construct the momenta 
    abs1 = lamo/(2*wrkm(j0))
    p0 = dsqrt( wrks(j1) + abs1*abs1 )
    st = dsqrt(1d0-ct*ct)
    ct = abs1*ct
    st = abs1*st
    p1 = st*dsin(phi)
    p2 = st*dcos(phi)
    p3 = ct
    call boost( p0,p1,p2,p3 ,&
                wrkp(0,j0),wrkp(1,j0),wrkp(2,j0),wrkp(3,j0),wrkm(j0) )
    wrkp(0,obj%i1) = p0
    wrkp(1,obj%i1) = p1
    wrkp(2,obj%i1) = p2
    wrkp(3,obj%i1) = p3
    wrkp(0,obj%i2) = wrkp(0,obj%i0) - p0
    wrkp(1,obj%i2) = wrkp(1,obj%i0) - p1
    wrkp(2,obj%i2) = wrkp(2,obj%i0) - p2
    wrkp(3,obj%i2) = wrkp(3,obj%i0) - p3
!
  endif
  end subroutine
!
  subroutine vertex_wght( mdl ,kin ,obj ,weight )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: mdl
  type(vertex_type) ,intent(inout) :: obj
  type(kinem_type)  ,intent(inout) :: kin
  real(kind(1d0))   ,intent(out)   :: weight
  integer :: j0,j1,j2,iz,ityp,jj
  real(kind(1d0)) :: lami,lamo,ww,phi,ct,st,cp,sp,p0,p1,p2,p3, &
                     tmin,tmax,absz,abs1,h1,hz,E1,Ez
  logical :: above
!
  j0 = obj%j0
  j1 = obj%j1
  j2 = obj%j2
  ityp = obj%typ
  iz = obj%iz                              
!
  if (ityp.eq.0) then
! No generation, no weight factor
    weight = 1d0
!
  elseif (iz.gt.0) then
    call wght_s1s2( mdl ,kin ,obj ,weight )
    if (weight.eq.0d0) return
    hz = wrkp(0,j0)*wrkp(0,iz) &
       - wrkp(1,j0)*wrkp(1,iz) &
       - wrkp(2,j0)*wrkp(2,iz) &
       - wrkp(3,j0)*wrkp(3,iz)
    hz = 4*( hz*hz - wrks(iz)*wrks(j0) )
    h1 = (wrkm(j0)-wrkm(j1)-wrkm(j2))*(wrkm(j0)+wrkm(j1)+wrkm(j2)) &
        *(wrkm(j0)+wrkm(j1)-wrkm(j2))*(wrkm(j0)-wrkm(j1)+wrkm(j2))
    if (hz.lt.0d0.or.h1.lt.0d0) then
      if (unitw.gt.0) write(unitw,*) &
        'WARNING from vertex_wght T-chan:' ,&
        ' hz =',hz,', h1 =',h1,', returning weight=0'
      weight = 0d0
      return
    endif
    absz = dsqrt(hz)/( 2*wrkm(j0) )
    Ez = absz*dsqrt( 1d0 + 4*wrks(j0)*wrks(iz)/hz )
    abs1 = dsqrt(h1)/( 2*wrkm(j0) )
    E1 = abs1*dsqrt( 1d0 + 4*wrks(j0)*wrks(j1)/h1 )
    if (iz.eq.1.or.iz.eq.kin%b(kin%n1)-1) then
      ct = wrks(j1) + wrks(iz) - 2*E1*Ez
    else
      ct = wrks(j1) + wrks(iz) + 2*E1*Ez ! never tested for |wrks(j1/iz)|>0
    endif
    st = 2*abs1*absz
    tmin = ct-st
    tmax = ct+st
    jj = j1+iz
    if (jj.gt.kin%b(kin%n1)) jj = iz-j1
    call wght_t( mdl ,obj%rvt ,kin,jj ,tmin,tmax ,ww )
    weight = weight * ww * twopi / (8d0*wrkm(j0)*absz)
!
  else
! S-channel generation (without correlated cos(theta))
    call wght_s1s2( mdl ,kin ,obj ,weight )
    if (weight.eq.0d0) return
    lamo = (wrkm(j0)-wrkm(j1)-wrkm(j2))*(wrkm(j0)+wrkm(j1)+wrkm(j2))  &
          *(wrkm(j0)+wrkm(j1)-wrkm(j2))*(wrkm(j0)-wrkm(j1)+wrkm(j2))
    if (lamo.lt.0d0) then
      if (unitw.gt.0) write(unitw,*) &
        'WARNING from vertex_wght S-chan: lamo =',lamo,', returning weight=0'
      weight = 0d0
      return
    endif
    lamo = dsqrt(lamo)
    weight = weight * twopi * lamo / (4*wrks(j0))
!
  endif
  end subroutine
!
!
  subroutine vertex_collect( obj ,weight ,nbatch,nsteps )
!********************************************************************
!********************************************************************
  implicit none
  type(vertex_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(in)    :: weight
  integer           ,intent(in)    :: nbatch,nsteps
  real(kind(1d0)) :: hh
  if (obj%sch1.and.obj%sch2) then
    call ranvar_collect( obj%rv1 ,weight )
    call ranvar_collect( obj%rv2 ,weight )
    if (obj%istp.ge.nsteps) return
    if (weight.gt.0d0) then
      hh = obj%dn1*obj%prob + obj%dn2*(1d0-obj%prob)
      if (hh.ne.0d0) hh = weight*weight/hh
      obj%av1 = obj%av1 + hh*obj%dn1
      obj%av2 = obj%av2 + hh*obj%dn2
      obj%dn1 = 0d0
      obj%dn2 = 0d0
      obj%idat = obj%idat + 1
    endif
    if (obj%idat.eq.nbatch) then
      if     (obj%av1.eq.0d0) then
        if   (obj%av2.eq.0d0) return
        obj%prob = 1d0/dble(nbatch)
      elseif (obj%av2.eq.0d0) then
        obj%prob = 1d0-1d0/dble(nbatch)
      else
        hh = (1d0-obj%prob)/obj%prob * dsqrt( obj%av2/obj%av1 )
        obj%prob = 1d0/( 1d0 + hh )
      endif
      obj%av1 = 0d0
      obj%av2 = 0d0
      obj%idat = 0
      obj%istp = obj%istp + 1
    endif
  elseif (obj%sch1) then
    call ranvar_collect( obj%rv1 ,weight )
  elseif (obj%sch2) then
    call ranvar_collect( obj%rv2 ,weight )
  endif
  if (obj%iz.gt.0) call ranvar_collect( obj%rvt ,weight )
  end subroutine
!
!
  subroutine vertex_plotgrids( obj ,iunit )
!********************************************************************
!********************************************************************
  implicit none
  type(vertex_type) ,intent(in) :: obj
  integer           ,intent(in) :: iunit
  if (obj%sch1)    call ranvar_plot( obj%rv1 ,iunit )
  if (obj%sch2)    call ranvar_plot( obj%rv2 ,iunit )
  if (obj%iz.gt.0) call ranvar_plot( obj%rvt ,iunit )
  end subroutine
!
!
  subroutine gnrt_s1s2( mdl ,kin ,obj ,discard )
!********************************************************************
!* Given ss(i0), generate ss(i1) and ss(i2).
!* i0,i1,i2  should be even and correspond to positive invariants.
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: mdl
  type(vertex_type) ,intent(inout) :: obj
  type(kinem_type)  ,intent(inout) :: kin
  logical           ,intent(out)   :: discard
  real(kind(1d0)) :: xx
  discard = .false.
  if (obj%sch1.and.obj%sch2) then
    call avh_random( xx )
    if (xx.lt.obj%prob) then
      call gnrt_s( mdl ,obj%rv1 ,kin,obj%j0,obj%j1,0      ,discard )
      if (discard) return
      call gnrt_s( mdl ,obj%rv2 ,kin,obj%j0,obj%j2,obj%j1 ,discard )
    else
      call gnrt_s( mdl ,obj%rv2 ,kin,obj%j0,obj%j2,0      ,discard )
      if (discard) return
      call gnrt_s( mdl ,obj%rv1 ,kin,obj%j0,obj%j1,obj%j2 ,discard )
    endif
  elseif (obj%sch1) then
    call gnrt_s( mdl ,obj%rv1 ,kin,obj%j0,obj%j1,0 ,discard )
  elseif (obj%sch2) then
    call gnrt_s( mdl ,obj%rv2 ,kin,obj%j0,obj%j2,0 ,discard )
  endif
  end subroutine
!
  subroutine wght_s1s2( mdl ,kin ,obj ,weight )
!********************************************************************
!* Given ss(i0), generate ss(i1) and ss(i2).
!* i0,i1,i2  should be even and correspond to positive invariants.
!********************************************************************
  implicit none
  type(model_type)  ,intent(in)    :: mdl
  type(vertex_type) ,intent(inout) :: obj
  type(kinem_type)  ,intent(inout) :: kin
  real(kind(1d0))   ,intent(out)   :: weight
  real(kind(1d0)) :: w1,w2,ww
  weight = 0d0
  if (obj%sch1.and.obj%sch2) then
    obj%dn1 = 0d0
    call wght_s( mdl ,obj%rv1 ,kin,obj%j0,obj%j1,0      ,w1 )
    call wght_s( mdl ,obj%rv2 ,kin,obj%j0,obj%j2,obj%j1 ,ww )
    w1 = ww*w1
    if (w1.ne.0d0) obj%dn1 = 1d0/w1
    obj%dn2 = 0d0
    call wght_s( mdl ,obj%rv2 ,kin,obj%j0,obj%j2,0      ,w2 )
    call wght_s( mdl ,obj%rv1 ,kin,obj%j0,obj%j1,obj%j2 ,ww )
    w2 = ww*w2
    if (w2.ne.0d0) obj%dn2 = 1d0/w2
    weight = obj%prob*obj%dn1 + (1d0-obj%prob)*obj%dn2
    if (weight.ne.0d0) weight = 1d0/weight
  elseif (obj%sch1) then
    call wght_s( mdl ,obj%rv1 ,kin,obj%j0,obj%j1,0 ,weight )
  elseif (obj%sch2) then
    call wght_s( mdl ,obj%rv2 ,kin,obj%j0,obj%j2,0 ,weight )
  else
    weight = 1d0
  endif
  end subroutine
!
!
  subroutine boost( p0,p1,p2,p3 ,q0,q1,q2,q3,mq )
!********************************************************************
!* apply on (p0,p1,p2,p3) the boost that
!* boosts (mq,0,0,0) to (q0,q1,q2,q3)
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(inout) :: p0,p1,p2,p3
  real(kind(1d0)) ,intent(in)    :: mq,q0,q1,q2,q3 
  real(kind(1d0)) :: aa,bb
!      mq = dsqrt(q0*q0 - q1*q1 - q2*q2 - q3*q3)
  aa = (p0*q0 + p1*q1 + p2*q2 + p3*q3)/mq
  bb = (p0 + aa)/(q0 + mq)
  p0 = aa
  p1 = p1 + bb*q1
  p2 = p2 + bb*q2
  p3 = p3 + bb*q3
  end subroutine
!
  subroutine tsoob( p0,p1,p2,p3 ,q0,q1,q2,q3,mq )
!********************************************************************
!* apply on (p0,p1,p2,p3) the boost that
!* boosts (q0,q1,q2,q3) to (mq,0,0,0)
!********************************************************************
  implicit none
  real(kind(1d0)) :: p0,p1,p2,p3 ,mq,q0,q1,q2,q3 ,aa,bb
!      mq = dsqrt(q0*q0 - q1*q1 - q2*q2 - q3*q3)
  aa = (p0*q0 - p1*q1 - p2*q2 - p3*q3)/mq
  bb = (p0 + aa)/(q0 + mq)
  p0 = aa
  p1 = p1 - bb*q1
  p2 = p2 - bb*q2
  p3 = p3 - bb*q3
  end subroutine
!
  subroutine rot3( v1,v2,v3 ,x1,x2,x3,lx )
!********************************************************************
!* apply on (v1,v2,v3) the inverse of the rotation
!* that rotates (x1,x2,x3) to the 3-axis
!* lx is the length of (x1,x2,x3)
!********************************************************************
  implicit none
  real(kind(1d0)) :: v1,v2,v3 ,lx,x1,x2,x3  ,h1,h2,h3 ,&
                     stsp,stcp,ct,st,sp,cp,ctsp,ctcp
!  lx = dsqrt( x1*x1 + x2*x2 + x3*x3 )
  stsp = x1/lx
  stcp = x2/lx
  ct   = x3/lx
  st   = dsqrt(stsp*stsp + stcp*stcp)
  if (st.ne.0d0) then
    sp = stsp/st
    cp = stcp/st
  else
    sp = 0d0
    cp = 1d0
  endif
  ctsp = ct*sp
  ctcp = ct*cp
  h1 = v1
  h2 = v2
  h3 = v3
  v1 =  cp*h1 + ctsp*h2 + stsp*h3
  v2 = -sp*h1 + ctcp*h2 + stcp*h3
  v3 =        - st  *h2 + ct  *h3
  end subroutine
!
  subroutine tor3( v1,v2,v3 ,x1,x2,x3,lx )
!********************************************************************
!* apply on (v1,v2,v3) the rotation
!* that rotates (x1,x2,x3) to the 3-axis
!* lx is the length of (x1,x2,x3)
!********************************************************************
  implicit none
  real(kind(1d0)) :: v1,v2,v3 ,lx,x1,x2,x3  ,h1,h2,h3 ,&
                     stsp,stcp,ct,st,sp,cp,ctsp,ctcp
!  lx = dsqrt( x1*x1 + x2*x2 + x3*x3 )
  stsp = x1/lx
  stcp = x2/lx
  ct   = x3/lx
  st   = dsqrt(stsp*stsp + stcp*stcp)
  if (st.ne.0d0) then
    sp = stsp/st
    cp = stcp/st
  else
    sp = 0d0
    cp = 1d0
  endif
  ctsp = ct*sp
  ctcp = ct*cp
  h1 = v1
  h2 = v2
  h3 = v3
  v1 =   cp*h1 -   sp*h2
  v2 = ctsp*h1 + ctcp*h2 - st*h3
  v3 = stsp*h1 + stcp*h2 + ct*h3
  end subroutine
!
!
end module
