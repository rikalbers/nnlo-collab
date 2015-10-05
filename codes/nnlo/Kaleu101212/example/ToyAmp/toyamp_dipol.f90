module avh_toyamp_dipol
  use avh_toyamp
  use avh_toyamp_model
!
  private ! Everything is private except the following list
  public :: dipol_type,dipol_close,dipol_init,dipol_mom,dipol_calc  &
           ,maxndip
!
! Maximal number of dipole terms
  integer ,parameter :: maxndip = 252
! Unit messages are send to
  integer ,parameter :: nunit = 6
! The number 2*pi
  real(kind(1d0)) ,parameter :: twopi = 6.2831853071795864769252867665590d0
! Hard cut-off
  real(kind(1d0)) ,parameter :: cutoff = 1d-8
!
!
  type :: dipol_type
    private
    integer :: ndip=0,npair=0
    integer           ,allocatable :: i(:),j(:),k(:),pair(:,:),riap(:,:)
    type(toyamp_type) ,allocatable :: toy(:)
    type(toyamp_type)              :: dip
    real(kind(1d0))   ,allocatable :: s(:)
    real(kind(1d0)) :: aff,afi,aif,aii
  end type
!
!
contains
!
!
  subroutine dipol_close( obj )
!**********************************************************************
!**********************************************************************
  implicit none
  type(dipol_type) ,intent(inout) :: obj
  integer :: ii
  do ii=1,obj%npair
    call toyamp_close( obj%toy(ii) )
  enddo
  call toyamp_close( obj%dip )
  deallocate( obj%toy )
  deallocate( obj%s   )
  deallocate( obj%i )
  deallocate( obj%j )
  deallocate( obj%k )
  deallocate( obj%pair )
  obj%ndip  = 0
  obj%npair = 0
  end subroutine
!
!
  subroutine dipol_init( mdl,mlv,obj ,process,nfinst ,cancel &
                        ,aff,afi,aif,aii                     &
                        ,dip,ndip )
!**********************************************************************
!* dip(1,l)=i, dip(2,l)=j, dip(3,l)=k
!**********************************************************************
  use particles
  implicit none
  type(model_type) ,intent(inout) :: mdl
  type(vertx_type) ,intent(inout) :: mlv
  type(dipol_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: process(-2:17),nfinst
  real(kind(1d0))  ,intent(in)    :: aff,afi,aif,aii
  integer          ,intent(out)   :: dip(3,maxndip),ndip
  integer :: tmpproc(-2:17),tmpdip(3,(nfinst+2)**3),fi,fj,fk,ii,jj,kk
  logical :: cancel
!
  obj%npair = 0
  do ii=-2,nfinst
    if (ii.eq.0) cycle
    fi = sign(1,ii)*process(ii)
    if (.not.is_qcd(abs(fi))) cycle
    do jj=ii+1,nfinst
      if (jj.le.0) cycle
      fj = process(jj)
      if (.not.is_qcd(abs(fj))) cycle
      if (abs(fi).ne.gluon.and.abs(fj).ne.gluon) cycle
      obj%npair = obj%npair+1
      if (abs(fj).eq.gluon) then
        tmpdip(1,obj%npair)=ii; tmpdip(2,obj%npair)=jj
      else
        tmpdip(2,obj%npair)=ii; tmpdip(1,obj%npair)=jj
      endif
    enddo
  enddo
!
  if (obj%npair.lt.2) then
    write(*,*) 'ERROR in avh_toyamp_dipol: not enough partons in process.'
    stop
  endif
!
  allocate( obj%pair(1:2,1:obj%npair) )
  obj%pair(1:2,1:obj%npair) = tmpdip(1:2,1:obj%npair)
!
  allocate( obj%riap(-2:nfinst,-2:nfinst) )
  obj%riap = 0
  do ii=1,obj%npair
    obj%riap( obj%pair(1,ii) ,obj%pair(2,ii) ) = ii
    obj%riap( obj%pair(2,ii) ,obj%pair(1,ii) ) = ii
  enddo
!
  allocate( obj%toy(1:obj%npair) )
  allocate( obj%s(  1:obj%npair) )
  do kk=1,obj%npair
    ii = obj%pair(1,kk)
    jj = obj%pair(2,kk)
    call toyamp_put_process( mdl,mlv,obj%toy(kk) ,process,nfinst ,ii,jj ,cancel )
    if (cancel) then 
      write(*,*) 'ERROR in avh_toyamp_dipol: no graphs with' &
                ,process(ii),process(jj),' in one vertex.'
      stop
    endif
  enddo
!
  obj%ndip = 0
  do jj=1,nfinst
    fj = process(jj)
    if (abs(fj).ne.gluon) cycle
    do ii=-2,nfinst
      if (ii.eq.0.or.ii.eq.jj) cycle
      fi = sign(1,ii)*process(ii)
      if (.not.is_qcd(abs(fi))) cycle
      if (abs(fi).eq.gluon.and.ii.gt.jj) cycle ! final-state gluons are equival
      do kk=-2,nfinst
        if (kk.eq.0.or.kk.eq.jj.or.kk.eq.ii) cycle
        fk = sign(1,kk)*process(kk)
        if (.not.is_qcd(abs(fk))) cycle
        obj%ndip = obj%ndip+1
        if (obj%ndip.gt.maxndip) then
          write(*,*) 'ERROR in avh_toyamp_dipol_init: increase the parameter ' &
                    ,'maxndip to at least',obj%ndip
          stop
        endif
        tmpdip(1,obj%ndip) = ii
        tmpdip(2,obj%ndip) = jj
        tmpdip(3,obj%ndip) = kk
      enddo
    enddo
  enddo
!
  allocate( obj%i(1:obj%ndip) )
  allocate( obj%j(1:obj%ndip) )
  allocate( obj%k(1:obj%ndip) )
  obj%i(1:obj%ndip) = tmpdip(1,1:obj%ndip)
  obj%j(1:obj%ndip) = tmpdip(2,1:obj%ndip)
  obj%k(1:obj%ndip) = tmpdip(3,1:obj%ndip)
  tmpproc(-2:nfinst) = process(-2:nfinst)
  do ii=1,nfinst
    if (abs(process(ii)).eq.gluon) then
      tmpproc(ii:nfinst-1) = tmpproc(ii+1:nfinst)
      call toyamp_put_process( mdl,mlv,obj%dip ,tmpproc,nfinst-1 ,cancel=cancel )
      if (cancel) then 
        write(*,*) 'ERROR in avh_toyamp_dipol: ' &
                  ,'underlying process not possible.'
        stop
      endif
      exit
    endif
  enddo
!
  ndip = obj%ndip
  dip(1,1:obj%ndip) = obj%i(1:obj%ndip)
  dip(2,1:obj%ndip) = obj%j(1:obj%ndip)
  dip(3,1:obj%ndip) = obj%k(1:obj%ndip)
!
!
! Put alpha_max
  obj%aff = aff 
  obj%afi = afi 
  obj%aif = aif 
  obj%aii = aii 
!  write(6,*) 'obj%aff', obj%aff !DEBUG
!  write(6,*) 'obj%afi', obj%afi !DEBUG
!  write(6,*) 'obj%aif', obj%aif !DEBUG
!  write(6,*) 'obj%aii', obj%aii !DEBUG
! 
  end subroutine
!
!
  subroutine dipol_mom( obj ,pkaleu ,idip ,pdip ,notzero )
!*********************************************************************
!*********************************************************************
  implicit none
  type(dipol_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(in)  :: pkaleu(0:3,-2:17)
  integer          ,intent(in)  :: idip
  real(kind(1d0))  ,intent(out) :: pdip(0:3,-2:17)
  logical          ,intent(out) :: notzero
  real(kind(1d0)) :: mijsq,si,sj,sk,mi,mj,mk,xx,zz,pipj,pkpj,w1
  integer         :: ii,jj,kk,ll,nn,ij,aa
!
  nn = toyamp_get_nfinst( obj%toy(1) )
!
  notzero = .false.
  pdip(0:3,-2:nn) = pkaleu(0:3,-2:nn)
!
  if (idip.eq.0) then
    do kk=1,obj%npair
      ii = obj%pair(1,kk)
      jj = obj%pair(2,kk)
      obj%s(kk) = pkaleu(0,ii)*pkaleu(0,jj) - pkaleu(1,ii)*pkaleu(1,jj) &
                - pkaleu(2,ii)*pkaleu(2,jj) - pkaleu(3,ii)*pkaleu(3,jj)
      if (2*abs(obj%s(kk)).lt.cutoff*(pkaleu(0,-1)+pkaleu(0,-2))**2) then
        return
      endif
    enddo
  else
    ii = obj%i(idip)
    jj = obj%j(idip)
    kk = obj%k(idip)
    if (ii.gt.0.and.kk.gt.0) then
! Final-final
      ij = ii
      if (ij.gt.jj) ij = ij-1
      si    = toyamp_get_s( obj%toy(1) ,ii )
      sj    = toyamp_get_s( obj%toy(1) ,jj )
      sk    = toyamp_get_s( obj%toy(1) ,kk )
      mi    = toyamp_get_m( obj%toy(1) ,ii )
      mj    = toyamp_get_m( obj%toy(1) ,jj )
      mk    = toyamp_get_m( obj%toy(1) ,kk )
      mijsq = toyamp_get_s( obj%dip    ,ij )
      call wght_ff( w1 ,pdip ,xx,zz ,ii,jj,kk ,mijsq,si,sj,sk,mi,mj,mk )
      if (w1.eq.0d0) return
      xx = xx/obj%aff
      if (xx.le.0d0.or.1d0.le.xx) return
      if (zz.le.0d0.or.1d0.le.zz) return
    elseif ( ii.gt.0.and.kk.lt.0 .or. ii.lt.0.and.kk.gt.0 ) then
! Final-initial or initial-final
      aa = min(ii,kk)
      ii = max(ii,kk) ! ii possibly changed
      si = toyamp_get_s( obj%toy(1) ,ii )
      call wght_fi( w1 ,pdip ,xx,zz ,ii,jj,aa ,si )
      if (w1.eq.0d0) return
      if (kk.lt.0) then
! Final-initial
        xx = xx/obj%afi
        if (xx.le.0d0.or.1d0.le.xx) return
        if (zz.le.0d0.or.1d0.le.zz) return
      else
! Initial-final
        zz = zz/obj%aif
        if (xx.le.0d0.or.1d0.le.xx) return
        if (zz.le.0d0.or.1d0.le.zz) return
      endif
    else
! Initial-initial
      call wght_ii( w1 ,pdip ,xx,zz ,ii,jj,kk,nn )
      if (w1.eq.0d0) return
      xx = 1d0-xx
      if (xx.le.0d0.or.1d0.le.xx) return
      zz = zz/(xx*obj%aii)
      if (zz.le.0d0.or.1d0.le.zz) return
    endif
    pdip(0:3,jj:nn-1) = pdip(0:3,jj+1:nn)
  endif
!
  notzero = .true.
!
  end subroutine 
!
!
  subroutine dipol_calc( mdl,obj ,idip ,pkaleu ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(dipol_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: idip
  real(kind(1d0))  ,intent(in)    :: pkaleu(0:3,-2:17)
  real(kind(1d0))  ,intent(out)   :: weight
  real(kind(1d0))    :: sij,skj,amp(1:obj%npair)
  complex(kind(1d0)) :: ztmp
  integer            :: ii,jj,kk,ll,ipair,jpair
  real(kind(1d0)),parameter :: eps=1d-32
!
  if (idip.eq.0) then
!
    do ipair=1,obj%npair
      call toyamp_put_mom( obj%toy(ipair) ,pkaleu )
      call toyamp_calc( mdl,obj%toy(ipair) ,ztmp ) 
      ii = obj%pair(1,ipair)
      jj = obj%pair(2,ipair)
      amp(ipair) = abs( ztmp /sqrt( cmplx(2*obj%s(ipair),eps) ) )
    enddo
    weight = 0d0
    do ipair=1,obj%npair-1
      ii = obj%pair(1,ipair)
      jj = obj%pair(2,ipair)
      do jpair=ipair+1,obj%npair
        kk = obj%pair(1,jpair)
        ll = obj%pair(2,jpair)
        if (ii.ne.kk.and.ii.ne.ll.and.jj.ne.kk.and.jj.ne.ll) cycle
        weight = weight + amp(ipair)*amp(jpair)
      enddo
    enddo
!
  else
!
    ii = obj%i(idip)
    jj = obj%j(idip)
    kk = obj%k(idip)
    call toyamp_put_mom( obj%dip ,pkaleu )
    call toyamp_calc( mdl,obj%dip ,ztmp )
    weight = real( ztmp*conjg(ztmp) )
    sij = abs( obj%s(obj%riap(ii,jj)) )
    skj = abs( obj%s(obj%riap(kk,jj)) )
    weight = weight/( 4*sij*(sij+skj) )
!
  endif
!
  end subroutine 
!
!
!*********************************************************************
!* From here private routines
!*********************************************************************
!
!
  subroutine wght_ff( weight ,pp ,xy,xz ,ii,jj,kk ,mijsq,si,sj,sk,mi,mj,mk )
!*********************************************************************
!*********************************************************************
  implicit none
  real(kind(1d0)) ,intent(out)   :: weight ,xy,xz
  real(kind(1d0)) ,intent(inout) :: pp(0:3,-2:17)
  integer         ,intent(in) :: ii,jj,kk
  real(kind(1d0)) ,intent(in) :: mijsq,si,sj,sk,mi,mj,mk
  real(kind(1d0)) :: pij(0:3),pik(0:3),QQ(0:3) &
                    ,sij,sik,sQ,mQ,hh,yy,ymin,ymax,lam,vk,vi,rr
!
  weight = 0d0
!
  pij(0:3) = pp(0:3,ii) + pp(0:3,jj)
  pik(0:3) = pp(0:3,ii) + pp(0:3,kk)
  QQ(0:3)  = pij(0:3) + pp(0:3,kk) 
  sij = (pij(0)+pij(3))*(pij(0)-pij(3)) - pij(1)**2 - pij(2)**2
  sik = (pik(0)+pik(3))*(pik(0)-pik(3)) - pik(1)**2 - pik(2)**2
  sQ  = (QQ(0)+QQ(3))*(QQ(0)-QQ(3)) - QQ(1)**2 - QQ(2)**2
  if (sQ.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in toyamp_dipol in wght_ff: sQ<0, putting weight=0'
    return
  endif
  mQ = dsqrt(sQ)
  hh = sQ-si-sj-sk
  yy = (sij-si-sj) !/hh ! yy=hh*y_usual, has dimension mass^2
  ymin = 2*mi*mj
  ymax = hh - 2*mk*(mQ-mk)
  xy = (yy-ymin)/(ymax-ymin)
!
  lam = (sQ-mijsq-sk)**2 - 4*mijsq*sk
  if (lam.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in toyamp_dipol in wght_ff: lam<0, putting weight=0'
    return
  endif
  vk = (sQ-sij-sk)**2 - 4*sij*sk
  if (vk.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in toyamp_dipol in wght_ff: vk<0, putting weight=0'
    return
  endif
  vi = yy**2 - 4*si*sj
  if (vi.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in toyamp_dipol in wght_ff: vi<0, putting weight=0'
    return
  endif
!
!  zz = (sik-si-sk)/(hh-yy)
!  xz = ( zz*sij*2 - (sij+si-sj) )*(hh-yy)/dsqrt( 4*vk*vi ) + 0.5d0
  xz = ( (sik-si-si)*sij*2 - (sij+si-sj)*(hh-yy) )/dsqrt( 4*vk*vi ) + 0.5d0
!
  rr = dsqrt(lam/vk)
  weight = twopi*dsqrt(vi)*(ymax-ymin)/(4*rr*sij) ! extra factor (2pi)^3
!
  pp(0:3,kk) = rr * pp(0:3,kk) &
             + (sQ-mijsq+sk-rr*(2*sk+hh-yy))/(2*sQ) * QQ(0:3)
  pp(0:3,ii) = QQ(0:3) - pp(0:3,kk)
!
  end subroutine


  subroutine wght_fi( weight ,pp ,xx,zz ,ii,jj,aa ,si )
!********************************************************************
!* Assumes pp(0:3,aa) and pp(0:3,jj)  massless
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(out)   :: weight ,xx,zz
  real(kind(1d0)) ,intent(inout) :: pp(0:3,-2:17)
  integer         ,intent(in)    :: ii,jj,aa
  real(kind(1d0)) ,intent(in)    :: si
  real(kind(1d0)) :: pij(0:3),pia(0:3),QQ(0:3) &
                    ,sij,sia,sQ,pipj2,hh,zmin
  weight = 0d0
!
  pij(0:3) = pp(0:3,ii) + pp(0:3,jj)
  pia(0:3) = pp(0:3,ii) + pp(0:3,aa)
  QQ(0:3) = pij(0:3)    + pp(0:3,aa)
  sij = ( pij(0)+pij(3) )*( pij(0)-pij(3) ) - pij(1)**2 - pij(2)**2
  sia = ( pia(0)+pia(3) )*( pia(0)-pia(3) ) - pia(1)**2 - pia(2)**2
  sQ  = (  QQ(0)+ QQ(3) )*(  QQ(0)- QQ(3) ) -  QQ(1)**2 -  QQ(2)**2
  pipj2 = sij-si
  hh = sQ-si
  xx = -pipj2/(hh-pipj2)
  zz = (sia-si)/(hh-pipj2)
  zmin = (1d0-xx)*si/dabs(hh)
  zmin = zmin/(xx+zmin)
  zz = (zz-zmin)/(1d0-zmin)
  weight = twopi*dabs(hh)*(1d0-zmin)/(4*(1d0-xx)) ! extra factor (2pi)^3
  pp(0:3,aa) = (1d0-xx)*pp(0:3,aa)
  pp(0:3,ii) = QQ(0:3) - pp(0:3,aa) 
  end subroutine


  subroutine wght_ii( weight ,pp ,xx,vv ,aa,jj,bb,nn )
!********************************************************************
!* Assumes pp(0:3,aa,jj,bb) is massless
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(out)   :: weight ,xx,vv
  real(kind(1d0)) ,intent(inout) :: pp(0:3,-2:17)
  integer         ,intent(in)    :: aa,jj,bb,nn
  real(kind(1d0)) :: KK(0:3),Kt(0:3),Ks(0:3),sKK,sKt,sKs,cs,cK
  integer :: ll
!
  weight = 0d0
!
  Kt(0:3) = -pp(0:3,aa) - pp(0:3,bb)
  Ks(0:3) = -pp(0:3,aa) + pp(0:3,jj)
  KK(0:3) = Kt(0:3) - pp(0:3,jj)
  sKt = (Kt(0)+Kt(3))*(Kt(0)-Kt(3)) - Kt(1)**2 - Kt(2)**2
  sKs = (Ks(0)+Ks(3))*(Ks(0)-Ks(3)) - Ks(1)**2 - Ks(2)**2
  sKK = (KK(0)+KK(3))*(KK(0)-KK(3)) - KK(1)**2 - KK(2)**2
  xx = sKK/sKt
  vv = sKs/sKt
!
  pp(0:3,aa) = xx*pp(0:3,aa)
  Kt(0:3) = -pp(0:3,aa) - pp(0:3,bb)
  Ks(0:3) = KK(0:3) + Kt(0:3)
  sKt = (Kt(0)+Kt(3))*(Kt(0)-Kt(3)) - Kt(1)**2 - Kt(2)**2
  sKs = (Ks(0)+Ks(3))*(Ks(0)-Ks(3)) - Ks(1)**2 - Ks(2)**2
!
  weight = twopi*sKt*(1d0-xx)/(4*xx) ! extra factor (2pi)^3
!
  do ll=1,nn
    if (ll.eq.jj) cycle
    cs = Ks(0)*pp(0,ll)-Ks(1)*pp(1,ll)-Ks(2)*pp(2,ll)-Ks(3)*pp(3,ll)
    cK = KK(0)*pp(0,ll)-KK(1)*pp(1,ll)-KK(2)*pp(2,ll)-KK(3)*pp(3,ll)
    cs = 2*cs/sKs
    cK = 2*cK/sKK
    pp(0:3,ll) = pp(0:3,ll) - cs*Ks(0:3) + cK*Kt(0:3)
  enddo
!
  end subroutine


  subroutine rot3( vv ,xx ,calc )
!********************************************************************
!* apply on (vv(1),vv(2),vv(3)) the inverse of the rotation
!* that rotates (xx(1),xx(2),xx(3)) to the 3-axis
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(inout) :: vv(0:3)
  real(kind(1d0)) ,intent(in)    :: xx(0:3)
  logical         ,intent(in)    :: calc
  real(kind(1d0)) :: lx  ,h1,h2,h3
  real(kind(1d0)) :: stsp,stcp,ct,st,sp,cp,ctsp,ctcp
  save :: stsp,stcp,ct,st,sp,cp,ctsp,ctcp
  if (calc) then
    lx = dsqrt( xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3) )
    stsp = xx(1)/lx
    stcp = xx(2)/lx
    ct   = xx(3)/lx
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
  endif
  h1 = vv(1)
  h2 = vv(2)
  h3 = vv(3)
  vv(1) =  cp*h1 + ctsp*h2 + stsp*h3
  vv(2) = -sp*h1 + ctcp*h2 + stcp*h3
  vv(3) =        - st  *h2 + ct  *h3
  end subroutine
!
  subroutine boost( pp ,qq,mq )
!********************************************************************
!* apply on pp the boost that boosts (mq,0,0,0) to qq
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(inout) :: pp(0:3)
  real(kind(1d0)) ,intent(in)    :: qq(0:3),mq
  real(kind(1d0)) :: aa,bb
!      mq = dsqrt(qq(0)*qq(0) - qq(1)*qq(1) - qq(2)*qq(2) - qq(3)*qq(3))
  aa = (pp(0)*qq(0) + pp(1)*qq(1) + pp(2)*qq(2) + pp(3)*qq(3))/mq
  bb = (pp(0) + aa)/(qq(0) + mq)
  pp(0)   = aa
  pp(1:3) = pp(1:3) + bb*qq(1:3)
  end subroutine
!
  subroutine tsoob( pp ,qq,mq )
!********************************************************************
!* apply on pp the boost that boosts qq to (mq,0,0,0)
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(inout) :: pp(0:3)
  real(kind(1d0)) ,intent(in)    :: qq(0:3),mq
  real(kind(1d0)) :: aa,bb
!      mq = dsqrt(qq(0)*qq(0) - qq(1)*qq(1) - qq(2)*qq(2) - qq(3)*qq(3))
  aa = (pp(0)*qq(0) - pp(1)*qq(1) - pp(2)*qq(2) - pp(3)*qq(3))/mq
  bb = (pp(0) + aa)/(qq(0) + mq)
  pp(0)   = aa
  pp(1:3) = pp(1:3) - bb*qq(1:3)
  end subroutine
!
!
!
!
  subroutine printmom( pp ,nn )
!********************************************************************
!********************************************************************
  implicit none
  real(kind(1d0)) ,intent(in) :: pp(0:3,-2:17)
  integer         ,intent(in) :: nn
  integer ,parameter :: iunit=6
  real(kind(1d0)) :: ss,ptot(0:3)
  integer         :: ii
  ptot(0:3) = 0d0
  do ii=-2,nn
    if (ii.eq.0) cycle
    ss = pp(0,ii)**2 - pp(1,ii)**2 - pp(2,ii)**2 - pp(3,ii)**2
    write(iunit,'(i3,4d16.8,2x,d16.8)') ii,pp(0:3,ii),ss
    ptot(0:3) = ptot(0:3) + pp(0:3,ii)
  enddo
  ss = ptot(0)**2 - ptot(1)**2 - ptot(2)**2 - ptot(3)**2
  write(iunit,'(a3,4d16.8,2x,d16.8)') 'sum',ptot(0:3),ss
  end subroutine
!
!
end module
