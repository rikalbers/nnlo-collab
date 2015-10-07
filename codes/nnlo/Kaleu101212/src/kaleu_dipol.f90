module avh_kaleu_dipol
  use avh_kaleu
  use avh_kaleu_model
  use avh_kaleu_splvar
  use avh_kaleu_grid
  use avh_kaleu_strf
!
  private ! Everything is private except the following list
  public :: dipol_close,dipol_type,dipol_init,dipol_init_strf &
           ,dipol_gnrt_strf,dipol_inst,dipol_gnrt,dipol_wght  &
           ,dipol_collect,dipol_plotgrids,dipol_printdata
!
!
! For the 1-dim grids
  integer         ,parameter :: def_nbatch_g = 1000
  integer         ,parameter :: def_nchmax_g = 200
  real(kind(1d0)) ,parameter :: def_fractn_g = 2d0/(2*def_nchmax_g-1)
! Number of events for optimizing splvars before optimizing grids
  integer ,parameter :: def_noffset = 400000
! Unit messages are send to
  integer ,parameter :: nunit = 6
! The number 2*pi
  real(kind(1d0)) ,parameter :: twopi = 6.2831853071795864769252867665590d0
!
!
  type :: mch_type
    logical :: yes = .false.
    integer :: i0,i1
    real(kind(1d0)) ,allocatable :: wch(:),ave(:),dns(:)
    integer         :: idat,ndat,ntot,istp,nstp
    real(kind(1d0)) :: thrs
  end type
!
!
  type :: dipol_type
    private
    integer :: ndip=0,nproc=0,nfinst
    integer          ,allocatable :: i(:),j(:),k(:),proc(:)
    include 'dipol_splvar_type.h'
    type(kaleu_type) ,allocatable :: kal(:)
    type(mch_type) :: mch
    real(kind(1d0)) :: aff,afi,aif,aii
    real(kind(1d0)) :: p(0:3,-2:17)
    integer :: idip=-1 ,idat=0 ,noffset=def_noffset
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
  call mch_close( obj%mch )
  do ii=0,obj%nproc
    call kaleu_close( obj%kal(ii) )
  enddo
  deallocate( obj%kal )
  deallocate( obj%i )
  deallocate( obj%j )
  deallocate( obj%k )
  deallocate( obj%proc )
  include 'dipol_splvar_close.h'
  obj%ndip = 0
  obj%nproc = 0
  end subroutine
!
!
  subroutine dipol_init( mdl,vtx,obj ,process,nfinst ,ecm &
                        ,smin                             &
                        ,nbatch,nstep,thrs                &
                        ,noffset,nbatch_grid,nbatch_cut   &
                        ,dip,ndip,gluon                   &
                        ,aff,afi,aif,aii                  )
!**********************************************************************
!* dip(1,l)=i, dip(2,l)=j, dip(3,l)=k
!**********************************************************************
  use avh_print
  implicit none
  type(model_type) ,intent(inout) :: mdl
  type(vertx_type) ,intent(inout) :: vtx
  type(dipol_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: process(-2:17),nfinst,ndip,gluon
  integer          ,intent(in)    :: nbatch,nstep
  integer,optional ,intent(in)    :: noffset,nbatch_grid,nbatch_cut
  integer          ,intent(in)    :: dip(3,ndip)
  real(kind(1d0))  ,intent(in)    :: ecm,smin(-2:17,-2:17),thrs
  real(kind(1d0))  ,intent(in)    :: aff,afi,aif,aii
  real(kind(1d0)) :: smin1(-2:17,-2:17)
  integer         :: proclist(-2:17,ndip),tmpproc(-2:17)
  integer         :: l1,l2,l3,iproc,j_proc(ndip),jj,nbatch_g
  logical         :: next,cancel,notexists(0:ndip)
  character(8)    :: lbl
  character(2) ,parameter :: symb(-2:9) = &
      (/'-2','-1','+0','+1','+2','+3','+4','+5','+6','+7','+8','+9'/)
  character(160) :: line
!
  call kaleu_hello
!
  obj%ndip = ndip
  obj%nproc = 0
!
  allocate( obj%i(   0:obj%ndip) )
  allocate( obj%j(   0:obj%ndip) )
  allocate( obj%k(   0:obj%ndip) )
  allocate( obj%proc(0:obj%ndip) )
  obj%i = 0
  obj%j = 0
  obj%k = 0
  obj%proc = 0
!
! Put dipoles and find processes
  do l1=1,ndip
    obj%i(l1) = dip(1,l1)
    obj%j(l1) = dip(2,l1)
    obj%k(l1) = dip(3,l1)
    tmpproc(-2:nfinst) = process(-2:nfinst)
    if     (process(obj%i(l1)).eq.gluon) then
      tmpproc(obj%i(l1)) = process(obj%j(l1))
    elseif (process(obj%j(l1)).eq.gluon) then
      tmpproc(obj%i(l1)) = process(obj%i(l1))
    else
      tmpproc(obj%i(l1)) = gluon
    endif
    tmpproc(obj%j(l1):nfinst-1) = tmpproc(obj%j(l1)+1:nfinst)
    tmpproc(0) = 0
    if (nunit.gt.0) then
      line = 'MESSAGE from kaleu_dipol'//trim(printsint(l1))//':' &
      //' '//trim(printsint(obj%i(l1)))//trim(sparticle(mdl,process(obj%i(l1)))) &
      //' '//trim(printsint(obj%j(l1)))//trim(sparticle(mdl,process(obj%j(l1)))) &
      //' '//trim(printsint(obj%k(l1)))//trim(sparticle(mdl,process(obj%k(l1)))) 
      line = trim(line)//' , '//trim(sparticle(mdl,process(-2)))        &
                         //' '//trim(sparticle(mdl,process(-1)))//' ->'
      do l2=1,nfinst-1
        line = trim(line)//' '//trim(sparticle(mdl,tmpproc(l2)))
      enddo
      write(nunit,*) trim(line)
    endif
! ! check if tmpproc is already in the list, add if not
    next = .true.
    l2 = 0
    do while (next.and.l2.lt.obj%nproc)
      l2 = l2+1
      do l3=-2,nfinst-1
!Each dipole its own        next = (tmpproc(l3).ne.proclist(l3,l2))
        if (next) exit
      enddo
    enddo
    if (next) then
      obj%nproc = obj%nproc+1
      proclist(-2:nfinst-1,obj%nproc) = tmpproc(-2:nfinst-1)
      j_proc(obj%nproc) = obj%j(l1)
      l2 = l2+1
    endif
    obj%proc(l1) = l2
  enddo
!
! Put processes
  allocate( obj%kal(0:obj%nproc) )
  call kaleu_put_process( mdl,vtx,obj%kal(0) ,process,nfinst ,ecm ,cancel )
  if (cancel) stop
  if (present(nbatch_cut)) then
    call kaleu_init_adapt( mdl,obj%kal(0) ,nbatch,nstep,thrs &
                                          ,nbatch_cut=nbatch_cut )
  else
    call kaleu_init_adapt( mdl,obj%kal(0) ,nbatch,nstep,thrs )
  endif
!  call printsmin( mdl,process ,smin,nfinst ) !DEBUG
  call kaleu_updt_smin( obj%kal(0) ,smin )
  call kaleu_get_p( obj%kal(0) ,obj%p )
  do iproc=1,obj%nproc
    tmpproc(-2:nfinst-1) = proclist(-2:nfinst-1,iproc)
    call kaleu_put_process( mdl,vtx,obj%kal(iproc) ,tmpproc,nfinst-1 ,ecm ,cancel )
    notexists(iproc) = cancel
    if (cancel) cycle
    if (present(nbatch_cut)) then
      call kaleu_init_adapt( mdl,obj%kal(iproc) ,nbatch,nstep,thrs &
                                                ,nbatch_cut=nbatch_cut )
    else
      call kaleu_init_adapt( mdl,obj%kal(iproc) ,nbatch,nstep,thrs )
    endif
    jj = j_proc(iproc)
    smin1( -2:jj-1     ,-2:jj-1     )=smin(   -2:jj-1   ,  -2:jj-1   )
    smin1( jj:nfinst-1 ,-2:jj-1     )=smin( jj+1:nfinst ,  -2:jj-1   )
    smin1( -2:jj-1     ,jj:nfinst-1 )=smin(   -2:jj-1   ,jj+1:nfinst )
    smin1( jj:nfinst-1 ,jj:nfinst-1 )=smin( jj+1:nfinst ,jj+1:nfinst )
!    call printsmin( mdl,tmpproc ,smin1,nfinst-1 ) !DEBUG
    call kaleu_updt_smin( obj%kal(iproc) ,smin1 )
  enddo

! Remove dipoles with non-existing underlying processes
  iproc = obj%nproc
  do while (iproc.ge.1)
    if (notexists(iproc)) then
      l1 = obj%ndip
      do while (l1.ge.1)
        if (obj%proc(l1).eq.iproc) then
          obj%ndip = obj%ndip-1
          obj%i(   l1:obj%ndip) = obj%i(   l1+1:obj%ndip+1)
          obj%j(   l1:obj%ndip) = obj%j(   l1+1:obj%ndip+1)
          obj%k(   l1:obj%ndip) = obj%k(   l1+1:obj%ndip+1)
          obj%proc(l1:obj%ndip) = obj%proc(l1+1:obj%ndip+1)
        endif
        l1 = l1-1
      enddo
      obj%nproc = obj%nproc-1
      obj%kal(iproc:obj%nproc) = obj%kal(iproc+1:obj%nproc+1)
      do l1=1,obj%ndip
        if (obj%proc(l1).gt.iproc) obj%proc(l1) = obj%proc(l1)-1
      enddo
    endif
    iproc = iproc-1
  enddo
!
! Initialize multi-channeling
  call mch_init( obj%mch ,0,obj%ndip ,nbatch,nstep,thrs )
!
! Initialize splvars and grids
  nbatch_g = def_nbatch_g
  if (present(nbatch_grid)) nbatch_g = nbatch_grid
  if (present(noffset)) obj%noffset = noffset
  if (obj%noffset.eq.0) obj%noffset = -1
  include 'dipol_splvar_init.h'
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
  subroutine dipol_init_strf( str,obj )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type)  ,intent(inout) :: str
  type(dipol_type) ,intent(in)    :: obj
  real(kind(1d0)) :: xmin2,xmin2min
  integer         :: ii
  xmin2min = kaleu_get_xmin2( obj%kal(0) )
  do ii=1,obj%nproc
    xmin2 = kaleu_get_xmin2( obj%kal(ii) )
    if (xmin2.lt.xmin2min) xmin2min = xmin2
  enddo
  call strf_init( str ,dsqrt(xmin2min) )
  end subroutine
!
!
  subroutine dipol_gnrt_strf( str,obj ,discard ,x1kaleu,x2kaleu )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type)  ,intent(inout) :: str
  type(dipol_type) ,intent(inout) :: obj
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(out)   :: x1kaleu,x2kaleu
  real(kind(1d0)) :: ehat,stot,s1,s2
!
  call strf_gnrt( str ,x1kaleu,x2kaleu )
  ehat = dsqrt(x1kaleu*x2kaleu)*kaleu_get_ecm( obj%kal(0) )
! Construct initial-state momenta
  stot = ehat*ehat
  s1 = kaleu_get_s( obj%kal(0) ,-1 )
  s2 = kaleu_get_s( obj%kal(0) ,-2 )
  obj%p(0,-1) = -( stot + s1 - s2 ) / (2*ehat)
  obj%p(3,-1) = -dsqrt( obj%p(0,-1)**2 - s1 )
  obj%p(2,-1) =  0d0
  obj%p(1,-1) =  0d0
  obj%p(0,-2) = -( stot + s2 - s1 ) / (2*ehat)
  obj%p(3,-2) = -obj%p(3,-1)
  obj%p(2,-2) = -obj%p(2,-1)
  obj%p(1,-2) = -obj%p(1,-1)
  end subroutine
!
!
  subroutine dipol_inst( obj ,pkaleu )
!*********************************************************************
!*********************************************************************
  implicit none
  type(dipol_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(in)    :: pkaleu(0:3,-2:17)
  obj%p(0:3,-2:-1) = pkaleu(0:3,-2:-1)
  end subroutine
!
!
  subroutine dipol_gnrt( mdl,obj ,discard ,pkaleu )
!*********************************************************************
!*********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(dipol_type) ,intent(inout) :: obj
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(inout) :: pkaleu(0:3,-2:17)
  real(kind(1d0)) :: si,sj,sk,mi,mj,mk,xx,zz
  real(kind(1d0)) :: gridx,gridz
  integer         :: idip,ii,jj,kk,nn,aa,bb,iproc
!
  call mch_gnrt( obj%mch ,idip )
  obj%idip = idip
!
  nn = kaleu_get_nfinst( obj%kal(0) )
!
  if (idip.eq.0) then
    call kaleu_inst( obj%kal(0) ,discard ,obj%p )
    if (discard) return
    call kaleu_gnrt( mdl,obj%kal(0) ,discard ,pkaleu )
    if (discard) return
  else
    iproc = obj%proc(idip)
    ii = obj%i(idip)
    jj = obj%j(idip)
    kk = obj%k(idip)
    include 'dipol_splvar_gnrt.h'
    if (ii.gt.0.and.kk.gt.0) then
! Final-final
      call kaleu_inst( obj%kal(iproc) ,discard ,obj%p )
      if (discard) return
      call kaleu_gnrt( mdl,obj%kal(iproc) ,discard ,pkaleu )
      if (discard) return
      pkaleu(0:3,jj+1:nn) = pkaleu(0:3,jj:nn-1)
      si = kaleu_get_s( obj%kal(0) ,ii )
      sj = kaleu_get_s( obj%kal(0) ,jj )
      sk = kaleu_get_s( obj%kal(0) ,kk )
      mi = kaleu_get_m( obj%kal(0) ,ii )
      mj = kaleu_get_m( obj%kal(0) ,jj )
      mk = kaleu_get_m( obj%kal(0) ,kk )
      xx = xx*obj%aff
      call gnrt_ff( discard ,pkaleu ,xx,zz ,ii,jj,kk ,si,sj,sk,mi,mj,mk )
!      call printmom( pkaleu ,nn ) !DEBUG
      if (discard) return
    elseif ( ii.gt.0.and.kk.lt.0 .or. ii.lt.0.and.kk.gt.0 ) then
! Final-initial or initial-final
      aa = min(ii,kk)
      ii = max(ii,kk) ! ii possibly changed
      if (kk.lt.0) then
! ! Final-initial
        xx = xx*obj%afi
      else
! ! Initial-final
        zz = zz*obj%aif
      endif
      pkaleu(0:3,-2:-1) = obj%p(0:3,-2:-1)
      pkaleu(0:3,aa) = (1d0-xx)*pkaleu(0:3,aa)
      call kaleu_inst( obj%kal(iproc) ,discard ,pkaleu )
      if (discard) return
      call kaleu_gnrt( mdl,obj%kal(iproc) ,discard ,pkaleu )
      if (discard) return
      pkaleu(0:3,aa) = obj%p(0:3,aa)
      pkaleu(0:3,jj+1:nn) = pkaleu(0:3,jj:nn-1)
      si = kaleu_get_s( obj%kal(0) ,ii )
      call gnrt_fi( discard ,pkaleu ,xx,zz ,ii,jj,aa ,si )
      if (discard) return
    else
! Initial-initial
      aa = ii
      bb = kk
      zz = zz * xx*obj%aii
      xx = 1d0-xx
      pkaleu(0:3,-2:-1) = obj%p(0:3,-2:-1)
      pkaleu(0:3,aa) = xx*pkaleu(0:3,aa)
      call kaleu_inst( obj%kal(iproc) ,discard ,pkaleu )
      if (discard) return
      call kaleu_gnrt( mdl,obj%kal(iproc) ,discard ,pkaleu )
      if (discard) return
      pkaleu(0:3,aa) = obj%p(0:3,aa)
      pkaleu(0:3,jj+1:nn) = pkaleu(0:3,jj:nn-1)
      call gnrt_ii( discard ,pkaleu ,xx,zz ,aa,jj,bb,nn )
      if (discard) return
    endif
  endif
  obj%p(0:3,-2:nn) = pkaleu(0:3,-2:nn)
!  write(6,*) 'dipol_gnrt b',idip !DEBUG
!  call printmom( obj%p ,nn ) !DEBUG
  end subroutine 
!
!
  subroutine dipol_wght( mdl,obj ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(dipol_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(out)   :: weight
  real(kind(1d0)) :: mijsq,si,sj,sk,mi,mj,mk,xx,zz,dnstot,w1,w2,w3
  real(kind(1d0)) :: gridx,gridz,wgx,wgz,wsx,wsz,wx,wz
  real(kind(1d0)) :: pp(0:3,-2:17)
  integer         :: idip,ii,jj,kk,nn,ij,aa,iproc
!
  weight = 0d0
  dnstot = 0d0
!
  nn = kaleu_get_nfinst( obj%kal(0) )
!  write(6,*) 'dipol_wght',obj%idip !DEBUG
!
  do idip=0,obj%ndip
!
    if (obj%mch%wch(idip).eq.0d0) cycle
!  
    pp(0:3,-2:nn) = obj%p(0:3,-2:nn)
    obj%mch%dns(idip) = 0d0
!  
    if (idip.eq.0) then
      call kaleu_put_mom( obj%kal(0) ,pp )
      call kaleu_wght( mdl,obj%kal(0) ,w1 )
      if (w1.eq.0d0) cycle
      obj%mch%dns(idip) = 1d0/w1
    else
!      call printmom( pp ,nn ) !DEBUG
      iproc = obj%proc(idip)
      ii = obj%i(idip)
      jj = obj%j(idip)
      kk = obj%k(idip)
      if (ii.gt.0.and.kk.gt.0) then
!   Final-final
        ij = ii
        if (ij.gt.jj) ij = ij-1
        si    = kaleu_get_s( obj%kal(0) ,ii )
        sj    = kaleu_get_s( obj%kal(0) ,jj )
        sk    = kaleu_get_s( obj%kal(0) ,kk )
        mi    = kaleu_get_m( obj%kal(0) ,ii )
        mj    = kaleu_get_m( obj%kal(0) ,jj )
        mk    = kaleu_get_m( obj%kal(0) ,kk )
        mijsq = kaleu_get_s( obj%kal(iproc) ,ij )
        call wght_ff( w1 ,pp ,xx,zz ,ii,jj,kk ,mijsq,si,sj,sk,mi,mj,mk )
        if (w1.eq.0d0) cycle
        xx = xx/obj%aff
        if (xx.le.0d0.or.1d0.le.xx) cycle
        if (zz.le.0d0.or.1d0.le.zz) cycle
        w2 = obj%aff
      elseif ( ii.gt.0.and.kk.lt.0 .or. ii.lt.0.and.kk.gt.0 ) then
!   Final-initial or initial-final
        aa = min(ii,kk)
        ii = max(ii,kk) ! ii possibly changed
        si = kaleu_get_s( obj%kal(0) ,ii )
        call wght_fi( w1 ,pp ,xx,zz ,ii,jj,aa ,si )
        if (w1.eq.0d0) cycle
        if (kk.lt.0) then
! !  Final-initial
          xx = xx/obj%afi
          if (xx.le.0d0.or.1d0.le.xx) cycle
          if (zz.le.0d0.or.1d0.le.zz) cycle
          w2 = obj%afi
        else
! !  Initial-final
          zz = zz/obj%aif
          if (xx.le.0d0.or.1d0.le.xx) cycle
          if (zz.le.0d0.or.1d0.le.zz) cycle
          w2 = obj%aif
        endif
      else
!   Initial-initial
        call wght_ii( w1 ,pp ,xx,zz ,ii,jj,kk,nn )
        if (w1.eq.0d0) cycle
        xx = 1d0-xx
        if (xx.le.0d0.or.1d0.le.xx) cycle
        zz = zz/(xx*obj%aii)
        if (zz.le.0d0.or.1d0.le.zz) cycle
        w2 = obj%aii
      endif
      pp(0:3,jj:nn-1) = pp(0:3,jj+1:nn)
      call kaleu_put_mom( obj%kal(iproc) ,pp )
      call kaleu_wght( mdl,obj%kal(iproc) ,w3 )
      if (w3.eq.0d0) cycle
      include 'dipol_splvar_wght.h'
      obj%mch%dns(idip) = 1d0/(w1*w2*w3*wx*wz)
    endif
!  
    dnstot = dnstot + obj%mch%wch(idip)*obj%mch%dns(idip)
!
  enddo
!
  if (dnstot.eq.0d0) return
  weight = 1d0/dnstot
!  if (.not.weight.le.0d0.and..not.weight.ge.0d0) then !DEBUG
!    write(6,*) 'dipol_wght',dnstot,obj%mch%wch(0)*obj%mch%dns(0) !DEBUG
!    stop !DEBUG
!  endif !DEBUG
!  write(6,*) 'tot',weight,dnstot !DEBUG
!  write(6,*) 'tot',weight/1d16,obj%mch%dns(0),obj%mch%dns(34) !DEBUG
  end subroutine 
!
!
  subroutine dipol_collect( obj ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(dipol_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(in)    :: weight
  integer :: ii
  logical :: cleanup ,active(0:obj%nproc)
!
  call mch_collect( obj%mch ,weight ,cleanup )
!
  if (weight.gt.0d0) obj%idat = obj%idat+1
  include 'dipol_splvar_collect.h'
!
  active(0:obj%nproc) = .false.
  do ii=0,obj%ndip
    active(obj%proc(ii)) = (active(obj%proc(ii)).or.obj%mch%wch(ii).gt.0d0)
  enddo
!
  do ii=0,obj%nproc
    if (active(ii)) call kaleu_collect( obj%kal(ii) ,weight )
  enddo
!
  end subroutine
!
!
  subroutine dipol_plotgrids( obj ,iunit )
!*********************************************************************
!*********************************************************************
  implicit none
  type(dipol_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: iunit
  integer :: ii
  do ii=1,obj%ndip
    include 'dipol_splvar_plot.h'
  enddo
  end subroutine
!
!
  subroutine dipol_printdata( obj ,iunit )
!*********************************************************************
!*********************************************************************
  implicit none
  type(dipol_type) ,intent(in) :: obj
  integer          ,intent(in) :: iunit
  if (iunit.le.0) return
  write(iunit,*) 'MESSAGE from Kaleu_dipol: idip was',obj%idip
  end subroutine
!
!
!*********************************************************************
!* From here private routines
!*********************************************************************
!
!
  subroutine gnrt_ff( discard ,pp ,xy,xz ,ii,jj,kk ,si,sj,sk,mi,mj,mk )
!*********************************************************************
!* In the following,  yy = (sijk-si-sj-sk)*y_usual
!*********************************************************************
  implicit none
  logical         ,intent(inout) :: discard
  real(kind(1d0)) ,intent(inout) :: pp(0:3,-2:17)
  integer         ,intent(in) :: ii,jj,kk
  real(kind(1d0)) ,intent(in) :: xy,xz,si,sj,sk,mi,mj,mk
  real(kind(1d0)) :: pkt(0:3),QQ(0:3),pi(0:3),pj(0:3),pk(0:3) &
                    ,sQ,mQ,hh,sij,ymin,ymax,yy,vk,vi,zmin,zdif,zz,phi &
                    ,pipj2,pipk2,pjpk2,abspj,abspk,cosjk,sinjk
!
  discard = .false.
!
  pkt(0:3) = pp(0:3,kk)
  QQ(0:3)  = pp(0:3,ii) + pkt(0:3)
!
  sQ = (QQ(0)+QQ(3))*(QQ(0)-QQ(3)) - QQ(1)**2 - QQ(2)**2
  if (sQ.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ff: sQ<0, discard event'
    discard = .true.
    return
  endif
  mQ = dsqrt(sQ)
  hh = sQ-si-sj-sk
  ymin = 2*mi*mj
  ymax = hh - 2*mk*(mQ-mk)
  yy = (ymax-ymin)*xy + ymin
  sij = yy+si+sj
  vk = (sQ-sij-sk)**2 - 4*sij*sk
  if (vk.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ff: vk<0, discard event'
    discard = .true.
    return
  endif
  vk = dsqrt( vk )/( hh-yy )
  vi = yy**2 - 4*si*sj
  if (vi.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ff: vi<0, discard event'
    discard = .true.
    return
  endif
  vi = dsqrt( vi )/( 2*si+yy )
  zmin = (2*si+yy)/(2*(si+sj+yy))
  zdif = zmin*vi*vk
  zmin = zmin-zdif
  zdif = 2*zdif
  zz = zdif*xz + zmin
  call avh_random( phi )
  phi = twopi*phi
!
  pipj2 = yy
  pipk2 = zz*(hh-yy)
  pjpk2 = (1d0-zz)*(hh-yy)
  pi(0) = (sQ-pjpk2+si-sj-sk)/(2*mQ)
  pj(0) = (sQ-pipk2-si+sj-sk)/(2*mQ)
  pk(0) = (sQ-pipj2-si-sj+sk)/(2*mQ)
  abspj = pj(0)*pj(0)-sj
  abspk = pk(0)*pk(0)-sk
  if (abspj.lt.0d0.or.abspk.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ff: abspj or abspk <0, discard event'
    discard = .true.
    return
  endif
  abspj = dsqrt(abspj)
  abspk = dsqrt(abspk)
  cosjk = (2*pj(0)*pk(0)-pjpk2)/(2*abspk)
  sinjk = (abspj+cosjk)*(abspj-cosjk)
  if (sinjk.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ff: |cosjk| >1, discard event'
    discard = .true.
    return
  endif
  sinjk = dsqrt(sinjk)
  pj(1) = sinjk*dcos(phi)
  pj(2) = sinjk*dsin(phi)
  pj(3) = cosjk
  pk(1) = 0d0
  pk(2) = 0d0
  pk(3) = abspk
  pi(1:3) = -pk(1:3) - pj(1:3)
  call tsoob( pkt ,QQ,mQ )
  call rot3( pi ,pkt ,.true. )
  call rot3( pj ,pkt ,.false.)
  call rot3( pk ,pkt ,.false.)
  call boost( pi ,QQ,mQ )
  call boost( pj ,QQ,mQ )
  call boost( pk ,QQ,mQ )
!
  pp(0:3,ii) = pi(0:3)
  pp(0:3,jj) = pj(0:3)
  pp(0:3,kk) = pk(0:3)
!
  end subroutine

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
      'ERROR in kaleu_dipol in wght_ff: sQ<0, putting weight=0'
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
      'ERROR in kaleu_dipol in wght_ff: lam<0, putting weight=0'
    return
  endif
  vk = (sQ-sij-sk)**2 - 4*sij*sk
  if (vk.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in wght_ff: vk<0, putting weight=0'
    return
  endif
  vi = yy**2 - 4*si*sj
  if (vi.le.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in wght_ff: vi<0, putting weight=0'
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


  subroutine gnrt_fi( discard ,pp ,xx,xz ,ii,jj,aa ,si )
!********************************************************************
!* Assumes pp(0:3,aa) and pp(0:3,jj) are massless
!********************************************************************
  implicit none
  logical         ,intent(inout) :: discard
  real(kind(1d0)) ,intent(inout) :: pp(0:3,-2:17)
  integer         ,intent(in) :: ii,jj,aa
  real(kind(1d0)) ,intent(in) :: xx,xz ,si
  real(kind(1d0)) :: pa(0:3),QQ(0:3),pi(0:3),pj(0:3) &
            ,zz,pijpa2,sQ,mQ,zmin,phi,abspi,pipa2,cosia,sinia
!
  discard = .false.
!
  pa(0:3) = -pp(0:3,aa)
  pijpa2 = pp(0,ii)*pa(0) &   ! positive pijpa2, pipa2
         - pp(1,ii)*pa(1) &
         - pp(2,ii)*pa(2) &
         - pp(3,ii)*pa(3)
  pijpa2 = pijpa2*2
  QQ(0:3) = pp(0:3,ii)+xx*pa(0:3)  ! Q = Q_usual - pa
  sQ = (QQ(0)+QQ(3))*(QQ(0)-QQ(3)) - QQ(1)**2 - QQ(2)**2
  if (sQ.lt.0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_fi: sQ<0, discard event'
    discard = .true.
!    write(nunit,*) ii,jj,aa    !DEBUG
!    write(nunit,*) xx,xz       !DEBUG
!    write(nunit,*) pa          !DEBUG
!    write(nunit,*) pp(0:3,ii)  !DEBUG
!    write(nunit,*) QQ          !DEBUG
    return
  endif
  mQ = dsqrt(sQ)
!
  zmin = si/pijpa2 ! positive pijpa2, pipa2
  zmin = zmin/(xx+zmin)
!  if (zz.le.zmin.or.1d0.le.zz) then
!    discard = .true.
!    return
!  endif
  zz = (1d0-zmin)*xz + zmin
  call avh_random( phi )
  phi = twopi*phi
!
  call tsoob( pa ,QQ,mQ )
  pipa2 = zz*pijpa2 ! positive pijpa2, pipa2
  pi(0) = (sQ+si)/(2*mQ)
  abspi = pi(0)*pi(0) - si
  if (abspi.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_fi: abspi<0, discard event'
    discard = .true.
    return
  endif
  abspi = dsqrt(abspi)
  cosia = (2*pi(0)*pa(0)-pipa2)/(2*pa(0)) ! positive pijpa2, pipa2
  sinia = (abspi+cosia)*(abspi-cosia)
  if (sinia.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_fi: |cosia| >1, discard event'
    discard = .true.
    return
  endif
  sinia = dsqrt(sinia)
  pi(1) = sinia*dcos(phi)
  pi(2) = sinia*dsin(phi)
  pi(3) = cosia
  pj(0) = abspi
  pj(1:3) = -pi(1:3)
  call rot3( pi ,pa ,.true. )
  call rot3( pj ,pa ,.false.)
  call boost( pi ,QQ,mQ )
  call boost( pj ,QQ,mQ )
!
  pp(0:3,ii) = pi(0:3)
  pp(0:3,jj) = pj(0:3)
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


  subroutine gnrt_ii( discard ,pp ,xx,vv ,aa,jj,bb,nn )
!********************************************************************
!* Assumes pp(0:3,aa,jj,bb) is massless
!********************************************************************
  implicit none
  logical         ,intent(inout) :: discard
  real(kind(1d0)) ,intent(inout) :: pp(0:3,-2:17)
  integer         ,intent(in) :: aa,jj,bb,nn
  real(kind(1d0)) ,intent(in) :: xx,vv
  real(kind(1d0)) :: QQ(0:3),pa(0:3),pi(0:3),KK(0:3),Kt(0:3),Ks(0:3) &
                    ,phi,sQ,mQ,pipa2,Ea,Ei,cosia,sinia,sKt,sKs,cs,ct
  integer :: ll
!
  discard = .false.
!
  call avh_random( phi )
  phi = twopi*phi
!
  QQ(0:3) = -pp(0:3,aa) - pp(0:3,bb)
  sQ = (QQ(0)+QQ(3))*(QQ(0)-QQ(3)) - QQ(1)**2 - QQ(2)**2
  if (sQ.lt.0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ii: sQ<0, discard event'
    discard = .true.
    return
  endif
  mQ = dsqrt(sQ)
  pipa2 = sQ*vv    ! positive pipa, pa(0)
  pa(0) = mQ/2 ! positive pipa, pa(0)
  pi(0) = (1d0-xx)*sQ/(4*pa(0)) ! ( (xpa+pb)^2-2papb )/(4pa(0))
  cosia = pi(0) - pipa2/(2*pa(0))
  sinia = (pi(0)+cosia)*(pi(0)-cosia)
  if (sinia.lt.0d0) then
    if (nunit.gt.0) write(nunit,*) &
      'ERROR in kaleu_dipol in gnrt_ii: |cosia| >1, discard event'
    discard = .true.
    return
  endif
  sinia = dsqrt(sinia)
  pi(1) = sinia*dcos(phi)
  pi(2) = sinia*dsin(phi)
  pi(3) = cosia
! "pi(0:3)" is momentum in restframe of "Q" with "Lpa=L(Q->rest)pa"
! positive along the z-axis. Notice that "Q" and "pa" are already along
! the z-axis, so the rotating "Lpa" out of the z-axis is 1 or -1.
  if (pp(3,aa).gt.0d0) pi(1:3) = -pi(1:3) ! pa = -pp(,aa) positive energy
  call boost( pi ,QQ,mQ )
!
  KK(0:3) = QQ(0:3) - pi(0:3)
  Kt(0:3) = -xx*pp(0:3,aa) - pp(0:3,bb)
  Ks(0:3) = KK(0:3) + Kt(0:3)
  sKt = xx*sQ
  sKs = (Ks(0)+Ks(3))*(Ks(0)-Ks(3)) - Ks(1)**2 - Ks(2)**2
!
  do ll=1,nn
    if (ll.eq.jj) then
      pp(0:3,jj) = pi(0:3)
    else
      cs = Ks(0)*pp(0,ll)-Ks(1)*pp(1,ll)-Ks(2)*pp(2,ll)-Ks(3)*pp(3,ll)
      ct = Kt(0)*pp(0,ll)-Kt(1)*pp(1,ll)-Kt(2)*pp(2,ll)-Kt(3)*pp(3,ll)
      cs = 2*cs/sKs
      ct = 2*ct/sKt
      pp(0:3,ll) = pp(0:3,ll) - cs*Ks(0:3) + ct*KK(0:3)
    endif
  enddo
!
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
  subroutine mch_init( obj ,i0,i1 ,ndat,nstp,thrs )
!********************************************************************
!********************************************************************
  implicit none
  type(mch_type)  ,intent(inout) :: obj
  integer         ,intent(in)    :: i0,i1,ndat,nstp
  real(kind(1d0)) ,intent(in)    :: thrs
  integer ii
!
  obj%yes = (nstp.gt.0)
  obj%i0 = min(i0,i1)
  obj%i1 = max(i0,i1)
!
  if (.not.obj%yes) return
!
  allocate( obj%wch(obj%i0:obj%i1) )
  allocate( obj%ave(obj%i0:obj%i1) )
  allocate( obj%dns(obj%i0:obj%i1) )
!
  obj%idat = 0
!DEBUG  obj%ndat = 1000000 !DEBUG
  obj%ndat = ndat
  obj%ntot = 0
  obj%istp = 0
  obj%nstp = nstp
  obj%thrs = thrs
  obj%ave = 0d0
  obj%wch = 1d0
  obj%dns = 0d0
  obj%wch = obj%wch/sum( obj%wch )
!  obj%wch = 0d0     !DEBUG
!  obj%wch( 0) = 1d0 !DEBUG                g g  -> t t b b g 
!  obj%wch( 1) = 1d0 !DEBUG -1g  3b  -2g   g b  -> t t b g 
!  obj%wch( 2) = 1d0 !DEBUG -1g  3b   1t   g b  -> t t b g 
!  obj%wch( 3) = 1d0 !DEBUG -1g  3b   2t   g b  -> t t b g 
!  obj%wch( 4) = 1d0 !DEBUG -1g  3b   4b   g b  -> t t b g 
!  obj%wch( 5) = 1d0 !DEBUG -1g  3b   5g   g b  -> t t b g 
!!!  obj%wch( 6) = 1d0 !DEBUG -1g  4b  -2g   g b  -> t t b g 
!!!  obj%wch( 7) = 1d0 !DEBUG -1g  4b   1t   g b  -> t t b g 
!!!  obj%wch( 8) = 1d0 !DEBUG -1g  4b   2t   g b  -> t t b g 
!!!  obj%wch( 9) = 1d0 !DEBUG -1g  4b   3b   g b  -> t t b g 
!!!  obj%wch(10) = 1d0 !DEBUG -1g  4b   5g   g b  -> t t b g 
!!  obj%wch(11) = 1d0 !DEBUG -1g  5g  -2g   g g  -> t t b b 
!!  obj%wch(12) = 1d0 !DEBUG -1g  5g   1t   g g  -> t t b b 
!!  obj%wch(13) = 1d0 !DEBUG -1g  5g   2t   g g  -> t t b b 
!!  obj%wch(14) = 1d0 !DEBUG -1g  5g   3b   g g  -> t t b b 
!!  obj%wch(15) = 1d0 !DEBUG -1g  5g   4b   g g  -> t t b b 
!  obj%wch(16) = 1d0 !DEBUG -2g  3b  -1g   b g  -> t t b g 
!  obj%wch(17) = 1d0 !DEBUG -2g  3b   1t   b g  -> t t b g 
!  obj%wch(18) = 1d0 !DEBUG -2g  3b   2t   b g  -> t t b g 
!  obj%wch(19) = 1d0 !DEBUG -2g  3b   4b   b g  -> t t b g 
!  obj%wch(20) = 1d0 !DEBUG -2g  3b   5g   b g  -> t t b g 
!!!  obj%wch(21) = 1d0 !DEBUG -2g  4b  -1g   b g  -> t t b g 
!!!  obj%wch(22) = 1d0 !DEBUG -2g  4b   1t   b g  -> t t b g 
!!!  obj%wch(23) = 1d0 !DEBUG -2g  4b   2t   b g  -> t t b g 
!!!  obj%wch(24) = 1d0 !DEBUG -2g  4b   3b   b g  -> t t b g 
!!!  obj%wch(25) = 1d0 !DEBUG -2g  4b   5g   b g  -> t t b g 
!!  obj%wch(26) = 1d0 !DEBUG -2g  5g  -1g   g g  -> t t b b 
!!  obj%wch(27) = 1d0 !DEBUG -2g  5g   1t   g g  -> t t b b 
!!  obj%wch(28) = 1d0 !DEBUG -2g  5g   2t   g g  -> t t b b 
!!  obj%wch(29) = 1d0 !DEBUG -2g  5g   3b   g g  -> t t b b 
!!  obj%wch(30) = 1d0 !DEBUG -2g  5g   4b   g g  -> t t b b 
!!  obj%wch(31) = 1d0 !DEBUG  1t  5g  -1g   g g  -> t t b b 
!!  obj%wch(32) = 1d0 !DEBUG  1t  5g  -2g   g g  -> t t b b 
!!  obj%wch(33) = 1d0 !DEBUG  1t  5g   2t   g g  -> t t b b 
!!  obj%wch(34) = 1d0 !DEBUG  1t  5g   3b   g g  -> t t b b 
!!  obj%wch(35) = 1d0 !DEBUG  1t  5g   4b   g g  -> t t b b 
!!  obj%wch(36) = 1d0 !DEBUG  2t  5g  -1g   g g  -> t t b b 
!!  obj%wch(37) = 1d0 !DEBUG  2t  5g  -2g   g g  -> t t b b 
!!  obj%wch(38) = 1d0 !DEBUG  2t  5g   1t   g g  -> t t b b 
!!  obj%wch(39) = 1d0 !DEBUG  2t  5g   3b   g g  -> t t b b 
!!  obj%wch(40) = 1d0 !DEBUG  2t  5g   4b   g g  -> t t b b 
!!!  obj%wch(41) = 1d0 !DEBUG  3b  4b  -1g   g g  -> t t g g 
!!!  obj%wch(42) = 1d0 !DEBUG  3b  4b  -2g   g g  -> t t g g 
!!!  obj%wch(43) = 1d0 !DEBUG  3b  4b   1t   g g  -> t t g g 
!!!  obj%wch(44) = 1d0 !DEBUG  3b  4b   2t   g g  -> t t g g 
!!!  obj%wch(45) = 1d0 !DEBUG  3b  4b   5g   g g  -> t t g g 
!!  obj%wch(46) = 1d0 !DEBUG  3b  5g  -1g   g g  -> t t b b 
!!  obj%wch(47) = 1d0 !DEBUG  3b  5g  -2g   g g  -> t t b b 
!!  obj%wch(48) = 1d0 !DEBUG  3b  5g   1t   g g  -> t t b b 
!!  obj%wch(49) = 1d0 !DEBUG  3b  5g   2t   g g  -> t t b b 
!!  obj%wch(50) = 1d0 !DEBUG  3b  5g   4b   g g  -> t t b b 
!!  obj%wch(51) = 1d0 !DEBUG  4b  5g  -1g   g g  -> t t b b 
!!  obj%wch(52) = 1d0 !DEBUG  4b  5g  -2g   g g  -> t t b b 
!!  obj%wch(53) = 1d0 !DEBUG  4b  5g   1t   g g  -> t t b b 
!!  obj%wch(54) = 1d0 !DEBUG  4b  5g   2t   g g  -> t t b b 
!!  obj%wch(55) = 1d0 !DEBUG  4b  5g   3b   g g  -> t t b b 
!  obj%wch = obj%wch/sum( obj%wch ) !DEBUG
!  do ii=obj%i0,obj%i1                            !DEBUG
!    if (obj%wch(ii).eq.0d0) cycle                !DEBUG
!    write(6,*) 'wch_dipol(',ii,') =',obj%wch(ii) !DEBUG
!  enddo                                          !DEBUG
  end subroutine
!
!
  subroutine mch_close( obj )
!********************************************************************
!********************************************************************
  implicit none
  type(mch_type) ,intent(inout) :: obj
  obj%yes = .false.
  obj%idat = 0
  obj%ndat = 0 
  obj%ntot = 0
  obj%istp = 0
  obj%nstp = 0
  obj%thrs = 0d0
  if (allocated(obj%wch)) deallocate( obj%wch )
  if (allocated(obj%ave)) deallocate( obj%ave )
  if (allocated(obj%dns)) deallocate( obj%dns )
  end subroutine
!
!
  subroutine mch_gnrt( obj ,ii )
!********************************************************************
!********************************************************************
  implicit none
  type(mch_type) ,intent(inout) :: obj
  integer        ,intent(out)   :: ii
  real(kind(1d0)) :: xx,sumw
  call avh_random( xx )
  sumw = 0d0
  ii = obj%i0-1
  if (obj%yes) then
    do while ( xx.gt.sumw .and. ii.lt.obj%i1 )
      ii = ii+1
      sumw = sumw + obj%wch(ii)
    enddo
    if (xx.gt.sumw) then
      write(*,*) 'ERROR in dipol_mch_gnrt: no channel chosen'
      stop
    endif
  else
    ii = 1 + int( (obj%i1-obj%i0+1)*xx )
  endif  
  end subroutine
!
!
  subroutine mch_collect( obj ,wght ,cleanup )
!********************************************************************
!********************************************************************
  implicit none
  type(mch_type)  ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)    :: wght
  logical         ,intent(out)   :: cleanup
  integer :: ii
  integer :: nrem,jrem(obj%i1-obj%i0+1)
  real(kind(1d0)) :: factor,hsum,thrs
!
  cleanup = .false.
  if (.not.obj%yes) return
  if (obj%istp.ge.obj%nstp) return
! gather data
  if (wght.ne.0d0) then
    factor = dot_product(obj%wch,obj%dns)    
    if (factor.ne.0d0) factor = wght*wght/factor ! variance
!   if (factor.ne.0d0) factor = wght/factor      ! entropy
    obj%ave = obj%ave + factor*obj%dns
    obj%dns = 0d0
    obj%idat = obj%idat + 1
  endif
  obj%ntot = obj%ntot + 1
  if (obj%idat.eq.obj%ndat) then
    obj%ave = dsqrt( obj%ave/dble(obj%ntot) ) ! variance
!   obj%ave =        obj%ave/dble(obj%ntot)   ! entropy
    if (sum(obj%ave).eq.0d0) return
    obj%wch = obj%wch*obj%ave
    obj%ave = 0d0
    hsum = sum(obj%wch)
    obj%wch = obj%wch/hsum
    if (obj%istp.eq.obj%nstp-1) then
      cleanup = .true.
      thrs = obj%thrs/dble(obj%i1-obj%i0+1)
      where (obj%wch.lt.thrs) obj%wch = 0d0
      hsum = sum(obj%wch)
      obj%wch = obj%wch/hsum
    endif
    obj%idat = 0
    obj%ntot = 0
    obj%istp = obj%istp+1
!    do ii=obj%i0,obj%i1                                   !DEBUG
!      if (obj%wch(ii).eq.0d0) cycle                       !DEBUG
!      write(6,'(a37,i3,a3,d16.8)') &                      !DEBUG
!       ' MESSAGE from kaleu_dipol: wch_dipol(',ii,') =' & !DEBUG
!                                  ,obj%wch(ii)            !DEBUG
!    enddo                                                 !DEBUG
  endif !(obj%idat.eq.obj%ndat)
  end subroutine
!
!
  subroutine printsmin( mdl,process ,smin,nn )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type) ,intent(in) :: mdl
  integer          ,intent(in) :: process(-2:17)
  real(kind(1d0))  ,intent(in) :: smin(-2:17,-2:17)
  integer          ,intent(in) :: nn
  character(2) :: ai,aj
  integer      :: ii,jj
  do ii=-2,nn-1                         
    if (ii.eq.0) cycle                      
    do jj=ii+1,nn                       
      if (jj.le.0) cycle
      ai = sparticle( mdl,process(ii) )
      aj = sparticle( mdl,process(jj) )
      write(nunit,'(a7,i2,a2,1x,i1,a2,1x,f13.4)') &
        ' smin: ',ii,ai,jj,aj,smin(ii,jj)
    enddo                                   
  enddo                                     
  end subroutine
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
end module
