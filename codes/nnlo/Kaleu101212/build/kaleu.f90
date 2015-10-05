module avh_kaleu
  use avh_kaleu_model
  use avh_kaleu_tree
  use avh_kaleu_ranvar
  use avh_kaleu_mulcha
  use avh_kaleu_kinem
  use avh_kaleu_vertex
  use avh_kaleu_strf
!
  private ! everything is private, except the following list
  public :: kaleu_hello,kaleu_type &
           ,kaleu_put_process &
           ,kaleu_updt_smin,kaleu_print_smin &
           ,kaleu_init_adapt &
           ,kaleu_init_strf &
           ,kaleu_updt_strf &
           ,kaleu_updt_strf_nlo &
           ,kaleu_gnrt_strf &
           ,kaleu_gnrt_strf_lab &
           ,kaleu_inst &
           ,kaleu_gnrt &
           ,kaleu_put_mom &
           ,kaleu_wght &
           ,kaleu_collect &
           ,kaleu_close &
           ,kaleu_closeall &
           ,kaleu_plotgrids &
           ,kaleu_get_ecm &
           ,kaleu_get_m &
           ,kaleu_get_s &
           ,kaleu_get_nfinst &
           ,kaleu_get_xmin2 &
           ,kaleu_get_p
!
! The unit to which messages are send
  integer ,private ,parameter :: nunit = 6
!
  type :: kaleu_type
    private
    type(kinem_type)  :: kin
    type(mulcha_type) :: mch
    type(vertex_type) ,allocatable :: vtx(:)
  end type
!
contains
!
!
  subroutine kaleu_hello
!**********************************************************************
!**********************************************************************
  implicit none
  logical :: init = .true.
  if (init) then
    init = .false.
    write(*,'(a72)') '########################################################################'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '#                         You are using Kaleu                          #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '#                     for phase space generation                       #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
    write(*,'(a72)') '#   date: 12-12-2010                                                   #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '# Please cite                                                          #'
    write(*,'(a72)') '#    A. van Hameren, arXiv:1003.4953 [hep-ph]                          #'
    write(*,'(a72)') '# in publications with results obtained with the help of this program. #'
    write(*,'(a72)') '#                                                                      #'
    write(*,'(a72)') '########################################################################'
  endif
  end subroutine
!
!
!  subroutine kaleu_init_adapt
  include 'kaleu_init_adapt.h'
!  end subroutine
!
!
  subroutine kaleu_put_process( mdl,mlv,obj ,process,nfinst ,ecm ,cancel )
!**********************************************************************
!**********************************************************************
  use avh_bint
  implicit none
  type(model_type) ,intent(inout) :: mdl
  type(vertx_type) ,intent(inout) :: mlv ! this is not vertex_type !
  type(kaleu_type) ,intent(inout) :: obj
  integer          ,intent(in)    :: process(-2:17),nfinst
  real(kind(1d0))  ,intent(in)    :: ecm
  logical ,optional,intent(out)   :: cancel
  real(kind(1d0)) :: pkaleu(0:3,-2:17),stot,m1,m2,s1,s2
  type(tree_type) :: tree
  integer :: ii
  logical :: discard,cncl
!
  call kaleu_hello
!
! Initialize some sizes
  call kinem_init( nfinst )
!
! Construct tree
  call tree_checkinput( mdl ,process, nfinst )
  call tree_dress( mdl,mlv ,tree ,process,nfinst ,0 )
  call tree_checkoutput( mdl ,tree ,process, nfinst ,cncl )
  if (present(cancel)) cancel = cncl
  if (cncl) return
!
! Initialize multi-channeling
  call mulcha_init( tree ,obj%mch )
!
! Allocate and intialize vertices
  allocate( obj%vtx(1:obj%mch%nv) )
  do ii=1,obj%mch%nv
    call vertex_init( tree,ii ,obj%vtx(ii) )
  enddo
!
! Put on-shell invariants
  call kinem_put_nfinst( obj%kin ,nfinst ,ecm )
  do ii=-2,nfinst
    if (ii.eq.0) cycle
    call kinem_put_mass( obj%kin ,mparticle(mdl,process(ii)) ,ii )
  enddo
!
! Put default initial-state momenta
  stot = ecm**2
  if (stot.lt.(obj%kin%m(-1)+obj%kin%m(-2))**2) then
    write(*,*) 'ERROR in avh_kaleu: ecm too small w.r.t. the ' ,&
               'initial-state masses'
    stop
  endif
  pkaleu(0,-1) = -( stot + obj%kin%s(-1) - obj%kin%s(-2) ) / (2*ecm)
  pkaleu(3,-1) = -dsqrt( pkaleu(0,-1)**2 - obj%kin%s(-1) )
  pkaleu(2,-1) =  0d0
  pkaleu(1,-1) =  0d0
  pkaleu(0,-2) = -( stot + obj%kin%s(-2) - obj%kin%s(-1) ) / (2*ecm)
  pkaleu(3,-2) = -pkaleu(3,-1)
  pkaleu(2,-2) = -pkaleu(2,-1)
  pkaleu(1,-2) = -pkaleu(1,-1)
  call kaleu_inst( obj ,discard ,pkaleu )
  if (discard) then
    write(*,*) 'ERROR in avh_kaleu: ' ,&
               'something wrong with initial-state momenta'
    stop
  endif
!
! Put kinematical minima on positive invariants
  call kinem_init_smin( obj%kin )
!
  end subroutine
!
!
  subroutine kaleu_updt_smin( obj ,smin )
!********************************************************************
!* Update lower kinematical limits on positive invariants, given a
!* matrix of lower limits on 2-particle invariants
!********************************************************************
  implicit none
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(in)    :: smin(-2:17,-2:17)
  call kinem_updt_smin( obj%kin ,smin )
  end subroutine
!
  subroutine kaleu_print_smin( mdl,obj ,process,nfinst ,ecm,smin ,nunit )
!********************************************************************
!********************************************************************
  use avh_print
  implicit none
  type(model_type) ,intent(in) :: mdl
  type(kaleu_type) ,intent(in) :: obj
  integer          ,intent(in) :: process(-2:17),nfinst,nunit
  real(kind(1d0))  ,intent(in) :: ecm,smin(-2:17,-2:17)
  integer :: ii,jj
  if (nunit.le.0) return
  write(nunit,*) 'MESSAGE from avh_kaleu: ecm=',printdbl(0,8,ecm)
  do ii=-2,nfinst
    if (ii.eq.0) cycle
    do jj=ii+1,nfinst
      if (jj.le.0) cycle
      write(nunit,*) 'MESSAGE from avh_kaleu: smin(' &
                    ,trim(printint(ii)),trim(sparticle(mdl,process(ii))) &
                ,',',trim(printint(jj)),trim(sparticle(mdl,process(jj))) &
                ,')=',printdbl(0,8,smin(ii,jj))
    enddo
  enddo
  end subroutine
!
!
  subroutine kaleu_init_strf( str,obj ,xmin )
!*********************************************************************
!* Let Kaleu do the integration of the structure functions
!* Input: xmin = sqrt(tau0)
!*             = the minimum value for sqrt(x1*x2), ie, the minimum
!*               value for  ehat/ecm = sqrt(shat/s)
!* The routine eventually chooses the largest of this  xmin  and the
!* one it determines using the final-state masses.
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout)  :: str
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)     :: xmin
  real(kind(1d0)) :: xx,smin(-2:17,-2:17),mass(obj%kin%b(obj%kin%n1))
! Check if the xmin due to final-state masses isn't actually larger
!  xx = obj%kin%smin( obj%kin%b(obj%kin%n1)-2 )/obj%kin%ecm**2
!  xx = max( xmin ,dsqrt(xx) )
  smin = 0d0
  call kinem_mass_smin( obj%kin ,mass ,smin )
  xx = mass( obj%kin%b(obj%kin%n1)-2 )/obj%kin%ecm
  xx = max( xmin ,xx )
  call strf_init( str ,xx )
  end subroutine
!
  subroutine kaleu_updt_strf( str,obj ,smin )
!*********************************************************************
!* Update xmin given a matrix of minima on 2-particle invariants
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout)  :: str
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)     :: smin(-2:17,-2:17)
  real(kind(1d0)) :: xx ,mass( obj%kin%b(obj%kin%n1) )
  call kinem_mass_smin( obj%kin ,mass ,smin )
  xx = mass( obj%kin%b(obj%kin%n1)-2 )/obj%kin%ecm
  call strf_updt_xmin( str ,xx )
  end subroutine
!
  subroutine kaleu_updt_strf_nlo( str,obj ,smin )
!*********************************************************************
!* Update xmin given a matrix of minima on 2-particle invariants.
!* Take the minimum xmin obtained by going through all cases where
!* one final-state row/collumn in smin is put to zero
!*********************************************************************
  implicit none
  type(strf_type) ,intent(inout)  :: str
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(in)     :: smin(-2:17,-2:17)
  real(kind(1d0)) :: xx ,mass( obj%kin%b(obj%kin%n1) )
  real(kind(1d0))  :: loc_smin(-2:17,-2:17),loc_xx
  integer          :: nfin,ii
!
  call kinem_mass_smin( obj%kin ,mass ,smin )
  xx = mass( obj%kin%b(obj%kin%n1)-2 )/obj%kin%ecm
!
  nfin = obj%kin%n1-1
  do ii=1,nfin
    loc_smin( -2:nfin ,-2:nfin ) = smin( -2:nfin ,-2:nfin )
    loc_smin(   ii    ,-2:nfin ) = 0d0
    loc_smin( -2:nfin ,  ii    ) = 0d0
    call kinem_mass_smin( obj%kin ,mass ,loc_smin )
    loc_xx = mass( obj%kin%b(obj%kin%n1)-2 )/obj%kin%ecm
    if (loc_xx.lt.xx) xx = loc_xx
  enddo
!
  call strf_updt_xmin( str ,xx )
  end subroutine
!
!
  subroutine kaleu_gnrt_strf( str,obj ,discard ,x1kaleu,x2kaleu )
!*********************************************************************
!*********************************************************************
  implicit none
  type(strf_type)  ,intent(inout) :: str
  type(kaleu_type) ,intent(inout) :: obj
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(out)   :: x1kaleu,x2kaleu
  real(kind(1d0)) :: ehat,stot,s1,s2
!
  discard = .false.
  call strf_gnrt( str ,x1kaleu,x2kaleu )
  ehat = dsqrt( x1kaleu*x2kaleu )*obj%kin%ecm
! Construct initial-state momenta
  stot = ehat*ehat
  s1 = obj%kin%s(-1)
  s2 = obj%kin%s(-2)
  obj%kin%p(0,-1) = -( stot + s1 - s2 ) / (2*ehat)
  obj%kin%p(3,-1) = -dsqrt( obj%kin%p(0,-1)**2 - s1 )
  obj%kin%p(2,-1) =  0d0
  obj%kin%p(1,-1) =  0d0
  obj%kin%p(0,-2) = -( stot + s2 - s1 ) / (2*ehat)
  obj%kin%p(3,-2) = -obj%kin%p(3,-1)
  obj%kin%p(2,-2) = -obj%kin%p(2,-1)
  obj%kin%p(1,-2) = -obj%kin%p(1,-1)
  end subroutine
!
!
  subroutine kaleu_gnrt_strf_lab( str,obj ,discard ,x1kaleu,x2kaleu )
!*********************************************************************
!* Only for massless initial-state momenta
!*********************************************************************
  implicit none
  type(strf_type)  ,intent(inout) :: str
  type(kaleu_type) ,intent(inout) :: obj
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(out)   :: x1kaleu,x2kaleu
  real(kind(1d0)) :: hecm,ee
!
  discard = .false.
  call strf_gnrt( str ,x1kaleu,x2kaleu )
  hecm = obj%kin%ecm/2d0
! Construct initial-state momenta
  ee = x1kaleu*hecm
  obj%kin%p(0,-1) = -ee
  obj%kin%p(3,-1) = -ee
  obj%kin%p(2,-1) =  0d0
  obj%kin%p(1,-1) =  0d0
  ee = x2kaleu*hecm
  obj%kin%p(0,-2) = -ee
  obj%kin%p(3,-2) =  ee
  obj%kin%p(2,-2) =  0d0
  obj%kin%p(1,-2) =  0d0
  end subroutine
!
!
  subroutine kaleu_inst( obj ,discard ,pkaleu )
!********************************************************************
!* Put initial-state momenta
!********************************************************************
  implicit none
  type(kaleu_type) ,intent(inout) :: obj
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(in)    :: pkaleu(0:3,-2:17)
  real(kind(1d0)) :: h1,h2
  discard = .false.
  h1 = dsign(1d0,pkaleu(0,-1))
  h2 = dsign(1d0,pkaleu(0,-2))
  if (h1.ne.h2) then
    if (nunit.gt.0) write(nunit,*) &
       'ERROR in kaleu_inst: energies of incoming momenta' ,&
       ' have different signs, discard event.'
    discard = .true.
    return
  endif
  if (h1.eq.1d0) then
    obj%kin%p(0:3,-2:-1) =-pkaleu(0:3,-2:-1)
  else 
    obj%kin%p(0:3,-2:-1) = pkaleu(0:3,-2:-1)
  endif
  end subroutine
!
!
  subroutine kaleu_gnrt( mdl,obj ,discard ,pkaleu )
!*********************************************************************
!*********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(kaleu_type) ,intent(inout) :: obj
  logical          ,intent(out)   :: discard
  real(kind(1d0))  ,intent(out)   :: pkaleu(0:3,-2:17)
  integer :: nn,ii,vrtx(obj%mch%n1)
!
  discard = .false.
!
! Put inst-momenta to work-arrays
  call kinem_inst( obj%kin ,discard )
  if (discard) return
!
! Generate a graph, ie, a list of vertices
  nn = 0
  ii = obj%mch%ov(0,obj%mch%nv)
  call mulcha_vrtx( obj%mch ,vrtx,nn ,ii ) ! recursive routine!
!
! Go through the vertices of the graph
  do ii=1,nn
    call vertex_gnrt( mdl ,obj%kin ,obj%vtx(vrtx(ii)) ,discard )
    if (discard) return
  enddo
!
! External momenta from work-arrays for output
  call kinem_result( obj%kin )
  pkaleu(0:3,-2:obj%kin%n1-1) = obj%kin%p(0:3,-2:obj%kin%n1-1)
!  do ii=-2,obj%kin%n1-1                      !DEBUG
!    if (ii.eq.0) cycle                       !DEBUG
!    write(6,'(i3,4e16.8)') ii,pkaleu(0:3,ii) !DEBUG
!  enddo                                      !DEBUG
!
  end subroutine
!
!
  subroutine kaleu_put_mom( obj ,pkaleu )
!*********************************************************************
!* Put final state momenta (optional)
!*********************************************************************
  implicit none
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(in)    :: pkaleu(0:3,-2:17)
  integer :: nfinst
  nfinst = obj%kin%n1-1
  obj%kin%p(0:3,-2:-1    ) = pkaleu(0:3,-2:-1)
  obj%kin%p(0:3, 1:nfinst) = pkaleu(0:3, 1:nfinst)
  end subroutine
!
!
  subroutine kaleu_wght( mdl,obj ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: mdl
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0)) ,intent(out)    :: weight
  real(kind(1d0)) :: dd(obj%mch%no),ww
  integer :: ii,i0,i1,i2,i3
  logical :: discard
  weight = 0d0
  call kinem_calcall( obj%kin ,discard )
  if (discard) return
  dd( 1:obj%mch%n1 ) = 1d0            ! obj%mch%no may be too large,
  dd( obj%mch%n1+1:obj%mch%no ) = 0d0 ! but not too small
  if (obj%mch%yes) then
    do ii=1,obj%mch%nv
      if (obj%mch%wtv(ii).eq.0d0) cycle
      obj%mch%dns(ii) = 0d0
      i0 = obj%mch%ov(0,ii)
      i1 = obj%mch%ov(1,ii)
      i2 = obj%mch%ov(2,ii)
      i3 = obj%mch%ov(3,ii)
      if (i3.ne.0) i1 = i3 ! specific T-channel
      call vertex_wght( mdl ,obj%kin ,obj%vtx(ii) ,ww )
      if (ww.eq.0d0) cycle !return
      obj%mch%dns(ii) = dd(i1)*dd(i2)/ww
      dd(i0) = dd(i0) + obj%mch%wtv(ii)*obj%mch%dns(ii)
    enddo
  else
    do ii=1,obj%mch%nv
      i0 = obj%mch%ov(0,ii)
      i1 = obj%mch%ov(1,ii)
      i2 = obj%mch%ov(2,ii)
      i3 = obj%mch%ov(3,ii)
      if (i3.ne.0) i1 = i3 ! specific T-channel
      call vertex_wght( mdl ,obj%kin ,obj%vtx(ii) ,ww )
      if (ww.eq.0d0) cycle !return
      dd(i0) = dd(i0) + dd(i1)*dd(i2)/ww/dble(obj%mch%vo(i0,0))
    enddo
  endif
  weight = dd( obj%mch%ov(0,obj%mch%nv) )
  if (weight.ne.0d0) weight = 1d0/weight
  end subroutine
!
!
  subroutine kaleu_collect( obj ,weight )
!*********************************************************************
!*********************************************************************
  implicit none
  type(kaleu_type) ,intent(inout) :: obj
  real(kind(1d0))  ,intent(in)    :: weight
  integer :: ii,nrem,jrem(obj%mch%nv)
!
  call mulcha_collect( obj%mch ,weight ,jrem,nrem )
!
  do ii=1,nrem
    call vertex_close( obj%vtx(jrem(ii)) )
  enddo
!
  if (weight.gt.0d0) then
    do ii=1,obj%mch%nv
      call vertex_collect( obj%vtx(ii) ,weight ,obj%mch%ndat,obj%mch%nstp )
    enddo
  endif
!
  end subroutine
!
!
  subroutine kaleu_close( obj )
!*********************************************************************
!*********************************************************************
  implicit none
  type(kaleu_type) ,intent(inout) :: obj
  integer :: ii
  do ii=1,obj%mch%nv
    call vertex_close( obj%vtx(ii) )
  enddo
  deallocate( obj%vtx )
  call kinem_close( obj%kin )
  call mulcha_close( obj%mch )
  end subroutine
!
  subroutine kaleu_closeall
!*********************************************************************
!*********************************************************************
  implicit none
  call kinem_closeall
  end subroutine
!
!
  subroutine kaleu_plotgrids( obj ,iunit )
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in) :: obj
  integer          ,intent(in) :: iunit
  integer :: ii
  do ii=1,obj%mch%nv
    call vertex_plotgrids( obj%vtx(ii) ,iunit )
  enddo
  end subroutine
!
!
!
  function kaleu_get_ecm( obj ) result(value)
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in) :: obj
  real(kind(1d0)) :: value
  value = obj%kin%ecm
  end function
!
  function kaleu_get_m( obj, ii ) result(value)
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in) :: obj
  integer          ,intent(in) :: ii
  real(kind(1d0)) :: value
  value = obj%kin%m( ii )
  end function
!
  function kaleu_get_s( obj, ii ) result(value)
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in) :: obj
  integer          ,intent(in) :: ii
  real(kind(1d0)) :: value
  value = obj%kin%s( ii )
  end function
!
  function kaleu_get_nfinst( obj ) result(value)
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in) :: obj
  integer :: value
  value = obj%kin%n1-1
  end function
!
  function kaleu_get_xmin2( obj ) result(value)
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in) :: obj
  real(kind(1d0)) :: value
  value = obj%kin%smin( obj%kin%b(obj%kin%n1)-2 )/obj%kin%ecm**2
  end function
!
  subroutine kaleu_get_p( obj ,pkaleu )
!*********************************************************************
!*********************************************************************
  type(kaleu_type) ,intent(in)  :: obj
  real(kind(1d0))  ,intent(out) :: pkaleu(0:3,-2:17)
  integer :: nfinst
  nfinst = obj%kin%n1-1
  pkaleu(0:3,-2:-1    ) = obj%kin%p(0:3,-2:-1)
  pkaleu(0:3, 1:nfinst) = obj%kin%p(0:3, 1:nfinst)
  end subroutine
!
!
end module
