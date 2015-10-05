module avh_kaleu_tree
  private ! Everything is private except the following list:
  public :: tree_type,tree_checkinput,tree_dress,tree_checkoutput ,&
            tree_nv,tree_copy
  public :: tree_check,tree_order,printtree
!
! Realize that the object with the sizes defined below only exists
! at initialization-stage, and is not stored. So don't be afraid to
! put large numbers.
!
! Maximal number of external particles minus 1 (= #final-state + 1)
  integer ,parameter :: size_n1 = 11
! Maximal number of vertices in a tree
  integer ,parameter :: size_nvx = 30000
! Maximal number of "off-shell currents"
  integer ,parameter :: size_noc = 2000
! Maximal number of vertices with the same 0-leg
  integer ,parameter :: size_nvo = 2000
! Unit to which messages are send
  integer ,parameter :: nunit = 6
! Seperate unit to print the whole tree
  integer ,parameter :: treeunit = 0
! Seperate unit for messages regarding array sizes
  integer ,parameter :: sizeunit = 0
!
  type :: tree_type
    integer :: n1 ! number of external particles minus 1
    integer :: nv ! number of vertices
    integer :: no ! number of o(ff)s(hell)c(urrents)     
    integer :: ov(0:3,size_nvx)           ! osc as function of vertex   
    integer :: po(0:size_noc)             ! momentum as function of osc 
    integer :: fo(0:size_noc)             ! flavor as function of osc   
    integer :: vo(0:size_noc,0:size_nvo)  ! vertices as function of osc 
    integer :: typ(size_nvx)              ! type of vertex              
    integer :: zax(size_nvx)              ! momentum defining z-axis    
    logical :: frst1(size_nvx)            ! related to possible ordering
  end type
!
!
contains
!
!
  function tree_nv( tree ) result(value)
!********************************************************************
!********************************************************************
  type(tree_type) ,intent(in) :: tree
  integer :: value
  value = tree%nv
  end function
!
!
  subroutine tree_checkinput( model ,process, nfinst )
!********************************************************************
!********************************************************************
  use avh_bint
  use avh_kaleu_model ,only : model_type,nparticle
  implicit none
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: process(-2:17),nfinst
  integer :: ii,jj,nn
  if (nfinst.lt.2) then
    write(*,*) 'ERROR in avh_kaleu_tree: nfinst =',nfinst ,&
               ', but should be at least 2'
    stop
  endif
  nn = nparticle( model )
  do ii=-2,nfinst
  if (ii.ne.0) then
    jj = process(ii)
    if (jj.lt.1) then
      write(*,*) 'ERROR in avh_kaleu_tree: process(',ii,') =' ,&
                 jj,', but should be larger than 0'
      stop
    endif
    if (jj.gt.nn) then
      write(*,*) 'ERROR in avh_kaleu_tree: process(',ii,') =' ,&
                 jj,', while nparticle =',nn
      stop
    endif
  endif
  enddo
  end subroutine
!
!
  subroutine tree_checkoutput( model ,tree ,process, nfinst ,cancel )
!********************************************************************
!********************************************************************
  use avh_bint
  use avh_kaleu_model ,only : model_type,sparticle
  implicit none
  type(model_type) ,intent(in) :: model
  type(tree_type)  ,intent(in) :: tree
  integer          ,intent(in) :: process(-2:17),nfinst
  logical          ,intent(out) :: cancel
  integer :: ii
  character(160) :: line
!
  cancel = .false.
!
  if (tree_check(tree).ne.0) then
    line = 'WARNING from avh_kaleu_tree: the process '    &
         //' '//trim(sparticle(model,process(-2)))        &
         //' '//trim(sparticle(model,process(-1)))//' ->'
    do ii=1,nfinst
      line = trim(line)//' '//trim(sparticle(model,process(ii)))
    enddo
    line = trim(line)//'  is not possible'
    write(*,*) trim(line)
    cancel = .true.
  else
    if (nunit.gt.0) then
      line = 'MESSAGE from avh_kaleu_tree: initializing for '  &
           //' '//trim(sparticle(model,process(-2)))           &
           //' '//trim(sparticle(model,process(-1)))//' ->'
      do ii=1,nfinst
        line = trim(line)//' '//trim(sparticle(model,process(ii)))
      enddo
      write(nunit,*) trim(line)
    endif
    call printtree( tree ,model ) !DEBUG
  endif
  end subroutine
!
  function tree_check( tree ) result(value)
!********************************************************************
!********************************************************************
  use avh_bint
  implicit none
  type(tree_type) ,intent(in) :: tree
  integer :: value
  value = 0
  if (tree%nv.eq.0) then
    value = 1
  elseif (tree%po(tree%ov(0,tree%nv)).ne.avh_bint_b(tree%n1)-1) then
    value = 2
  endif
  end function
!
!
  subroutine tree_copy( copy ,tree )
!********************************************************************
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: copy
  type(tree_type) ,intent(in)    :: tree
  integer :: nn
  nn = maxval( tree%vo(:,0) )
  copy%n1 = tree%n1
  copy%nv = tree%nv
  copy%no = tree%no
  copy%ov(0:3,1:copy%nv)  = tree%ov(0:3,1:copy%nv)
  copy%po(0:copy%no)      = tree%po(0:copy%no)
  copy%fo(0:copy%no)      = tree%fo(0:copy%no)
  copy%vo(0:copy%no,0:nn) = tree%vo(0:copy%no,0:nn)
  copy%typ(1:copy%nv)     = tree%typ(1:copy%nv)
  copy%zax(1:copy%nv)     = tree%zax(1:copy%nv)
  copy%frst1(1:copy%nv)   = tree%frst1(1:copy%nv)
  end subroutine
!
!
  subroutine tree_dress( model,vertx ,tree ,process,nfinst ,ioption )
!********************************************************************
!* Dress up the  nkinv  kinematical vertices  kinv
!* to obtain  tree%nv  dressed vertices  tree%ov ,
!* i.e. vertices refering to flavored oc's following the rules
!* coming from  fusion .
!********************************************************************
  use avh_bint
  use avh_kaleu_model ,only : model_type,vertx_type
  implicit none
  type(model_type) ,intent(in)    :: model
  type(vertx_type) ,intent(inout) :: vertx
  type(tree_type)  ,intent(out)   :: tree
  integer          ,intent(in)    :: process(-2:17),nfinst ,ioption
!
  call dress0( model,vertx ,tree ,process,nfinst )
!
  if (ioption.eq.0) then
! Remove vertices with leg having momentum 1, except for the last oc
    call dress_rm_1( tree )
    if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
  endif
! Re-find for each off-shell current the vertices it appears as 0-leg
  call dress_vo( tree )
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
! Put the type of vertices: no_generation=0, t-channel=1, s-channel=2
  call dress_inittype( tree )
  103 format(' MESSAGE from avh_kaleu_tree:  nvx =',i9,',  noc =',i9)
!
  end subroutine
!
!
  subroutine dress0( model,vertx ,tree ,process,nfinst )
!********************************************************************
!* Dress up the  nkinv  kinematical vertices  kinv
!* to obtain  tree%nv  dressed vertices  tree%ov ,
!* i.e. vertices refering to flavored oc's following the rules
!* coming from  fusion .
!********************************************************************
  use avh_bint
  use avh_kaleu_model ,only : model_type,vertx_type,nparticle,fillfusion
  implicit none
  type(model_type) ,intent(in)    :: model
  type(vertx_type) ,intent(inout) :: vertx
  type(tree_type)  ,intent(inout) :: tree
  integer          ,intent(in)    :: process(-2:17),nfinst
  integer ,allocatable :: fusion(:,:,:)
  integer ,allocatable :: o_p(:,:)
  integer :: nkinv,kinv(0:2,avh_bint_maxv)
  integer :: ii,p0,p1,p2,l1,l2,oc1,oc2,f1,f2,nfl,ifl,f0,jfl,oc0
  logical :: next
!
  tree%n1 = 0
  do p0=0,nfinst
    if (p0.eq.0) then
      ii = -1
    else
      ii = p0
    endif
    tree%n1 = tree%n1+1
    tree%po(tree%n1) = avh_bint_b(p0)
    tree%fo(tree%n1) = process(ii)
  enddo
!
  nfl = nparticle( model )
!
  allocate( fusion(1:nfl,1:nfl,0:nfl) )
  fusion(1:nfl,1:nfl,0:nfl) = 0
  call fillfusion( model,vertx ,fusion,nfl )
!
  allocate( o_p(1:2**tree%n1,0:nfl) )
  o_p(1:2**tree%n1,0:nfl) = 0
!
! Initialize tree%no
  tree%no = tree%n1
!
! Clean all oc's of level higher than 1
  tree%po(0) = 0
  tree%fo(0) = 0
  do ii=tree%no+1,size_noc
    tree%po(ii) = 0
    tree%fo(ii) = 0
  enddo
!
! Initialize list of vertices
  do ii=1,size_nvx
    tree%ov(0,ii) = 0
    tree%ov(1,ii) = 0
    tree%ov(2,ii) = 0
    tree%ov(3,ii) = 0
  enddo
! Initialize o_p
  do ii=1,tree%no
    o_p(tree%po(ii),0) = 1
    o_p(tree%po(ii),1) = ii
  enddo
!
  tree%nv = 0
!
  call avh_bint_kinv( kinv,nkinv ,tree%n1 )
!
  do ii=1,nkinv
    p0 = kinv(0,ii)
    p1 = kinv(1,ii)
    p2 = kinv(2,ii)
!
    do l1=1,o_p(p1,0)
    do l2=1,o_p(p2,0)
      oc1 = o_p(p1,l1)
      oc2 = o_p(p2,l2)
      f1 = tree%fo(oc1)
      f2 = tree%fo(oc2)
      nfl = fusion( f1 ,f2 ,0 )
      if (nfl.ne.0) then
        do ifl=1,nfl
          f0 = fusion( f1 ,f2 ,ifl )
          jfl = 0
          next = .true.
          do while (next)
            jfl = jfl+1
            oc0 = o_p(p0,jfl)
            if ( f0.eq.tree%fo(oc0) ) then
              tree%nv = tree%nv+1
              if (tree%nv.gt.size_nvx) then
                write(*,101) tree%nv
                stop
              endif
              tree%ov(0,tree%nv) = oc0
              tree%ov(1,tree%nv) = oc1
              tree%ov(2,tree%nv) = oc2
              next = .false.
            endif
            if (next.and.(jfl.ge.o_p(p0,0))) then
              tree%no = tree%no+1
              if (tree%no.gt.size_noc) then
                write(*,102) tree%no
                stop
              endif
              tree%po(tree%no) = p0
              tree%fo(tree%no) = f0
              o_p(p0,0) = o_p(p0,0) + 1
              o_p(p0,o_p(p0,0)) = tree%no
              tree%nv = tree%nv+1
              if (tree%nv.gt.size_nvx) then
                write(*,101) tree%nv
                stop
              endif
              tree%ov(0,tree%nv) = tree%no
              tree%ov(1,tree%nv) = oc1
              tree%ov(2,tree%nv) = oc2
              next = .false.
            endif
          enddo
        enddo
      endif
    enddo
    enddo
!
  enddo
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
!
  101 format(' ERROR in avh_kaleu_tree: increase the parameter' ,&
             ' "size_nvx" to at least',i9)
  102 format(' ERROR in avh_kaleu_tree: increase the parameter' ,&
             ' "size_noc" to at least',i9)
  103 format(' MESSAGE from avh_kaleu_tree:  nvx =',i9,',  noc =',i9)
!
  deallocate( fusion )
  deallocate( o_p )
!
! Remove vertices with final momentum having the wrong flavor
  call dress_rm_last( tree ,process(-2) )
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
!
! Remove dead ends
  call dress_rm_dead( tree )
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
!
! Until here we have a list of vertices which would be sufficient for
! amplitude calculation. Now adjust for phase space generation
!
! Find for each off-shell current the vertices it appears as 0-leg
  call dress_vo( tree )
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
!
! Make sure the 1-leg of each vertex has the odd momentum
  call dress_odd( tree )
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
!
! Identify T-channel vertices, and increase the list of vertices
! corresponding to the different s-invariants
  call dress_tchan( tree )
  if (sizeunit.gt.0) write(sizeunit,103) tree%nv,tree%no
!
  end subroutine
!
  subroutine dress_shup( tree ,ii )
!********************************************************************
!* Remove vertex ii and shift the ones below one up
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer ,intent(in) :: ii
  integer :: jj,kk
  tree%nv = tree%nv-1
  do jj=ii,tree%nv
    kk = jj+1
    tree%ov(0,jj) = tree%ov(0,kk)
    tree%ov(1,jj) = tree%ov(1,kk)
    tree%ov(2,jj) = tree%ov(2,kk)
    tree%ov(3,jj) = tree%ov(3,kk)
  enddo
  end subroutine
!
  subroutine dress_shdown( tree ,ii )
!********************************************************************
!* Shift all vertices below ii one down, and copy vertex ii one down
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer ,intent(in) :: ii
  integer :: jj,kk
  tree%nv = tree%nv+1
  if (tree%nv.gt.size_nvx) then
    write(*,*) &
      'ERROR in avh_kaleu_tree: increase the parameter' ,&
      ' "size_nvx" to at least',tree%nv
    stop
  endif
  do jj=tree%nv,ii+1,-1
    kk = jj-1
    tree%ov(0,jj) = tree%ov(0,kk)
    tree%ov(1,jj) = tree%ov(1,kk)
    tree%ov(2,jj) = tree%ov(2,kk)
    tree%ov(3,jj) = tree%ov(3,kk)
  enddo
  end subroutine
!
  subroutine dress_rm_last( tree ,flast )
!********************************************************************
!* Remove vertices with final momentum having the wrong flavor
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer ,intent(in) :: flast
  integer :: p0,ii,jj
  p0 = tree%po(tree%ov(0,tree%nv))
  ii = tree%nv
  do while ( tree%po(tree%ov(0,ii)).eq.p0 )
    if ( tree%fo(tree%ov(0,ii)).ne.flast ) call dress_shup( tree ,ii )
    ii = ii-1
  enddo
  end subroutine
!
  subroutine dress_rm_dead( tree )
!********************************************************************
!* Remove dead ends
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  logical :: isonlhs(0:tree%no),isonrhs(0:tree%no)
  integer :: ii,jj,oc0,oc1,oc2,oc3
  jj = 0
  do while (jj.ne.tree%nv)
    isonlhs(0:tree%no) = .false.
    isonrhs(0:tree%no) = .false.
    isonlhs(0:tree%n1) = .true.
    isonrhs(0:tree%n1) = .true.
    isonlhs(tree%ov(0,tree%nv)) = .true.
    isonrhs(tree%ov(0,tree%nv)) = .true.
    do ii=1,tree%nv
      isonlhs(tree%ov(0,ii)) = .true.
      isonrhs(tree%ov(1,ii)) = .true.
      isonrhs(tree%ov(2,ii)) = .true.
      isonrhs(tree%ov(3,ii)) = .true.
    enddo
    jj = tree%nv
    tree%nv = 0
    do ii=1,jj
      oc0 = tree%ov(0,ii)
      oc1 = tree%ov(1,ii)
      oc2 = tree%ov(2,ii)
      oc3 = tree%ov(3,ii)
      if (oc3.ne.0) oc1 = oc3
      if (     isonlhs(oc0).and.isonrhs(oc0) &
          .and.isonlhs(oc1).and.isonrhs(oc1) &
          .and.isonlhs(oc2).and.isonrhs(oc2) ) then
        tree%nv = tree%nv+1
        tree%ov(0:3,tree%nv) = tree%ov(0:3,ii)
        tree%typ(tree%nv)    = tree%typ(ii)
        tree%zax(tree%nv)    = tree%zax(ii)
        tree%frst1(tree%nv)  = tree%frst1(ii)
      endif
    enddo
  enddo
!  p0 = tree%po(tree%ov(0,tree%nv))
!  ii = tree%nv
!  do while (tree%po(tree%ov(0,ii)).eq.p0) 
!    ii = ii-1
!  enddo
!  do while( ii.ge.1 )
!    oc0 = tree%ov(0,ii)
!    jj = ii+1
!    do while ( jj.le.tree%nv &
!              .and.tree%ov(1,jj).ne.oc0 &
!              .and.tree%ov(2,jj).ne.oc0 &
!              .and.tree%ov(3,jj).ne.oc0 )
!      jj = jj+1
!    enddo
!    if (jj.gt.tree%nv) call dress_shup( tree ,ii )
!    ii = ii-1
!  enddo
  end subroutine
!
  subroutine dress_vo( tree )
!********************************************************************
!* Find for each off-shell current the vertices it appears as 0-leg
!********************************************************************
  use avh_bint
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer :: ii,oc0
  tree%vo(1:tree%no,0) = 0
  do ii=1,tree%nv
    oc0 = tree%ov(0,ii)
    tree%vo(oc0,0) = tree%vo(oc0,0)+1 ! the number of vertices
    if (tree%vo(oc0,0).gt.size_nvo) then
      write(*,*) &
        'ERROR in avh_kaleu_tree: increase the parameter' ,&
        ' "size_nvo" to at least',tree%vo(oc0,0)
      stop
    endif
    tree%vo(oc0,tree%vo(oc0,0)) = ii ! the vertex
  enddo
!
! Only external particles may occure in exactly 0 vertices 
  tree%vo(0,0) = -1
  do ii=1,tree%no
    if ( avh_bint_l(tree%po(ii)).gt.1 &
             .and. tree%vo(ii,0).eq.0 ) tree%vo(ii,0) = -1
  enddo
  end subroutine
!
  subroutine dress_tchan( tree )
!********************************************************************
!* Identify T-channel vertices, and increase the list of vertices
!* corresponding to the different s-invariants
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer :: ii,oc0,p0,oc1,nn,l0,l1,p1,jj
  ii=tree%nv+1
!  do while (ii.gt.0)
! Corrected by AK 25.08.2014
  do while (ii.gt.1)
    ii = ii-1
    oc0 = tree%ov(0,ii)
    p0 = tree%po(oc0)
    if (p0.ne.p0/2*2) then
      oc1 = tree%ov(1,ii)
      nn = tree%vo(oc1,0)
      l0 = 0
      l1 = ii-1
      do jj=1,nn
        p1 = tree%po(tree%ov(1,tree%vo(oc1,jj)))
        if (p1.eq.1.or.l0.eq.0) then
          l1 = l1+1
          if (l1.gt.ii) call dress_shdown( tree ,ii )
          if (p1.ne.1) then
            l0 = 1
            tree%ov(3,ii) = 0
          else
            tree%ov(3,ii) = tree%ov(2,tree%vo(oc1,jj))
          endif
        endif
      enddo
    endif
  enddo
  end subroutine
!
  subroutine dress_odd( tree )
!********************************************************************
!* Make sure the 1-leg of each vertex has the odd momentum,
!* or else the smallest momentum
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer :: ii,p0,oc1,p1,oc2,p2
  do ii=1,tree%nv
    p0 = tree%po(tree%ov(0,ii))
    oc1 = tree%ov(1,ii)
    p1 = tree%po(oc1)
    if (p0.ne.p0/2*2) then
      if (p1.eq.p1/2*2) then
        oc1 = tree%ov(2,ii)
        tree%ov(2,ii) = tree%ov(1,ii)
        tree%ov(1,ii) = oc1
      endif
    else
      oc2 = tree%ov(2,ii)
      p2 = tree%po(oc2)
      if (p1.gt.p2) then
        oc1 = tree%ov(2,ii)
        tree%ov(2,ii) = tree%ov(1,ii)
        tree%ov(1,ii) = oc1
      endif
    endif
  enddo
  end subroutine
!
  subroutine dress_rm_1( tree )
!********************************************************************
!* Remove vertices with leg having momentum 1,
!* except for the last oc
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer :: ii,jj,p0
  p0 = tree%po(tree%ov(0,tree%nv))
  jj = 0
  do ii=1,tree%nv
    if (ii.eq.jj) cycle
    if (tree%po(tree%ov(0,ii)).eq.p0.or.tree%po(tree%ov(1,ii)).ne.1) then
      jj = jj+1
      tree%ov(0:3,jj) = tree%ov(0:3,ii)
      tree%typ(jj)    = tree%typ(ii)  
      tree%zax(jj)    = tree%zax(ii)  
      tree%frst1(jj)  = tree%frst1(ii)
    endif
  enddo
  tree%nv = jj
  end subroutine
!  
  subroutine dress_inittype ( tree )
!********************************************************************
!* put the type of vertex, zax, and frst1
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  integer :: ii,p0,p1,p2
  do ii=1,tree%nv
    p0 = tree%po(tree%ov(0,ii))
    p1 = tree%po(tree%ov(1,ii))
    p2 = tree%po(tree%ov(2,ii))
    tree%frst1(ii) = .true.
    if     (p1.eq.1.or.p2.eq.1) then
      tree%typ(ii) = 0
      tree%zax(ii) = 0
    elseif (p0.ne.p0/2*2) then
      tree%typ(ii) = 1
      tree%zax(ii) = 1
    else
      tree%typ(ii) = 2
      tree%zax(ii) = 0
    endif
  enddo
  end subroutine
!  
!  
  subroutine printtree( tree ,model )
!********************************************************************
!********************************************************************
  use avh_kaleu_model ,only : model_type,sparticle
  implicit none
  type(tree_type)  ,intent(in) :: tree
  type(model_type) ,intent(in) :: model
  integer :: ii,o0,o1,o2,o3
  character(12) :: f0,f1,f2,f3
  if (treeunit.gt.0) then
  do ii=1,tree%nv
    if (tree%typ(ii).eq.1.or.tree%typ(ii).eq.3.or.tree%typ(ii).eq.5) then
      o0 = tree%ov(0,ii)
      o1 = tree%ov(1,ii)
      o2 = tree%ov(2,ii)
      o3 = tree%ov(3,ii)
      f0 = sparticle(model,tree%fo(o0))
      f1 = sparticle(model,tree%fo(o1))
      f2 = sparticle(model,tree%fo(o2))
      f3 = sparticle(model,tree%fo(o3))
      write(treeunit,101) &
        'vrtx(',ii,'):  ','[',o0,tree%po(o0),f0,']' ,&
                '  <--  ','[',o1,tree%po(o1),f1,']' ,&
                     '  ','[',o2,tree%po(o2),f2,']' ,&
                     '  ','(',o3,tree%po(o3),f3,')' ,&
                     '  ','(',tree%zax(ii),')',tree%typ(ii),tree%vo(o0,0)
  101   format(a5,i5,a4 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     a7 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     a2 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     a2 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     a2 ,a1,i3,1x,a1,i2,i3)
    else
      o0 = tree%ov(0,ii)
      o1 = tree%ov(1,ii)
      o2 = tree%ov(2,ii)
      f0 = sparticle(model,tree%fo(o0))
      f1 = sparticle(model,tree%fo(o1))
      f2 = sparticle(model,tree%fo(o2))
      write(treeunit,102) &
        'vrtx(',ii,'):  ','[',o0,tree%po(o0),f0,']' ,&
                '  <--  ','[',o1,tree%po(o1),f1,']' ,&
                     '  ','[',o2,tree%po(o2),f2,']' ,&
                     '  ','(',tree%zax(ii),')',tree%typ(ii),tree%vo(o0,0)
  102   format(a5,i5,a4 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     a7 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     a2 ,a1,i3,i4,1x,a2,1x,a1 ,&
                     15x                      ,&
                     a2 ,a1,i3,1x,a1,i2,i3)
    endif
  enddo
  endif
  end subroutine
!  
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
!
  subroutine tree_order( ordt ,tree ,order,nparton )
!***********************************************************************
!***********************************************************************
  use avh_bint
  implicit none
  type(tree_type) ,intent(out) :: ordt
  type(tree_type) ,intent(in)  :: tree
  integer         ,intent(in)  :: nparton
  integer         ,intent(in)  :: order(nparton)
  logical :: isparton(0:2**tree%n1),keep,remove(1:tree%nv)
  logical :: isonlhs(0:tree%no),isonrhs(0:tree%no)
  integer :: frstparton(1:2**tree%n1),lastparton(1:2**tree%n1)
  integer :: ii,jj,p0,p1,p2,p3,oc0,order0(0:nparton)
!
  isparton( 0:avh_bint_b(tree%n1) ) = .false.
  do ii=1,nparton-1
    jj = order(ii)
    if (jj.eq.tree%n1) then
      write(*,*) 'ERROR: initial-state parton -2 should be the last ',&
                 'one in any ordering'
      stop
    endif
    order0(ii) = avh_bint_b(jj)
    isparton( order0(ii) ) = .true.
  enddo
  jj = order(nparton)
  if (jj.eq.tree%n1) then
    order0(nparton) = avh_bint_b(jj)-1
    order0(0) = order0(nparton)
    isparton( order0(nparton) ) = .true.
    call ordered( frstparton,lastparton ,order,nparton-1 ,tree%n1-1 )
  else
    if (isparton(1)) then
      write(*,*) 'ERROR: if there is only 1 parton in the initial ',&
                 'state, it should be the -2'
      stop
    endif
    order0(nparton) = avh_bint_b(jj)
    order0(0) = 0
    isparton( order0(nparton) ) = .true.
    call ordered( frstparton,lastparton ,order,nparton   ,tree%n1-1 )
  endif
!
! Remove vertices violating the ordering
  ordt%nv = 0
  do ii=1,tree%nv
    p0 = tree%po(tree%ov(0,ii))
    p1 = tree%po(tree%ov(1,ii))
    p2 = tree%po(tree%ov(2,ii))
    p3 = tree%po(tree%ov(3,ii))
!
    keep = .false.
    if (frstparton(p1).gt.-1.and.frstparton(p2).gt.-1) then
      if     (frstparton(p1).eq.0.and.frstparton(p2).eq.0) then
        keep = .true.
      elseif (frstparton(p1).eq.0) then
        keep = .true.
      elseif (frstparton(p2).eq.0) then
        if ( p3.eq.0.or.lastparton(p3)+1.eq.frstparton( 1) &
                    .or.lastparton( 1)+1.eq.frstparton(p3) ) then
          keep = .true.
        endif
      elseif (lastparton(p1)+1.eq.frstparton(p2)) then
        if     (p3.eq.0) then ! includes even p0
          if (avh_bint_l(p1).eq.1.or.avh_bint_l(p2).gt.1) then
            keep = .true.
          endif
        elseif (    lastparton(p3)+1.eq.frstparton( 1) &
                .or.lastparton( 1)+1.eq.frstparton(p3) ) then
          if (avh_bint_l(p3).eq.1.or.avh_bint_l(p2).gt.1) then
            keep = .true.
          endif
        endif
      elseif (lastparton(p2)+1.eq.frstparton(p1)) then
        if     (p3.eq.0) then !includes even p0
          if (avh_bint_l(p2).eq.1.or.avh_bint_l(p1).gt.1) then
            keep = .true.
          endif
        elseif (    lastparton(p3)+1.eq.frstparton( 1) &
                .or.lastparton( 1)+1.eq.frstparton(p3) ) then
          if (avh_bint_l(p2).eq.1.or.avh_bint_l(p3).gt.1) then
            keep = .true.
          endif
        endif
      endif
    endif
!        
!    if     (frstparton(p1).eq.-1.or.frstparton(p2).eq.-1) then
!      keep = .false.
!    elseif (frstparton(p1).eq.0.and.frstparton(p2).eq.0) then
!      keep = .true.
!    elseif (    (frstparton(p1).eq.0) &
!            .or.(frstparton(p2).eq.0) &
!            .or.( (lastparton(p1)+1.eq.frstparton(p2)) &
!                 .and.(avh_bint_l(p1).eq.1.or.avh_bint_l(p2).gt.1) ) &
!            .or.( (lastparton(p2)+1.eq.frstparton(p1)) &
!                 .and.(avh_bint_l(p2).eq.1.or.avh_bint_l(p1).gt.2.or.(avh_bint_l(p1).eq.2.and.p1.eq.p1/2*2)) ) &
!           ) then
!      if     (p0.eq.p0/2*2) then
!        keep = .true.
!      elseif (p3.eq.0) then
!        keep = .true.
!      elseif (frstparton(1).eq.0) then
!        keep = .true.
!      elseif (    lastparton(p3)+1.eq.frstparton( 1) &
!              .or.lastparton( 1)+1.eq.frstparton(p3) ) then
!        keep = .true.
!      else
!        keep = .false.
!      endif
!    else
!      keep = .false.
!    endif
!
    if (keep) then
      ordt%nv = ordt%nv+1
      ordt%ov(0:3,ordt%nv) = tree%ov(0:3,ii)
      ordt%typ(ordt%nv) = tree%typ(ii)
      ordt%zax(ordt%nv) = tree%zax(ii)
      ordt%frst1(ordt%nv) = (frstparton(p1).le.frstparton(p2)) 
! Start special cases involving on-shell partons
      if (p0.eq.p0/2*2) then ! s-channel
        if (isparton(p1).and.isparton(p2)) then
          if (frstparton(p1).lt.frstparton(p2)) then
            jj = frstparton(p1)-1
            if (jj.ne.0.or.order0(jj).ne.0) then
!            if (order0(jj).eq.1.or.order0(jj).eq.avh_bint_b(tree%n1)-1) then !DEBUG
              ordt%typ(ordt%nv) = 4
              ordt%zax(ordt%nv) = order0(jj)
!            endif !DEBUG
            endif
          else
            jj = frstparton(p2)-1
            if (jj.ne.0.or.order0(jj).ne.0) then
!            if (order0(jj).eq.1.or.order0(jj).eq.avh_bint_b(tree%n1)-1) then !DEBUG
              ordt%typ(ordt%nv) = 4
              ordt%zax(ordt%nv) = order0(jj)
              ordt%ov(1,ordt%nv) = tree%ov(2,ii)
              ordt%ov(2,ordt%nv) = tree%ov(1,ii)
!            endif !DEBUG
            endif
          endif
        elseif (isparton(p1)) then
          jj = frstparton(p1)-1
          if (jj.ne.0.or.order0(jj).ne.0) then
!          if (order0(jj).eq.1.or.order0(jj).eq.avh_bint_b(tree%n1)-1) then !DEBUG
            ordt%typ(ordt%nv) = 4
            ordt%zax(ordt%nv) = order0(jj)
!          endif !DEBUG
          endif
        elseif (isparton(p2)) then
          jj = frstparton(p2)-1
          if (jj.ne.0.or.order0(jj).ne.0) then
!          if (order0(jj).eq.1.or.order0(jj).eq.avh_bint_b(tree%n1)-1) then !DEBUG
            ordt%typ(ordt%nv) = 4
            ordt%zax(ordt%nv) = order0(jj)
            ordt%ov(1,ordt%nv) = tree%ov(2,ii)
            ordt%ov(2,ordt%nv) = tree%ov(1,ii)
!          endif !DEBUG
          endif
        endif
      else
        jj = frstparton(p3)-1
        if (isparton(p3).and.isparton(p2)) then
          if (frstparton(p3).lt.frstparton(p2)) then
            jj = frstparton(p3)-1
            if (jj.ne.0.or.order0(jj).ne.0) then
!            if (order0(jj).ne.1.and.order0(jj).ne.avh_bint_b(tree%n1)-1) then !DEBUG
              ordt%typ(ordt%nv) = 3
              ordt%zax(ordt%nv) = order0(jj)
!            endif !DEBUG
            endif
          else
            jj = frstparton(p2)-1
            if (jj.ne.0.or.order0(jj).ne.0) then
!            if (order0(jj).ne.1.and.order0(jj).ne.avh_bint_b(tree%n1)-1) then !DEBUG
              ordt%typ(ordt%nv) = 5
              ordt%zax(ordt%nv) = order0(jj)
!            endif !DEBUG
            endif
          endif
        elseif (isparton(p3)) then
          jj = frstparton(p3)-1
          if (jj.ne.0.or.order0(jj).ne.0) then
!          if (order0(jj).ne.1.and.order0(jj).ne.avh_bint_b(tree%n1)-1) then !DEBUG
            ordt%typ(ordt%nv) = 3
            ordt%zax(ordt%nv) = order0(jj)
!          endif !DEBUG
          endif
        elseif (isparton(p2)) then
          jj = frstparton(p2)-1
          if (jj.ne.0.or.order0(jj).ne.0) then
!          if (order0(jj).ne.1.and.order0(jj).ne.avh_bint_b(tree%n1)-1) then !DEBUG
            ordt%typ(ordt%nv) = 5
            ordt%zax(ordt%nv) = order0(jj)
!          endif !DEBUG
          endif
        endif
      endif
! End special cases involving on-shell partons
    endif
  enddo
!
  ordt%n1 = tree%n1
  ordt%no = tree%no
  ordt%po(0:ordt%no) = tree%po(0:tree%no)
  ordt%fo(0:ordt%no) = tree%fo(0:tree%no)
!
  call dress_rm_1( ordt )
  call dress_rm_dead( ordt )
  call dress_vo( ordt )
!
  end subroutine
!
  subroutine ordered( frstparton,lastparton ,order,nparton ,nfinst )
!***********************************************************************
!* Find all numbers 1,2,..,2^(nfinst+1)-1 compatible with order, eg
!*                                   2^(anything_not_in_order) --> 0,0
!*  2^order(1)+2^order(2)           +2^(anything_not_in_order) --> 1,2
!*             2^order(2)+2^order(3)+2^(anything_not_in_order) --> 2,3
!*  2^order(1)           +2^order(3)+2^(anything_not_in_order) --> -1,-1
!*  2^order(1)+2^order(2)+2^order(3)+2^(anything_not_in_order) --> 1,3
!***********************************************************************
  use avh_bint
  implicit none
  integer ,intent(in)  :: nparton,nfinst
  integer ,intent(in)  :: order(nparton)
  integer ,intent(out) :: frstparton(1:2**(nfinst+1)),&
                          lastparton(1:2**(nfinst+1))
  integer :: ifrst,ilast,nn,kk,ii,jj,list(nfinst+2),aa(nfinst+2)
  integer :: nb,bb(0:2**(nfinst+1))
  logical :: next,isparton(0:nfinst)
!
  frstparton(1:avh_bint_b(nfinst+1)) = -1
  lastparton(1:avh_bint_b(nfinst+1)) = -1
!
! find all numbers from 0...nfinst
! that are not equal to any number in order(1:nparton)
  isparton(0:nfinst) = .false.
  do ii=1,nparton
    isparton(order(ii)) = .true.
  enddo
  nn = 0
  do ii=0,nfinst
  if (.not.isparton(ii)) then
    nn = nn+1
    list(nn) = ii
  endif
  enddo
!
! construct all possible numbers from 2^list(1), 2^list(2),...
  nb = 0
  bb(nb) = 0
  do ii=1,nn
  call avh_bint_antl_i( next,aa ,ii,nn )
  do while (next)
    jj = 0
    do kk=1,ii
      jj = jj + avh_bint_b( list(aa(kk)) )
    enddo
    frstparton(jj) = 0
    lastparton(jj) = 0
    nb = nb+1
    bb(nb) = jj
  call avh_bint_antl_n( next,aa ,ii,nn )
  enddo!nexta
  enddo!ii=1,nn 
!
  do ii=0,nparton-1
  do ifrst=1,nparton-ii
    ilast = ifrst+ii
! calculate jj = 2^order(ifrst)+2^order(ifrst+1)+...+2^order(ilast)
    jj = 0
    do kk=ifrst,ilast
      jj = jj + avh_bint_b(order(kk))
    enddo
! all numbers that can be written as jj+other_powers_of_2, where
! other_powers_of_2 does not contain any of the 2^order(1:nparton)
    do kk=0,nb
      frstparton(jj+bb(kk)) = ifrst
      lastparton(jj+bb(kk)) = ilast
    enddo
  enddo!ifrst=1,nparton-ii
  enddo!ii=0,nparton-1
!
  end subroutine
!

end module
