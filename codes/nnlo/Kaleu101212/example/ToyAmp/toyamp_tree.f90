module avh_toyamp_tree
  private ! Everything is private except the following list:
  public :: tree_type,tree_checkinput,tree_dress,tree_checkoutput ,&
            tree_nv,tree_copy,tree_ij
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
  use avh_toyamp_model ,only : model_type,nparticle
  implicit none
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: process(-2:17),nfinst
  integer :: ii,jj,nn
  if (nfinst.lt.2) then
    write(*,*) 'ERROR in avh_toyamp_tree: nfinst =',nfinst ,&
               ', but should be at least 2'
    stop
  endif
  nn = nparticle( model )
  do ii=-2,nfinst
  if (ii.ne.0) then
    jj = process(ii)
    if (jj.lt.1) then
      write(*,*) 'ERROR in avh_toyamp_tree: process(',ii,') =' ,&
                 jj,', but should be larger than 0'
      stop
    endif
    if (jj.gt.nn) then
      write(*,*) 'ERROR in avh_toyamp_tree: process(',ii,') =' ,&
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
  use avh_toyamp_model ,only : model_type,sparticle
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
    line = 'WARNING from avh_toyamp_tree: the process '   &
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
      line = 'MESSAGE from avh_toyamp_tree: initializing for ' &
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
  end subroutine
!
!
  subroutine tree_dress( model,vertx ,tree ,process,nfinst )
!********************************************************************
!* Dress up the  nkinv  kinematical vertices  kinv
!* to obtain  tree%nv  dressed vertices  tree%ov ,
!* i.e. vertices refering to flavored oc's following the rules
!* coming from  fusion .
!********************************************************************
  use avh_bint
  use avh_toyamp_model ,only : model_type,vertx_type,nparticle,fillfusion
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
  101 format(' ERROR in avh_toyamp_tree: increase the parameter' ,&
             ' "size_nvx" to at least',i9)
  102 format(' ERROR in avh_toyamp_tree: increase the parameter' ,&
             ' "size_noc" to at least',i9)
  103 format(' MESSAGE from avh_toyamp_tree:  nvx =',i9,',  noc =',i9)
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
! Find for each off-shell current the vertices it appears as 0-leg
  call dress_vo( tree )
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
      'ERROR in avh_toyamp_tree: increase the parameter' ,&
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
  tree%no = tree%ov(0,tree%nv)
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
      endif
    enddo
  enddo
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
        'ERROR in avh_toyamp_tree: increase the parameter' ,&
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
  subroutine printtree( tree ,model )
!********************************************************************
!********************************************************************
  use avh_toyamp_model ,only : model_type,sparticle
  implicit none
  type(tree_type)  ,intent(in) :: tree
  type(model_type) ,intent(in) :: model
  integer :: ii,o0,o1,o2,o3
  character(2) :: f0,f1,f2,f3
  if (treeunit.gt.0) then
  do ii=1,tree%nv
    o0 = tree%ov(0,ii)
    o1 = tree%ov(1,ii)
    o2 = tree%ov(2,ii)
    f0 = sparticle(model,tree%fo(o0))
    f1 = sparticle(model,tree%fo(o1))
    f2 = sparticle(model,tree%fo(o2))
    write(treeunit,102) &
      'vrtx(',ii,'):  ','[',o0,tree%po(o0),f0,']' ,&
              '  <--  ','[',o1,tree%po(o1),f1,']' ,&
                   '  ','[',o2,tree%po(o2),f2,']'
  102 format(a5,i5,a4 ,a1,i3,i4,1x,a2,1x,a1 ,&
                   a7 ,a1,i3,i4,1x,a2,1x,a1 ,&
                   a2 ,a1,i3,i4,1x,a2,1x,a1 ) 
  enddo
  endif
  end subroutine
!  
!
  subroutine tree_ij( tree ,cancel ,ii,jj )
!********************************************************************
!* Reduce tree such that external particles ii and jj are only
!* connected to the same vertex.
!********************************************************************
  implicit none
  type(tree_type) ,intent(inout) :: tree
  logical         ,intent(out)   :: cancel
  integer         ,intent(in)    :: ii,jj
  integer :: iv,i1,j1
  if (ii.eq.jj) then; write(*,*) 'ERROR in tree_ij: i=j'; stop; endif
  if (ii.lt.-2) then; write(*,*) 'ERROR in tree_ij: i<-2'; stop; endif
  if (jj.lt.-2) then; write(*,*) 'ERROR in tree_ij: j<-2'; stop; endif
  if (ii.eq. 0) then; write(*,*) 'ERROR in tree_ij: i=0'; stop; endif
  if (jj.eq. 0) then; write(*,*) 'ERROR in tree_ij: j=0'; stop; endif
  if (ii.ge.tree%n1) then; write(*,*) 'ERROR in tree_ij: i>n'; stop; endif
  if (jj.ge.tree%n1) then; write(*,*) 'ERROR in tree_ij: j>n'; stop; endif
!
  if (ii.ne.-2.and.jj.ne.-2) then
    if (ii.eq.-1) then; i1=1; else; i1=ii+1; endif
    if (jj.eq.-1) then; j1=1; else; j1=jj+1; endif
    iv = tree%nv
    do while (iv.gt.0)
      if ((tree%ov(1,iv).eq.i1.and.tree%ov(2,iv).ne.j1).or. &
          (tree%ov(2,iv).eq.i1.and.tree%ov(1,iv).ne.j1) ) then
        tree%nv = tree%nv-1
        tree%ov(0:3,iv:tree%nv) = tree%ov(0:3,iv+1:tree%nv+1)
      endif
      iv = iv-1
    enddo
  else
    if   (ii.eq.-2) then; if (jj.eq.-1) then; i1=1; else; i1=jj+1; endif
    else;                 if (ii.eq.-1) then; i1=1; else; i1=ii+1; endif
    endif
    j1 = tree%ov(0,tree%nv) ! last off-shell current
    iv = tree%nv
    do while (iv.gt.0)
      if ((tree%ov(1,iv).eq.i1.and.tree%ov(0,iv).ne.j1).or. &
          (tree%ov(2,iv).eq.i1.and.tree%ov(0,iv).ne.j1) ) then
        tree%nv = tree%nv-1
        tree%ov(0:3,iv:tree%nv) = tree%ov(0:3,iv+1:tree%nv+1)
      endif
      iv = iv-1
    enddo
  endif
! Remove dead ends
  call dress_rm_dead( tree )
  cancel = (tree_check(tree).ne.0)
! Find for each off-shell current the vertices it appears as 0-leg
  if (.not.cancel) call dress_vo( tree )
  end subroutine
!
!
end module
