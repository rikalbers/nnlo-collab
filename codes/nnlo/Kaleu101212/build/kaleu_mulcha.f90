module avh_kaleu_mulcha
  use avh_kaleu_tree ,only : tree_type
!
! Unit to which messages are send
  integer ,parameter ,private :: nunit = 6
!
  type :: mulcha_type 
    logical :: yes = .false.
    integer ,allocatable :: ov(:,:),vo(:,:)
    logical ,allocatable :: frst1(:)
    integer              :: n1,nv,no
    real(kind(1d0)) ,allocatable :: wtv(:),ave(:),dns(:)
    integer                      :: idat,ndat,ntot,istp,nstp
    real(kind(1d0)) :: small
  end type
!
!
contains
!
!
  subroutine mulcha_init( tree ,obj )
!********************************************************************
!********************************************************************
  implicit none
  type(tree_type)   ,intent(in)  :: tree
  type(mulcha_type) ,intent(out) :: obj
  integer nn
!
  obj%n1 = tree%n1
  obj%nv = tree%nv
  obj%no = tree%no
!
  allocate( obj%ov(0:3,1:obj%nv) )
  nn = maxval( tree%vo(:,0) )
  allocate( obj%vo(0:obj%no,0:nn) )
  allocate( obj%frst1(1:obj%nv) )
  obj%ov(0:3,1:obj%nv)  = tree%ov(0:3,1:obj%nv)
  obj%vo(0:obj%no,0:nn) = tree%vo(0:obj%no,0:nn)
  obj%frst1(1:obj%nv)   = tree%frst1(1:obj%nv)
  end subroutine
!
!
  subroutine mulcha_init_adapt( obj ,ndat,nstp,small )
!********************************************************************
!********************************************************************
  use avh_print
  implicit none
  type(mulcha_type) ,intent(inout) :: obj
  integer           ,intent(in)    :: ndat,nstp
  real(kind(1d0))   ,intent(in)    :: small
  integer :: ii,nn,iv,jj
  real(kind(1d0)) :: wtv0
  character(10) :: chndat,chnstp,chsmall
  character(160) :: line
!
  obj%yes = (nstp.gt.0)
  if (.not.obj%yes) return
!
  allocate( obj%wtv(1:obj%nv) )
  allocate( obj%ave(1:obj%nv) )
  allocate( obj%dns(1:obj%nv) )
!
  obj%idat = 0
  obj%ndat = ndat
  obj%ntot = 0
  obj%istp = 0
  obj%nstp = nstp
  obj%small = small
!  if (nunit.gt.0) then
!    write(chndat,'(i10)') ndat
!    write(chnstp,'(i10)') nstp
!    write(chsmall,'(es8.2)') small
!    line = 'MESSAGE from avh_kaleu_mulcha:' &
!         //' nbatch='//trim(adjustl(chndat)) &
!         //' nstep='//trim(adjustl(chnstp)) &
!         //' thrs='//trim(adjustl(chsmall))
!    write(nunit,*) trim(line)
!  endif
  if (nunit.gt.0) write(nunit,*) 'MESSAGE from avh_kaleu_mulcha:' &
                 ,' nbatch=',trim(printint(ndat)) &
                 ,' nstep=',trim(printint(nstp)) &
                 ,' thrs=',trim(printdbl(1,2,small))
  do ii=1,obj%no
    nn = obj%vo(ii,0)
    if (nn.ge.1) then
      wtv0 = 1d0/dble(nn)
      do jj=1,nn
        iv = obj%vo(ii,jj)
        obj%ave(iv) = 0d0
        obj%wtv(iv) = wtv0
        obj%dns(iv) = 0d0
      enddo
    endif
  enddo
  end subroutine
!
!
  subroutine mulcha_close( obj )
!********************************************************************
!********************************************************************
  implicit none
  type(mulcha_type) ,intent(inout) :: obj
  obj%yes = .false.
  obj%nv = 0
  obj%no = 0 
  obj%idat = 0
  obj%ndat = 0 
  obj%ntot = 0
  obj%istp = 0
  obj%nstp = 0
  obj%small = 0d0
  if (allocated(obj%ov   )) deallocate( obj%ov    )
  if (allocated(obj%vo   )) deallocate( obj%vo    )
  if (allocated(obj%frst1)) deallocate( obj%frst1 )
  if (allocated(obj%wtv)) deallocate( obj%wtv )
  if (allocated(obj%ave)) deallocate( obj%ave )
  if (allocated(obj%dns)) deallocate( obj%dns )
  end subroutine
!
!
  subroutine mulcha_collect( obj ,wght ,jremtot,nremtot )
!********************************************************************
!********************************************************************
  implicit none
  type(mulcha_type) ,intent(inout) :: obj
  real(kind(1d0))   ,intent(in)    :: wght
  integer           ,intent(out)   :: nremtot,jremtot(obj%nv)
  integer :: nn,kk,jj,ii,nrem
  integer :: jrem(obj%nv)
  real(kind(1d0)) :: factor,hsum,thrs
!
  nremtot = 0
  if (.not.obj%yes) return
  if (obj%istp.ge.obj%nstp) return
! gather data
  if (wght.ne.0d0) then
    do ii=1,obj%no
      nn = obj%vo(ii,0)
      if (nn.gt.1) then
        factor = 0d0
        do jj=1,nn
          kk = obj%vo(ii,jj)
          factor = factor + obj%wtv(kk)*obj%dns(kk)
        enddo
        if (factor.ne.0d0) factor = wght*wght/factor ! variance
!       if (factor.ne.0d0) factor = wght/factor      ! entropy
        do jj=1,nn
          kk = obj%vo(ii,jj)
          obj%ave(kk) = obj%ave(kk) + factor*obj%dns(kk)
          obj%dns(kk) = 0d0
        enddo
      endif !(nn.gt.1)
    enddo !ii=1,obj%no
    obj%idat = obj%idat + 1
  endif
  obj%ntot = obj%ntot + 1
  if (obj%idat.eq.obj%ndat) then
! multi-channel step
    do ii=1,obj%no
      nn = obj%vo(ii,0)
      if (nn.gt.1) then
! Update channel weights
        hsum = 0d0
        do jj=1,nn
          kk = obj%vo(ii,jj)
          obj%ave(kk) = dsqrt( obj%ave(kk)/dble(obj%ntot) ) ! variance
!         obj%ave(kk) = obj%ave(kk)/dble(obj%ntot)          ! entropy
          hsum = hsum + dabs( obj%ave(kk) )
        enddo
        if (hsum.eq.0d0) cycle
        hsum = 0d0
        do jj=1,nn
          kk = obj%vo(ii,jj)
          obj%wtv(kk) = obj%wtv(kk)*obj%ave(kk)
          obj%ave(kk) = 0d0
          hsum = hsum + obj%wtv(kk)
        enddo
        thrs = obj%small/dble(nn)
        nrem = 0
        if (hsum.eq.0d0) hsum = 1d0
! Normalized adapted channels
        do jj=1,nn
          kk = obj%vo(ii,jj)
          obj%wtv(kk) = obj%wtv(kk)/hsum
          if (obj%wtv(kk).lt.thrs) then
! Tag channel for removal
            nrem = nrem+1
            jrem(nrem) = jj
          endif
        enddo
! Remove tagged channels
        if (obj%istp.eq.obj%nstp-1) then
        if (nrem.gt.0) then
          do jj=nrem,1,-1
            kk = obj%vo(ii,jrem(jj))
            nremtot = nremtot+1
            jremtot(nremtot) = kk
            obj%wtv(kk) = 0d0
            nn = nn-1
            do kk=jrem(jj),nn
              obj%vo(ii,kk) = obj%vo(ii,kk+1)
            enddo
          enddo
          obj%vo(ii,0) = nn 
! Re-normalize channels
          hsum = 0d0
          do jj=1,nn
            kk = obj%vo(ii,jj)
            hsum = hsum + obj%wtv(kk)
          enddo
          if (hsum.eq.0d0) hsum = 1d0
          do jj=1,nn
            kk = obj%vo(ii,jj)
            obj%wtv(kk) = obj%wtv(kk)/hsum
          enddo
        endif !(nrem.gt.0)
        endif !(obj%istp.eq.obj%nstp-1)
      endif !(nn.gt.1)
    enddo !ii=1,obj%no
    obj%idat = 0
    obj%ntot = 0
    obj%istp = obj%istp+1
  endif !(obj%idat.eq.obj%ndat)
  end subroutine
!
!
  subroutine mulcha_print( obj )
!********************************************************************
!********************************************************************
  implicit none 
  type(mulcha_type) ,intent(inout) :: obj
  integer :: i0,jj,nn,ii
!
  if (nunit.gt.0) then
    do i0=1,obj%no
      nn = obj%vo(i0,0)
      if (nn.gt.1) then
        write(nunit,'(a17,i4)') 'avh_kaleu_mulcha:',i0
        do jj=1,nn
          ii = obj%vo(i0,jj)
          write(nunit,'(a17,4x,i3,i4,d24.16)') &
            'avh_kaleu_mulcha:',jj,ii,obj%wtv(ii)
        enddo
      endif
    enddo
  endif
  end subroutine
!
!
  recursive subroutine mulcha_vrtx( obj ,vrtx,nn ,i0 ) 
!********************************************************************
!* Recursively create a list of splitting vertices
!********************************************************************
  implicit none
  type(mulcha_type) ,intent(in) :: obj
  integer :: nn,i0,vrtx(obj%n1)
  integer :: iv,i1,i2,i3
  call mulcha_split( obj ,iv ,i0 )
  if (iv.ne.0) then
    nn = nn+1
    vrtx(nn) = iv
    i1 = obj%ov(1,iv)
    i2 = obj%ov(2,iv)
    i3 = obj%ov(3,iv)
    if (i3.ne.0) i1 = i3 ! specific T-channel
    if (obj%frst1(iv)) then
      call mulcha_vrtx( obj ,vrtx,nn ,i1 ) 
      call mulcha_vrtx( obj ,vrtx,nn ,i2 ) 
    else
      call mulcha_vrtx( obj ,vrtx,nn ,i2 ) 
      call mulcha_vrtx( obj ,vrtx,nn ,i1 ) 
    endif
  endif
  end subroutine
!
  subroutine mulcha_split( obj ,iv ,i0 )
!*********************************************************************
!* Given  i0  choose random splitting vertex
!*********************************************************************
  implicit none
  type(mulcha_type) ,intent(in)  :: obj
  integer           ,intent(out) :: iv 
  integer           ,intent(in)  :: i0
  integer :: nn,ii
  real(kind(1d0)) :: xx,sumw
!
  nn = obj%vo(i0,0)
  if (nn.eq.0) then
! No splitting
    iv = 0
  elseif (nn.eq.1) then
! Only one splitting
    iv = obj%vo(i0,1)
  else
! Choose splitting
    call avh_random( xx )
    sumw = 0d0
    ii = 0
    if (obj%yes) then
      do while ( xx.gt.sumw .and. ii.lt.nn )
        ii = ii+1
        iv = obj%vo(i0,ii)
        sumw = sumw + obj%wtv(iv)
      enddo
      if (xx.gt.sumw) then
        write(*,*) 'ERROR in mulcha_split: no splitting chosen'
        stop
      endif
    else
      ii = 1 + int(nn*xx)
      iv = obj%vo(i0,ii)
    endif  
  endif
  end subroutine
!
end module
