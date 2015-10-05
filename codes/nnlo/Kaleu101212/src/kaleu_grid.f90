module avh_kaleu_grid
!
  private! Everything is private, except the following list
  public :: grid_type &
           ,grid_init &
           ,grid_close &
           ,grid_gnrt &
           ,grid_wght &
           ,grid_collect &
           ,grid_plot
!
! Character lenght label
  integer ,parameter :: nlbl = 8
! Message unit
  integer ,parameter :: nunit = 0
!
!
  type :: grid_type
    private
    real(kind(1d0)) ,allocatable :: xri(:),wch(:),sch(:)
    integer                      :: idat,ndat,nch,nchmax
    real(kind(1d0))              :: frac
    integer                      :: iw,ic
    character(nlbl)              :: lbl
  end type
!
!
contains
!
!
  subroutine grid_init( grid ,nbatch,nchmax,frac ,lbl )
!********************************************************************
!* nbatch = number of events to be collected before adaptation step
!* nchmax = the maximal number of bins/channels
!*   frac = fraction of the number of channels/bins to be merged in
!*          an adaptation step
!*    lbl = a character label for the grid
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  integer         ,intent(in)    :: nbatch,nchmax
  real(kind(1d0)) ,intent(in)    :: frac
  character(nlbl) ,intent(in)    :: lbl
  integer :: ii
  grid%lbl  = lbl
  grid%ndat = nbatch
  grid%idat = 0
  grid%nchmax = nchmax
  grid%frac   = frac
  grid%nch = 1
  grid%iw = 0
  grid%ic = 0
  if (grid%nchmax.le.1) then
    write(*,*) 'ERROR in grid_init: nchmax should be at least 2'
    stop
  endif
  allocate( grid%xri(0:nchmax) )
  allocate( grid%wch(1:nchmax) )
  allocate( grid%sch(0:nchmax) )
  grid%xri(0) = 0d0
  grid%sch(0) = 0d0
  do ii=1,grid%nch
    grid%xri(ii) = dble(ii)/dble(grid%nch)
    grid%wch(ii) = 0d0
    grid%sch(ii) = dble(ii)/dble(grid%nch)
  enddo
  end subroutine
!
  subroutine grid_close( grid )
!********************************************************************
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  if (allocated(grid%xri)) deallocate( grid%xri )
  if (allocated(grid%wch)) deallocate( grid%wch )
  if (allocated(grid%sch)) deallocate( grid%sch )
  end subroutine
!
  subroutine grid_gnrt( grid ,xx )
!********************************************************************
!* xx is output
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  real(kind(1d0)) ,intent(out)   :: xx
  integer         :: i0,i1,ii
  call avh_random( xx )
  i0 = 0
  i1 = grid%nch
  do while (i1-i0.gt.1)
    ii = (i0+i1)/2
    if (xx.le.grid%sch(ii)) then
      i1 = ii
    else
      i0 = ii
    endif
  enddo
  ii = i1
  grid%iw = ii
  grid%ic = ii
  call avh_random( xx )
  xx = ( grid%xri(ii)-grid%xri(ii-1) )*xx + grid%xri(ii-1)
  end subroutine
!
  subroutine grid_wght( grid ,weight ,xx )
!********************************************************************
!* Put xx=-1d0 to use the latest generated, and stored, channel
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  real(kind(1d0)) ,intent(out)   :: weight
  real(kind(1d0)) ,intent(in)    :: xx
  integer :: ii
  if (xx.lt.0d0) then
    ii = grid%iw 
  else
    call findch( grid ,ii ,xx )
    if (ii.eq.0) then
      if (nunit.gt.0) write(nunit,*) &
        'ERROR in grid_wght: no channel found, putting weight=0'
      weight = 0d0
      return
    endif
    grid%ic = ii
  endif
  weight = ( grid%xri(ii) - grid%xri(ii-1) ) &
         / ( grid%sch(ii) - grid%sch(ii-1) )
  end subroutine
!
  subroutine grid_collect( grid ,weight ,xx )
!********************************************************************
!* Put xx=-1d0 to update the latest active, and stored, channel
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  real(kind(1d0)) ,intent(in)    :: weight,xx
  integer :: ii,jj
  logical :: update
!
! Find the bin into which xx belongs.
! If xx<0, use the latest active (by gnrt or wght) bin.
  if (xx.lt.0d0) then
    ii = grid%ic
  else
    call findch( grid ,ii ,xx )
  endif
!
! Make sure the weight coming with a channel is used only once
  if (ii.eq.0) return
  grid%ic = 0
!
! Update the weights of the bins
  if (weight.ne.0d0) then
    grid%wch(ii) = grid%wch(ii) + weight    ! entropy
!   grid%wch(ii) = grid%wch(ii) + weight**2 ! variance
    grid%idat = grid%idat+1
  endif
!
! Adapt the bins
  update = .false.
  if     (grid%idat.eq.grid%ndat) then
    call splitch( grid )
    grid%idat = 0
    update = .true.
  elseif (grid%idat.eq.grid%ndat/2) then
    call mergech( grid )
    update = .true.
  endif
!
! Activate the updated bins for gnrt/wght
  if (update) then
    do jj=1,grid%nch
      grid%sch(jj) = grid%sch(jj-1) + grid%wch(jj)        ! entropy
!     grid%sch(jj) = grid%sch(jj-1) + dsqrt(grid%wch(jj)) ! variance
    enddo
    do jj=1,grid%nch
      grid%sch(jj) = grid%sch(jj)/grid%sch(grid%nch)
    enddo
    if (nunit.gt.0) write(nunit,*) &
      'MESSAGE from grid ',grid%lbl,': nch =',grid%nch
  endif
  end subroutine
!
!
  subroutine splitch( grid )
!********************************************************************
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  real(kind(1d0)) :: ineff,ineff_old,wchtmp,wchmax,xx
  integer         :: nch,ii,lch(grid%nchmax),nsame
  real(kind(1d0)) ,parameter :: ex=0.0d0
!
  ineff_old = 1d300
  do while (grid%nch.lt.grid%nchmax)
    nch = grid%nch
    wchmax = -1d0
    do ii=1,nch
      wchtmp = grid%wch(ii) !**(1d0-ex) * (grid%xri(ii)-grid%xri(ii-1))**ex
      if     (wchtmp.gt.wchmax) then
        wchmax = wchtmp
        lch(1) = ii
        nsame = 1
      elseif (wchtmp.eq.wchmax) then
        nsame = nsame+1
        lch(nsame) = ii
      endif
    enddo
!  
    ineff = wchmax*nch
    if (ineff.ge.ineff_old) then
      exit
    else
      ineff_old = ineff
    endif
!  
    if (nsame.gt.1) then
      call avh_random( xx )
      ii = lch(1+int(nsame*xx))
    else
      ii = lch(1)
    endif
!
    if (grid%xri(ii)-grid%xri(ii-1).le.epsilon(1d0)) exit
!  
    grid%xri(ii+1:nch+1) = grid%xri(ii:nch)
    grid%wch(ii+1:nch+1) = grid%wch(ii:nch)
    grid%nch = nch+1
    grid%xri(ii)   = ( grid%xri(ii-1)+grid%xri(ii) )/2
    grid%wch(ii)   = grid%wch(ii)/2
    grid%wch(ii+1) = grid%wch(ii)
  enddo
!
  end subroutine
!
!
  subroutine mergech( grid )
!********************************************************************
!********************************************************************
  implicit none
  type(grid_type) ,intent(inout) :: grid
  real(kind(1d0)) :: wchtmp,wchmin,xx
  integer         :: nch,ii,lch(grid%nchmax),nsame,newnch
  real(kind(1d0)) ,parameter :: ex=0.0d0
!
  newnch = 1+int( (1d0-grid%frac)*dble(grid%nch) )
  do while (grid%nch.gt.newnch)
    nch = grid%nch
    wchmin = 1d300
    do ii=1,nch-1
!      wchtmp = grid%wch(ii)+grid%wch(ii+1)
      wchtmp = grid%xri(ii+1)-grid%xri(ii-1)
      if     (wchtmp.lt.wchmin) then
        wchmin = wchtmp
        lch(1) = ii
        nsame = 1
      elseif (wchtmp.eq.wchmin) then
        nsame = nsame+1
        lch(nsame) = ii
      endif
    enddo
!  
    if (nsame.gt.1) then
      call avh_random( xx )
      ii = lch(1+int(nsame*xx))
    else
      ii = lch(1)
    endif
!
    grid%xri(ii) = grid%xri(ii+1)
    grid%wch(ii) = grid%wch(ii)+grid%wch(ii+1)
    grid%xri(ii+1:nch-1) = grid%xri(ii+2:nch)
    grid%wch(ii+1:nch-1) = grid%wch(ii+2:nch)
    grid%nch = nch-1
  enddo
!
  end subroutine
!
!
  subroutine findch( grid ,ii ,xx )
!********************************************************************
!********************************************************************
  implicit none
  type(grid_type) ,intent(in)   :: grid
  integer         ,intent(out)  :: ii
  real(kind(1d0)) ,intent(in)   :: xx
  integer i0,i1
!
  if (xx.lt.0d0.or.xx.gt.1d0) then
    ii = 0
    return
  endif
!
  i0 = 0
  i1 = grid%nch
  do while (i1-i0.gt.1)
    ii = (i0+i1)/2
    if (xx.le.grid%xri(ii)) then
      i1 = ii
    else
      i0 = ii
    endif
  enddo
  ii = i1
  end subroutine
!
!
  subroutine grid_plot( grid  ,iunit )
!********************************************************************
!********************************************************************
  implicit none
  type(grid_type) ,intent(in) :: grid
  integer         ,intent(in) :: iunit
  character(4+nlbl) :: filename
  real(kind(1d0)) ,parameter :: o0 = 0d0
  real(kind(1d0)) :: h1,h2
  integer :: ii
  if (iunit.le.0) return
  if (grid%nchmax.le.1) return
  filename = 'grid'//grid%lbl
  open( unit=iunit, file=filename, status="replace" )
  do ii=1,grid%nch
    h1 = ( grid%sch(ii) - grid%sch(ii-1) ) &
       / ( grid%xri(ii) - grid%xri(ii-1) )
    h2 = grid%sch(ii)
    write(iunit,'(99e16.8)') grid%xri(ii-1),o0,o0
    write(iunit,'(99e16.8)') grid%xri(ii-1),h1,h2
    write(iunit,'(99e16.8)') grid%xri(ii  ),h1,h2
    write(iunit,'(99e16.8)') grid%xri(ii  ),o0,o0
  enddo
  close( unit=iunit )
  end subroutine
!
!
  function grid_get_nch( grid ) result(value)
!********************************************************************
!********************************************************************
  implicit none
  type(grid_type) ,intent(in) :: grid
  integer :: value
  value = grid%nch
  end function
!
!
end module
