module avh_qgrid
  use avh_kaleu_grid
  type(grid_type) ,allocatable ,save :: grid(:)
  logical         ,allocatable ,save :: actv(:)
  integer                      ,save :: ngrid=0
  integer         ,parameter :: nbatch=1000 ,nchmax=1000
  real(kind(1d0)) ,parameter :: frac=2d0/(2d0*nchmax-1d0)
end module

  subroutine avh_qgrid_close
!**********************************************************************
!**********************************************************************
  use avh_qgrid
  implicit none
  if (allocated(grid)) deallocate( grid )
  if (allocated(actv)) deallocate( actv )
  ngrid = 0
  end subroutine

  subroutine avh_qgrid_collect( ID ,xx ,ww )
!**********************************************************************
!**********************************************************************
  use avh_qgrid
  implicit none
  integer         ,intent(in) :: ID
  real(kind(1d0)) ,intent(in) :: xx,ww
  type(grid_type) :: tmpg(1:ID)
  logical         :: tmpa(1:ID)
  character(8) :: label
  if (ID.le.0) then
    write(*,*) 'ERROR in qgrid_collect: ID=',ID,', but should be positive'
    return
  endif
  if (ID.gt.ngrid) then
    if (.not.allocated(grid)) then
      allocate( grid(ID) )
      allocate( actv(ID) )
      actv(1:ID) = .false.
    else
      tmpg(1:ngrid) = grid(1:ngrid)
      tmpa(1:ngrid) = actv(1:ngrid)
      deallocate( grid )
      deallocate( actv )
      allocate( grid(ID) )
      allocate( actv(ID) )
      actv(1:ID) = .false.
      grid(1:ngrid) = tmpg(1:ngrid)
      actv(1:ngrid) = tmpa(1:ngrid)
    endif
    ngrid = ID
  endif
  if (.not.actv(ID)) then
    label(1:1) = 'Q'
    write(label(2:8),'(i7.7)') ID
    call grid_init( grid(ID) ,nbatch,nchmax,frac ,label )
    actv(ID) = .true.
  endif
  call grid_collect( grid(ID) ,ww ,xx )
  end subroutine

  subroutine avh_qgrid_plot( iunit )
!**********************************************************************
!**********************************************************************
  use avh_qgrid
  implicit none
  integer ,intent(in) :: iunit
  integer ::ii
  if (.not.allocated(grid)) return
  do ii=1,ngrid
    if (actv(ii)) call grid_plot( grid(ii) ,iunit )
  enddo
  end subroutine
