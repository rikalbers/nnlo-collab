  do ii=1,obj%ndip
    if (obj%mch%wch(ii).eq.0d0) cycle
    call grid_collect( obj%grid(2*ii-1) ,weight,-1d0 )
    call grid_collect( obj%grid(2*ii  ) ,weight,-1d0 )
  enddo
