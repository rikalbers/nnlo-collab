  do ii=1,obj%ndip
    if (obj%mch%wch(ii).eq.0d0) cycle
    call splvar_collect( obj%var(2*ii-1) ,weight )
    call splvar_collect( obj%var(2*ii  ) ,weight )
  enddo
