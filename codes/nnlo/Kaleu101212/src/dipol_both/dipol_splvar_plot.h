    call splvar_plot( obj%var(2*ii-1) ,iunit )
    call splvar_plot( obj%var(2*ii  ) ,iunit )
    call grid_plot( obj%grid(2*ii-1) ,iunit )
    call grid_plot( obj%grid(2*ii  ) ,iunit )
