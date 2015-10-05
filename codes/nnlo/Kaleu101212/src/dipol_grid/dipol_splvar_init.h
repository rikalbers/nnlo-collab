  allocate( obj%grid(1:2*obj%ndip) )
  lbl = 'D       '
  do l1=1,obj%ndip
    lbl(2:3) = symb( obj%i(l1) )
    lbl(4:5) = symb( obj%j(l1) )
    lbl(6:7) = symb( obj%k(l1) )
    lbl(8:8) = 'x'
    call grid_init( obj%grid(2*l1-1) ,nbatch_g,def_nchmax_g,def_fractn_g ,lbl )
    lbl(8:8) = 'z'
    call grid_init( obj%grid(2*l1  ) ,nbatch_g,def_nchmax_g,def_fractn_g ,lbl )
  enddo
