  allocate( obj%var(1:2*obj%ndip) )
  lbl = 'D       '
  do l1=1,obj%ndip
    lbl(2:3) = symb( obj%i(l1) )
    lbl(4:5) = symb( obj%j(l1) )
    lbl(6:7) = symb( obj%k(l1) )
    lbl(8:8) = 'x'
    call splvar_init( obj%var(2*l1-1) ,0.5d0 ,lbl )
    lbl(8:8) = 'z'
    call splvar_init( obj%var(2*l1  ) ,0.5d0 ,lbl )
  enddo
