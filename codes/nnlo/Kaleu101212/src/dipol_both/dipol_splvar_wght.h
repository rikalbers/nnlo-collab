      call splvar_wght( obj%var(2*idip-1) ,wsx,gridx ,xx )
      call splvar_wght( obj%var(2*idip  ) ,wsz,gridz ,zz )
      call grid_wght( obj%grid(2*idip-1) ,wgx ,gridx )
      call grid_wght( obj%grid(2*idip  ) ,wgz ,gridz )
      wx = wsx*wgx
      wz = wsz*wgz
