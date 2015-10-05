    call grid_gnrt( obj%grid(2*idip-1) ,gridx )
    call grid_gnrt( obj%grid(2*idip  ) ,gridz )
    call splvar_gnrt( obj%var(2*idip-1) ,xx ,gridx )
    call splvar_gnrt( obj%var(2*idip  ) ,zz ,gridz )
