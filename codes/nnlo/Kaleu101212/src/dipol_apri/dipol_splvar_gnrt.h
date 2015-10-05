    call avh_random( gridx )
    call avh_random( gridz )
    call splvar_gnrt( obj%var(2*idip-1) ,xx ,gridx )
    call splvar_gnrt( obj%var(2*idip  ) ,zz ,gridz )
