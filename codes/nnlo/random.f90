! TODO There is a module called random and now there is a function too
function myrandom() result(randnum)
use random
implicit none
!
!
  real(kind(1d0)) :: randnum
!
  randnum = gen_rand()
!
end function myrandom

