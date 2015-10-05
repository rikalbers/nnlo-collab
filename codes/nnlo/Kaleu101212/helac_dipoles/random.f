      subroutine avh_random( rho )
!* ********************************************************************
!* * return a double precision random number
!* ********************************************************************
      implicit none
      double precision rho
!      write(*,*) 'ERROR in avh_random: please implement the random '
!                ,'number generator in the file avh_random.f'
!      stop
! Out-comment the 3 lines above and put a line like for example:
!      call rangen(rho,1)
      double precision rnmy
      rho = rnmy(0)
      end
