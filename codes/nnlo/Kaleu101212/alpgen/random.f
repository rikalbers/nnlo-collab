      subroutine avh_random( rho )
!* ********************************************************************
!* * return a double precision random number
!* ********************************************************************
      implicit none
      double precision rho,rangen
!      write(*,*) 'ERROR in avh_random: please implement the random '
!     &          ,'number generator in the file random.f'
!      stop
! Out-comment the 3 lines above and put a line like for example:
      rho = rangen(1)
      end
