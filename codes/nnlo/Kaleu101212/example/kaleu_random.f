      subroutine avh_random( rho )
!* ********************************************************************
!* * return a double precision random number
!* ********************************************************************
      implicit none
      double precision rho
!      write(*,*) 'ERROR in avh_random: please implement the random '
!     &          ,'number generator in the file random.f'
!      stop
! Out-comment the 3 lines above and put something like for example:
      real rvec
      call ranlux(rvec,1)
      rho = dble(rvec)
      end
