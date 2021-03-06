      subroutine avh_random( rho )
!* ********************************************************************
!* * return a double precision random number
!* ********************************************************************
      implicit none
      real*8 rho
      real*8 myrandom
      external myrandom
!      write(*,*) 'ERROR in avh_random: please implement the random '
!     &          ,'number generator in the file random.f'
!      stop
! Out-comment the 3 lines above and put something like for example:
!      real rvec
!      call ranlux(rvec,1)
!      rho = dble(rvec)
      rho = myrandom()
      end
!
      subroutine avh_randarr( rho , n )
!* ********************************************************************
!* * return an array of double precision random numbers
!* ********************************************************************
      implicit none
      integer n
      integer irand
      real*8 rho(n)
      real*8 myrandom
      external myrandom
!      write(*,*) 'ERROR in avh_random: please implement the random '
!     &          ,'number generator in the file random.f'
!      stop
! Out-comment the 3 lines above and put something like for example:
!      real rvec
!      call ranlux(rvec,1)
!      rho = dble(rvec)
      do irand=1,n
        rho(irand) = myrandom()
      end do
      end
