!      integer kinv(0:2,avh_bint_maxv),nkinv,nn,idummy(10),ii
!      integer iev,nev1,nev2
!      parameter( nev1=20 000 000 )
!      parameter( nev2=1 )
!      call avh_bint_base_init( 10 )
!      nn = 4
!      do iev=1,nev1
!        call avh_bint_kinv_exp( kinv,nkinv ,idummy,nn )
!      enddo
!      do ii=1,nkinv
!        write(6,*) ii,':',kinv(0,ii),' =',kinv(1,ii),' +',kinv(2,ii)
!      enddo
!      write(6,*) ' '
!      do ii=1,nn
!        idummy(ii) = avh_bint_b(ii-1)
!      enddo
!      do iev=1,nev2
!        nkinv = 0
!        call avh_bint_kinv( kinv,nkinv ,idummy,nn )
!      enddo
!      do ii=1,nkinv
!        write(6,*) ii,':',kinv(0,ii),' =',kinv(1,ii),' +',kinv(2,ii)
!      enddo
!      end

module avh_bint
  private
  public :: avh_bint_b ,&
            avh_bint_l ,&
            avh_bint_base_init ,&
            avh_bint_kinv ,&
            avh_bint_antl_i ,&
            avh_bint_antl_n ,&
            avh_bint_esab ,&
            avh_bint_maxn ,&
            avh_bint_maxp ,&
            avh_bint_maxv
!
! maxn: 2 ,3  ,4  ,5   ,6   ,7    ,8    ,9     ,10    ,11     ,12     
! maxv: 6 ,25 ,90 ,301 ,966 ,3025 ,9330 ,28501 ,86526 ,261625 ,788970
  integer ,parameter :: avh_bint_maxn = 10
  integer ,parameter :: avh_bint_maxv = 86526
  integer ,parameter :: avh_bint_maxp = 2**(avh_bint_maxn+1)
!
  integer ,parameter :: nunit = 6 ! unit messages are send to
!
  integer ,dimension(0:avh_bint_maxn+1) :: avh_bint_b
  integer ,dimension(0:avh_bint_maxp  ) :: avh_bint_l
!
contains
!
      subroutine avh_bint_kinv( kinv,nkinv ,nn )
!  ********************************************************************
!  * Create list of "kinematical vertices"
!  * nn=4 : 5.5 times slower
!  * nn=9 : 4.3 times slower
!  ********************************************************************
      implicit none
      integer :: kinv(0:2,avh_bint_maxv),nkinv,nn
      integer :: aa(nn+1),bb(nn+1),ii,jj,kk,kinv0,kinv1,kinv2
      logical :: nexta,nextb
      nkinv = 0
      do ii=2,nn
        call avh_bint_antl_i(nexta,aa ,ii,nn)
        do while (nexta)
          kinv0 = 0
          do kk=1,ii
            kinv0 = kinv0 + avh_bint_b(aa(kk)-1)
          enddo
          do jj=1,ii/2
            call avh_bint_antl_i(nextb,bb ,jj,ii)
            do while (nextb)
              kinv1 = 0
              do kk=1,jj
                kinv1 = kinv1 + avh_bint_b(aa(bb(kk))-1)
              enddo
              kinv2 = kinv0-kinv1
              nkinv = nkinv+1
              kinv(0,nkinv) = kinv0
              if (kinv1.lt.kinv2) then
                kinv(1,nkinv) = kinv1
                kinv(2,nkinv) = kinv2
              else
                kinv(1,nkinv) = kinv2
                kinv(2,nkinv) = kinv1
              endif
            call avh_bint_antl_n(nextb,bb ,jj,ii)
            if (2*jj.eq.ii) nextb = (nextb.and.bb(1).eq.1)
            enddo
          enddo
        call avh_bint_antl_n(nexta,aa ,ii,nn)
        enddo
      enddo
      end subroutine

      
      subroutine avh_bint_base_init( nmax )
!  ********************************************************************
!  ********************************************************************
      implicit none
      integer :: nmax
      integer :: ilevel,lbl(nmax+2),ii,jj
      logical :: next
!
      if (avh_bint_maxn.lt.nmax) then
        write(*,*) &
          'ERROR in avh_bint_base_init:' ,&
          ' parameter avh_bint_maxn to small.' ,&
          ' Increase it to at least',nmax
        stop
      endif
!
      do ii=0,nmax+1
        avh_bint_b(ii) = 2**ii
      enddo
!
      avh_bint_l(0) = 0
      do ilevel=1,nmax+1
        call avh_bint_antl_i(next,lbl ,ilevel,nmax+1)
        do while (next)
          ii = 0
          do jj=1,ilevel
!c            ii = ii + avh_bint_b(lbl(jj))
            ii = ii + avh_bint_b(lbl(jj)-1)
          enddo
          avh_bint_l(ii) = ilevel
        call avh_bint_antl_n(next,lbl ,ilevel,nmax+1)
        enddo
      enddo
!
      end subroutine


      subroutine avh_bint_antl_i(next,lbl ,ndim,ntrm)
!  ********************************************************************
!  * Loop over all integer arrays lbl = (l_1,l_2,...,l_ndim) such that
!  *               l_1 < l_2 < ... < l_ndim
!  * and where each entry runs from minimally 1 to maximally ntrm.
!  * usage:
!  *
!  *       call avh_bint_antl_i(next,lbl ,ndim,ntrm)
!  *       do while (next)
!  *         write(6,*) (lbl(kdim),kdim=1,ndim)
!  *       call avh_bint_antl_n(next,lbl ,ndim,ntrm)
!  *       enddo
!  *
!  * "next" is a logical
!  * DO NOT CHANGE lbl inside the loop.
!  ********************************************************************
      implicit none
      integer :: ndim
      integer :: lbl(1+ndim),ntrm ,jdim
      logical :: next
      do jdim=1,ndim
        lbl(jdim) = jdim
      enddo
      if (ndim.eq.0) then
        next = .true.
      else
        next = (ndim.gt.0.and.ntrm.gt.0)
      endif
      end subroutine
!
      subroutine avh_bint_antl_n(next,lbl ,ndim,ntrm)
      implicit none
      integer :: ndim
      integer :: lbl(1+ndim),ntrm ,jdim,kdim
      logical :: next
      if (ndim.eq.0) then
        next = .false.
        return
      endif
      if (lbl(ndim).lt.ntrm) then
        lbl(ndim) = lbl(ndim)+1
        next = .true.
        return
      endif
      do jdim=ndim-1,1,-1
      if (lbl(jdim).lt.lbl(jdim+1)-1) then
        lbl(jdim) = lbl(jdim)+1
        do kdim=jdim+1,ndim
          lbl(kdim) = lbl(jdim)+(kdim-jdim)
        enddo
        next = .true.
        return
      endif
      enddo
      next = .false.
      end subroutine


      subroutine avh_bint_esab( l1,n1,l2,n2 ,ii ,nmax )
!  ********************************************************************
!  * Return the expansion in powers of 2 of  ii  and  2^nmax-1-ii
!  ********************************************************************
      implicit none
      integer :: nmax
      integer :: l1(nmax+1),n1,l2(nmax+1),n2,ii ,&
                 tmp,jj,dif,l1o(nmax+1),l2o(nmax+1)
!
      tmp = ii
      n1 = 0
      n2 = 0
      do jj=nmax,0,-1
        dif = tmp-avh_bint_b(jj)
        if (dif.ge.0) then
          n1 = n1+1
          l1o(n1) = jj
          tmp = dif
        else
          n2 = n2+1
          l2o(n2) = jj
        endif
      enddo
      do jj=1,n1
        l1(jj) = l1o(n1-jj+1)
      enddo
      do jj=1,n2
        l2(jj) = l2o(n2-jj+1)
      enddo
      end subroutine
!
end module
