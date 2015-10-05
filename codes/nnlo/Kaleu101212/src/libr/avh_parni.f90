!***********************************************************************
!                                                                      *
!                          This is PARNI95                             *
!                                                                      *
!          for importance sampling and/or denstiy estimation.          *
!                                                                      *
! author: Andreas van Hameren <Andre.HamerenREMOVETHIS@ifj.edu.pl>     *
!   date: 18-12-2009                                                   *
!                                                                      *
! Please cite                                                          *
!    A. van Hameren, 
!    Acta Phys.Polon.B40:259-272(2009), arXiv:0710.2448 [hep-ph]       *
! in publications with results obtained with the help of this program. *
!***********************************************************************

module avh_parni
  private ! Everything is private, except the following list
  public :: avh_parni_type ,&
            avh_parni_init ,&
            avh_parni_generate ,&
            avh_parni_calcwght ,&
            avh_parni_density ,&
            avh_parni_weight ,&
            avh_parni_adapt ,&
            avh_parni_adaptone ,&
            avh_parni_delete ,&
            avh_parni_charnum ,&
            avh_parni_result ,&
            avh_parni_marg ,&
            avh_parni_plot ,&
            avh_parni_dim
!
! The maximal number of dimension:
  integer ,parameter :: avh_parni_dim = 2
! The minimal allowed volume for a box, as 2**(-size_ncut):
  integer ,parameter :: size_ncut = 52
! If the following (overflow) parameter is 0, channel creation will 
! stop if the maximum is reached. If it is >0, creation will continue,
! but combined with merging:
  integer ,parameter :: size_over = 8
! The maximum number of channels:
!  integer ,parameter :: size_nch  = 200                    !FIXED
!  integer ,parameter :: size_tot  = (size_nch+size_over)*2 !FIXED
!
  type :: avh_parni_type
    private
    real(kind(1d0)) ,allocatable :: weight(:) ,& !ALLOC
                                    sumwt1(:)    !ALLOC
    integer         ,allocatable :: mother(:) ,& !ALLOC
                                    daught(:) ,& !ALLOC
                                    cutdim(:) ,& !ALLOC
                                    channl(:) ,& !ALLOC
                                    nodech(:) ,& !ALLOC
                                    pointr(:)    !ALLOC
!    real(kind(1d0)) :: weight(size_tot) ,& !FIXED
!                       sumwt1(size_tot)    !FIXED
!    integer         :: mother(size_tot) ,& !FIXED
!                       daught(size_tot) ,& !FIXED
!                       cutdim(size_tot) ,& !FIXED
!                       channl(size_tot) ,& !FIXED
!                       nodech(size_tot) ,& !FIXED
!                       pointr(size_tot)    !FIXED
    real(kind(1d0)) :: xstr(avh_parni_dim),wstr ,ineff_old
    integer :: nchann ,&
               nchmax ,&
               ndimen ,&
               nbatch ,&
               nsofar ,&
               iitask
    logical :: do_adapt
  end type
!
  real(kind(1d0)) ,parameter :: verylarge = 1d300
  real(kind(1d0)) :: halfpwr(0:size_ncut)
!
contains
!
      subroutine avh_parni_init(obj,task,ndimen_in,nbatch_in,nchmax_in)
!* *********************************************************************
!* * All arguments are input.
!* * "obj" is the label of the 'copy' of the algorithm.
!* * "task" sets the kind of optimizations used internally. For the
!* *    moment the options are 11,12 ,16,17.
!* *    See extended explanation below. In case of importance sampling
!* *    task=11 seems to be the best choice. In case of density
!* *    estimation, task MUST be set to 11 or 16.
!* * "ndimen_in" is the number of dimensions of the hypercube.
!* * "nbatch_in" is the size of the batches of data after which an
!* *    optimization step is performed. The square-root of the number
!* *    of data points that is going to be used for the optimization
!* *    seems to be a good choice. The larger this number, the more
!* *    conservative the algorithm, i.e. the smaller the number of
!* *    channels.
!* * "nchmax_in" is the maximal number of channels. If put to 0, the
!* *    number of channels keeps growing. If larger than 0, boxes
!* *    originating from the same parent with the smallest weight are
!* *    merged if nchmax is exeeded after splitting boxes. This process
!* *    continues: boxes are splitted and others merged and the number 
!* *    of channels is kept around nchmax.
!* *
!* * Let f be the integrand and denote
!* *
!* *   < f > = integral( f(x), whole hypercube )
!* *
!* * Furthermore, the density built by parni is given by
!* *
!* *   p(x) = sum_{i=1}^{nchannel} w_i*t_i(x)/v_i
!* *
!* * where t_i is the characteristic function of the hyperbox
!* * corresponding to channel i, and v_i is the volume of that box.
!* * For the different tasks, the weights w_i are optimized as follows:
!* *
!* * task=11,16: w_i = < f*t_i >                  (histogram-like)
!* * task=12,17: w_i = sqrt( < f^2*t_i > * v_i )  (minimal variance)
!* *
!* * divided by the normalisation such that their sum is equal to 1.
!* * Channels are dissected if they have the largest w_i or w_i*v_i
!* *
!* * task=11,12: dissect channel with largest w_i
!* * task=16,17: dissect channel with largest w_i*v_i
!* *
!* * The integrals are estimated as follows: (up to an overall
!* * normalization)
!* *
!* * task<10:  < f*t_i > = sum( f(xj)*t_i(xj) )/sum( t_i(xj) )
!* * task>10:  < f*t_i > = sum( f(xj)*t_i(xj) )*v_i/w_i
!* *
!* * and the same for f^2. The sum is over all events. For the
!* * adaptation, this means the following. Let W(xj) be the numbers
!* * going into  avh_parni_adapt( ID, W(xj), xj ) , then
!* * 
!* * task<1 :  w_i <-- sum( W(xj)*t_i(xj)*w_i/v_i )/sum( t_i(xj) )
!* * task>11:  w_i <-- sum( W(xj)*t_i(xj) )
!* *
!* * summary
!* *
!* * task=11 : histogram-like  , largest w_i    
!* * task=12 : minimal variance, largest w_i    
!* * task=16 : histogram-like  , largest w_i*v_i
!* * task=17 : minimal variance, largest w_i*v_i
!* *
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      integer ,intent(in) :: task,ndimen_in,nbatch_in,nchmax_in
      integer :: node,ii,isize
      logical :: firstcall = .true.
!
      if (firstcall) then
        firstcall = .false.
        write(*,'(a72)') '########################################################################'
        write(*,'(a72)') '#                                                                      #'
        write(*,'(a72)') '#                        You are using PARNI95                         #'
        write(*,'(a72)') '#                                                                      #'
        write(*,'(a72)') '#          for importance sampling and/or density estimation.          #'
        write(*,'(a72)') '#                                                                      #'
        write(*,'(a72)') '# author: Andreas van Hameren <hamerenREMOVETHIS@ifj.edu.pl>           #'
        write(*,'(a72)') '#   date: 05-02-2010                                                   #'
        write(*,'(a72)') '#                                                                      #'
        write(*,'(a72)') '# Please cite                                                          #'
        write(*,'(a72)') '#    A. van Hameren,                                                   #'
        write(*,'(a72)') '#    Acta Phys.Polon.B40:259-272(2009), arXiv:0710.2448 [hep-ph]       #'
        write(*,'(a72)') '# in publications with results obtained with the help of this program. #'
        write(*,'(a72)') '#                                                                      #'
        write(*,'(a72)') '########################################################################'
!
        halfpwr(0) = 1d0
        do ii=1,size_ncut
          halfpwr(ii) = halfpwr(ii-1)/2d0
        enddo
!
      endif
!
!      if (nchmax_in.gt.size_nch) then                          !FIXED
!        write(*,*) 'ERROR in avh_parni_init: increase the ' ,& !FIXED
!                   'parameter size_nch to at least',nchmax_in  !FIXED
!        stop                                                   !FIXED
!      endif                                                    !FIXED
!
      obj%ndimen = ndimen_in
      obj%nbatch = nbatch_in
      obj%nchmax = nchmax_in
!
      isize = 2*(nchmax_in+size_over) ! possibility to overflow !ALLOC
      allocate( obj%weight(isize) )                             !ALLOC
      allocate( obj%sumwt1(isize) )                             !ALLOC
      allocate( obj%mother(isize) )                             !ALLOC
      allocate( obj%daught(isize) )                             !ALLOC
      allocate( obj%cutdim(isize) )                             !ALLOC
      allocate( obj%channl(isize) )                             !ALLOC
      allocate( obj%nodech(isize) )                             !ALLOC
      allocate( obj%pointr(isize) )                             !ALLOC
!      isize = size_tot                                         !FIXED
!
      do node=1,isize
        obj%cutdim(node) = 0
        obj%daught(node) = 0
        obj%sumwt1(node) = 0d0
        obj%pointr(node) = node
      enddo
!
      obj%nchann = 1
      obj%mother(1) = 0
      obj%daught(1) = 0
      obj%channl(1) = 1
      obj%nodech(1) = 1
      obj%pointr(1) = 1
      obj%weight(1) = 1d0
!
      obj%nsofar = 0
      obj%iitask = task
      obj%ineff_old = verylarge
      obj%do_adapt = .true.
      end subroutine


      subroutine avh_parni_delete( obj )                         !ALLOC
!* ********************************************************************
!* ********************************************************************
      implicit none                                              !ALLOC
      type(avh_parni_type) ,intent(inout) :: obj                 !ALLOC
      if (allocated(obj%weight)) deallocate( obj%weight )        !ALLOC
      if (allocated(obj%sumwt1)) deallocate( obj%sumwt1 )        !ALLOC
      if (allocated(obj%mother)) deallocate( obj%mother )        !ALLOC
      if (allocated(obj%daught)) deallocate( obj%daught )        !ALLOC
      if (allocated(obj%cutdim)) deallocate( obj%cutdim )        !ALLOC
      if (allocated(obj%channl)) deallocate( obj%channl )        !ALLOC
      if (allocated(obj%nodech)) deallocate( obj%nodech )        !ALLOC
      if (allocated(obj%pointr)) deallocate( obj%pointr )        !ALLOC
      end subroutine                                             !ALLOC

!      subroutine avh_parni_delete( obj )                        !FIXED
!**********************************************************************
!**********************************************************************
!      implicit none                                             !FIXED
!      type(avh_parni_type) ,intent(inout) :: obj                !FIXED
!      end subroutine                                            !FIXED


      subroutine avh_parni_weight( obj ,ww,xx )
!* ********************************************************************
!* * Return the latest generated point "xx" and corresponding weight
!* * "ww".
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(out) :: xx(avh_parni_dim),ww
      integer :: idm
!
      do idm=1,obj%ndimen
        xx(idm) = obj%xstr(idm)
      enddo
      ww = obj%wstr
      end subroutine


      subroutine avh_parni_calcwght( obj ,ww ,xx )
!* *********************************************************************
!* * Returns the inverse value of of the density built by the
!* * algorithm at the point "xx" in the hypercube [0,1]^ndimen
!* * This routine also stores (overwrites) both "ww" and "xx"
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(in) :: xx(avh_parni_dim)
      real(kind(1d0)) ,intent(out) :: ww
!      real(kind(1d0)) :: halfpwr
      integer :: ich,nstp,idm
!
      call channel( obj ,ich,nstp ,xx )
!
      do idm=1,obj%ndimen
        obj%xstr(idm) = xx(idm)
      enddo
      obj%wstr = halfpwr(nstp)/obj%weight(obj%pointr(ich))
!
      ww = obj%wstr
      end subroutine


      real(kind(1d0)) function avh_parni_density( obj ,xx )
!* *********************************************************************
!* * Returns the value of the density built by the algorithm at
!* * the point "xx" in the hypercube [0,1]^ndimen
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(in) :: xx(avh_parni_dim)
!      real(kind(1d0)) :: halfpwr
      integer :: ich,nstp
!
      call channel( obj ,ich,nstp ,xx )
      avh_parni_density = obj%weight(obj%pointr(ich)) / halfpwr( nstp )
      end function 


      subroutine channel( obj ,ich,nstp ,xx )
!* *********************************************************************
!* * Returns the number "ich" of the channel containing "xx".
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(in) :: xx(avh_parni_dim)
      integer ,intent(out) :: ich,nstp
      real(kind(1d0)) :: x0(avh_parni_dim),x1(avh_parni_dim),xh
      integer :: node,dght,idm
!
      node = 1
      dght = obj%daught(node)
      do idm=1,obj%ndimen
        x0(idm) = 0d0
        x1(idm) = 1d0
      enddo
      nstp = 0
      do while(dght.ne.0)
        idm = obj%cutdim(node)
        xh  = (x0(idm)+x1(idm))/2
        if (xx(idm).lt.xh) then
          node = dght
          x1(idm) = xh
        else
          node = dght+1
          x0(idm) = xh
        endif
        nstp = nstp+1
        dght = obj%daught(node)
      enddo
      ich = obj%channl(node)
      end subroutine


      subroutine avh_parni_generate( obj, xx )
!* *********************************************************************
!* *  "xx" is a data point, distributed following the density
!* *  built by the algorithm.
!* * The value of xx and the corresponding weight are stored.
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(out) :: xx(avh_parni_dim)
      real(kind(1d0)) x0(avh_parni_dim),x1(avh_parni_dim),rho,rch,wght
      integer :: ich,idm,nstp,node,jch,nch
!
      nch = obj%nchann
!
      call avh_random( rch )
      ich = 2*nch-1
      do while(ich.gt.nch)
        jch = 2*(ich-nch)-1
        wght = obj%weight(obj%pointr(jch))
        if (rch.ge.wght) then
          rch = rch - wght
          ich = jch+1
        else
          ich = jch
        endif
      enddo
!
      node = obj%nodech(ich)
      call get_box(obj ,x0,x1,nstp ,node)
!
      do idm=1,obj%ndimen
        call avh_random( rho )
        xx(idm) = x0(idm) + (x1(idm)-x0(idm))*rho
      enddo
!
      do idm=1,obj%ndimen
        obj%xstr(idm) = xx(idm)
      enddo
      obj%wstr = halfpwr(nstp)/obj%weight(obj%pointr(ich))
      end subroutine

      subroutine get_box( obj ,x0,x1,nstp ,node )
!* ********************************************************************
!* * Return the box that corresponds to "node"
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(in) :: obj
      integer ,intent(in) :: node
      real(kind(1d0)) ,intent(out) :: x0(avh_parni_dim),x1(avh_parni_dim)
      integer ,intent(out) :: nstp
      integer :: idm,istp,dimnsn(size_ncut)
      logical :: upperh(size_ncut)
!
      call history(obj ,dimnsn,upperh,nstp ,node)
      do idm=1,obj%ndimen
        x0(idm) = 0d0
        x1(idm) = 1d0
      enddo
      do istp=nstp,1,-1
        idm = dimnsn(istp)
        if (upperh(istp)) then
          x0(idm) = (x0(idm)+x1(idm))/2
        else
          x1(idm) = (x0(idm)+x1(idm))/2
        endif
      enddo
      end subroutine

      subroutine history( obj ,dimnsn,upperh,nstp ,node_in )
!* ********************************************************************
!* * Return the history of cutting dimensions and choosing halfs
!* * that led to "node_in"
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(in) :: obj
      integer ,intent(in) :: node_in
      integer ,intent(out) :: dimnsn(size_ncut),nstp
      logical ,intent(out) :: upperh(size_ncut)
      integer :: mthr,node
!
      nstp = 0
      node = node_in
      mthr = obj%mother(node)
      do while (mthr.ne.0)
        nstp = nstp+1
        dimnsn(nstp) = obj%cutdim(mthr)
        upperh(nstp) = (node.ne.2*(node/2))
        node = mthr
        mthr = obj%mother(node)
      enddo
      end subroutine

      real(kind(1d0)) function volume( obj ,node_in )
!* ********************************************************************
!* * Volume of box at node_in
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      integer ,intent(in) :: node_in
      integer :: nstp,mthr,node
!      real(kind(1d0)) :: halfpwr
!
      nstp = 0
      node = node_in
      mthr = obj%mother(node)
      do while (mthr.ne.0)
        nstp = nstp+1
        node = mthr
        mthr = obj%mother(node)
      enddo
      volume = halfpwr(nstp)
      end function


      subroutine avh_parni_adaptone( obj ,wdat )
!* ********************************************************************
!* * Collect "wdat" for adaptation of the instance ID of PARNI
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(in) :: wdat
      real(kind(1d0)) :: xx(avh_parni_dim)
      integer :: idm
!
      do idm=1,obj%ndimen
        xx(idm) = obj%xstr(idm)
      enddo
      call avh_parni_adapt( obj ,wdat ,xx )
      end subroutine


      subroutine avh_parni_adapt( obj ,wdat,xdat )
!* *********************************************************************
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) ,intent(in) :: wdat,xdat(avh_parni_dim)
      real(kind(1d0)) cnst
      integer :: node,ich,jch,kch,nstp,nnew
      logical :: cont
!
      call channel( obj ,ich,nstp ,xdat )
      kch = obj%pointr(obj%nchann+ich)
      if     (obj%iitask.eq.11.or.obj%iitask.eq.16) then
        obj%sumwt1(kch) = obj%sumwt1(kch) + wdat
      elseif (obj%iitask.eq.12.or.obj%iitask.eq.17) then
        obj%sumwt1(kch) = obj%sumwt1(kch) + wdat*wdat
      endif
      obj%nsofar = obj%nsofar + 1
!
      if (obj%nsofar.lt.obj%nbatch) return
      obj%nsofar = 0
!
      do ich=1,obj%nchann
        jch = obj%pointr(ich)
        kch = obj%pointr(obj%nchann+ich)
        obj%sumwt1(jch) = obj%sumwt1(jch) + obj%sumwt1(kch)
      enddo
!
      if     (obj%iitask.eq.11.or.obj%iitask.eq.16) then
        do ich=1,obj%nchann
          jch = obj%pointr(ich)
          obj%weight(jch) = obj%sumwt1(jch)
        enddo
      elseif (obj%iitask.eq.12.or.obj%iitask.eq.17) then
        do ich=1,obj%nchann
          jch = obj%pointr(ich)
          obj%weight(jch) = dsqrt( obj%sumwt1(jch) )
        enddo
      endif
!
      if (obj%do_adapt) then
        cont = (obj%nchann.lt.obj%nchmax+size_over)
        do while (cont)
          call findmax( obj ,cont )
          if (cont) call dissect( obj )
          cont = (cont.and.(obj%nchann.lt.obj%nchmax+size_over))
        enddo
!
        if (obj%nchmax.gt.0) then
        do while (obj%nchann.gt.obj%nchmax)
          call mergebox( obj )
        enddo
        endif
      endif
!
      do ich=1,obj%nchann
        kch = obj%pointr(obj%nchann+ich)
        obj%sumwt1(kch) = 0d0
        obj%daught(obj%nodech(ich)) = 0
      enddo
!
! Add the weights together in a tree-structure, and normalize them
      nnew = obj%nchann
      do ich=1,2*obj%nchann-3,2
        nnew = nnew+1
        obj%weight(obj%pointr(nnew)) &
        = obj%weight(obj%pointr(ich)) + obj%weight(obj%pointr(ich+1))
      enddo
      cnst = 1d0/obj%weight(obj%pointr(nnew))
      do ich=1,2*obj%nchann-1
        jch = obj%pointr(ich)
        obj%weight(jch) = obj%weight(jch) * cnst
      enddo
!
      end subroutine


      subroutine findmax( obj ,cont )
!* *********************************************************************
!* * determine the channel with the largest weight and put it last
!* *********************************************************************
      implicit none 
      type(avh_parni_type) ,intent(inout) :: obj
      logical ,intent(inout) :: cont
      integer :: ichmax,node,ich,lch(2*obj%nchann),nsame,nch
      real(kind(1d0)) :: wchmax,ineff,rho,vlm,weightmax,wch
      logical :: case1
      real(kind(1d0)) ,parameter :: par1=0.5d0
      real(kind(1d0)) ,parameter :: par2=1d0-par1
!
      nch = obj%nchann
!
      nsame = 1
      wchmax = -1d0
      weightmax = -1d0
      case1 = (obj%iitask.le.15)
!
      if (case1) then
        do ich=1,nch
          wch = obj%weight(obj%pointr(ich))
          if     (wch.gt.wchmax) then
            wchmax = wch
            lch(1) = ich
            nsame = 1
          elseif (wch.eq.wchmax) then
            nsame = nsame+1
            lch(nsame) = ich
          endif
          weightmax = wchmax
        enddo
      else
        do ich=1,nch
          wch = obj%weight(obj%pointr(ich))
          vlm = volume(obj,obj%nodech(ich))
          wch = wch**par1 * vlm**par2
          if     (wch.gt.wchmax) then
            wchmax = wch
            lch(1) = ich
            nsame = 1
          elseif (wch.eq.wchmax) then
            nsame = nsame+1
            lch(nsame) = ich
          endif
          weightmax = wchmax
        enddo
      endif
!
      ineff = weightmax*nch
!
      if (ineff.ge.obj%ineff_old) then
        obj%ineff_old = verylarge
        cont = .false.
      else
        obj%ineff_old = ineff
        cont = .true.
      endif
!
      if (nsame.gt.1) then
        call avh_random( rho )
        ichmax = lch(1+int(nsame*rho))
      else
        ichmax = lch(1)
      endif
!
      call switch(obj ,ichmax,nch)
      end subroutine


      subroutine dissect( obj )
!* *********************************************************************
!* * dissect the last channel into two new ones.
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) :: rho,stdev,w1old,weold
      integer :: node,idm,istp,nstp,ncut(avh_parni_dim),ncutmin,maxidm ,&
                 nch,dimnsn(size_ncut),nsame,ldm(avh_parni_dim),node1 ,&
                 node2,jch1,jch2
      logical :: init,upperh(size_ncut)
!
      nch = obj%nchann
      node = obj%nodech(nch)
!
! determine the dimension in which to dissect box 1:
! determine largerst dimension in box 1,
! which is the one with the fewest number of cuts
      call history(obj ,dimnsn,upperh,nstp ,node)
      if (nstp.ge.size_ncut) return ! restrict minimal bin-size
      do idm=1,obj%ndimen
        ncut(idm) = 0
      enddo
      do istp=1,nstp
        ncut(dimnsn(istp)) = ncut(dimnsn(istp)) + 1
      enddo
      maxidm = 1
      do idm=1,obj%ndimen
        if (ncut(idm).lt.nstp) then
          nstp = ncut(idm)
          maxidm = idm
        endif
      enddo
! if several dimensions have this size, choose one at random
      nsame = 0
      do idm=1,obj%ndimen
        if (ncut(idm).eq.nstp) then
          nsame = nsame+1
          ldm(nsame) = idm
        endif
      enddo
      if (nsame.gt.1) then
        call avh_random( rho )
        maxidm = ldm( 1+int(nsame*rho) )
      endif
!
      jch1 = obj%pointr(nch)
      weold = obj%weight(jch1)
      w1old = obj%sumwt1(jch1)
!
      node1 = 2*nch
      node2 = 2*nch+1 
      obj%cutdim(node ) = maxidm
      obj%daught(node ) = node1
      obj%mother(node1) = node
      obj%daught(node1) = 0
      obj%cutdim(node1) = 0
      obj%mother(node2) = node
      obj%daught(node2) = 0
      obj%cutdim(node2) = 0
!
      obj%nodech( 2*nch)   = node1
      obj%nodech( 2*nch+1) = node2
      obj%channl( node1 )  = 2*nch
      obj%channl( node2 )  = 2*nch+1
      call switch(obj ,nch  ,2*nch  )
      call switch(obj ,nch+1,2*nch+1)
!
      jch1 = obj%pointr(nch)
      jch2 = obj%pointr(nch+1)
      obj%sumwt1(jch1) = w1old/2d0
      obj%sumwt1(jch2) = w1old/2d0
      obj%weight(jch1) = weold/2d0
      obj%weight(jch2) = weold/2d0
!
      obj%nchann = obj%nchann+1
      end subroutine


      subroutine mergebox( obj )
!* *********************************************************************
!* * determine the pair of nodes with the smallest weights
!* * and merge the boxes to one box
!* *********************************************************************
      implicit none 
      type(avh_parni_type) ,intent(inout) :: obj
      real(kind(1d0)) :: wchmin,rho,wchtmp
      integer :: ich1,ich2,node,ich,nsame,jch,nch,sstr,mthr ,&
                 jch1,jch2,lch(2*obj%nchann)
!
      nch = obj%nchann
!
      nsame = 1
      lch(1) = 1
      wchmin = verylarge
      do ich=1,nch
        node = obj%nodech(ich)
        sstr = (node/2)*2
        if (node.eq.sstr) sstr = sstr+1
        if (obj%daught(sstr).eq.0) then
          jch = obj%channl(sstr)
          wchtmp = obj%weight(obj%pointr(ich)) &
                 + obj%weight(obj%pointr(jch))
          if     (wchtmp.lt.wchmin) then
            wchmin = wchtmp
            lch(1) = ich
            nsame = 1
          elseif (wchtmp.eq.wchmin) then
            nsame = nsame+1
            lch(nsame) = ich
          endif
        endif
      enddo
!
      if (nsame.gt.2) then
        call avh_random( rho )
        ich1 = lch(1+int(nsame*rho))
      else
        ich1 = lch(1)
      endif
      node = obj%nodech(ich1)
      sstr = (node/2)*2
      if (node.ne.sstr) node = sstr
      sstr = sstr+1
!
      ich1 = obj%channl(node)
      ich2 = obj%channl(sstr)
!
      mthr = obj%mother(node)
      ich  = obj%channl(mthr)
      jch  = obj%pointr(ich)
      jch1 = obj%pointr(ich1)
      jch2 = obj%pointr(ich2)
      obj%weight(jch) = obj%weight(jch1) + obj%weight(jch2)
      obj%sumwt1(jch) = obj%sumwt1(jch1) + obj%sumwt1(jch2)
      obj%cutdim(mthr) = 0
      obj%daught(mthr) = 0
!
      call switch(obj ,ich ,2*nch-1)
      call switch(obj ,ich1,2*nch-1)
      call switch(obj ,ich2,  nch  )
      call switch(obj ,nch ,2*nch-2)
!
      call rename(obj ,node,2*nch-2)
!
      obj%nchann = obj%nchann-1
      end subroutine


      subroutine switch( obj ,ich,jch )
!* ********************************************************************
!* * Switch number for channel ich and jch
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      integer ,intent(in) :: ich,jch
      integer :: itmp,ii,jj
!
      ii = ich
      jj = jch
!
      itmp       = obj%pointr(ii)
      obj%pointr(ii) = obj%pointr(jj)
      obj%pointr(jj) = itmp
      itmp       = obj%nodech(ii)
      obj%nodech(ii) = obj%nodech(jj)
      obj%nodech(jj) = itmp
      obj%channl(obj%nodech(ii)) = ich
      obj%channl(obj%nodech(jj)) = jch
      end subroutine


      subroutine rename( obj ,node1_in,node2_in )
!* ********************************************************************
!* * Switch the numbers of node1 and her sister with those of node2
!* * and her sister
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      integer ,intent(in) :: node1_in,node2_in
      integer :: node1,node2,sstr1,sstr2,n1,n2,s1,s2,itmp
!
      node1 = node1_in
      node2 = node2_in
!
      sstr1 = (node1/2)*2
      if (node1.ne.sstr1) node1 = sstr1
      sstr1 = sstr1+1
      sstr2 = (node2/2)*2
      if (node2.ne.sstr2) node2 = sstr2
      sstr2 = sstr2+1
!
      if (node1.eq.node2) return
      n1 = node1
      s1 = sstr1
      n2 = node2
      s2 = sstr2
      itmp       = obj%mother(n1)
      obj%mother(n1) = obj%mother(n2)
      obj%mother(n2) = itmp
      itmp       = obj%daught(n1)
      obj%daught(n1) = obj%daught(n2)
      obj%daught(n2) = itmp
      itmp       = obj%cutdim(n1)
      obj%cutdim(n1) = obj%cutdim(n2)
      obj%cutdim(n2) = itmp
      itmp       = obj%channl(n1)
      obj%channl(n1) = obj%channl(n2)
      obj%channl(n2) = itmp
      obj%nodech(obj%channl(n1)) = node1
      if (obj%daught(n1).ne.0) then
        obj%mother(obj%daught(n1)  ) = node1
        obj%mother(obj%daught(n1)+1) = node1
      endif
      if (obj%mother(n1).ne.0) obj%daught(obj%mother(n1)) = node1
      obj%nodech(obj%channl(n2)) = node2
      if (obj%daught(n2).ne.0) then
        obj%mother(obj%daught(n2)  ) = node2
        obj%mother(obj%daught(n2)+1) = node2
      endif
      if (obj%mother(n2).ne.0) obj%daught(obj%mother(n2)) = node2
      itmp       = obj%mother(s1)
      obj%mother(s1) = obj%mother(s2)
      obj%mother(s2) = itmp
      itmp       = obj%daught(s1)
      obj%daught(s1) = obj%daught(s2)
      obj%daught(s2) = itmp
      itmp       = obj%cutdim(s1)
      obj%cutdim(s1) = obj%cutdim(s2)
      obj%cutdim(s2) = itmp
      itmp       = obj%channl(s1)
      obj%channl(s1) = obj%channl(s2)
      obj%channl(s2) = itmp
      obj%nodech(obj%channl(s1)) = sstr1
      if (obj%daught(s1).ne.0) then
        obj%mother(obj%daught(s1)  ) = sstr1
        obj%mother(obj%daught(s1)+1) = sstr1
      endif
      obj%nodech(obj%channl(s2)) = sstr2
      if (obj%daught(s2).ne.0) then
        obj%mother(obj%daught(s2)  ) = sstr2
        obj%mother(obj%daught(s2)+1) = sstr2
      endif
      end subroutine


      subroutine avh_parni_marg( obj ,nunit ,label )
!* ********************************************************************
!* * The marginal distribution for each dimension, that is the
!* * distribution obtained by integrating over all other dimensions,
!* * is written into files. The input "nunit" is the unit to write to,
!* * and "label" is a label for those files. To view, for
!* * example, the distribution for label=37 and dimension=2, use
!* *   gnuplot> plot 'PARNI37d2' w l
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      integer ,intent(in) :: nunit,label
      integer :: nch,idm,length,ndig,ich,nstp,node,nx,ix,jx
      real(kind(1d0)) :: x0(avh_parni_dim),x1(avh_parni_dim),xx(2*obj%nchann) ,&
                         height(2*obj%nchann),rr
      character(24) :: filename,chii
      character(24) ,parameter :: spaces = '                        '
!
      nch = obj%nchann
!
      call avh_parni_charnum(chii,ndig ,label)
      write(6,'(a28,a24)') 'Writing marginals for PARNI ',chii
!
      do idm=1,obj%ndimen
!
      filename = spaces
      filename(1:5) = 'PARNI'
      length = 5
      call avh_parni_charnum(chii,ndig ,label)
      filename(length+1:length+ndig) = chii(1:ndig)
        length = length+ndig
      filename(length+1:length+1) = 'd'
        length = length+1
      call avh_parni_charnum(chii,ndig ,idm)
      filename(length+1:length+ndig) = chii(1:ndig)
        length = length+ndig
      open(nunit,file=filename)
!
      do ich=1,nch
        node = obj%nodech(ich)
        call get_box(obj ,x0,x1,nstp ,node)
        xx(2*ich-1) = x0(idm)
        xx(2*ich  ) = x1(idm)
      enddo
      nx = 2*nch
      call sort(xx,nx)
      nx = nx+1
      ix = 0
      do while (ix.lt.nx)
        ix = ix+1
        if (xx(ix).eq.xx(ix+1)) then
          do jx=ix+1,nx-1
            xx(jx) = xx(jx+1)
          enddo
          ix = ix-1
          nx = nx-1
        endif
      enddo
      do ix=1,nx-1
        height(ix) = 0d0
        do ich=1,nch
          node = obj%nodech(ich)
          call get_box(obj ,x0,x1,nstp ,node)
          rr = min(xx(ix+1),x1(idm)) - max(xx(ix),x0(idm))
          if (rr.lt.0) rr = 0d0
          rr = rr/(x1(idm)-x0(idm))
          height(ix) = height(ix) + rr * obj%weight(obj%pointr(ich))
        enddo
        height(ix) = height(ix)/(xx(ix+1)-xx(ix))
      enddo
      do ix=1,nx-1
        write(nunit,101) xx(ix  ),'0d0'
        write(nunit,102) xx(ix  ),height(ix)
        write(nunit,102) xx(ix+1),height(ix)
        write(nunit,101) xx(ix+1),'0d0'
      enddo
  101 format(e16.8,1x,a3)
  102 format(2e16.8)
!
      close(nunit)
      call avh_parni_charnum(chii,ndig ,idm)
      write(6,*) 'To view marginal for dim ',chii(1:ndig) ,&
                 ', use gnuplot> plot "',filename(1:length),'" w l'
      enddo
      end subroutine

      subroutine sort( xx,nn )
!* ********************************************************************
!* ********************************************************************
      implicit none
      integer ,intent(in) :: nn
      real(kind(1d0)) ,intent(inout) :: xx(nn)
      integer :: ll,ii,jj,ir
      real(kind(1d0)) :: xxb 
      if (nn.le.1) goto 30
!     
      ll = nn/2 + 1
      ir = nn
  10  continue
        if (ll.gt.1) then 
          ll = ll - 1
          xxb = xx(ll)
        else
          xxb = xx(ir)
          xx(ir) = xx(1)
          ir = ir - 1
          if (ir.eq.1) then
            xx(1) = xxb
            goto 30
          endif
        endif
        ii = ll
        jj = ll + ll
  20    continue
        if (jj.le.ir) then
          if (jj.lt.ir) then
            if (xx(jj).lt.xx(jj+1)) jj = jj + 1
          endif
          if (xxb.lt.xx(jj)) then
            xx(ii) = xx(jj)
            ii = jj
            jj = jj + jj
          else
            jj = ir + 1
          endif
        goto 20
        endif
        xx(ii) = xxb
      goto 10
  30  continue
!
      end subroutine


      subroutine avh_parni_plot( obj ,nunit ,label )
!* *********************************************************************
!* * Only works for dimension=1 and dimension=2.
!* * Consider the example that label=37 . A plot of the built
!* * density is written in the file with the name PARNI37p. To view,
!* * open gnuplot and put 
!* *   gnuplot> plot 'PARNIp37' w l
!* * for the 1-dim case. For the the 2-dim case, the cell structure
!* * will be depicted this way. To get the full 3-dim plot for the
!* * 2-dim case, use
!* *   gnuplot> splot 'PARNIp37' w l
!* * Input: "nunit" is a number for the unit to write to.
!* *        "label" is a label to the filename
!* *********************************************************************
      implicit none
      type(avh_parni_type) ,intent(in) :: obj
      integer ,intent(in) :: nunit,label
      real(kind(1d0)) :: height,x0(avh_parni_dim),x1(avh_parni_dim)
      integer :: ich,node,nstp,ndig,length
      character(24) :: filename,chii
      character(24) ,parameter :: spaces = '                        '
!
      call avh_parni_charnum(chii,ndig ,label)
!      write(6,'(a31,a24)') 'Writing distribution for PARNI ',chii
!
      filename = spaces
      filename(1:5) = 'PARNI'
        length = 5
      call avh_parni_charnum(chii,ndig ,label)
      filename(length+1:length+ndig) = chii(1:ndig)
        length = length+ndig
      filename(length+1:length+1) = 'p'
        length = length+1
      open(nunit,file=filename)
!
      if (obj%ndimen.eq.1) then
!
  101 format(e16.8,1x,a3)
  102 format(e16.8,e16.8)
      do ich=1,obj%nchann
        node = obj%nodech(ich)
        call get_box(obj ,x0,x1,nstp ,node)
        height = obj%weight(obj%pointr(ich))/halfpwr(nstp)
        write(nunit,101) x0(1),'0d0' 
        write(nunit,102) x0(1),height
        write(nunit,102) x1(1),height
        write(nunit,101) x1(1),'0d0' 
      enddo
!      write(6,*) 'To view distribution, please use '
!     &          ,'gnuplot> plot "',filename(1:length),'" w l'
!
      elseif (obj%ndimen.eq.2) then
!
  201 format(3e16.8)
  202 format()
      do ich=1,obj%nchann
        node = obj%nodech(ich)
        call get_box(obj ,x0,x1,nstp ,node)
        height = 0d0
        write(nunit,201) x0(1),x0(2),height
        write(nunit,201) x0(1),x1(2),height
        write(nunit,201) x1(1),x1(2),height
        write(nunit,201) x1(1),x0(2),height
        write(nunit,201) x0(1),x0(2),height
        write(nunit,202)
        height = obj%weight(obj%pointr(ich))/halfpwr(nstp)
        write(nunit,201) x0(1),x0(2),height
        write(nunit,201) x0(1),x1(2),height
        write(nunit,201) x1(1),x1(2),height
        write(nunit,201) x1(1),x0(2),height
        write(nunit,201) x0(1),x0(2),height
        write(nunit,202)
        write(nunit,202)
      enddo
      write(6,*) 'To view cell-structure, please use ' ,&
                 'gnuplot> plot "',filename(1:length),'" w l'
      write(6,*) 'To view distribution, please use ' ,&
                 'gnuplot> splot "',filename(1:length),'" w l'
!
      endif
!
      close(nunit)
!
      end subroutine

  
      subroutine avh_parni_result( obj ,label )
!* ********************************************************************
!* ********************************************************************
      implicit none
      type(avh_parni_type) ,intent(inout) :: obj
      integer ,intent(in) :: label
      integer :: length,ndig ,ich,jch,node,nstp,dimnsn(size_ncut)
      logical :: upperh(size_ncut)
      character(24) :: filename,chii
      integer ,parameter :: nunit = 6
!
      filename = '                        '
      filename(1:6) = 'PARNIr'
      length = 6
      call avh_parni_charnum(chii,ndig ,label)
      filename(length+1:length+ndig) = chii(1:ndig)
      length = length+ndig
!
!      open(nunit,file=filename)
      write(nunit,'(a36,a36)') '####################################' ,&
                               '####################################'
      write(nunit,'(a1,a42,a1,i2,a1,a25)') &
       '#','Diagnostics for the optimization of PARNI',' ',label,':','#'
      write(nunit,'(a1,a7,a1,i3 ,a1,a15,a1,i7 ,a1,a20,a1,i6 ,a8)') &
       '#','task',':',obj%iitask ,&
       ',','batch-size',':',obj%nbatch ,&
       ',','number of channels',':',obj%nchann,'#'
      write(nunit,'(a36,a36)') '####################################' ,&
                               '####################################'
!      close(nunit)
      end subroutine
      

      subroutine avh_parni_charnum(chii,ndig ,ii)
!* ********************************************************************
!* * Determines the character representing "ii".
!* ********************************************************************
      implicit none
      integer ,intent(in) :: ii
      integer ,intent(out) :: ndig
      character(24) ,intent(out) :: chii
      integer :: istp,iden,jj,ipos,idig
      character(24) :: nmbr
      integer ,parameter :: ndigmax = 9
      character(1) :: arab(0:9)
      data arab/'0','1','2','3','4','5','6','7','8','9'/
!
      chii = '                        '
      istp = ii
      iden = 10**(ndigmax-1)
      if (ii.ge.10*iden) then
        write(6,*) 'WARNING from avh_parni_charnum: ii =',ii ,&
                   ' is too large, returning "_x".'
        chii(1:2) = '_x'
        ndig = 2
        return
      endif
      do jj=ndigmax,1,-1
        ipos = ndigmax-jj+1
        idig = istp/iden
        nmbr(ipos:ipos) = arab(idig)
        istp = istp - iden*idig
        iden = iden/10
      enddo
      ipos=1
      do while (nmbr(ipos:ipos).eq.arab(0))
        ipos = ipos+1
      enddo
      ndig = ndigmax-ipos+1
      chii(1:ndig) = nmbr(ipos:ndigmax)
      end subroutine

end module avh_parni
