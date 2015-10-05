module statistics
implicit none
!
  integer :: numconts
  character(len=2) , dimension(:) , allocatable :: conts_str
!
  real(kind(1d0)) , dimension(:) , allocatable :: totevents
  real(kind(1d0)) , dimension(:,:) , allocatable :: acceptedevents
  real(kind(1d0)) , dimension(:,:) , allocatable :: sumweight
  real(kind(1d0)) , dimension(:,:) , allocatable :: sumweightsqr
  real(kind(1d0)) , dimension(:,:) , allocatable :: ave,sig,rel
!
contains
!
subroutine init_stat
use flags
use scales
implicit none
!
!
  integer :: istat,icont
!
! We allocate the arrays which will hold the statistics for
! all the contributions.
! We have to investigate how many contributions we have:
  numconts = 0
  if (flg_LO) numconts = numconts + 1
  if (flg_NLO_V) numconts = numconts + 1
  if (flg_NLO_R) numconts = numconts + 1
  if (flg_NNLO_VV) numconts = numconts + 1
  if (flg_NNLO_RV) numconts = numconts + 1
  if (flg_NNLO_RR) numconts = numconts + 1
!
  allocate(conts_str(numconts), &
           totevents(numconts), &
           acceptedevents(nscales,numconts), &
           sumweight(nscales,numconts), &
           sumweightsqr(nscales,numconts), &
           ave(nscales,numconts), &
           sig(nscales,numconts), &
           rel(nscales,numconts), &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating stat arrays..."
    stop
  end if
  icont = 0
  if (flg_LO) then
    icont = icont + 1
    conts_str(icont) = 'B '
  end if
  if (flg_NLO_V) then
    icont = icont + 1
    conts_str(icont) = 'V '
  end if
  if (flg_NLO_R) then
    icont = icont + 1
    conts_str(icont) = 'R '
  end if
  if (flg_NNLO_VV) then
    icont = icont + 1
    conts_str(icont) = 'VV'
  end if
  if (flg_NNLO_RV) then
    icont = icont + 1
    conts_str(icont) = 'RV'
  end if
  if (flg_NNLO_RR) then
    icont = icont + 1
    conts_str(icont) = 'RR'
  end if
!
  print *,"Allocations are made for the following contributions: "
  do icont=1,numconts
    write(*,'(A,1x)',advance='no') conts_str(icont)
  end do
  write(*,*)
!
  totevents      = 0d0
  acceptedevents = 0d0
  sumweight      = 0d0
  sumweightsqr   = 0d0
!
end subroutine init_stat
!
subroutine accu_stat(cont,weight)
use scales
implicit none
!
  character (len=2) , intent(in) :: cont
  real(kind(1d0)) , dimension(:) , intent(in) :: weight
!
  integer :: icont,iscale
!
!
  do icont=1,numconts
    if (conts_str(icont).eq.cont) then
      totevents(icont) = totevents(icont) + 1d0
      do iscale=1,nscales
        if (weight(iscale).ne.0d0) acceptedevents(iscale,icont) = &
                                   acceptedevents(iscale,icont) + 1d0
      end do
      sumweight(:,icont) = sumweight(:,icont) + weight(:)
      sumweightsqr(:,icont) = sumweightsqr(:,icont) + weight(:)**2
    end if
  end do
!
end subroutine accu_stat
!
subroutine null_stat
implicit none
!
!
  print *,"Nullifying the stat arrays..."
!
! We nullify all statistics:
  totevents      = 0d0
  acceptedevents = 0d0
  sumweight      = 0d0
  sumweightsqr   = 0d0
!
end subroutine null_stat
!
subroutine null_stati(cont)
implicit none
!
  character (len=5) , intent(in) :: cont
!
  integer :: icont
!
  do icont=1,numconts
    if ((cont.eq.'Bkin ').and. &
        ((conts_str(icont).eq.'B ').or. &
         (conts_str(icont).eq.'V ').or. &
         (conts_str(icont).eq.'VV'))) then
      totevents(icont)        = 0d0
      acceptedevents(:,icont) = 0d0
      sumweight(:,icont)      = 0d0
      sumweightsqr(:,icont)   = 0d0
    elseif ((cont.eq.'Rkin ').and. &
            ((conts_str(icont).eq.'R ').or. &
            (conts_str(icont).eq.'RV'))) then
      totevents(icont)        = 0d0
      acceptedevents(:,icont) = 0d0
      sumweight(:,icont)      = 0d0
      sumweightsqr(:,icont)   = 0d0
    elseif ((cont.eq.'RRkin').and. &
            (conts_str(icont).eq.'RR')) then
      totevents(icont)        = 0d0
      acceptedevents(:,icont) = 0d0
      sumweight(:,icont)      = 0d0
      sumweightsqr(:,icont)   = 0d0
    end if
  end do
!
end subroutine null_stati
!
subroutine finalize_stat(iun)
use flags
use scales
use misc
implicit none
!
  integer , intent(in) :: iun
!
  integer :: icont,iscale
  real(kind(1d0)) :: tot,toterr,totrel, &
                     totLO,totLOerr,totLOrel, &
                     totNLO,totNLOerr,totNLOrel, &
                     totNNLO,totNNLOerr,totNNLOrel
 
!
!
  do icont=1,numconts
    if (totevents(icont).eq.0) cycle
    ave(:,icont) = sumweight(:,icont)/totevents(icont)
    sig(:,icont) = sqrt(abs(sumweightsqr(:,icont)/totevents(icont) &
               - ave(:,icont)**2)/(totevents(icont) - 1d0))
    rel(:,icont) = sig(:,icont)/abs(ave(:,icont))
  end do
!
! We print out the results: 
  do icont=1,numconts
    if (totevents(icont).eq.0) cycle
    if (nscales.eq.1) then
!      write(iun,'(a,1x,e16.8,1x,e12.4,1x,f8.4)') conts_str(icont)//"  ", &
!        ave(1,icont),sig(1,icont),rel(1,icont)
      write(iun,'(a,1x,g0,1x,g0,1x,f8.4)') conts_str(icont)//"  ", &
        ave(1,icont),sig(1,icont),rel(1,icont)
    else
      write(iun,'(a,I0,a)') "=================== tot events: ", &
                       int(totevents(icont))," ========================"
      do iscale=1,nscales
        write(iun,'(2(a,1x,F5.3,1x),a,a,1x,e16.8,1x,e12.4,1x,f8.4)') &
          "xiR:",xirf(1,iscale),"xiF:",xirf(2,iscale),", ", &
          conts_str(icont)//"  ",ave(iscale,icont),sig(iscale,icont), &
          rel(iscale,icont)
      end do
    end if
  end do
! Final stat in each order:
! We can have several different scales:
  do iscale=1,nscales
    tot = 0d0 ; toterr = 0d0 ; totrel = 0d0
    totLO = 0d0 ; totLOerr = 0d0 ; totLOrel = 0d0
    totNLO = 0d0 ; totNLOerr = 0d0 ; totNLOrel = 0d0
    totNNLO = 0d0 ; totNNLOerr = 0d0 ; totNNLOrel = 0d0
! We print out \xi_R and \xi_F if we have more than one scales:
    if (nscales.ne.1) then
      write(iun,'(10x,a,1x,F5.3,4x,a,4x,a,1x,F5.3)') &
        "xiR:",xirf(1,iscale)," , ","xiF:",xirf(2,iscale)
    end if
! Loop over all contributions:
    do icont=1,numconts
! If the calculation of the given contribution is not started yet
! we skip it:
      if (totevents(icont).eq.0) cycle
      tot = tot + ave(iscale,icont)
      toterr = toterr + sig(iscale,icont)**2
      if (flg_LO.and.conts_str(icont).eq.'B ') then
        totLO = totLO + ave(iscale,icont)
        totLOerr = totLOerr + sig(iscale,icont)**2
      end if
      if (flg_NLO.and. &
          ((conts_str(icont).eq.'V ').or. &
           (conts_str(icont).eq.'R '))) then
        totNLO = totNLO + ave(iscale,icont)
        totNLOerr = totNLOerr + sig(iscale,icont)**2
      end if
      if (flg_NNLO.and. &
          ((conts_str(icont).eq.'VV').or. &
           (conts_str(icont).eq.'RV').or. &
           (conts_str(icont).eq.'RR'))) then
        totNNLO = totNNLO + ave(iscale,icont)
        totNNLOerr = totNNLOerr + sig(iscale,icont)**2
      end if
    end do
    if ((flg_LO.and.totLO.ne.0d0).or. &
        (flg_NLO.and.totNLO.ne.0d0).or. &
        (flg_NNLO.and.totNNLO.ne.0d0)) then
!      write(iun,"(a)") "++++++++++++++++++++++++++++++++++++++++++++"// &
!                     "++++++++++++++++++++++++++++++++++++"
    end if
    toterr = sqrt(toterr)
    totrel = toterr / tot
!    write(iun,'(a,1x,e16.8,1x,e12.4,1x,f8.4)') 'tot ', &
!      tot,toterr,totrel
    write(iun,'(a,1x,g0,1x,g0,1x,f8.4)') 'tot ', &
      tot,toterr,totrel
    if (flg_LO.and.totLO.ne.0d0) then
      totLOerr = sqrt(totLOerr)
      totLOrel = totLOerr / totLO
!      write(iun,'(a,1x,e16.8,1x,e12.4,1x,f8.4)') 'LO  ', &
!        totLO,totLOerr,totLOrel
      write(iun,'(a,1x,g0,1x,g0,1x,f8.4)') 'LO  ', &
        totLO,totLOerr,totLOrel
    end if
    if (flg_NLO.and.totNLO.ne.0d0) then
      totNLOerr = sqrt(totNLOerr)
      totNLOrel = totNLOerr / totNLO
!      write(iun,'(a,1x,e16.8,1x,e12.4,1x,f8.4)') 'NLO ', &
!        totNLO,totNLOerr,totNLOrel
      write(iun,'(a,1x,g0,1x,g0,1x,f8.4)') 'NLO ', &
        totNLO,totNLOerr,totNLOrel
    end if
    if (flg_NNLO.and.totNNLO.ne.0d0) then
      totNNLOerr = sqrt(totNNLOerr)
      totNNLOrel = totNNLOerr / totNNLO
!      write(iun,'(a,1x,e16.8,1x,e12.4,1x,f8.4)') 'NNLO', &
!        totNNLO,totNNLOerr,totNNLOrel
      write(iun,'(a,1x,g0,1x,g0,1x,f8.4)') 'NNLO', &
        totNNLO,totNNLOerr,totNNLOrel
    end if
! Only print if the contribution is present and we already have
! numbers for it:
    if ((flg_LO.and.totLO.ne.0d0).or. &
        (flg_NLO.and.totNLO.ne.0d0).or. &
        (flg_NNLO.and.totNNLO.ne.0d0)) then
! This calculation does make sense in the extreme when
! one contribution is turned off provided by initialization to
! zero:
      tot = totLO + totNLO + totNNLO
      toterr = sqrt(totLOerr**2 + totNLOerr**2 + totNNLOerr**2)
      totrel = toterr / tot
!      write(iun,'(a,1x,e16.8,1x,e12.4,1x,f8.4)') 'tot:', &
!        tot,toterr,totrel
      write(iun,'(a,1x,g0,1x,g0,1x,f8.4)') 'tot:', &
        tot,toterr,totrel
    end if
    if ((flg_LO.and.totLO.ne.0d0).or. &
        (flg_NLO.and.totNLO.ne.0d0).or. &
        (flg_NNLO.and.totNNLO.ne.0d0)) then
!      write(iun,"(a)") "++++++++++++++++++++++++++++++++++++++++++++"// &
!                     "++++++++++++++++++++++++++++++++++++"
    end if
  end do
!
end subroutine finalize_stat
!
! We write out the statistics to a simple text file for later
! purposes:
subroutine output_stat
use process
use random
implicit none
!
!
  integer :: filestat
  integer , parameter :: iun = 21
  character (len=72) :: fname
!
!
  fname = trim(process_name)//'-stat'
! If multiseed mode is active an additional number is needed to 
! distinguish between different runs:
  if (flg_manyseeds) then
    write(fname,'(a,a,I0.4)') trim(fname),'-',iseed_manyseeds
  end if
! File extension:
  fname = trim(fname)//'.dat'
!
  open(file=fname,unit=iun,status='unknown',iostat=filestat)
  if (filestat.ne.0) then
    print *,"Problem ocurred during creation of the stat output..."
    stop
  end if
!
  call finalize_stat(iun)
  close(iun)
!
end subroutine output_stat
!
end module statistics
!
module histo
implicit none
!
  integer :: maxhist = 50
  integer :: numhist = 0
  integer :: nweight
  integer , parameter :: maxbin    = 100
!
  integer , dimension(:) , allocatable :: numbins
  integer , dimension(:,:) , allocatable :: nhits
  integer , dimension(:) , allocatable :: nevnts
  real(kind(1d0)) , dimension(:,:) , allocatable :: xhist
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: & 
    yhist,yhist_out,yhist_tmp
  real(kind(1d0)) , dimension(:,:,:) , allocatable :: & 
    errhist,errhist_out
  character (len=100) , dimension(:) , allocatable :: histtitle
!
interface bookup_hist
  module procedure bookup_eqv_hist
  module procedure bookup_var_hist
end interface bookup_hist
!
contains
!
subroutine init_hist(nhist)
use scales
implicit none
!
  integer , optional , intent(in) :: nhist
!
  integer :: istat
!
! We have to set nweight to coincide with nscales:
  nweight = nscales
! If nhist is present we use it as the maximal number of histos:
  if (present(nhist)) then
    maxhist = nhist
    print *,"The maximum number of bookable histos is changed to: ", &
            maxhist
  end if
!
! We allocate all the arrays needed to store the histos:
  allocate(numbins(maxhist),                        &
           nevnts(maxhist),                         &
           nhits(0:maxbin+1,maxhist),               &
           xhist(maxbin+1,maxhist),                 &
           yhist(nweight,0:maxbin+1,maxhist),       &
           yhist_out(nweight,0:maxbin+1,maxhist),   &
           yhist_tmp(nweight,0:maxbin+1,maxhist),   &
           errhist(nweight,0:maxbin+1,maxhist),     &
           errhist_out(nweight,0:maxbin+1,maxhist), &
           histtitle(maxhist),                      &
           stat=istat)
  if (istat.ne.0) then
    print *,"Problem occured during allocating arrays for histos..."
    stop
  end if
!
! We nullify everything:
  yhist = 0d0 ; yhist_out = 0d0 ; yhist_tmp = 0d0 ; nhits = 0
  errhist = 0d0; errhist_out = 0d0 ; nevnts = 0d0
!
end subroutine init_hist
!
! This routine manages the bookup of the histos:
! We have two possibilities: equal binwidth or variable:
! To have equal bins the user should specify the 
! binwidth (dx), the lower (xmin) and upper (xmax) y limits,
! To have variable bin widths the user should specify:
! the number of bins (nbins) and the boundary for all the bins
! (bins):
subroutine bookup_eqv_hist(title,dx,xmin,xmax)
implicit none
!
  character (len=*) , intent(in) :: title
  real(kind(1d0)) , intent(in) :: dx,xmin,xmax
!
  integer :: ihist,ibin
  real(kind(1d0)) :: xpos
  real(kind(1d0)) , parameter :: eps = 1d-8
!
! We fill in the title of the histo to histtitle, but we also check 
! whether a histo with the same name is already booked or not:
  do ihist=1,numhist
    if (trim(histtitle(ihist)).eq.trim(title)) then
      write(*,'(a,a,a)') "A histo with name ",trim(title), &
                         " is already booked"
      stop
    end if
  end do
! Increase the number of histos:
  numhist = numhist + 1
! Fill in the title:
  histtitle(numhist) = trim(title)
! Check whether we have more than bookable:
  if (numhist.gt.maxhist) then
    print *,"The number of histos exceeds maxhist, change it..."
    print *,"maxhist: ",maxhist
    stop
  end if
!  print *,"Booking histo with equal bins..."
! We fill in the x values:
  ibin = 0
  xpos = xmin
  do while(.true.)
    ibin = ibin + 1
    xhist(ibin,numhist) = xpos
! We increment the position:
    xpos = xpos + dx
    if (xpos.gt.(xmax+eps)) exit
  end do
  numbins(numhist) = ibin - 1
!
end subroutine bookup_eqv_hist
!
! This routine manages the bookup of the histos:
! To have variable bin widths the user should specify:
! the number of bins (nbins) and the boundary for all the bins
! (bins):
subroutine bookup_var_hist(title,nbins,bins)
implicit none
!
  character (len=*) , intent(in) :: title
  integer , intent(in) :: nbins
  real(kind(1d0)) , dimension(:) , intent(in) :: bins
!
  integer :: ihist,ibin
  real(kind(1d0)) :: xpos
!
! We fill in the title of the histo to histtitle, but we also check 
! whether a histo with the same name is already booked or not:
  do ihist=1,numhist
    if (trim(histtitle(ihist)).eq.trim(title)) then
      write(*,'(a,a,a)') "A histo with name ",trim(title), &
                         " is already booked"
      stop
    end if
  end do
! Increase the number of histos:
  numhist = numhist + 1
! Fill in the title:
  histtitle(numhist) = trim(title)
! Check whether we have more than bookable:
  if (numhist.gt.maxhist) then
    print *,"The number of histos exceeds maxhist, change it..."
    print *,"maxhist: ",maxhist
    stop
  end if
!  print *,"Booking variable bins..."
! We don not have to count down the number of bins it is in nbins:
  numbins(numhist) = nbins
! The bounds of each bin is in bins:
  xhist(1:nbins+1,numhist) = bins(1:nbins+1)
!
end subroutine bookup_var_hist
!
subroutine fill_hist(title,xpos,weights)
implicit none
!
  character (len=*), intent(in) :: title
  real(kind(1d0)) , intent(in) :: xpos
  real(kind(1d0)) , dimension(:) , intent(in) :: weights
!
  integer :: ihist,ibin
!
! We look for the histo in quetion:
  ihist = find_hist(title)
!
! We fill in to the temporary histo only:
! we can have an under or overflow:
! underflow:
  if (xpos.lt.xhist(1,ihist)) then
    nhits(0,ihist) = nhits(0,ihist) + 1
    yhist_tmp(:,0,ihist) = yhist_tmp(:,0,ihist) + weights
! overflow:
  elseif (xpos.gt.xhist(numbins(ihist)+1,ihist)) then
    nhits(numbins(ihist)+1,ihist) = nhits(numbins(ihist)+1,ihist) + 1
    yhist_tmp(:,numbins(ihist)+1,ihist) = &
      yhist_tmp(:,numbins(ihist)+1,ihist) + weights
! fill in an ordinary bin:
  else
! We have to find the bin:
    do ibin=1,numbins(ihist)
      if ((xhist(ibin,ihist).lt.xpos).and. &
          (xhist(ibin+1,ihist).gt.xpos)) then
        nhits(ibin,ihist) = nhits(ibin,ihist) + 1
        yhist_tmp(:,ibin,ihist) = &
          yhist_tmp(:,ibin,ihist) &
          + weights/(xhist(ibin+1,ihist) - xhist(ibin,ihist))
        exit
      end if
    end do
  end if
!
end subroutine fill_hist
!
subroutine accu_hist
implicit none
!
!
  integer :: ihist,ibin
!
! We have to go through all the histos and copy everything:
  do ihist=1,numhist
    do ibin=0,numbins(ihist)+1
      yhist(:,ibin,ihist) = yhist(:,ibin,ihist) &
                          + yhist_tmp(:,ibin,ihist)
      errhist(:,ibin,ihist) = errhist(:,ibin,ihist) &
                            + yhist_tmp(:,ibin,ihist)**2
      yhist_tmp(:,ibin,ihist) = 0d0
    end do
    nevnts(ihist) = nevnts(ihist) + 1
  end do
!
end subroutine accu_hist
!
subroutine store_hist
use process
implicit none
!
!
  integer :: ihist,ibin,iwgt
  integer :: n
  real(kind(1d0)) , dimension(:) , allocatable :: wgt
!
! We take yhist and errhist and process their content and
! add it to yhist_out and errhist_out this is needed because
! we can have several different contributions with different
! statistics:
  allocate(wgt(nweight))
!
  do ihist=1,numhist
    n = nevnts(ihist)
    do ibin=1,numbins(ihist)+1
      wgt = yhist(:,ibin,ihist)
      yhist_out(:,ibin,ihist) = yhist_out(:,ibin,ihist) + wgt/n
! The square root is not taken here, because more contributions
! can further added, the square root is taken only during 
! writting out the result:
      errhist_out(:,ibin,ihist) = errhist_out(:,ibin,ihist) &
                                + (errhist(:,ibin,ihist)/n &
                                - wgt**2/n/n)/(n - 1)
    end do
  end do
!
! Deallocating:
  deallocate(wgt)
!
! If we stored everything we have to nullify the arrays:
  nevnts = 0 ; yhist = 0d0 ; errhist = 0d0
!
end subroutine store_hist
!
subroutine output_hist
use process
use flags
use random
implicit none
!
!
  integer :: ihist,ibin,iwgt
  integer :: filestat
  integer , parameter :: iun = 99
  character (len=72) :: fname
!
! We dump the histos into files:
  do iwgt=1,nweight
    if (.not.flg_manyseeds) then
      fname = trim(process_name)//'-output-W'
    else
      write(fname,'(a,I0.4,a)') trim(process_name)//'-output-', &
        iseed_manyseeds,'-W'
    end if
    write(fname,'(a,I0,a)') trim(fname),iwgt,'.dat'
! We open the file, if it exists we simply overwrite:
    open(file=fname,unit=iun,status='unknown',iostat=filestat)
    if (filestat.ne.0) then
      print *,"Problem during creation of ",trim(fname)
      stop
    end if
! We dump the histos in:
    do ihist=1,numhist
      write(iun,'(a,I3)') '# '//trim(histtitle(ihist))//' index ',ihist - 1 
      do ibin=1,numbins(ihist)
! Note that square root has to be taken for the uncertainty:
        write(iun,'(4(1x,e15.8))') xhist(ibin,ihist), &
                                   xhist(ibin+1,ihist), &
                                   yhist_out(iwgt,ibin,ihist), &
                                   sqrt(errhist_out(iwgt,ibin,ihist))
      end do
      write(iun,*)
      write(iun,*)
    end do
    close(iun)
  end do
!
end subroutine output_hist
!
subroutine null_hist
!
!
!
  print *,"Nullifying the histos..."
  yhist = 0d0 ; yhist_out = 0d0 ; yhist_tmp = 0d0 ; nhits = 0
  errhist = 0d0; errhist_out = 0d0 ; nevnts = 0d0
!
!
end subroutine null_hist
!
function find_hist(title) result(ihist)
implicit none
!
  integer :: ihist
!
  character (len=*) , intent(in) :: title
!
! We run through all the booked histos and try to find the
! needed one:
  do ihist=1,numhist
    if (trim(histtitle(ihist)).eq.trim(title)) then
      return
    end if
  end do
  ! If we reach this point the histo was not booked...
  write(*,'(a,a,a)') "The histo ",title," was not booked..."
  stop
!
end function find_hist
!
subroutine print_hist
implicit none
!
!
  integer :: ibin,ihist,iweight
!
  write(*,'(a,1x,I3,1x,a)') "We booked",numhist,"histos"
  do ihist=1,numhist
    write(*,'(a,1x,I3,a,a)') "The",ihist,"th histo is: ", &
                             histtitle(ihist)
    write(*,'(a,1x,I3,1x,a)') "The histo holds",numbins(ihist),"bins"
    write(*,*) "The x positions are: "
    do ibin=1,numbins(ihist) + 1
      write(*,advance='no',fmt='(F8.4,1x)') xhist(ibin,ihist)
    end do
    write(*,*)
  end do
!
end subroutine print_hist
!
end module histo
