module avh_print
!
contains
!
  function printdbl( ibefor,iafter,xx ) result(value)
  implicit none
  integer         ,intent(in) :: ibefor,iafter
  real(kind(1d0)) ,intent(in) :: xx
  character(6+abs(ibefor)+iafter) :: value
  character(10) :: aa,bb
  character(72) :: cc,dd
  real(kind(1d0)) :: yy
  integer :: ibef,idec
  ibef = max(0,ibefor)
  idec = abs(ibef)+iafter
  write(aa,'(i10)') 7+idec;  aa=adjustl(aa)
  write(bb,'(i10)') idec;    bb=adjustl(bb)
  aa = '(e'//trim(aa)//'.'//trim(bb)//')'
  write(cc,aa) xx;
  yy = xx/10**ibef
  write(dd,aa) yy;
  value(1:1) = cc(1:1)
  if     (ibef.gt.0) then
    value(2:ibef+1) = cc(4:ibef+3)
    value(ibef+2:ibef+2) = '.'
    value(ibef+3:idec+2) = cc(ibef+4:idec+3)
    value(idec+3:idec+6) = dd(idec+4:idec+7)
  elseif (ibef.eq.0) then
    value(2:2) = '.'
    value(3:idec+2) = cc(4:idec+3)
    value(idec+3:idec+6) = cc(idec+4:idec+7)
  endif
  if (value(idec+5:idec+6).eq.'00'.and.value(idec+3:idec+3).eq.'E') &
    value(idec+3:idec+6) = '    '
  value = adjustl(value)
  end function
!
  function printdble( ibefor,iafter,xx ) result(value)
  implicit none
  integer         ,intent(in) :: ibefor,iafter
  real(kind(1d0)) ,intent(in) :: xx
  character(6+abs(ibefor)+iafter) :: value
  character(10) :: aa,bb
  character(72) :: cc,dd
  real(kind(1d0)) :: yy
  integer :: ibef,idec
  ibef = max(0,ibefor)
  idec = abs(ibef)+iafter
  write(aa,'(i10)') 7+idec;  aa=adjustl(aa)
  write(bb,'(i10)') idec;    bb=adjustl(bb)
  aa = '(e'//trim(aa)//'.'//trim(bb)//')'
  write(cc,aa) xx;
  yy = xx/10**ibef
  write(dd,aa) yy;
  value(1:1) = cc(1:1)
  if     (ibef.gt.0) then
    value(2:ibef+1) = cc(4:ibef+3)
    value(ibef+2:ibef+2) = '.'
    value(ibef+3:idec+2) = cc(ibef+4:idec+3)
    value(idec+3:idec+6) = dd(idec+4:idec+7)
  elseif (ibef.eq.0) then
    value(2:2) = '.'
    value(3:idec+2) = cc(4:idec+3)
    value(idec+3:idec+6) = cc(idec+4:idec+7)
  endif
  value = adjustl(value)
  end function
!
  function printsdbl( ibefor,iafter,xx ) result(value)
  implicit none
  integer         ,intent(in) :: ibefor,iafter
  real(kind(1d0)) ,intent(in) :: xx
  character(7+abs(ibefor)+iafter) :: value
  value = printdbl( ibefor,iafter,xx )
  if (value(1:1).ne.'-') value = ' '//value
  end function
!
  function printsdble( ibefor,iafter,xx ) result(value)
  implicit none
  integer         ,intent(in) :: ibefor,iafter
  real(kind(1d0)) ,intent(in) :: xx
  character(7+abs(ibefor)+iafter) :: value
  value = printdble( ibefor,iafter,xx )
  if (value(1:1).ne.'-') value = ' '//value
  end function
!
  function printint( ii ) result(value)
  implicit none
  integer ,intent(in) :: ii
  character(11) :: value
  write(value,'(i11)') ii
  value = adjustl(value)
  end function
!
  function printsint( ii ) result(value)
  implicit none
  integer ,intent(in) :: ii
  character(11) :: value
  write(value,'(i11)') ii
  value = adjustl(value)
  if (value(1:1).ne.'-') value = ' '//value
  end function
!
end module


!program test
!use avh_print
!implicit none
!real(kind(1d0)) xx
!xx = 0.1234567d-8
!write(*,*) 'Format:',printdbl(11,4, xx)
!write(*,*) 'Format:',printdbl(11,4,-xx)
!write(*,*) 'Format:',printsdbl(11,4, xx)
!write(*,*) 'Format:',printsdbl(11,4,-xx)
!xx = 3.1234567d-0
!write(*,*) 'Format:',printdbl(0,4, xx)
!write(*,*) 'Format:',printdbl(1,4, xx)
!write(*,*) 'Format:',printdble(1,4, xx)
!write(*,*) 'Format:',printsdble(1,4, xx)
!end program
