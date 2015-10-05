module avh_kaleu_model
  use avh_print
  private ! Everything is private except the following list
  public :: model_type ,&
            vertx_type ,&
            addparticle ,&
            setexponent,&
            makeparton,&
            addvertex ,&
            model_gnrt_above ,&
            model_wght_above ,&
            model_gnrt_below ,&
            model_wght_below ,&
            nparticle ,&
            sparticle ,&
            mparticle ,&
            wparticle ,&
            haswidth ,&
            nowidth ,&
            isparton ,&
            fillfusion
!
! Maximum number of particles in a model
  integer ,parameter :: size_npa = 17
! Maximum number of vertices in a model
  integer ,parameter :: size_nvx = 100
! Unit to which messages are send
  integer ,parameter :: nunit = 6
! Default exponent for widthless propagators
  real(kind(1d0)) ,parameter :: expdef = 0.05d0
! Parameter for the distribution for low values of invariants
  real(kind(1d0)) ,parameter :: frac = 0.1d0
!
! The particle_type is private
  type :: particle_type
    real(kind(1d0)) :: mass=0d0 ,width=0d0 ,expon=0d0
    logical         :: haswidth=.false. ,nowidth =.false.
    character(12)   :: symb='0           '
    logical         :: isparton=.false.
  end type
!
! The model_type and vertx_type are public with private components,
! ie, the user is only allowed to declare them and pass them to
! subroutines from this module
  type :: model_type
    private
    type(particle_type) ,dimension(0:size_npa) :: pa ! particle
    integer :: npa = 0 ! number of particles
  end type
!
  type :: vertx_type
    private
    integer ,dimension(3,size_nvx) :: vx ! vertex
    integer :: nvx = 0 ! number of vertices
    logical :: printed=.false.
  end type
!
!
contains
!
!
  subroutine addparticle( model ,label ,symb,mass,width )
!**********************************************************************
!**********************************************************************
  implicit none
  type(model_type) ,intent(inout) :: model
  integer          ,intent(in)    :: label
  character(*)     ,intent(in)    :: symb
  real(kind(1d0))  ,intent(in)    :: mass,width
  integer :: nn
!
  model%npa = model%npa+1
  nn = model%npa
  if (nn.gt.size_npa) then
    write(*,*) 'ERROR in avh_kaleu_model: ' ,&
               'increase the parameter size_npa to at least',nn
    stop
  endif
!
  if (label.ne.nn) then
    write(*,*) 'ERROR in avh_kaleu_model: ' ,&
               'the,',nn,'-th call of addparticle has label =',label ,&
               ' as input, which however should also be',nn
    write(*,*) 'This rule may be annoying, but is for your own good.'
    stop
  endif
!
  model%pa(nn)%symb  = adjustl(symb)
  model%pa(nn)%mass  = mass
  model%pa(nn)%width = width
!
  if     (mass.gt.0d0.and.width.gt.0d0) then
    model%pa(nn)%haswidth = .true.
  elseif (mass.ge.0d0) then
    model%pa(nn)%nowidth = .true.
    model%pa(nn)%expon = expdef
  else
    write(*,*) 'ERROR in avh_kaleu_model: particle',label ,&
               ' has mass =',mass,' < 0'
    stop
  endif
  if (nunit.gt.0) write(nunit,*) 'MESSAGE from avh_kaleu_model:' &
                 ,' particle=',trim(model%pa(nn)%symb) &
                 ,', mass=',trim(printdbl(0,7,model%pa(nn)%mass)) &
                 ,' width=',trim(printdbl(0,7,model%pa(nn)%width))
  end subroutine
!
  subroutine setexponent( model ,label,expon )
!**********************************************************************
!**********************************************************************
  implicit none
  type(model_type) ,intent(inout) :: model
  integer          ,intent(in)    :: label
  real(kind(1d0))  ,intent(in)    :: expon
  if (label.lt.1.or.label.gt.model%npa) then
    write(*,*) 'ERROR in avh_kaleu_model_setexponent: label =',label
    stop
  endif
  if (expon.ge.1d0) then
    write(*,*) 'ERROR in avh_kaleu_model_setexponent: '&
              ,'invalid exponent(=',expon,') for label =',label
    stop
  endif
  if (nunit.gt.0) write(nunit,*) 'MESSAGE from avh_kaleu_model:' &
                    ,' setting expon=',trim(printdbl(0,6,expon)) &
                    ,' for particle=',model%pa(label)%symb
  model%pa(label)%expon = 1d0-expon
  end subroutine
!
  subroutine makeparton( model ,label )
!**********************************************************************
!**********************************************************************
  implicit none
  type(model_type) ,intent(inout) :: model
  integer          ,intent(in)    :: label
  if (label.lt.1.or.label.gt.model%npa) then
    write(*,*) 'ERROR in avh_kaleu_model makeparton: label =',label
    stop
  endif
  model%pa(label)%isparton = .true.
  end subroutine
!
  subroutine addvertex( vertx ,i1,i2,i3 )
!**********************************************************************
!**********************************************************************
  implicit none
  type(vertx_type) ,intent(inout) :: vertx
  integer          ,intent(in)    :: i1,i2,i3
!    
  vertx%nvx = vertx%nvx+1
  if (vertx%nvx.gt.size_nvx) then
    write(*,*) 'ERROR in avh_kaleu_model: ' ,&
               'increase the parameter size_nvx to at least',vertx%nvx
    stop
  endif
  vertx%vx(1,vertx%nvx) = i1
  vertx%vx(2,vertx%nvx) = i2
  vertx%vx(3,vertx%nvx) = i3
  end subroutine
!
!
!********************************************************************
  include 'kaleu_model.h'
!********************************************************************
!
!
  function nparticle( model ) result(value)
!********************************************************************
!* The number of particles in model
!********************************************************************
  implicit none
  integer :: value
  type(model_type) ,intent(in) :: model
  value = model%npa
  end function
!
  function sparticle( model ,ii ) result(value)
!********************************************************************
!* The symbol of particle number ii in model
!********************************************************************
  implicit none
  character(12) :: value
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  value = model%pa(ii)%symb
  end function
!
  function mparticle( model ,ii ) result(value)
!********************************************************************
!* The mass of particle number ii in model
!********************************************************************
  implicit none
  real(kind(1d0)) :: value
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  value = model%pa(ii)%mass
  end function
!
  function wparticle( model ,ii ) result(value)
!********************************************************************
!* The width of particle number ii in model
!********************************************************************
  implicit none
  real(kind(1d0)) :: value
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  value = model%pa(ii)%width
  end function
!
  function nowidth( model ,ii ) result(value)
!********************************************************************
!* The symbol of particle number ii in model
!********************************************************************
  implicit none
  logical :: value
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  value = model%pa(ii)%nowidth
  end function
!
  function haswidth( model ,ii ) result(value)
!********************************************************************
!* The symbol of particle number ii in model
!********************************************************************
  implicit none
  logical :: value
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  value = model%pa(ii)%haswidth
  end function
!
  function isparton( model ,ii ) result(value)
!********************************************************************
!* The symbol of particle number ii in model
!********************************************************************
  implicit none
  logical :: value
  type(model_type) ,intent(in) :: model
  integer          ,intent(in) :: ii
  value = model%pa(ii)%isparton
  end function
!
!
!
  subroutine fillfusion( model,vertx ,fusion,nsize )
!********************************************************************
!********************************************************************
  implicit none
  type(model_type) ,intent(in)    :: model
  type(vertx_type) ,intent(inout) :: vertx
  integer          ,intent(in)    :: nsize
  integer          ,intent(out)   :: fusion(1:nsize,1:nsize,0:nsize)
  integer :: i1,i2,i3,j1,j2,j3,ii,jj,kk,nn
!
! Remove equal vertices
  ii = vertx%nvx
  do while (ii.gt.0)
    i1 = vertx%vx(1,ii)
    i2 = vertx%vx(2,ii)
    i3 = vertx%vx(3,ii)
    jj = ii-1
    do while (jj.gt.0)
      j1 = vertx%vx(1,jj)
      j2 = vertx%vx(2,jj)
      j3 = vertx%vx(3,jj)
      if (   i1.eq.j1.and.i2.eq.j2.and.i3.eq.j3 &
         .or.i1.eq.j2.and.i2.eq.j3.and.i3.eq.j1 &
         .or.i1.eq.j3.and.i2.eq.j1.and.i3.eq.j2 &
         .or.i1.eq.j3.and.i2.eq.j2.and.i3.eq.j1 &
         .or.i1.eq.j2.and.i2.eq.j1.and.i3.eq.j3 &
         .or.i1.eq.j1.and.i2.eq.j3.and.i3.eq.j2 ) then
        if (nunit.gt.0) write(nunit,*) 'WARNING from avh_kaleu_model:' &
          ,' removing multiple-defined (',trim(model%pa(i1)%symb) &
                                     ,',',trim(model%pa(i2)%symb) &
                                     ,',',trim(model%pa(i3)%symb),')-vertices'
        vertx%nvx = vertx%nvx-1
        do kk=ii,vertx%nvx
          vertx%vx(1,kk) = vertx%vx(1,kk+1)
          vertx%vx(2,kk) = vertx%vx(2,kk+1)
          vertx%vx(3,kk) = vertx%vx(3,kk+1)
        enddo
        jj = 1
      endif
      jj = jj-1
    enddo
    ii = ii-1
  enddo
!
  if (nunit.gt.0.and..not.vertx%printed) then
    vertx%printed = .true.
    write(nunit,*) 'MESSAGE from avh_kaleu_model:' &
                  ,' here follows the list of vertices'
    do ii=1,vertx%nvx/2
      jj = ii+vertx%nvx/2
      i1=vertx%vx(1,ii); i2=vertx%vx(2,ii); i3=vertx%vx(3,ii)
      j1=vertx%vx(1,jj); j2=vertx%vx(2,jj); j3=vertx%vx(3,jj)
      write(nunit,*) 'MESSAGE from avh_kaleu_model:   ' &
                    ,'(',trim(model%pa(i1)%symb) &
                    ,',',trim(model%pa(i2)%symb) &
                    ,',',trim(model%pa(i3)%symb),') ' &
                    ,'(',trim(model%pa(j1)%symb) &
                    ,',',trim(model%pa(j2)%symb) &
                    ,',',trim(model%pa(j3)%symb),')'
    enddo
    ii = vertx%nvx
    if (ii.ne.ii/2*2) then
      i1=vertx%vx(1,ii); i2=vertx%vx(2,ii); i3=vertx%vx(3,ii)
      write(nunit,*) 'MESSAGE from avh_kaleu_model:   ' &
                    ,'       ' &
                    ,'(',trim(model%pa(i1)%symb) &
                    ,',',trim(model%pa(i2)%symb) &
                    ,',',trim(model%pa(i3)%symb),')'
    endif
  endif
!
! Fill fusion
  do ii=1,vertx%nvx
    i1 = vertx%vx(1,ii)
    i2 = vertx%vx(2,ii)
    i3 = vertx%vx(3,ii)
    if     (i1.eq.i2.and.i2.eq.i3) then
      nn = fusion( i1,i2 ,0 ) + 1
      fusion( i1,i2 ,0  ) = nn
      fusion( i1,i2 ,nn ) = i3
    elseif (i1.eq.i2.and.i2.ne.i3) then
      nn = fusion( i1,i2 ,0 ) + 1
      fusion( i1,i2 ,0  ) = nn
      fusion( i1,i2 ,nn ) = i3
      nn = fusion( i1,i3 ,0 ) + 1
      fusion( i1,i3 ,0  ) = nn
      fusion( i1,i3 ,nn ) = i2
      fusion( i3,i1 ,0  ) = nn
      fusion( i3,i1 ,nn ) = i2
    elseif (i3.eq.i1.and.i1.ne.i2) then
      nn = fusion( i3,i1 ,0 ) + 1
      fusion( i3,i1 ,0  ) = nn
      fusion( i3,i1 ,nn ) = i2
      nn = fusion( i3,i2 ,0 ) + 1
      fusion( i3,i2 ,0  ) = nn
      fusion( i3,i2 ,nn ) = i1
      fusion( i2,i3 ,0  ) = nn
      fusion( i2,i3 ,nn ) = i1
    elseif (i2.eq.i3.and.i3.ne.i1) then
      nn = fusion( i2,i3 ,0 ) + 1
      fusion( i2,i3 ,0  ) = nn
      fusion( i2,i3 ,nn ) = i1
      nn = fusion( i2,i1 ,0 ) + 1
      fusion( i2,i1 ,0  ) = nn
      fusion( i2,i1 ,nn ) = i3
      fusion( i1,i2 ,0  ) = nn
      fusion( i1,i2 ,nn ) = i3
    else
      nn = fusion( i1,i2 ,0 ) + 1
      fusion( i1,i2 ,0  ) = nn
      fusion( i1,i2 ,nn ) = i3
      fusion( i2,i1 ,0  ) = nn
      fusion( i2,i1 ,nn ) = i3
      nn = fusion( i3,i1 ,0 ) + 1
      fusion( i3,i1 ,0  ) = nn
      fusion( i3,i1 ,nn ) = i2
      fusion( i1,i3 ,0  ) = nn
      fusion( i1,i3 ,nn ) = i2
      nn = fusion( i2,i3 ,0 ) + 1
      fusion( i2,i3 ,0  ) = nn
      fusion( i2,i3 ,nn ) = i1
      fusion( i3,i2 ,0  ) = nn
      fusion( i3,i2 ,nn ) = i1
    endif
  enddo
  end subroutine
!
end module
