MODULE ED_SECTOR
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private


  public :: build_sector
  public :: delete_sector
  public :: show_sector
  !
  public :: get_Nup
  public :: get_DimUp
  !
  public :: get_sector_dimension
  public :: binomial
  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate



contains


  !##################################################################
  !##################################################################
  !BUILD SECTORS
  !##################################################################
  !##################################################################
  subroutine build_sector(isector,self)
    integer,intent(in)                  :: isector
    type(sector)                        :: self
    integer                             :: i,iup,idw,ipup,ipdw
    integer                             :: nup_,ndw_
    integer                             :: dim,iel,iimp
    !
    if(self%status)call delete_sector(self)
    !
    self%index = isector
    !
    allocate(self%H(1))
    !
    call get_Nup(isector,self%Nup)
    call get_DimUp(isector,self%DimUp)
    self%Dim=self%DimUp
    !
    call map_allocate(self%H,[self%DimUp])
    !UP    
    dim=0
    do iup=0,2**Ns-1
       nup_ = popcnt(iup)
       if(nup_ /= self%Nup)cycle
       dim  = dim+1
       self%H(1)%map(dim) = iup
    enddo
    !
    self%Nlanc = min(self%Dim,lanc_nGFiter)
    self%status=.true.
    !
  end subroutine build_sector


  subroutine delete_sector(self)
    type(sector) :: self
    call map_deallocate(self%H)
    if(allocated(self%H))deallocate(self%H)
    self%index=0
    self%DimUp=0
    self%Dim=0
    self%Nup=0
    self%Nlanc=0
    self%status=.false.
  end subroutine delete_sector



  subroutine show_sector(isector)
    integer,intent(in)                  :: isector
    integer                             :: i,iup,idw,ipup
    integer                             :: nup_,Nup
    integer                             :: dim,DimUp
    !
    !
    call get_Nup(isector,Nup)
    call get_DimUp(isector,DimUp)
    !
    !UP
    do iup=0,2**Ns-1
       if(popcnt(iup) /= nup) cycle
       write(*,"(I6,A1)",advance='no')iup,"-"
       call print_conf(iup,Ns,.true.)
    enddo
    !
  end subroutine show_sector





  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector






  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector





  subroutine get_Nup(isector,nup)
    integer                :: isector,nup
    Nup = getNup(isector)
  end subroutine get_Nup

  subroutine  get_DimUp(isector,DimUp)
    integer                :: isector,DimUp,Nup
    call get_Nup(isector,Nup)
    DimUp = get_sector_dimension(Nup)
  end subroutine get_DimUp





  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################
  elemental function get_sector_dimension(Nup) result(dim)
    integer,intent(in)          :: Nup
    integer                     :: dim
    dim = binomial(Ns,Nup)
  end function get_sector_dimension


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2>n1) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial


  subroutine print_conf(i,Ntot,advance)
    integer :: dim,i,j,Ntot
    logical :: advance
    integer :: ivec(Ntot)
    ivec = bdecomp(i,Ntot)
    write(LOGfile,"(A1)",advance="no")"|"
    write(LOGfile,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    if(advance)then
       write(LOGfile,"(A1)",advance="yes")">"
    else
       write(LOGfile,"(A1)",advance="no")">"
    endif
  end subroutine print_conf




end MODULE ED_SECTOR
