module ED_GRAPH_MATRIX
  USE SF_IOTOOLS
  USE SF_TIMER,     only: start_timer,stop_timer
  USE SF_CONSTANTS, only: zero,one
  USE SF_LINALG,    only: eye,inv,diag
  USE SF_MISC,      only: assert_shape
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  type link
     integer            :: siteI,siteJ
     integer            :: orbI,orbJ
     complex(8)         :: Hval
     type(link),pointer :: next !link to next box (chain)
  end type link

  type Hij_structure
     type(link),pointer                    :: root !head/root of the list\== list itself
     character(len=:),allocatable          :: file     !Name of the W90 file 
     integer                               :: Size=0   !Number of hopping elements
     integer,allocatable,dimension(:)      :: Nsites    !Number of sites per orbital
     integer                               :: Ns=0     !total Number of lattice sites
     integer                               :: Norb=0   !Number of orbitals
     complex(8),allocatable,dimension(:,:) :: Hij    !The H(Ri,Rj)_ab^sr Hamiltonian [Ns,Ns]
     complex(8),allocatable,dimension(:,:) :: Hloc   !The local part of H(Ri,Ri)_ab^sr
     logical                               :: built  =.false. !Hij is built
     logical                               :: status =.false. !Allocated
  end type hij_structure

  interface Hij_get
     module procedure :: Hij_get_main
  end interface Hij_get

  interface Hij_local
     module procedure :: Hij_local_main
  end interface Hij_local


  public :: Hij_structure
  !
  public :: Hij_init
  public :: Hij_delete
  public :: Hij_info
  public :: Hij_add_link
  public :: Hij_read
  public :: Hij_build
  public :: Hij_get
  public :: Hij_local
  public :: Hij_write


  integer :: mpi_ierr
  integer :: mpi_rank
  integer :: mpi_size
  logical :: mpi_master


  type(Hij_structure),public :: Hmatrix


contains



  !< Setup/Delete the default structure
  subroutine Hij_init(Nsites)
    integer,dimension(:)      :: Nsites
    integer                   :: Ns
    integer                   :: Norb
    !
    if(Hmatrix%status)call Hij_delete()
    !
    Norb = size(Nsites)
    Ns   = sum(Nsites)
    !
    allocate(Hmatrix%root)
    Hmatrix%root%next => null()
    !
    allocate(Hmatrix%Nsites(Norb))
    Hmatrix%Nsites  = Nsites
    Hmatrix%file   = ""
    Hmatrix%size   = 0
    Hmatrix%Ns     = Ns
    Hmatrix%Norb   = Norb
    !
    allocate( Hmatrix%Hij(Ns,Ns) )
    Hmatrix%Hij    = zero
    allocate( Hmatrix%Hloc(Ns,Ns) )
    Hmatrix%Hloc   = zero
    Hmatrix%built  =.false.
    Hmatrix%status =.true.
    return
  end subroutine Hij_init
  !
  subroutine Hij_delete()
    type(link),pointer  :: p,c
    !
    if(.not.Hmatrix%status)return
    !
    do
       p => Hmatrix%root
       c => p%next    !current is the first node (root's next)
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next => null()
       call del_link(c)
       deallocate(c)
    end do
    deallocate(Hmatrix%root)
    deallocate(Hmatrix%Nsites)
    deallocate(Hmatrix%file)
    deallocate(Hmatrix%Hij)
    deallocate(Hmatrix%Hloc)
    Hmatrix%size   = 0
    Hmatrix%Ns     = 0
    Hmatrix%Norb   = 0
    Hmatrix%built  =.false.
    Hmatrix%status =.false.
  end subroutine Hij_delete

  subroutine Hij_info()
    !
    if(.not.Hmatrix%status)then
       write(*,"(A)")"Hij_info WARNING: Hmatrix not allocated"
       return
    endif
    !
    write(*,"(A,10I8)")"Nsites  =",Hmatrix%Nsites
    write(*,"(A,I8)")"Ns      =",Hmatrix%Ns
    write(*,"(A,I8)")"Norb    =",Hmatrix%Norb
    write(*,"(A,I8)")"# link  =",Hmatrix%size
    if(Hmatrix%file /= "")write(*,"(A,1x,A)")"file    =",Hmatrix%file
    write(*,"(A,L8)")"built H =",Hmatrix%built
    return
  end subroutine Hij_info





  !Add/Remove link to a given structure
  subroutine Hij_add_link(siteI,siteJ,orbI,orbJ,Hval)
    integer                           :: siteI,siteJ
    integer                           :: orbI,orbJ
    complex(8) ,intent(in)            :: Hval
    type(link),pointer                :: p,c
    integer                           :: k
    !
    !
    if(.not.Hmatrix%status)stop "Hij_add_link ERROR: Hmatrix not allocated"
    !
    !
    if(orbI>Hmatrix%Norb  .OR. orbJ>Hmatrix%Norb  .OR. orbI<=0 .OR. orbJ<=0)&
         stop "Hij_add_link ERROR: orbI or orbJ either > Norb OR < 0"
    if(siteI>Hmatrix%Nsites(orbI) .OR. siteJ>Hmatrix%Nsites(orbJ) .OR. siteI<=0 .OR. siteI<=0 )&
         stop "Hij_add_link ERROR: siteI or siteJ either > Nlat OR < 0"
    !
    p => Hmatrix%root
    c => p%next
    do k=1,Hmatrix%size                    !traverse the Hmatrix
       if(.not.associated(c))exit !beginning of the Hmatrix
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the Hmatrix
    !
    p%next%siteI = siteI
    p%next%siteJ = siteJ
    p%next%orbI  = orbI
    p%next%orbJ  = orbJ
    p%next%Hval  = Hval
    Hmatrix%size=Hmatrix%size+1
    !
    if(.not.associated(c))then !end of the Hmatrix special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine Hij_add_link




  !Set/Get/Delete link: INTERNAL USE
  subroutine set_link(self,siteI,siteJ,orbI,orbJ,Hval)
    type(link),intent(inout) :: self
    integer,intent(in)       :: siteI,siteJ
    integer,intent(in)       :: orbI,orbJ
    complex(8),intent(in)    :: Hval
    self%siteI = siteI
    self%siteJ = siteJ
    self%orbI  = orbI
    self%orbJ  = orbJ
    self%Hval  = Hval
  end subroutine set_link
  !

  subroutine get_link(self,siteI,siteJ,orbI,orbJ,Hval)
    type(link),intent(inout) :: self
    integer,intent(inout)    :: siteI,siteJ
    integer,intent(inout)    :: orbI,orbJ
    complex(8),intent(inout) :: Hval
    siteI = self%siteI
    siteJ = self%siteJ
    orbI  = self%orbI
    orbJ  = self%orbJ
    Hval  = self%Hval
  end subroutine get_link
  !
  subroutine del_link(self)
    type(link),intent(inout) :: self
    self%siteI = 0
    self%siteJ = 0
    self%orbI  = 0
    self%orbJ  = 0
    self%Hval  = zero
    self%next  => null()
  end subroutine del_link
  !
  subroutine print_link(self,unit)
    type(link),intent(in) :: self
    integer               :: unit
    write(unit,"(4I6,2F12.4)")&
         self%siteI,self%siteJ,&
         self%orbI,self%orbJ,&
         dreal(self%Hval),dimag(self%Hval)
  end subroutine print_link
  !
  function indices_link(siteI,orbI) result(Indx)
    integer,intent(in)          :: siteI,orbI
    integer                     :: Indx
    select case(orbI)
    case(1)
       Indx = siteI
    case default
       Indx = siteI + (orbI-1)*Hmatrix%Nsites(orbI-1)
    end select
  end function indices_link









  !< Read a given structure from file
  subroutine Hij_read(file)
    character(len=*)          :: file
    integer                   :: Nsize
    integer                   :: unitIO
    integer                   :: ih,i,j,a,b
    real(8)                   :: re,im
    !
    if(.not.Hmatrix%status)stop "Hij_read ERROR: Hmatrix not allocated"
    !
    !Read from file initial info OR reconstruct them is Header is missing
    open(free_unit(unitIO),file=reg(file),status="old",action="read")
    Nsize = file_length(reg(file))
    Hmatrix%file = reg(file)
    Hmatrix%Size = Nsize
    !
    do ih=1,Hmatrix%Size
       read(unitIO,*)i,j,a,b,re,im
       call Hij_add_link(i,j,a,b,dcmplx(re,im))
    enddo
    close(unitIO)
    return
  end subroutine Hij_read



  !< Build Hij from a given structure
  subroutine Hij_build()
    type(link),pointer :: c           
    integer            :: i,j,a,b
    integer            :: io,jo
    complex(8)         :: val
    !
    if(.not.Hmatrix%status)stop "Hij_build ERROR: Hmatrix not allocated"
    !
    c => Hmatrix%root%next
    do
       if(.not.associated(c))exit
       call get_link(c,i,j,a,b,val)
       io = indices_link(i,a)
       jo = indices_link(j,b)
       Hmatrix%Hij(io,jo) = val
       if(i==j)Hmatrix%Hloc(io,jo) = val
       c => c%next
    enddo
    if(.not.herm_check(Hmatrix%Hij(:,:),1d-6))stop "Hij_build ERROR: Hij is not Hermitian"
    Hmatrix%built=.true.
    return
  end subroutine Hij_build




  subroutine Hij_get_main(Hij)
    complex(8),dimension(:,:) :: Hij
    integer                   :: Nlso
    if(.not.Hmatrix%status)stop "Hij_get ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Ns
    call assert_shape(Hij,[Nlso,Nlso],"Hij_get","Hij")
    if(.not.Hmatrix%built)call Hij_build()
    Hij = Hmatrix%Hij
  end subroutine Hij_Get_Main




  subroutine Hij_local_main(Hloc)
    complex(8),dimension(:,:) :: Hloc
    integer                   :: Nlso
    if(.not.Hmatrix%status)stop "Hij_local ERROR: Hmatrix not allocated"
    !
    Nlso = Hmatrix%Ns
    call assert_shape(Hloc,[Nlso,Nlso],"Hij_local","Hloc")
    if(.not.Hmatrix%built)call Hij_build()
    Hloc = Hmatrix%Hloc
  end subroutine Hij_local_main
  





  subroutine Hij_write(unit,file)
    character(len=*),optional :: file
    integer,optional          :: unit
    integer                   :: unit_
    type(link),pointer        :: c           
    logical                   :: mpi_master
    !
    if(.not.Hmatrix%status)stop "Hij_write ERROR: Hmatrix not allocated"
    !
    mpi_master=.true.
#ifdef _MPI    
    if(check_MPI())mpi_master= get_master_MPI()
#endif
    !
    if(.not.Hmatrix%status)stop "write_Hmatrix ERROR: Hmatrix structure not allocated. Call setup_Hmatrix first."
    !
    if(mpi_master)then
       unit_ = 6 ; if(present(unit))unit_=unit
       if(present(file))open(free_unit(unit_),file=reg(file))
       c => Hmatrix%root%next
       do
          if(.not.associated(c))exit
          call print_link(c,unit_)
          c => c%next
       enddo
       if(present(file))close(unit_)
    endif
    return
  end subroutine Hij_write











  function herm_check(A,err) result(bool)
    complex(8),intent(in) :: A(:,:)
    real(8),optional      :: err
    real(8)               :: err_
    logical               :: bool
    err_ = 0d0;if(present(err))err_=err
    bool = .true.
    if( any( abs(A-conjg(transpose(A))) > err_ ) )bool=.false.
  end function herm_check






END MODULE ED_GRAPH_MATRIX
















