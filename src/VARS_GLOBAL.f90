MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only: str,free_unit
  USE ED_SPARSE_MATRIX
  USE ED_GRAPH_MATRIX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none





  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!

  type sector_map
     integer,dimension(:),allocatable :: map
     logical                          :: status=.false.
  end type sector_map

  type sector
     integer                                     :: index       !
     type(sector_map),dimension(:),allocatable   :: H
     integer                                     :: DimUp
     integer                                     :: Dim
     integer                                     :: Nup
     integer                                     :: Nlanc
     logical                                     :: status=.false.
  end type sector



  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
  end interface




  !-------------------------- ED  VARIABLES --------------------------!

  !SIZE OF THE PROBLEM |\s_1....\s_Ns>
  !=========================================================
  integer,save                                     :: Ns       !Number spin levels
  integer,save                                     :: Nsectors !Number of sectors


  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:)                 :: getDim             ! [Nsectors]
  integer,allocatable,dimension(:)                 :: GetSector
  integer,allocatable,dimension(:)                 :: getNup      ! [Nsectors]

  !internal Hmatrix storage
  !PRIVATE
  !=========================================================  
  real(8),allocatable,dimension(:,:,:,:,:,:)       :: impHij

  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix_csr)                          :: spH0d
  procedure(dd_sparse_HxV),pointer                 :: spHtimesV_p=>null()


  !Variables for DIAGONALIZATION
  !=========================================================  
  integer,allocatable,dimension(:)                 :: neigen_sector
  logical                                          :: trim_state_list=.false.

  !Partition function
  !=========================================================
  real(8)                                          :: zeta_function
  real(8)                                          :: gs_energy
  

  real(8),dimension(:),allocatable                 :: temperature_list


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                 :: wm,tau,wr,vm,vr




  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                :: ed_file_suffix=""       !suffix string attached to the output files.
  integer                                          :: site_indx_padding=4
  logical                                          :: finiteT             !flag for finite temperature calculation
  integer                                          :: lanc_nstates_total=1  !Max number of states hold in the finite T calculation


  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                          :: MpiComm_Global=MPI_COMM_NULL
  integer                                          :: MpiComm=MPI_COMM_NULL
  integer                                          :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                          :: MpiGroup=MPI_GROUP_NULL
#endif
  logical                                          :: MpiStatus=.false.
  logical                                          :: MpiMaster=.true.
  integer                                          :: MpiRank=0
  integer                                          :: MpiSize=1
  integer,allocatable,dimension(:)                 :: MpiMembers
  integer                                          :: mpiQup=0
  integer                                          :: mpiRup=0
  integer                                          :: mpiQdw=0
  integer                                          :: mpiRdw=0
  integer                                          :: mpiQ=0
  integer                                          :: mpiR=0
  integer                                          :: mpiIstart
  integer                                          :: mpiIend
  integer                                          :: mpiIshift
  logical                                          :: mpiAllThreads=.true.
  !

contains

  !IF code is compiled with MPI support
  !+  MPI is initialized:
  !THEN this routine setup the internal communicator
  !(inherited from MPI_COMM_WORLD) plus global variables
  !ELSE it does nothing
  !
  !
  subroutine ed_set_MpiComm()
#ifdef _MPI
    integer :: ierr,irank
    if(check_MPI())then
       MpiComm_Global = MPI_COMM_WORLD
       MpiComm        = MPI_COMM_WORLD
       call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
       MpiStatus      = .true.
       MpiSize        = get_Size_MPI(MpiComm_Global)
       MpiRank        = get_Rank_MPI(MpiComm_Global)
       MpiMaster      = get_Master_MPI(MpiComm_Global)
    endif
#endif
  end subroutine ed_set_MpiComm


  !IF code is compiled with MPI support
  !THEN this routine reset global variables to default values (SERIAL)
  subroutine ed_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_COMM_NULL
    MpiComm        = MPI_COMM_NULL
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#endif
  end subroutine ed_del_MpiComm



END MODULE ED_VARS_GLOBAL
