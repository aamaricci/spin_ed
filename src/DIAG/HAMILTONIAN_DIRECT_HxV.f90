! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec
#ifdef _MPI
  public  :: directMatVec_MPI
#endif



contains




  subroutine directMatVec(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    integer,dimension(Ns)               :: Nup
    real(8),dimension(Ns)               :: Sz
    complex(8),dimension(Ns,Ns)         :: Hij,Hloc
    real(8),dimension(Ns)               :: Hdiag
    integer                             :: i,j
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    Hdiag(:) = dreal(diagonal(Hloc(:,:)))
    !
    Hv=zero
    !
    !-----------------------------------------------!
    do j=MpiIstart,MpiIend
       m   = Hsector%H(1)%map(j)
       Nup = bdecomp(m,Ns)
       Sz  = Nup-0.5d0
       !
       !LOCAL HAMILTONIAN TERMS
       include "direct/HxV_diag.f90"
       !
       !HOPPING TERMS
       include "direct/HxV_hop.f90"
    enddo
    !-----------------------------------------------!
    !
    !
  end subroutine directMatVec




#ifdef _MPI
  subroutine directMatVec_MPI(Nloc,v,Hv)
    integer                             :: Nloc,N
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vin
    integer,dimension(Ns)              :: Nup
    real(8),dimension(Ns)              :: Sz
    complex(8),dimension(Ns,Ns)         :: Hij,Hloc
    real(8),dimension(Ns)               :: Hdiag
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index    
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    Hdiag(:) = dreal(diagonal(Hloc(:,:)))
    !
    !MPI part:
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    allocate(vin(N)) ; vin = zero
    call allgather_vector_MPI(MpiComm,v,vin)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    states: do j=MpiIstart,MpiIend
       m   = Hsector%H(1)%map(j)
       Nup = bdecomp(m,Ns)
       Sz  = Nup-0.5d0

       !LOCAL HAMILTONIAN TERMS
       include "direct/HxV_diag.f90"
       !
       !HOPPING TERMS
       include "direct/HxV_hop.f90"
       !
    enddo states
    !-----------------------------------------------!
    deallocate(vin)
    return
  end subroutine directMatVec_MPI
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
