! > BUILD SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_STORED_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Matric constructors
  public :: ed_buildh

  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec
#ifdef _MPI
  public  :: spMatVec_MPI
#endif


contains


  subroutine ed_buildh(Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    integer                               :: isector,ispin,i,j
    complex(8),dimension(:,:),allocatable :: Htmp_up,Hrdx,Hmat_tmp
    integer,dimension(1)                  :: QNs    ![1/Norb]
    integer,dimension(Ns)                 :: Nup
    real(8),dimension(Ns)                 :: Sz 
    complex(8),dimension(Ns,Ns)           :: Hij,Hloc
    real(8),dimension(Ns)                 :: Hdiag
    !
#ifdef _MPI
    if(Mpistatus .AND. MpiComm == MPI_COMM_NULL)return
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"ed_buildh_main","Hmat")
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    Hdiag(:) = dreal(diagonal(Hloc(:,:)))
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0d,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0d,DimUp)
    else
       call sp_init_matrix(spH0d,DimUp)
    endif
#else
    call sp_init_matrix(spH0d,DimUp)
#endif
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN TERMS
    include "stored/H_diag.f90"    
    !
    ! !NON-LOCAL HAMILTONIAN TERMS
    include "stored/H_hop.f90"
    !
    !-----------------------------------------------!
    if(present(Hmat))then
       Hmat = zero
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0d,Hmat)
       else
          call sp_dump_matrix(spH0d,Hmat)
       endif
#else
       call sp_dump_matrix(spH0d,Hmat)
#endif
    endif
    !
    return
    !
  end subroutine ed_buildh

















  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial
  ! - MPI
  !+------------------------------------------------------------------+
  subroutine spMatVec(Nloc,v,Hv)
    integer                         :: Nloc
    complex(8),dimension(Nloc)      :: v
    complex(8),dimension(Nloc)      :: Hv
    integer                         :: i,j
    Hv=zero
    do i=1,Nloc
       matmul: do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       end do matmul
    end do
  end subroutine spMatVec



#ifdef _MPI
  subroutine spMatVec_mpi(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i,j,mpiIerr
    integer                             :: N,MpiShift
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: Counts,Offset
    !
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    !
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    MpiShift = spH0d%Ishift
    Hv=0d0
    do i=1,Nloc
       local: do j=1,spH0d%loc(i)%Size
          Hv(i) = Hv(i) + spH0d%loc(i)%vals(j)*v(spH0d%loc(i)%cols(j)-MpiShift)
       end do local
    end do
    !
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)) ; vin = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin      ,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    do i=1,Nloc                 !==spH0d%Nrow
       matmul: do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*vin(spH0d%row(i)%cols(j))
       end do matmul
    end do
    !
  end subroutine spMatVec_mpi
#endif



end MODULE ED_HAMILTONIAN_STORED_HXV







