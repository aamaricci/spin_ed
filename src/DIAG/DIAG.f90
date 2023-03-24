module ED_DIAG
  USE SF_CONSTANTS
  USE SF_LINALG, only: eigh
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit,file_length
  USE SF_SP_LINALG
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN
  implicit none
  private


  public :: diagonalize_lattice
  public :: partition_function_lattice

  real(8),dimension(:),pointer       :: state_cvec


contains






  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine diagonalize_lattice()
    if(MPIMASTER)then
       write(LOGfile,"(A)")"Diagonalize H:"
       call start_timer()
    endif
    !
    call ed_diag_d
    !
    if(MPIMASTER)call stop_timer(unit=LOGfile)
  end subroutine diagonalize_lattice



  subroutine partition_function_lattice
    integer :: i,istate,isector
    integer :: Dim
    real(8) :: Egs,Ei
    integer :: Nup
    !
    zeta_function=0d0
    !
    if(MPIMASTER)then
       write(LOGfile,"(A)")"Get Z:"
       call start_timer()
    endif
    Egs = state_list%emin
    if(finiteT)then
       call es_trim_size(state_list,temp,cutoff)
       do istate=1,state_list%trimd_size
          Ei            = es_return_energy(state_list,istate)
          zeta_function = zeta_function + exp(-(Ei-Egs)/temp)
       enddo
    else
       zeta_function=dble(state_list%size)
    end if
    !
    if(MPIMASTER.AND.ed_verbose>=2)write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    !
    if(MPIMASTER)call stop_timer(unit=LOGfile)
    !
  end subroutine partition_function_lattice






  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine ed_diag_d
    integer                :: isector,Dim,istate
    integer                :: DimUp
    integer                :: Nup
    integer                :: i,j,iter,unit,vecDim,PvecDim
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin,Ei,Egs
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:),eig_basis_tmp(:,:)
    logical                :: lanc_solve,Tflag,lanc_verbose,bool
    integer                :: is,irank
    integer                :: nstates_below_cutoff
    real(8)                :: Ec
    !
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0    
    !
    lanc_verbose=.false.
    if(ed_verbose>2)lanc_verbose=.true.
    !
    iter=0
    sector: do isector=1,Nsectors
       call get_Nup(isector,nup)
       !
       if(ed_filling/=0 .AND. Nup/=ed_filling )cycle sector
       iter=iter+1
       !
       !
       Dim = getdim(isector)
       !
       select case (lanc_method)
       case default       !use P-ARPACK
          Neigen   = min(dim,neigen_sector(isector))
          Nitermax = min(dim,lanc_niter)
          Nblock   = min(dim,lanc_ncv_factor*max(Neigen,lanc_nstates_sector) + lanc_ncv_add)
       case ("lanczos")
          Neigen   = 1
          Nitermax = min(dim,lanc_niter)
          Nblock   = 1
       end select
       !
       lanc_solve  = .true.
       if(Neigen==dim)lanc_solve=.false.
       if(dim<=max(lanc_dim_threshold,MPISIZE))lanc_solve=.false.
       !
       if(MPIMASTER)then
          if(ed_verbose>1)then
             write(LOGfile,"(1X,I9,A,I9,A6,1I3,A7,I20,1X)",advance='no')&
                  iter,"-Solving sector:",isector,", nup:",nup,", dim=",getdim(isector)
             if(lanc_solve)write(LOGfile,"(A12,3I6,1X)",advance='no')"Lanc Info:",Neigen,Nitermax,Nblock
             write(LOGfile,*)""
          endif
       endif
       !
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
       if(ed_verbose>=3.AND.MPIMASTER)call start_timer()
       if(lanc_solve)then
          allocate(eig_values(Neigen)) ; eig_values=0d0 
          !
          call build_Hv_sector(isector) !For MPI: MpiComm==MpiComm_Global .OR. MpiComm subset of MpiComm_Global
          !
          vecDim = vecDim_Hv_sector(isector)
          allocate(eig_basis(vecDim,Neigen))
          eig_basis=zero
          !
          select case (lanc_method)
          case default       !use P-ARPACK
#ifdef _MPI
             if(MpiStatus)then
                call sp_eigh(MpiComm,spHtimesV_p,eig_values,eig_basis,&
                     Nblock,&
                     Nitermax,&
                     tol=lanc_tolerance,&
                     iverbose=(ed_verbose>3))
             else
                call sp_eigh(spHtimesV_p,eig_values,eig_basis,&
                     Nblock,&
                     Nitermax,&
                     tol=lanc_tolerance,&
                     iverbose=(ed_verbose>3))
             endif
#else
             call sp_eigh(spHtimesV_p,eig_values,eig_basis,&
                  Nblock,&
                  Nitermax,&
                  tol=lanc_tolerance,&
                  iverbose=(ed_verbose>3))
#endif             
             !
             !
          case ("lanczos")   !use Simple Lanczos
#ifdef _MPI
             if(MpiStatus)then
                call sp_lanc_eigh(MpiComm,spHtimesV_p,eig_values(1),eig_basis(:,1),Nitermax,&
                     iverbose=(ed_verbose>3),threshold=lanc_tolerance)
             else
                call sp_lanc_eigh(spHtimesV_p,eig_values(1),eig_basis(:,1),Nitermax,&
                     iverbose=(ed_verbose>3),threshold=lanc_tolerance)
             endif
#else
             call sp_lanc_eigh(spHtimesV_p,eig_values(1),eig_basis(:,1),Nitermax,&
                  iverbose=(ed_verbose>3),threshold=lanc_tolerance)
#endif
             !
             !          
          end select
          !
          !
          if(MpiMaster.AND.ed_verbose>3)write(LOGfile,*)""
          call delete_Hv_sector()
          call Bcast_MPI(MpiComm,eig_values)
          !
          !
       else                     !else LAPACK_SOLVE
          !
          !
          allocate(eig_values(Dim)) ; eig_values=0d0
          allocate(eig_basis_tmp(Dim,Dim)) ; eig_basis_tmp=0d0
          call build_Hv_sector(isector,eig_basis_tmp)
          if(MpiMaster)call eigh(eig_basis_tmp,eig_values)
          if(dim==1)eig_basis_tmp(dim,dim)=1d0
          !
          call delete_Hv_sector()
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,eig_values)
             vecDim = vecDim_Hv_sector(isector)
             allocate(eig_basis(vecDim,Neigen)) ; eig_basis=0d0
             call scatter_basis_MPI(MpiComm,eig_basis_tmp,eig_basis)
          else
             allocate(eig_basis(Dim,Neigen)) ; eig_basis=0d0
             eig_basis = eig_basis_tmp(:,1:Neigen)
          endif
#else
          allocate(eig_basis(Dim,Neigen)) ; eig_basis=0d0
          eig_basis = eig_basis_tmp(:,1:Neigen)
#endif
          !
       endif
       !
       if(ed_verbose>=3.AND.MPIMASTER)call stop_timer(unit=LOGfile)
       !
       if(ed_verbose>=4)then
          write(LOGfile,*)"EigValues: ",eig_values(:Neigen)
          write(LOGfile,*)""
          write(LOGfile,*)""
       endif
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(:,i),isector,size=lanc_nstates_total)
          enddo
       else
          do i=1,Neigen
             enemin = eig_values(i)
             if (enemin < oldzero-10.d0*gs_threshold)then
                oldzero=enemin
                call es_free_espace(state_list)
                call es_add_state(state_list,enemin,eig_basis(:,i),isector)
             elseif(abs(enemin-oldzero) <= gs_threshold)then
                oldzero=min(oldzero,enemin)
                call es_add_state(state_list,enemin,eig_basis(:,i),isector)
             endif
          enddo
       endif
       !
       !Print a list of the collected eigenvalues:
       if(MPIMASTER)then
          open(free_unit(unit),file="eigenvalues_list.ed",position='append',action='write')
          call print_eigenvalues_list(isector,eig_values(1:Neigen),unit,lanc_solve,mpiAllThreads)
          close(unit)
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis_tmp))deallocate(eig_basis_tmp)
       if(allocated(eig_basis))deallocate(eig_basis)
       !
    enddo sector
    !
    !
    !Print the obtained state_list, before any cut
    call save_state_list()
    if(ed_verbose>=2)call print_state_list(LOGfile)
    !
    !Print info about the ground state(s)
    if(MPIMASTER.AND.ed_verbose>=2)then
       do istate=1,es_return_gs_degeneracy(state_list,gs_threshold)
          isector = es_return_sector(state_list,istate)
          Egs     = es_return_energy(state_list,istate)
          call get_Nup(isector,Nup)
          write(LOGfile,"(A,F20.12,1I4)")'Egs =',Egs,nup
       enddo
    endif
    !
    !if finite T and loop over T we diagonalize first at the largest temperature.
    !to avoid keeping useless state we prune those above the cutoff at this T_max.
    if(finiteT)then
       !check if the number of states is enough to reach the required accuracy:
       !the condition to fullfill is:
       ! exp(-(Ec-Egs)/T) < \epsilon_c
       ! if this condition is violated then required number of states is increased
       ! if number of states is larger than those required to fullfill the cutoff: 
       ! trim the list and number of states.
       Egs  = state_list%emin
       Ec   = state_list%emax
       if(exp(-(Ec-Egs)/temp) > cutoff)stop "ED_DIAG:  exp(-(Ei-Egs)/T)>Cutoff condition not met, try increasing lanc_nstates_sector"
       write(LOGfile,*)
       Ei = es_return_energy(state_list,state_list%size)
       do while ( exp(-(Ei-Egs)/temp) <= cutoff )
          call es_pop_state(state_list)
          Ei = es_return_energy(state_list,state_list%size)
       enddo
       write(LOGfile,"(A,I4)")"Adjusting lanc_nstates_total to:",state_list%size
       !
    endif
  end subroutine ed_diag_d




















  !###################################################################################################
  !
  !    POST-PROCESSING ROUTINES
  !
  !###################################################################################################
  subroutine save_state_list()
    integer :: qn,isector
    integer :: istate
    integer :: unit
    if(MPIMASTER)then
       open(free_unit(unit),file="state_list.ed")
       do istate=1,state_list%size
          isector = es_return_sector(state_list,istate)
          call get_nup(isector,qn)
          write(unit,"(i8,i12,i8)")istate,isector,qn
       enddo
       close(unit)
    endif
  end subroutine save_state_list


  subroutine print_state_list(unit)
    integer :: qn,isector
    integer :: istate
    integer :: unit
    real(8) :: Estate
    if(MPIMASTER)then
       write(unit,"(A1,A6,A18,2x,A19,1x,2A10,A4)")"#","i","E_i","exp(-(E-E0)/T)","Sect","Dim","QN:"
       do istate=1,state_list%size
          Estate  = es_return_energy(state_list,istate)
          isector = es_return_sector(state_list,istate)
          call get_nup(isector,qn)
          write(unit,"(i6,f18.12,2x,ES19.12,1x,2I10,I4)")&
               istate,Estate,exp(-(Estate-state_list%emin)/temp),isector,getdim(isector),qn
       enddo
    endif
  end subroutine print_state_list



  subroutine print_eigenvalues_list(isector,eig_values,unit,lanc,allt)
    integer              :: isector
    real(8),dimension(:) :: eig_values
    integer              :: unit,i,qn
    logical              :: lanc,allt
    if(MPIMASTER)then
       if(lanc)then
          if(allt)then
             write(unit,"(A9,A15)")" # Sector","Indices"
          else
             write(unit,"(A10,A15)")" #T Sector","Indices"
          endif
       else
          write(unit,"(A10,A15)")" #X Sector","Indices"
       endif
       call get_Nup(isector,qn)
       write(unit,"(I9,I6)")isector,qn
       do i=1,size(eig_values)
          write(unit,*)eig_values(i)
       enddo
       write(unit,*)""
    endif
  end subroutine print_eigenvalues_list





end MODULE ED_DIAG













