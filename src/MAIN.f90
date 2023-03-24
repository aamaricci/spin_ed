module ED_MAIN
  USE SF_IOTOOLS, only: str,reg,free_unit,file_length
  USE SF_TIMER,only: start_timer,stop_timer  
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace,delete_eigenspace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_DIAG
  implicit none
  private
  !
  public :: ed_init_solver
  public :: ed_solve



contains





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate Memory and Initialize ED -+!
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver()
    integer             :: unit,itemp

    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm()
    !
    if(MpiMaster)write(LOGfile,"(A)")"INIT ED SOLVER"
    !
    !Init ED Structure & memory
    call init_ed_structure()
    !
    !
    !Check Hmatrix is allocated:
    if(.not.Hmatrix%status)stop "ED_INIT_SOLVER ERROR: Hmatrix is not allocated"
    !
    call setup_global
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve()
    integer :: itemp
    !
    !SET THE LOCAL MPI COMMUNICATOR
    call ed_set_MpiComm() 
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(.not.Hmatrix%status)stop "ED_INIT_SOLVER ERROR: Hmatrix is not allocated"
    call Hij_write(unit=LOGfile)
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_lattice            !-> store state_list, independent of TEMP
    !
    do itemp=1,size(temperature_list)
       temp = temperature_list(itemp)
       if(size(temperature_list)>1)then
          ed_file_suffix=str(temp,lead=6)
          write(LOGfile,"(A)")""
          write(LOGfile,"(A)")'Solving for temperature T='//str(temp,lead=6)
       endif
       call partition_function_lattice !-> get trimmed state_list
       ! call observables_lattice        !-> get static observables
       ed_file_suffix=""
    enddo
    !
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    nullify(spHtimesV_p)
    write(Logfile,"(A)")""
  end subroutine ed_solve

end module ED_MAIN







