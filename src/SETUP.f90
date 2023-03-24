MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SECTOR
  USE SF_TIMER
  USE SF_MISC, only: sort
  USE SF_IOTOOLS, only:free_unit,reg,file_length,save_array
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private



  public :: init_ed_structure
  public :: setup_global

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure()
    logical                          :: control
    integer                          :: i,iud,iorb,jorb,ispin,jspin,unit
    integer                          :: Nup,Ndw,Ntot,DimUp
    integer                          :: Tstep,Dim
    integer,allocatable              :: Tord(:)
    logical                          :: Tbool
    !
    !
    !>Setup Dimensions of the problem
    !
    Ns   = sum(Nsites(1:Norb))
    Ntot = Ns
    if(Norb>2)stop "ED ERROR: Norb > 2 is currently not supported"
    !
    !
    if(lanc_method=="lanczos")then
       if(ed_finite_temp)stop "ED ERROR: lanc_method==lanczos available only for T=0"
       if(lanc_nstates_sector>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_sector==1, T=0"
    endif
    !
    lanc_nstates_total=1
    !
    Nsectors = (Ns+1)
    !
    !
    if(MpiMaster)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I6)")'# of spin levels      = ',Ns
       write(LOGfile,"(A,I6)")'# of orbitals         = ',Norb
       write(LOGfile,"(A,3I6)")'# Nsites              = ',Nsites
       write(LOGfile,"(A,I6)")'# of  sectors         = ',Nsectors
    endif
    !
    Nup = Ntot/2
    DimUp = get_sector_dimension(Nup)
    if(MpiMaster)then
       write(LOGfile,"(A,1I8)")&
            'Largest Sector(s)     = ',DimUp
       write(LOGfile,"(A)")"--------------------------------------------"
    endif
    !
    !
    !Allocate indexing arrays
    allocate(getDim(Nsectors));getDim=0
    allocate(getNup(Nsectors));getNup=0
    allocate(getSector(0:Ntot));getSector=0
    allocate(neigen_sector(Nsectors))
    !
    !
    finiteT = ed_finite_temp
    !
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       !
       lanc_nstates_total=lanc_nstates_sector*Nsectors+10
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif
       write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
       write(LOGfile,"(A,I6)")"Nstates   Total  = ", lanc_nstates_total
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
    endif
    !
    !
    Tstep = 1
    allocate(temperature_list(Tstep))
    temperature_list = temp
    if(finiteT)then
       inquire(file=trim(Tfile)//".restart",exist=Tbool)
       if(Tbool)then
          deallocate(temperature_list)
          write(LOGfile,"(A)")'Reading temperature list from file '//trim(Tfile)//".restart"
          Tstep = file_length(trim(Tfile)//".restart")
          open(free_unit(unit),file=trim(Tfile)//".restart")
          allocate(temperature_list(Tstep),Tord(Tstep))
          do i=1,Tstep
             read(unit,*)temperature_list(i)
          enddo
          close(unit)
          call sort(temperature_list,Tord)                !sort from smallest to largest
          temperature_list = temperature_list(Tstep:1:-1) !invert order
          call save_array(trim(Tfile)//".used",temperature_list)
          temp             = temperature_list(1)          !set actual Temp to largest
       endif
    endif
    !
  end subroutine init_ed_structure





  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global
    integer                          :: DimUp
    integer                          :: i,iud,iorb,Nup
    integer                          :: isector,jsector,gsector,ksector,lsector
    integer                          :: unit,status,istate,ishift,isign,dim
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    type(sector) :: sectorI,sectorJ,sectorK,sectorG,sectorL
    !
    !Store full dimension of the sectors:
    isector=0
    do Nup=0,Ns
       isector=isector+1
       getNup(isector)  = Nup
       getDim(isector)  = get_sector_dimension(Nup)
       ! if(ed_verbose>4.AND.MpiMaster.AND.Ns<9)write(*,"(A,I6,A3,I2,A)")"--- sector",isector,"  (",Nup,")---- "       
       ! if(ed_verbose>4.AND.MpiMaster.AND.Ns<9)call show_sector(isector)
       ! if(ed_verbose>4.AND.MpiMaster.AND.Ns<9)write(*,"(A,2I6)")"--- DIM:",getDim(isector)
       ! if(ed_verbose>4.AND.MpiMaster.AND.Ns<9)write(*,*)""
    enddo
    !
    !
    do isector=1,Nsectors
       neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector)
    enddo
    !
  end subroutine setup_global

  




end MODULE ED_SETUP
