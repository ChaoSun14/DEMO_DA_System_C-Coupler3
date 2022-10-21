module da_ensgather_demo_mod
    use model_setting_mod_2
    use spmd_init_mod_2
    use da_ensgather_demo_ccpl_coupling_mod
    use CCPL_interface_mod

    implicit none
    integer                    :: start_time, end_time, ierr
    real                       :: total_time, max_total_time, max_total_time1, total_time2, max_total_time2, total_time3, total_time4, max_total_time3, max_total_time4
    public :: da_ensgather_demo_ccpl_init
    public :: da_ensgather_demo_ccpl_run
    public :: da_ensgather_demo_ccpl_finalize

contains

    subroutine da_ensgather_demo_ccpl_init(ccpl_procedure_inst_id) bind(c)
        use mpi
        implicit none
!include 'mpif.h'
        integer, intent(in) ::  ccpl_procedure_inst_id
        logical             ::  da_ensgather_demo_log
        character(len=512)  :: log_file_name

        call system_clock(start_time)
        da_ensgather_demo_instance_id = ccpl_procedure_inst_id 
        da_ensgather_demo_ccpl_comm = CCPL_external_procedures_get_local_comm(ccpl_procedure_inst_id, annotation="get local comm of da_ensgather_demo")
        
        call spmd_init_da(da_ensgather_demo_ccpl_comm)
        if (masterpe) then
            write(*,*) "++++++++++++++++++++++++++++++++++++++++"
            write(*,*) "[CCPL <da_ensgather_demo> PE ", mytask_id, " ] Start da_ensgather_demo initialization"
            write(*,*) "[CCPL <da_ensgather_demo> PE ", mytask_id, " ] da_ensgather_demo_ccpl_comm: ", da_ensgather_demo_ccpl_comm
        end if
        !da_ensgather_demo_log  = .false.
        !da_ensgather_demo_log  = CCPL_get_comp_log_file_name(da_ensgather_demo_instance_id, log_file_name, annotation="get the logfile name of da_ensgather_demo" )
        !open(da_stdout, file=trim(log_file_name), status="UNKNOWN")
       
        call da_algorithm_demo_init

        
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, da_ensgather_demo_ccpl_comm, ierr)
        if (masterpe) then
            write(*,*) "[CCPL <da_ensgather_demo>] total time: ", max_total_time  
            write(*,*) "[CCPL <da_ensgather_demo> PE ", mytask_id, " ] Finish da_ensgather_demo initialization"
            write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        end if
   
    end subroutine da_ensgather_demo_ccpl_init
 
 
    subroutine da_ensgather_demo_ccpl_run() bind(c)
        
        use mpi
        implicit none

        call system_clock(start_time)
        if (masterpe) then
            write(*,*) "++++++++++++++++++++++++++++++++++++++++"
            write(*,*) "[CCPL <da_ensgather_demo> ENSEMBLE ", ccpL_ens_id," PE ", mytask_id, " ] Start da_ensgather_demo run"
        end if
        !call da_algorithm_demo_step_on
        if (masterpe) then
            print *, "da_algorithm_demo finished analysis"
        end if
        !call finalize_da_algorithm_demo
        if (masterpe) then
            print *, "da_algorithm_demo has been finalized"
        end if
    
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, da_ensgather_demo_ccpl_comm, ierr)
        if (masterpe) then
            write(*,*) "[CCPL <da_ensgather_demo>] run total time: ", max_total_time
            write(*,*) "[CCPL <da_ensgather_demo> ENSEMBLE ", ccpL_ens_id," PE ", mytask_id, " ] Finish da_ensgather_demo run"
            write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        end if
 
    end subroutine da_ensgather_demo_ccpl_run
     
    subroutine da_ensgather_demo_ccpl_finalize() bind(c)
         
        implicit none
 
    end subroutine da_ensgather_demo_ccpl_finalize


end module da_ensgather_demo_mod
