module coupling_atm_model_mod

    use CCPL_interface_mod
    use spmd_init_mod, only     : mytask_id, npes, stdout, demo_log, masterproc
    use parse_namelist_mod,only : time_step, coupling_freq, num_3d_vars
    use grid_init_mod, only     : latlen, lonlen, levlen, lon, lat, lev
    use decomp_init_mod, only   : decomp_size, local_grid_cell_index
    use variable_mod, only      : sstm, sshm, saltm, saltm_1, vars_3d, TAUX, TAUY, PS
    !use mpi

    implicit none
include 'mpif.h'

    integer, public            :: ocn_demo_comp_id, ocn_comm
    character(len=512), public :: log_file_name
    integer, private           :: decomp_id, grid_h2d_id, ierr, grid_v1d_id, grid_m3d_id
    integer                    :: start_time, end_time
    real                       :: total_time, max_total_time, max_total_time1, total_time2, max_total_time2, total_time3, total_time4, max_total_time3, max_total_time4

#ifdef CCPL_DA
    integer, public            :: ensemble_comp_id, ensemble_comm
    integer, public            :: ocn_da_demo_individual_inst_id, ocn_da_demo_ensmean_inst_id, ocn_da_demo_ensgather_inst_id
    integer,  allocatable      :: local_cells_global_index_all(:)
    integer,  allocatable      :: da_field_ids(:)
    integer,  allocatable      :: control_vars(:)
    integer,  allocatable      :: comp_or_grid_ids(:)
    integer,  allocatable      :: decomp_ids(:)
    integer,  allocatable      :: field_inst_ids(:)
    integer,  allocatable      :: timer_ids(:)
    integer,  private          :: global_decomp_id
    integer                    :: year, month, day, hour, minute, second
    real(kind=8), allocatable  :: date
    character(len=3), public   :: num_str
#endif
    
contains

    subroutine register_ocn_demo_component(comm)
        implicit none
        integer, intent(inout) :: comm

        call system_clock(start_time)
#ifdef CCPL_DA
        ensemble_comp_id = CCPL_register_component(-1, "ensmbe_comp", "active_coupled_system", comm, change_dir=.false., annotation = "register model ensemble")
        ensemble_comm = comm
        ocn_comm = CCPL_NULL_COMM
        ocn_demo_comp_id = CCPL_register_component(ensemble_comp_id, "ocn_demo", "ocn", ocn_comm, change_dir=.true., annotation = "register ocn model ocn_demo")
        comm = ocn_comm
#else
        ocn_demo_comp_id = CCPL_register_component(-1, "ocn_demo", "ocn", comm, change_dir=.true., annotation = "register ocn model ocn_demo")
        ocn_comm = comm
#endif
        demo_log      = .false.
        demo_log      = CCPL_get_comp_log_file_name(ocn_demo_comp_id, log_file_name, annotation="get the logfile name of ocn_demo" )
        open(stdout, file=trim(log_file_name), status="UNKNOWN")
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, ocn_comm, ierr)
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] register_ocn_demo_component total time: ", max_total_time
        end if

    end subroutine register_ocn_demo_component

    subroutine register_component_coupling_configuration

        implicit none

        integer          :: export_interface_id, import_interface_id
        integer          :: timer_id, fields_id(5)
        integer          :: field_id_TAUX, field_id_TAUY, field_id_PS
        integer          :: field_id_sst, field_id_ssh, field_id_salt, field_id_salt_1
#ifdef CCPL_DA
        integer          :: i, j, k
        integer          :: field_id_lat, field_id_lon, field_id_date, field_id_lev
#endif
        call system_clock(start_time)
        call CCPL_set_normal_time_step(ocn_demo_comp_id, time_step, annotation="setting the time step for ocn_demo")
        if (mytask_id.eq.0) then
        
        endif
        !grid_h2d_id = CCPL_register_H2D_grid_via_global_data(ocn_demo_comp_id, "ocn_demo_H2D_grid", "LON_LAT", "degrees", "acyclic", lonlen, latlen, -999999., -999999.,-999999., -999999., real(lon(:,1),4), real(lat(1,:),4), mask, annotation="register ocn_demo H2D grid ")
        grid_h2d_id = CCPL_register_H2D_grid_via_global_data(ocn_demo_comp_id, "ocn_demo_H2D_grid", "LON_LAT", "degrees", "cyclic", lonlen, latlen, minval(real(lon,4)), maxval(real(lon,4)),minval(real(lat,4)), maxval(real(lat,4)), real(lon,4), real(lat,4), annotation="register ocn_demo H2D grid ")
        decomp_id = CCPL_register_normal_parallel_decomp("decomp_ocn_demo_grid", grid_H2D_id, decomp_size(mytask_id+1), local_grid_cell_index(1:decomp_size(mytask_id+1), mytask_id+1), annotation="allocate decomp for ocn_demo grid")

#ifdef CCPL_DA
    !    allocate(local_cells_global_index_all(lonlen*latlen))
    !    local_cells_global_index_all=CCPL_NULL_INT
    !    if(mytask_id<1) then
    !    k=0
    !    do j=1, latlen
    !        do i=1, lonlen
    !           k=k+1
    !           local_cells_global_index_all(k)=(j-1)*lonlen+i
    !        end do 
    !    end do
    !    end if
    !    global_decomp_id = CCPL_register_normal_parallel_decomp("global_decomp_atm_demo_grid", grid_h2d_id, lonlen*latlen, local_cells_global_index_all, annotation="allocate global decomp for atm_demo grid")
        !field_id_lat = CCPL_register_field_instance(lat, "LAT", global_decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="degree_north", annotation="register field instance of lat")
        !field_id_lon  = CCPL_register_field_instance(lon, "LON", global_decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="degree_east", annotation="register field instance of lon")
        grid_v1d_id = CCPL_register_V1D_Z_grid_via_model_data(ocn_demo_comp_id, "ocn_demo_vertical_grid", "unitless", lev(1:levlen), annotation="register vertical grid for ocn_demo")
        grid_m3d_id = CCPL_register_MD_grid_via_multi_grids(ocn_demo_comp_id, "ocn_demo_3d_grid", grid_H2D_id, grid_v1d_id, annotation="register 3D-grid for ocn_demo")
        field_id_lev = CCPL_register_field_instance(lev, "LEV", -1, grid_v1d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="unitless", annotation="register field instance of lev") 
        call CCPL_get_current_time(ocn_demo_comp_id, year, month, day, hour, minute, second, annotation = "get current time of atm_demo")
        allocate(date)
        date=real(year*10000000000+month*100000000+day*1000000+hour*3600+minute*60+second,8)
        
        field_id_date = CCPL_register_field_instance(date,"date", -1, ocn_demo_comp_id, 0, CCPL_TAG_CPL_REST, field_unit="seconds", annotation="register field instance of date")

#endif  
        !------------register field instances to C-Coupler2--------------

        field_id_TAUX = CCPL_register_field_instance(TAUX, "TAUX", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="N m-2", annotation="register field instance of TAUX") 
        field_id_TAUY = CCPL_register_field_instance(TAUY, "TAUY", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="N m-2", annotation="register field instance of TAUY")
        field_id_PS  = CCPL_register_field_instance(PS, "PS", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="Pa", annotation="register field instance of PS")
        field_id_sst  = CCPL_register_field_instance(sstm, "sst", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="C", annotation="register field instance of Sea surface temperature")
        field_id_ssh = CCPL_register_field_instance(sshm, "ssh", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of Sea surface height")
        field_id_salt = CCPL_register_field_instance(saltm, "salt", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="PSU", annotation="register field instance of Sea surface salinity")

        !--------register coupling frequency to C-Coupler2-------------
        timer_id = CCPL_define_single_timer(ocn_demo_comp_id, "seconds", coupling_freq, 0, 0, annotation="define a single timer for ocn_demo")

        !--------register export & import interface to C-Coupler2------
        fields_id(1) = field_id_sst
        fields_id(2) = field_id_ssh
        fields_id(3) = field_id_salt
        export_interface_id = CCPL_register_export_interface("send_data_to_atm", 3, fields_id, timer_id, annotation="register interface for sending data to atmosphere")

        fields_id(1) = field_id_TAUX
        fields_id(2) = field_id_TAUY
        fields_id(3) = field_id_PS
        import_interface_id = CCPL_register_import_interface("receive_data_from_atm", 3, fields_id, timer_id, 0, annotation="register interface for receiving data from atmosphere")
       
        !call CCPL_do_individual_coupling_generation(ocn_demo_comp_id, "Do coupling generation of ocn_demo")
        CALL CCPL_do_family_coupling_generation(ensemble_comp_id, ensemble_comp_id, "Do family coupling generation of the coupled model")
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, ocn_comm, ierr)
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] do_family_coupling_generation total time: ", max_total_time
        end if
#ifdef CCPL_DA
        field_id_salt_1 = CCPL_register_field_instance(saltm_1, "salt_3d", decomp_id, grid_m3d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="PSU", annotation="register field instance of 3d salinity")

        allocate(da_field_ids(num_3d_vars+5))
        da_field_ids(1) = field_id_date
        da_field_ids(2) = field_id_sst
        da_field_ids(3) = field_id_ssh
        da_field_ids(4) = field_id_salt_1
        da_field_ids(5) = field_id_lev

        do i = 1, num_3d_vars
            write(num_str, "(i3.3)") i
            da_field_ids(i+5) = CCPL_register_field_instance(vars_3d(:,:,i), "vars_3d_"//trim(adjustl(num_str)), decomp_id, grid_m3d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="PSU", annotation="register field instance of 3d vars")
        end do
        
        allocate(control_vars(1:4))
        control_vars(1) = levlen
        control_vars(2) = latlen
        control_vars(3) = lonlen
        control_vars(4) = npes
        
        if (masterproc) then
            write(*,*) "##############################################################################"
            write(*,*) "[CCPL <OCN_DEMO>] Start ensemble procedure ocn_da_demo_individual initialize"
        end if
        call mpi_barrier(ocn_comm, ierr)
        call system_clock(start_time)
        ocn_da_demo_individual_inst_id = CCPL_ensemble_procedures_inst_init(ocn_demo_comp_id, "ocn_da_demo_individual", da_field_ids, control_vars, annotation="do da procedure ocn_da_demo_individual initialize")
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, ocn_comm, ierr)
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] individual total time: ", max_total_time
            write(*,*) "[CCPL <OCN_DEMO>] Finish ensemble procedure ocn_da_demo_individual initialize"
            write(*,*) "##############################################################################"
        end if
        
        
        if (masterproc) then
            write(*,*) "##############################################################################"
            write(*,*) "[CCPL <OCN_DEMO>] Start ensemble procedure ocn_da_demo_ensmean initialize"
        end if
        call mpi_barrier(ocn_comm, ierr)
        call system_clock(start_time)
        ocn_da_demo_ensmean_inst_id = CCPL_ensemble_procedures_inst_init(ocn_demo_comp_id, "ocn_da_demo_ensmean", da_field_ids, control_vars, annotation="do da procedure ocn_da_demo_ensmean initialize")
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time1, 1, MPI_REAL, MPI_MAX, 0, ocn_comm, ierr)
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] ensmean total time: ", max_total_time1
            write(*,*) "[CCPL <OCN_DEMO>] Finish ensemble procedure ocn_da_demo_ensmean initialize"
            write(*,*) "##############################################################################"
        end if

        
        if (masterproc) then
            write(*,*) "##############################################################################"
            write(*,*) "[CCPL <OCN_DEMO>] Start ensemble procedure ocn_da_demo_ensgather initialize"
        end if
        call mpi_barrier(ocn_comm, ierr)
        call system_clock(start_time)
        ocn_da_demo_ensgather_inst_id = CCPL_ensemble_procedures_inst_init(ocn_demo_comp_id, "ocn_da_demo_ensgather", da_field_ids, control_vars, annotation="do da procedure ocn_da_demo_ensgather initialize")
        call system_clock(end_time)
        total_time = real(end_time-start_time)/10000.
        call mpi_reduce(total_time, max_total_time2, 1, MPI_REAL, MPI_MAX, 0, ocn_comm, ierr)
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] ensgather total time: ", max_total_time2
            write(*,*) "[CCPL <OCN_DEMO>] Finish ensemble procedure ocn_da_demo_ensgather initialize"
            write(*,*) "##############################################################################"
        end if

#endif      
        call CCPL_end_coupling_configuration(ocn_demo_comp_id, annotation = "component ocn_demo ends configuration")

    
    end subroutine register_component_coupling_configuration

    subroutine execute_ocn_coupling

        implicit none
        logical    :: interface_status
        
        interface_status = CCPL_execute_interface_using_name(ocn_demo_comp_id, "send_data_to_atm", .false., annotation = "execute interface for sending data to ocean")
          
        interface_status = CCPL_execute_interface_using_name(ocn_demo_comp_id, "receive_data_from_atm", .false., annotation = "execute interface for receiving data from ocean")
          
        !call CCPL_do_restart_write_IO(ocn_demo_comp_id,.false.)

#ifdef CCPL_DA
        call CCPL_get_current_time(ocn_demo_comp_id, year, month, day, hour, minute, second, annotation = "get current time of OCN_DEMO") 
        date=real(year*10000000000+month*100000000+day*1000000+hour*3600+minute*60+second,8)
        
        if (masterproc) then
            write(*,*) "OCN_DEMO date : ", year*10000+month*100+day,"-",hour,":",minute,":",second
            write(*,*) "OCN_DEMO date : ", date
            write(*,*) "OCN_DEMO SST : ", sstm(1:10)
        end if

        if (masterproc) then
            write(*,*) "##############################################################################"
            write(*,*) "[CCPL <OCN_DEMO>] Start ensemble procedure ocn_da_demo_individual run"
        end if
        call CCPL_ensemble_procedures_inst_run(ocn_da_demo_individual_inst_id, annotation="do ocn_da_demo_individual run")
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] Finish ensemble procedure ocn_da_demo_individual run"
            write(*,*) "##############################################################################"
            write(*,*) "##############################################################################"
            write(*,*) "[CCPL <OCN_DEMO>] Start ensemble procedure ocn_da_demo_ensmean run"
        end if
        call CCPL_ensemble_procedures_inst_run(ocn_da_demo_ensmean_inst_id, annotation="do ocn_da_demo_ensmean run") 
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] Finish ensemble procedure ocn_da_demo_ensmean run"
            write(*,*) "##############################################################################"
!
            write(*,*) "##############################################################################"
            write(*,*) "[CCPL <OCN_DEMO>] Start ensemble procedure ocn_da_demo_ensgather run"
        end if
        call CCPL_ensemble_procedures_inst_run(ocn_da_demo_ensgather_inst_id, annotation="do ocn_da_demo_ensgather run") 
        if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] Finish ensemble procedure ocn_da_demo_ensgather run"
            write(*,*) "##############################################################################"
        end if
        call CCPL_advance_time(ocn_demo_comp_id, "OCN_DEMO advances time for one step")
#endif

    end subroutine execute_ocn_coupling

    subroutine finalize_ocn_coupling

        implicit none

        CALL CCPL_finalize(.false., annotation="ocn_demo CALL CCPL_finalize") 

    end subroutine finalize_ocn_coupling

end module coupling_atm_model_mod
