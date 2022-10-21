module coupling_ocn_model_mod

    use CCPL_interface_mod
    use spmd_init_mod, only      : mytask_id, npes, stdout, demo_log, masterproc
    use parse_namelist_mod, only : time_step, coupling_freq
    use grid_init_mod, only      : latlen, lonlen, lon, lat, lev, levlen
    use decomp_init_mod, only    : decomp_size, local_grid_cell_index
    use variable_mod, only       : TAUX_m, TAUY_m, PS_m, T_m, sst, ssh, salt

    implicit none

    integer, public            :: atm_demo_comp_id, atm_comm
    integer, private           :: decomp_id, grid_h2d_id
    character(len=512), public :: log_file_name


    integer, public            :: ensemble_comp_id, ensemble_comm
#ifdef CCPL_DA
    integer, public            :: da_demo_individual_inst_id, da_demo_ensmean_inst_id, da_demo_ensgather_inst_id
    integer, private           :: grid_v1d_id, grid_m3d_id, global_decomp_id
    integer,  allocatable      :: local_cells_global_index_all(:)
    integer,  allocatable      :: da_field_ids(:)
    integer,  allocatable      :: control_vars(:)
    integer,  allocatable      :: comp_or_grid_ids(:)
    integer,  allocatable      :: decomp_ids(:)
    integer,  allocatable      :: field_inst_ids(:)
    integer,  allocatable      :: timer_ids(:)
    integer                    :: year, month, day, hour, minute, second
    real(kind=8),  allocatable :: date
#endif
    
contains

    subroutine register_atm_demo_component(comm)

        implicit none
        integer, intent(inout) :: comm

!#ifdef CCPL_DA
        ensemble_comp_id = CCPL_register_component(-1, "ensmbe_comp", "active_coupled_system", comm, change_dir=.false., annotation = "register model ensemble")
        ensemble_comm = comm
        atm_comm = CCPL_NULL_COMM
        atm_demo_comp_id = CCPL_register_component(ensemble_comp_id, "atm_demo", "atm", atm_comm, change_dir=.true., annotation = "register atm model atm_demo")
        comm = atm_comm
!#else
        !atm_demo_comp_id = CCPL_register_component(-1, "atm_demo", "atm", comm, change_dir=.true., annotation = "register atm model atm_demo")
        !atm_comm = comm
!#endif
        demo_log      = .false.
        demo_log      = CCPL_get_comp_log_file_name(atm_demo_comp_id, log_file_name, annotation="get the logfile name of atm_demo" )
        open(stdout, file=trim(log_file_name), status="UNKNOWN")

    end subroutine register_atm_demo_component

    subroutine register_component_coupling_configuration

        implicit none

        integer          :: export_interface_id, import_interface_id
        integer          :: timer_id, fields_id(5), landmask(latlen*lonlen)
        integer          :: field_id_TAUX, field_id_TAUY, field_id_PS, field_id_T
        integer          :: field_id_sst, field_id_ssh, field_id_salt
#ifdef CCPL_DA
        integer          :: i, j, k
        integer          :: field_id_lev, field_id_date
#endif

        call CCPL_set_normal_time_step(atm_demo_comp_id, time_step, annotation="setting the time step for atm_demo")
        
        grid_h2d_id = CCPL_register_H2D_grid_via_global_data(atm_demo_comp_id, "atm_demo_H2D_grid", "LON_LAT", "degrees", "cyclic", lonlen, latlen, -999999., -999999.,-999999., -999999., lon, lat, annotation="register atm_demo H2D grid ")
        decomp_id = CCPL_register_normal_parallel_decomp("decomp_atm_demo_grid", grid_H2D_id, decomp_size(mytask_id+1), local_grid_cell_index(1:decomp_size(mytask_id+1),mytask_id+1), annotation="allocate decomp for atm_demo grid")
#ifdef CCPL_DA
  !      allocate(local_cells_global_index_all(lonlen*latlen))
  !      local_cells_global_index_all=CCPL_NULL_INT
  !      if(mytask_id<1) then
  !      k=0
  !      do j=1, latlen
  !          do i=1, lonlen
  !             k=k+1
  !             local_cells_global_index_all(k)=(j-1)*lonlen+i
  !          end do 
  !      end do
  !      end if
        
  !      global_decomp_id = CCPL_register_normal_parallel_decomp("global_decomp_atm_demo_grid", grid_h2d_id_mask, lonlen*latlen, local_cells_global_index_all, annotation="allocate global decomp for atm_demo grid")
        grid_v1d_id = CCPL_register_V1D_Z_grid_via_model_data(atm_demo_comp_id, "atm_demo_vertical_grid", "unitless", lev(1:levlen), annotation="register vertical grid for atm_demo")
        grid_m3d_id = CCPL_register_MD_grid_via_multi_grids(atm_demo_comp_id, "atm_demo_3d_grid", grid_H2D_id, grid_v1d_id, annotation="register 3D-grid for atm_demo")
        field_id_lev = CCPL_register_field_instance(lev, "LEV", -1, grid_v1d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="unitless", annotation="register field instance of lev") 
        !field_id_lat = CCPL_register_field_instance(lat, "LAT", global_decomp_id, grid_h2d_id_mask, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="degree_north", annotation="register field instance of lat")
        !field_id_lon  = CCPL_register_field_instance(lon, "LON", global_decomp_id, grid_h2d_id_mask, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="degree_east", annotation="register field instance of lon")
        field_id_t  = CCPL_register_field_instance(t_m, "T", decomp_id, grid_m3d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="K", annotation="register field instance of temperature")
        call CCPL_get_current_time(atm_demo_comp_id, year, month, day, hour, minute, second, annotation = "get current time of atm_demo")
        allocate(date)
        date=real(year*10000000000+month*100000000+day*1000000+hour*3600+minute*60+second,8)
        
        field_id_date = CCPL_register_field_instance(date,"date", -1, atm_demo_comp_id, 0, CCPL_TAG_CPL_REST, field_unit="seconds", annotation="register field instance of date")

#endif        
        !------------register field instances to C-Coupler2--------------

        field_id_TAUX = CCPL_register_field_instance(TAUX_m, "TAUX", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="N m-2", annotation="register field instance of TAUX") 
        field_id_TAUY = CCPL_register_field_instance(TAUY_m, "TAUY", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="N m-2", annotation="register field instance of TAUY")
        field_id_PS  = CCPL_register_field_instance(PS_m, "PS", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="Pa", annotation="register field instance of PS")
        field_id_sst  = CCPL_register_field_instance(sst, "sst", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="C", annotation="register field instance of Sea surface temperature")
        field_id_ssh = CCPL_register_field_instance(ssh, "ssh", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of Sea surface height")
        field_id_salt = CCPL_register_field_instance(salt, "salt", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, field_unit="PSU", annotation="register field instance of Sea surface salinity")

        !--------register coupling frequency to C-Coupler2-------------
        timer_id = CCPL_define_single_timer(atm_demo_comp_id, "seconds", coupling_freq, 0, 0, annotation="define a single timer for atm_demo")

        !--------register export & import interface to C-Coupler2------
        fields_id(1) = field_id_sst
        fields_id(2) = field_id_ssh
        fields_id(3) = field_id_salt
        import_interface_id = CCPL_register_import_interface("receive_data_from_ocn", 3, fields_id, timer_id, 0, annotation="register interface for receiving data from atmosphere")

        fields_id(1) = field_id_TAUX
        fields_id(2) = field_id_TAUY
        fields_id(3) = field_id_PS
        export_interface_id = CCPL_register_export_interface("send_data_to_ocn", 3, fields_id, timer_id, annotation="register interface for sending data to atmosphere")
        
        !call CCPL_do_individual_coupling_generation(atm_demo_comp_id, "Do coupling generation of atm_demo")
        CALL CCPL_do_family_coupling_generation(ensemble_comp_id, ensemble_comp_id, "Do family coupling generation of the coupled model")
#ifdef CCPL_DA
       
        allocate(comp_or_grid_ids(1:3))
        comp_or_grid_ids(1) = grid_v1d_id
        comp_or_grid_ids(2) = grid_h2d_id
        comp_or_grid_ids(3) = grid_m3d_id
        allocate(decomp_ids(1:3))
        decomp_ids(1) = -1
        decomp_ids(2) = decomp_id
        decomp_ids(3) = -1
        allocate(da_field_ids(1:6))
        da_field_ids(1) = field_id_date
        da_field_ids(2) = field_id_lev
        da_field_ids(3) = field_id_TAUX
        da_field_ids(4) = field_id_TAUY
        da_field_ids(5) = field_id_PS
        da_field_ids(6) = field_id_t
        allocate(control_vars(1:5))
        control_vars(1) = levlen
        control_vars(2) = latlen
        control_vars(3) = lonlen
        control_vars(4) = npes
        control_vars(5) = decomp_size(mytask_id+1)
        allocate(timer_ids(1))
        timer_ids(1) = timer_id
        
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <ATM_DOME>] Start ensemble procedure da_demo_individual initialize"
        write(*,*) "[CCPL <ATM_DOME>] atm_demo_comm: ", atm_comm
        da_demo_individual_inst_id = CCPL_ensemble_procedures_inst_init(atm_demo_comp_id, "da_demo_individual", da_field_ids, control_vars, annotation="do da procedure da_demo_individual initialize")
        write(*,*) "[CCPL <ATM_DOME>] Finish ensemble procedure da_demo_individual initialize"
        write(*,*) "==============================================================================="

        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <ATM_DOME>] Start ensemble procedure da_demo_ensmean initialize"
        write(*,*) "[CCPL <ATM_DOME>] atm_demo_comm: ", atm_comm
        da_demo_ensmean_inst_id = CCPL_ensemble_procedures_inst_init(atm_demo_comp_id, "da_demo_ensmean", da_field_ids, control_vars, annotation="do da procedure da_demo_ensmean initialize")
        write(*,*) "[CCPL <ATM_DOME>] Finish ensemble procedure da_demo_ensmean initialize"
        write(*,*) "==============================================================================="

        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <ATM_DOME>] Start ensemble procedure da_demo_ensgather initialize"
        da_demo_ensgather_inst_id = CCPL_ensemble_procedures_inst_init(atm_demo_comp_id, "da_demo_ensgather", da_field_ids, control_vars, annotation="do da procedure da_demo_ensgather initialize")
        write(*,*) "[CCPL <ATM_DOME>] Finish ensemble procedure da_demo_ensgather initialize"
        write(*,*) "==============================================================================="
#endif      
        
        call CCPL_end_coupling_configuration(atm_demo_comp_id, annotation = "component atm_demo ends configuration")

    end subroutine register_component_coupling_configuration

    subroutine execute_atm_coupling

        implicit none
        logical    :: interface_status

        
        interface_status = CCPL_execute_interface_using_name(atm_demo_comp_id, "send_data_to_ocn", .false., annotation = "execute interface for sending data to atmosphere")
          
        interface_status = CCPL_execute_interface_using_name(atm_demo_comp_id, "receive_data_from_ocn", .false., annotation = "execute interface for receiving data from atmosphere")
        
        write(stdout,*) "TAUX : ", TAUX_m(1:20)

        !call CCPL_do_restart_write_IO(atm_demo_comp_id,.false.)

#ifdef CCPL_DA
        call CCPL_get_current_time(atm_demo_comp_id, year, month, day, hour, minute, second, annotation = "get current time of atm_demo")
        date=real(year*10000000000+month*100000000+day*1000000+hour*3600+minute*60+second,8)
        
        if (masterproc) then
            write(*,*) "ATM_DEMO date : ", year*10000000000+month*100000000+day*1000000+hour*10000+minute*100+second
            write(*,*) "ATM_DEMO date : ", date
            write(*,*) "ATM_DEMO PS : ", PS_m(1:10)
        end if
        if (masterproc) then
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <ATM_DOME>] Start ensemble procedure da_demo_individual run"
        end if
        call CCPL_ensemble_procedures_inst_run(da_demo_individual_inst_id, annotation="do da_demo_individual run")
        if (masterproc) then
        write(*,*) "[CCPL <ATM_DOME>] Finish ensemble procedure da_demo_individual run"
        write(*,*) "==============================================================================="
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <ATM_DOME>] Start ensemble procedure da_demo_ensmean run"
        end if
        call CCPL_ensemble_procedures_inst_run(da_demo_ensmean_inst_id, annotation="do da_demo_ensmean run") 
        if (masterproc) then
        !call CCPL_advance_time(atm_demo_comp_id, "atm_demo advances time for one step")
        write(*,*) "[CCPL <ATM_DOME>] Finish ensemble procedure da_demo_ensmean run"
        write(*,*) "==============================================================================="

        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <ATM_DOME>] Start ensemble procedure da_demo_ensgather run"
        end if
        call CCPL_ensemble_procedures_inst_run(da_demo_ensgather_inst_id, annotation="do da_demo_ensgather run") 
        if (masterproc) then
        !call CCPL_advance_time(atm_demo_comp_id, "atm_demo advances time for one step")
        write(*,*) "[CCPL <ATM_DOME>] Finish ensemble procedure da_demo_ensgather run"
        write(*,*) "==============================================================================="
        end if
#endif
        call CCPL_advance_time(atm_demo_comp_id, "atm_demo advances time for one step")

    end subroutine execute_atm_coupling

    subroutine finalize_atm_coupling

        implicit none

        CALL CCPL_finalize(.false., annotation="atm_demo CALL CCPL_finalize") 

    end subroutine finalize_atm_coupling
    
end module coupling_ocn_model_mod
