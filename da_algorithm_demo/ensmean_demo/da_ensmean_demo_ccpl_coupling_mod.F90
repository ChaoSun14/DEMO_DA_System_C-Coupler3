module da_ensmean_demo_ccpl_coupling_mod
    
    use CCPL_interface_mod
    use mpi
    use parse_namelist_mod_1, only: NANALS, NVARS

#ifdef ATM_DA

    implicit none
    integer, public  :: npe, mype, ierror
    integer, public  :: da_ensmean_demo_instance_id, da_ensmean_demo_ccpl_comm, da_stdout
    integer, public  :: da_ensmean_demo_v1d_grid_id_input, da_ensmean_demo_h2d_grid_id_input, da_ensmean_demo_m3d_grid_id_input, da_ensmean_demo_h2d_grid_id_mask_input
    integer, public  :: da_ensmean_demo_global_decomp_id_input, da_ensmean_demo_decomp_id_input, da_ensmean_demo_global_decomp_id, da_ensmean_demo_decomp_id
    integer, public  :: field_id_date, field_id_lev, field_id_lat, field_id_lon, field_id_TAUX, field_id_TAUY, field_id_t
    integer, public  :: ccpl_nlev, ccpl_nlat, ccpl_nlon, ccpl_npe_atm, ccpl_local_size, ccpl_ens_num, ccpl_ens_id
    real(kind=8), public, pointer :: ccpl_date=>null()
    real(kind=4), public, pointer :: ccpl_lev(:)=>null(), ccpl_lat(:)=>null(), ccpl_lon(:)=>null(), ccpl_TAUX(:)=>null(), ccpl_TAUY(:)=>null(), ccpl_t(:,:)=>null()
    integer, allocatable :: local_cells_index_all(:)
    real(kind=4), public, allocatable :: ccpl_grid_lat(:), ccpl_grid_lon(:)
  ! 
contains
    
    subroutine get_grids_info

        implicit none
         
        integer              :: i,j,k,ii,jj,kk

        call mpi_comm_size(da_ensmean_demo_ccpl_comm, npe, ierror)
        call mpi_comm_rank(da_ensmean_demo_ccpl_comm, mype, ierror)
 
        ccpl_nlev       = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 1, annotation="get the number of levels")
        ccpl_nlat       = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 2, annotation="get the number of latitudes")
        ccpl_nlon       = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 3, annotation="get the number of longitudes")
        ccpl_npe_atm    = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 4, annotation="get the number process of model")
        ccpl_local_size = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 5, annotation="get the number of local_size")

        ccpl_ens_num = CCPL_external_procedures_para_get_ensemble_size(da_ensmean_demo_instance_id, annotation="get ensemble number")
        ccpL_ens_id  = CCPL_external_procedures_para_get_ensemble_member_index(da_ensmean_demo_instance_id, annotation="get ensemble id")
        IF (mype .EQ. 0) THEN
        write(*,*) "+++++++++++++++++++++++++++"
        write(*,*) "[CCPL <atm_da_ensmean_demo>]  nlat, nlon, nlev: ",ccpl_nlat, ccpl_nlon, ccpl_nlev
        write(*,*) "[CCPL <atm_da_ensmean_demo>]  ccpl_npe_atm, ccpl_ens_num, ccpl_ens_id, ccpl_local_size: ", ccpl_npe_atm, ccpl_ens_num, ccpl_ens_id, ccpl_local_size
        END IF

    end subroutine get_grids_info


    subroutine register_grids_decomps_fields

        use decomp_init_mod_1, only : local_grid_index, local_grid_size

        implicit none
         
        integer              :: i,j,k,ii,jj,kk

        da_ensmean_demo_v1d_grid_id_input = CCPL_external_procedures_para_get_field_grid_ID(da_ensmean_demo_instance_id, "LEV", annotation = "get v1d_grid_id of model")
        da_ensmean_demo_h2d_grid_id_input = CCPL_external_procedures_para_get_field_grid_ID(da_ensmean_demo_instance_id, "TAUX", annotation = "get h2d_grid_id of model")
        da_ensmean_demo_m3d_grid_id_input = CCPL_external_procedures_para_get_field_grid_ID(da_ensmean_demo_instance_id, "T", annotation = "get m3d_grid_id of model")

    !    allocate(local_cells_index_all(1:ccpl_nlon*ccpl_nlat))
    !    local_cells_index_all=CCPL_NULL_INT
    !    k=0
    !    do j=1, ccpl_nlat
    !        do i=1, ccpl_nlon
    !           k=k+1
    !           local_cells_index_all(k)=(j-1)*ccpl_nlon+i
    !        end do 
    !    end do

        da_ensmean_demo_decomp_id = CCPL_register_normal_parallel_decomp("decomp_da_ensmean_demo_grid", da_ensmean_demo_h2d_grid_id_input, local_grid_size(mype+1), local_grid_index(:,mype+1), annotation="allocate decomp for da_ensmean_demo grid")

        da_ensmean_demo_decomp_id_input        = CCPL_external_procedures_para_get_field_decomp_ID(da_ensmean_demo_instance_id, "T", annotation = "get decomp_id of model")

       
        field_id_date = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_date, "date", CCPL_PARA_TYPE_IN, -1, -1, annotation="declare date")
        IF (mype .EQ. 0) write(*,*) "date: ", ccpl_date

        field_id_lev  = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_lev, "LEV", CCPL_PARA_TYPE_IN, -1, da_ensmean_demo_v1d_grid_id_input, (/ccpl_nlev/), annotation="declare lev values")
        field_id_TAUX  = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_TAUX, "TAUX", CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_h2d_grid_id_input, (/local_grid_size(mype+1)/), annotation="declare TAUX")
        field_id_TAUY  = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_TAUY, "TAUY", CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_h2d_grid_id_input, (/local_grid_size(mype+1)/), annotation="declare TAUY")
        field_id_t    = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_t, "T", CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_m3d_grid_id_input, (/local_grid_size(mype+1), ccpl_nlev/), annotation="declare temperature")

   end subroutine register_grids_decomps_fields
#endif

#ifdef OCN_DA

    implicit none
    integer, public  :: npe, mype, ierror
    integer, public  :: da_ensmean_demo_instance_id, da_ensmean_demo_ccpl_comm, da_stdout
    integer, public  :: da_ensmean_demo_v1d_grid_id_input, da_ensmean_demo_h2d_grid_id_input, da_ensmean_demo_m3d_grid_id_input
    integer, public  :: da_ensmean_demo_global_decomp_id_input, da_ensmean_demo_decomp_id_input, da_ensmean_demo_global_decomp_id, da_ensmean_demo_decomp_id
    integer, public  :: field_id_date, field_id_lat, field_id_lon, field_id_sst, field_id_ssh, field_id_salt, field_id_lev
    integer, public  :: ccpl_nlat, ccpl_nlon, ccpl_nlev, ccpl_npe_model, ccpl_ens_num, ccpl_ens_id
    real(kind=8), public, pointer :: ccpl_date=>null()
    real(kind=4), public, pointer :: ccpl_lev(:)=>null(), ccpl_lat(:)=>null(), ccpl_lon(:)=>null(), ccpl_sst(:)=>null(), ccpl_ssh(:)=>null(), ccpl_salt(:,:)=>null()
    integer, allocatable :: local_cells_index_all(:), field_id_3d_vars(:)
    real(kind=4), public, allocatable :: ccpl_grid_lat(:), ccpl_grid_lon(:)
    character(len=3)     :: num_str
    type, private :: vars_3d
        real(kind=4), pointer :: data_3d(:,:)=>null()
    end type vars_3d
    type(vars_3d), public, allocatable :: ccpl_vars_3d(:)
   
contains
    
    subroutine get_grids_info

        implicit none
         
        integer              :: i,j,k,ii,jj,kk

        call mpi_comm_size(da_ensmean_demo_ccpl_comm, npe, ierror)
        call mpi_comm_rank(da_ensmean_demo_ccpl_comm, mype, ierror)
 
        
        ccpl_nlev       = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 1, annotation="get the number of levels")
        ccpl_nlat       = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 2, annotation="get the number of latitudes")
        ccpl_nlon       = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 3, annotation="get the number of longitudes")
        ccpl_npe_model  = CCPL_external_procedures_para_get_control_var(da_ensmean_demo_instance_id, 4, annotation="get the number process of model")
        ccpl_ens_num = CCPL_external_procedures_para_get_ensemble_size(da_ensmean_demo_instance_id, annotation="get ensemble number")
        ccpL_ens_id  = CCPL_external_procedures_para_get_ensemble_member_index(da_ensmean_demo_instance_id, annotation="get ensemble id")
        IF (mype .EQ. 0) THEN
        write(*,*) "+++++++++++++++++++++++++++"
        write(*,*) "[CCPL <ocn_da_ensmean_demo>]  nlat, nlon: ",ccpl_nlat, ccpl_nlon
        write(*,*) "[CCPL <ocn_da_ensmean_demo>]  ccpl_npe_model, ccpl_ens_num, ccpl_ens_id: ", ccpl_npe_model, ccpl_ens_num, ccpl_ens_id
        END IF

    end subroutine get_grids_info


    subroutine register_grids_decomps_fields

        use decomp_init_mod_1, only : local_grid_index, local_grid_size

        implicit none
         
        integer              :: i,j,k,ii,jj,kk

        
        da_ensmean_demo_h2d_grid_id_input = CCPL_external_procedures_para_get_field_grid_ID(da_ensmean_demo_instance_id, "sst", annotation = "get h2d_grid_id of model")
        da_ensmean_demo_v1d_grid_id_input = CCPL_external_procedures_para_get_field_grid_ID(da_ensmean_demo_instance_id, "LEV", annotation = "get v1d_grid_id of model")
        da_ensmean_demo_m3d_grid_id_input = CCPL_external_procedures_para_get_field_grid_ID(da_ensmean_demo_instance_id, "salt_3d", annotation = "get m3d_grid_id of model")

        allocate(ccpl_grid_lat(ccpl_nlon*ccpl_nlat))
        allocate(ccpl_grid_lon(ccpl_nlon*ccpl_nlat))
        call CCPL_get_H2D_grid_data(da_ensmean_demo_h2d_grid_id_input, -1, "lat", ccpl_grid_lat, annotation="get latitudes of H2D grid")
        call CCPL_get_H2D_grid_data(da_ensmean_demo_h2d_grid_id_input, -1, "lon", ccpl_grid_lon, annotation="get longitudes of H2D grid")
        IF (mype .EQ. 0) THEN
            write(*,*) "[CCPL <da_ensmean_demo>] ccpl_grid_lat size, (1:10): ", size(ccpl_grid_lat), ccpl_grid_lat(1:10)
            write(*,*) "[CCPL <da_ensmean_demo>] ccpl_grid_lon size, (1:10): ", size(ccpl_grid_lon), ccpl_grid_lon(1:10)
        END IF

    !    allocate(local_cells_index_all(1:ccpl_nlon*ccpl_nlat))
    !    local_cells_index_all=CCPL_NULL_INT
    !    k=0
    !    do j=1, ccpl_nlat
    !        do i=1, ccpl_nlon
    !           k=k+1
    !           local_cells_index_all(k)=(j-1)*ccpl_nlon+i
    !        end do 
    !    end do

        da_ensmean_demo_decomp_id = CCPL_register_normal_parallel_decomp("decomp_da_ensmean_demo_grid", da_ensmean_demo_h2d_grid_id_input, local_grid_size(mype+1), local_grid_index(:,mype+1), annotation="allocate decomp for da_ensmean_demo grid")
       
        field_id_date = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_date, "date", CCPL_PARA_TYPE_IN, -1, -1, annotation="declare date")
        IF (mype .EQ. 0) write(*,*) "date: ", ccpl_date

        field_id_lev  = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_lev, "LEV", CCPL_PARA_TYPE_IN, -1, da_ensmean_demo_v1d_grid_id_input, (/ccpl_nlev/), annotation="declare lev values")
        field_id_sst  = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_sst, "sst", CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_h2d_grid_id_input, (/local_grid_size(mype+1)/), annotation="declare sea surface temperature")
        field_id_ssh  = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_ssh, "ssh", CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_h2d_grid_id_input, (/local_grid_size(mype+1)/), annotation="declare sea surface height")
        field_id_salt = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_salt, "salt_3d", CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_m3d_grid_id_input, (/local_grid_size(mype+1), ccpl_nlev/), annotation="declare 3d salinity")
        allocate(ccpl_vars_3d(NVARS), field_id_3d_vars(NVARS))
        do i=1, NVARS
            write(num_str, "(i3.3)") i
            field_id_3d_vars(i) = CCPL_external_procedures_para_declare_field(da_ensmean_demo_instance_id, ccpl_vars_3d(i)%data_3d, "vars_3d_"//trim(adjustl(num_str)), CCPL_PARA_TYPE_IN, da_ensmean_demo_decomp_id, da_ensmean_demo_m3d_grid_id_input, (/local_grid_size(mype+1), ccpl_nlev/), annotation="declare 3d vars")
        end do

   end subroutine register_grids_decomps_fields
#endif

   SUBROUTINE check_diff_1D(model_type, fields_name, fields, ccpl_fields_in)
     
        implicit none
        integer, intent(in)                  :: model_type
        character(len=*), intent(in)         :: fields_name
        real(kind=4), INTENT(IN), DIMENSION(:)   :: fields
        real(kind=4), INTENT(IN), DIMENSION(:)   :: ccpl_fields_in
        real(kind=4), allocatable, DIMENSION(:)  :: ccpl_fields
        real(kind=4), parameter    :: diff_ref=1.0e-4
        real(kind=4), parameter    :: Large = 1.0E+20
        real(kind=4), dimension(2) :: range
        real(kind=4), dimension(2) :: f1, f2
        character(len=37)      :: model_string(2)
        logical :: diff_pass
        integer :: numx,numy
        integer :: i,j,k,ierr
   
     20 FORMAT (a6,a10,a10)
        IF (model_type .EQ. 1) THEN
            write(model_string(2),20) 'DAFCC',fields_name,' Min/Max: '
            write(model_string(1),20) 'READ ',fields_name,' Min/Max: '
        END IF
          

        numx = SIZE (ccpl_fields_in,1)
        ALLOCATE (ccpl_fields(numx))
        ccpl_fields=ccpl_fields_in     
   
        range(1) = Large
        range(2) =-Large
   
        range(1) = MINVAL(fields)
        range(2) = MAXVAL(fields)
        CALL mpi_allreduce (range(1),f1(1),1,MPI_REAL,MPI_MIN,da_ensmean_demo_ccpl_comm,ierr)
        CALL mpi_allreduce (range(2),f1(2),1,MPI_REAL,MPI_MAX,da_ensmean_demo_ccpl_comm,ierr)
   
        range(1) = MINVAL(ccpl_fields)
        range(2) = MAXVAL(ccpl_fields)
        CALL mpi_allreduce (range(1),f2(1),1,MPI_REAL,MPI_MIN,da_ensmean_demo_ccpl_comm,ierr)
        CALL mpi_allreduce (range(2),f2(2),1,MPI_REAL,MPI_MAX,da_ensmean_demo_ccpl_comm,ierr)
   
     40 FORMAT (a37,2(1pe14.6))
   
        IF (MAXVAL(ABS(fields-ccpl_fields)) .LT. diff_ref) THEN
            diff_pass = .true.
            IF (mype .EQ. 0) THEN
              write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
              write(*,*)  "DA_ensmean_DEMO Fields ", fields_name, " Pass verification!"
              write(*,40) trim(model_string(1)),f1(1),f1(2)
              write(*,40) trim(model_string(2)),f2(1),f2(2)
              write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
            END IF
        ELSE
            diff_pass = .false.
            write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
            write(*,*) "DA_ensmean_DEMO Rank ", mype, " Fields ", fields_name, " Failed verification!"
            IF (mype .EQ. 0) THEN
              write(*,*) "NUMBERS:  ", SIZE (fields,1), SIZE (ccpl_fields,1)
              write(*,*) "MAXVAL(ABS(fields-ccpl_fields))= ", MAXVAL(ABS(fields-ccpl_fields))
              write(*,40) trim(model_string(1)),f1(1),f1(2)
              write(*,40) trim(model_string(2)),f2(1),f2(2)
              DO i = 1,numx
                IF (fields(i) .NE. ccpl_fields(i)) THEN
                   write(*,*) i,fields(i),ccpl_fields(i)
                END IF
              END DO
            END IF
            write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
        END IF
        DEALLOCATE (ccpl_fields)
     
    END SUBROUTINE check_diff_1D

       SUBROUTINE check_diff_2D(model_type, fields_name, fields_in, ccpl_fields_in)
     
        implicit none
        integer, intent(in)                  :: model_type
        character(len=*), intent(in)         :: fields_name
        real(kind=4), INTENT(IN), DIMENSION(:,:)   :: fields_in
        real(kind=4), INTENT(IN), DIMENSION(:,:)   :: ccpl_fields_in
        real(kind=4), allocatable, DIMENSION(:)  :: fields, ccpl_fields
        real(kind=4), parameter    :: diff_ref=1.0e-4
        real(kind=4), parameter    :: Large = 1.0E+20
        real(kind=4), dimension(2) :: range
        real(kind=4), dimension(2) :: f1, f2
        character(len=37)      :: model_string(2)
        logical :: diff_pass
        integer :: numx,numy,numall
        integer :: i,j,k,ierr
   
     20 FORMAT (a6,a10,a10)
        IF (model_type .EQ. 1) THEN
            write(model_string(2),20) 'DAFCC',fields_name,' Min/Max: '
            write(model_string(1),20) 'READ ',fields_name,' Min/Max: '
        END IF
          
        numx = SIZE (ccpl_fields_in,1)
        numy = SIZE (ccpl_fields_in,2)
        numall = numx*numy
        ALLOCATE (fields(numall))
        ALLOCATE (ccpl_fields(numall))
        k = 0
        DO j = 1, numy
         DO i = 1, numx
           k = k+1
           fields(k)=fields_in(i,j)
           ccpl_fields(k)=ccpl_fields_in(i,j)
         END DO
        END DO
   
        range(1) = Large
        range(2) =-Large
   
        range(1) = MINVAL(fields)
        range(2) = MAXVAL(fields)
        CALL mpi_allreduce (range(1),f1(1),1,MPI_REAL,MPI_MIN,da_ensmean_demo_ccpl_comm,ierr)
        CALL mpi_allreduce (range(2),f1(2),1,MPI_REAL,MPI_MAX,da_ensmean_demo_ccpl_comm,ierr)
   
        range(1) = MINVAL(ccpl_fields)
        range(2) = MAXVAL(ccpl_fields)
        CALL mpi_allreduce (range(1),f2(1),1,MPI_REAL,MPI_MIN,da_ensmean_demo_ccpl_comm,ierr)
        CALL mpi_allreduce (range(2),f2(2),1,MPI_REAL,MPI_MAX,da_ensmean_demo_ccpl_comm,ierr)
   
     40 FORMAT (a37,2(1pe14.6))
        
        IF (MAXVAL(ABS(fields-ccpl_fields)) .LT. diff_ref) THEN
            diff_pass = .true.
            IF (mype .EQ. 0) THEN
              write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
              write(*,*)  "DA_ensmean_DEMO Fields ", fields_name, " Pass verification!"
              write(*,40) trim(model_string(1)),f1(1),f1(2)
              write(*,40) trim(model_string(2)),f2(1),f2(2)
              write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
            END IF
        ELSE
            diff_pass = .false.
            write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
            write(*,*) "DA_ensmean_DEMO Rank ", mype, " Fields ", fields_name, " Failed verification!"
            IF (mype .EQ. 0) THEN
              write(*,*) "NUMBERS:  ", SIZE (fields,1), SIZE (ccpl_fields,1)
              write(*,*) "MAXVAL(ABS(fields-ccpl_fields))= ", MAXVAL(ABS(fields-ccpl_fields)), diff_ref
              write(*,40) trim(model_string(1)),f1(1),f1(2)
              write(*,40) trim(model_string(2)),f2(1),f2(2)
              DO i = 1,numall
                IF (fields(i) .NE. ccpl_fields(i)) THEN
                   write(*,*) i,fields(i),ccpl_fields(i)
                END IF
              END DO
            END IF
            write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@"
        END IF
        DEALLOCATE (ccpl_fields)
     
    END SUBROUTINE check_diff_2D
        
end module da_ensmean_demo_ccpl_coupling_mod


