module variable_mod_2
    use spmd_init_mod_2, only: masterpe, ier, mpi_comm, npes, mytask_id
    use decomp_init_mod_2, only : local_grid_index, local_grid_size
    use grid_init_mod_2, only: nlat, nlon, nlev
    use parse_namelist_mod_2, only: NANALS
    implicit none
    include 'mpif.h'
#ifdef ATM_DA
    real(kind=4), public, allocatable :: rTAUX(:,:), rTAUY(:,:), rPS(:,:)
    real(kind=4), public, allocatable :: rTAUX_l(:), rTAUY_l(:), rPS_l(:)
    real(kind=4), public, allocatable :: rTAUX_m(:,:), rTAUY_m(:,:), rPS_m(:,:)
    real(kind=4), public, allocatable :: rT(:,:,:), rT_l(:,:), rT_m(:,:,:)

#endif
#ifdef OCN_DA
    real, public, allocatable    :: rssh(:,:), rsst(:,:,:), rsalt(:,:,:)
    real, public, allocatable    :: rsst_l(:), rssh_l(:), rsalt_l(:,:)
    real, public, allocatable    :: rsst_m(:,:), rssh_m(:,:), rsalt_m(:,:,:)

#endif

contains

    subroutine variable_init_da

        call read_input_variables_da
        call scatter_fields_da

    end subroutine variable_init_da

#ifdef ATM_DA
    subroutine read_input_variables_da

        implicit none
        include "netcdf.inc"

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer :: ncid_input, ret
        integer :: rTAUXid, rTAUYid, rPSid, rTid
        integer :: i, j, k
    
        input_data_dir = ''
        input_file_name = "atm_demo.nc"
        input_file_dir_name = input_data_dir//input_file_name

        if (masterpe) then
            ret = nf_open (input_file_name, nf_nowrite, ncid_input)
            ret = nf_inq_varid (ncid_input, "TAUX", rTAUXid)
            ret = nf_inq_varid (ncid_input, "TAUY", rTAUYid)
            ret = nf_inq_varid (ncid_input, "PS", rPSid)
            ret = nf_inq_varid (ncid_input, "T", rTid)

            allocate(rTAUX(nlon, nlat))
            allocate(rTAUY(nlon, nlat))
            allocate(rPS(nlon, nlat))
            allocate(rT(nlon, nlat, nlev))

            ret = nf_get_var_real (ncid_input, rTAUXid, rTAUX)
            ret = nf_get_var_real (ncid_input, rTAUYid, rTAUY)
            ret = nf_get_var_real (ncid_input, rPSid, rPS)
            ret = nf_get_var_real (ncid_input, rTid, rT)
        else
            allocate(rTAUX(1,1),rTAUY(1,1),rPS(1,1),rT(1,1,1))
        end if
    end subroutine read_input_variables_da

    subroutine scatter_fields_da

        implicit none

        allocate(rTAUX_l(local_grid_size(mytask_id+1)))
        allocate(rTAUY_l(local_grid_size(mytask_id+1)))
        allocate(rPS_l(local_grid_size(mytask_id+1)))
        allocate(rT_l(local_grid_size(mytask_id+1), nlev))

        call scatter_field_2D_da(rTAUX, rTAUX_l)
        call scatter_field_2D_da(rTAUY, rTAUY_l)
        call scatter_field_2D_da(rPS, rPS_l)
        call scatter_field_3D_da(rT, rT_l)

        allocate(rTAUX_m(local_grid_size(mytask_id+1), NANALS))
        allocate(rTAUY_m(local_grid_size(mytask_id+1), NANALS))
        allocate(rPS_m(local_grid_size(mytask_id+1), NANALS))
        allocate(rT_m(local_grid_size(mytask_id+1), nlev, NANALS))

    end subroutine scatter_fields_da

    subroutine deallocate_fields_da
        implicit none

        if(allocated(rTAUX)) deallocate(rTAUX)
        if(allocated(rTAUY)) deallocate(rTAUY)
        if(allocated(rPS)) deallocate(rPS)
        if(allocated(rT)) deallocate(rT)
        if(allocated(rTAUX_l)) deallocate(rTAUX_l)
        if(allocated(rTAUY_l)) deallocate(rTAUY_l)
        if(allocated(rPS_l)) deallocate(rPS_l)
        if(allocated(rT_l)) deallocate(rT_l)
        if(allocated(rTAUX_m)) deallocate(rTAUX_m)
        if(allocated(rTAUY_m)) deallocate(rTAUY_m)
        if(allocated(rPS_m)) deallocate(rPS_m)
        if(allocated(rT_m)) deallocate(rT_m)

    end subroutine deallocate_fields_da
#endif

#ifdef OCN_DA
    subroutine read_input_variables_da
        
        use netcdf
        implicit none
        !include "netcdf.inc"

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer :: ncid_input, ret
        integer :: sshid, sstid, saltid
        integer :: i, j, k

        input_data_dir = ''
        input_file_name = "ocn_demo.nc"
        input_file_dir_name = input_data_dir//input_file_name

        if (masterpe) then
            ret = nf90_open (input_file_name, nf90_nowrite, ncid_input)
            ret = nf90_inq_varid (ncid_input, "z0", sshid)
            ret = nf90_inq_varid (ncid_input, "ts", sstid)
            ret = nf90_inq_varid (ncid_input, "ss", saltid)

            allocate(rsst(nlon, nlat, nlev))
            allocate(rsalt(nlon, nlat, nlev))
            allocate(rssh(nlon, nlat))

            ret = nf90_get_var (ncid_input, sshid, rssh)
            ret = nf90_get_var (ncid_input, sstid, rsst)
            ret = nf90_get_var (ncid_input, saltid, rsalt)
        else
            allocate(rssh(1,1),rsst(1,1,1),rsalt(1,1,1))
        end if
    end subroutine read_input_variables_da  
 
    subroutine scatter_fields_da

        implicit none

        allocate(rsst_l(local_grid_size(mytask_id+1)))
        allocate(rssh_l(local_grid_size(mytask_id+1)))
        allocate(rsalt_l(local_grid_size(mytask_id+1), nlev))

        !call scatter_field_2D_da(rsst(:,:,1), rsst_l)
        !call scatter_field_2D_da(rssh, rssh_l)
        !call scatter_field_3D_da(rsalt(:,:,:), rsalt_l)

        allocate(rsst_m(local_grid_size(mytask_id+1), NANALS))
        allocate(rssh_m(local_grid_size(mytask_id+1), NANALS))
        allocate(rsalt_m(local_grid_size(mytask_id+1), nlev, NANALS))

    end subroutine scatter_fields_da
    
    subroutine deallocate_fields_da
        implicit none

        if(allocated(rsst)) deallocate(rsst)
        if(allocated(rssh)) deallocate(rssh)
        if(allocated(rsalt)) deallocate(rsalt)

        if(allocated(rsst_l)) deallocate(rsst_l)
        if(allocated(rssh_l)) deallocate(rssh_l)
        if(allocated(rsalt_l)) deallocate(rsalt_l)
        
        if(allocated(rsst_m)) deallocate(rsst_m)
        if(allocated(rssh_m)) deallocate(rssh_m)
        if(allocated(rsalt_m)) deallocate(rsalt_m)
        
    end subroutine deallocate_fields_da
#endif 

   subroutine scatter_field_2D_da(global_field, local_field)

        implicit none
        real(kind=4), intent(in)  :: global_field(nlon, nlat)
        real(kind=4), intent(out) :: local_field(local_grid_size(mytask_id+1))
        !----------local variables-----------------------------------
        real(kind=4) gfield(max(local_grid_size(mytask_id+1),local_grid_size(npes)), npes)
        real(kind=4) lfield(max(local_grid_size(mytask_id+1),local_grid_size(npes)))
        integer :: p, i, j, m
        integer :: displs(1:npes)  !scatter displacements
        integer :: sndcnts(1:npes) !scatter send counts
        integer :: recvcnt(1:npes)  !scatter receive count
        
        gfield = -99999.
        lfield = -99999.
        sndcnts = local_grid_size(npes)
        displs(1) = 0
        do p=2, npes
            displs(p) = displs(p-1)+max(local_grid_size(mytask_id+1),local_grid_size(npes))
        end do
        recvcnt = local_grid_size(npes)
        if (masterpe) then
            j = 1
            do p=1,npes
                do i=1,local_grid_size(p)
                    m = local_grid_index(i,p)
                    gfield(i,p) = global_field(mod(m-1,nlon)+1,(m-1)/nlon+1)
                end do
            end do
        end if
        call mpi_scatterv(gfield, sndcnts, displs, mpi_real4, lfield, recvcnt, mpi_real4, 0, mpi_comm, ier)
        do i=1,local_grid_size(mytask_id+1)
            local_field(i) = lfield(i)
        end do
    end subroutine scatter_field_2D_da

    subroutine scatter_field_3D_da(global_field, local_field)

        implicit none
        real(kind=4), intent(in)  :: global_field(nlon, nlat, nlev)
        real(kind=4), intent(out) :: local_field(local_grid_size(mytask_id+1), nlev)
        !----------local variables-----------------------------------
        real(kind=4) gfield(max(local_grid_size(mytask_id+1),local_grid_size(npes))*nlev, npes)
        real(kind=4) lfield(max(local_grid_size(mytask_id+1),local_grid_size(npes))*nlev)
        integer :: p, i, j, k, m
        integer :: displs(1:npes)  !scatter displacements
        integer :: sndcnts(1:npes) !scatter send counts
        integer :: recvcnt(1:npes)  !scatter receive count
        
        gfield = -99999.
        lfield = -99999.
        sndcnts = local_grid_size*nlev
        displs(1) = 0
        do p=2, npes
            displs(p) = displs(p-1)+max(local_grid_size(mytask_id+1),local_grid_size(npes))*nlev
        end do
        recvcnt = local_grid_size*nlev
        if (masterpe) then
            do p=1,npes
              k = 0  
              do j=1, nlev
                do i=1,local_grid_size(p)
                    k = k+1
                    m = local_grid_index(i,p)
                    gfield(k,p) = global_field(mod(m-1,nlon)+1,(m-1)/nlon+1, j)
                end do
              end do
            end do
        end if
        call mpi_scatterv(gfield, sndcnts, displs, mpi_real4, lfield, recvcnt, mpi_real4, 0, mpi_comm, ier)
        k = 0
        do j=1,nlev
        do i=1,local_grid_size(mytask_id+1)
            k = k+1
            local_field(i,j) = lfield(k)
        end do
        end do
    end subroutine scatter_field_3D_da

end module
