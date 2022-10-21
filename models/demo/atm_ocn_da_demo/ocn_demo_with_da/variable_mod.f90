module variable_mod
    use netcdf
    !use mpi
    use spmd_init_mod, only   : masterproc, mytask_id, ier, mpicomm, npes, stdout
    use decomp_init_mod, only : local_grid_cell_index, decomp_size
    use grid_init_mod, only   : latlen, lonlen, levlen
    use parse_namelist_mod, only   : num_3d_vars

    implicit none
    include 'mpif.h'
    real, public, allocatable    :: ssh(:,:), sst(:,:,:), salt(:,:,:)
    real, public, allocatable    :: sst_l(:), ssh_l(:), salt_l(:), salt_l_1(:,:)
    real, public, allocatable    :: sstm(:), sshm(:), saltm(:), saltm_1(:,:), vars_3d(:,:,:)
    real, public, allocatable    :: TAUX(:), TAUY(:), PS(:)

contains
    subroutine variable_init
        call read_input_variables
        call scatter_fields
    end subroutine variable_init

    subroutine read_input_variables

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

        if (masterproc) then
            ret = nf90_open (input_file_name, nf90_nowrite, ncid_input)
            ret = nf90_inq_varid (ncid_input, "z0", sshid)
            ret = nf90_inq_varid (ncid_input, "ts", sstid)
            ret = nf90_inq_varid (ncid_input, "ss", saltid)

            allocate(sst(lonlen, latlen, levlen))
            allocate(salt(lonlen, latlen, levlen))
            allocate(ssh(lonlen, latlen))

            ret = nf90_get_var (ncid_input, sshid, ssh)
            ret = nf90_get_var (ncid_input, sstid, sst)
            ret = nf90_get_var (ncid_input, saltid, salt)

        else
            allocate(ssh(1,1),sst(1,1,1),salt(1,1,1))
        end if
    end subroutine read_input_variables

    subroutine scatter_fields

        implicit none

        allocate(sst_l(decomp_size(mytask_id+1)))
        allocate(ssh_l(decomp_size(mytask_id+1)))
        allocate(salt_l(decomp_size(mytask_id+1)))
        allocate(salt_l_1(decomp_size(mytask_id+1), levlen))

        call scatter_field_2D(sst(:,:,1), sst_l)
        call scatter_field_2D(ssh, ssh_l)
        call scatter_field_2D(salt(:,:,1), salt_l)
        call scatter_field_3D(salt, salt_l_1)
        
        if(allocated(sst)) deallocate(sst)
        if(allocated(ssh)) deallocate(ssh)
        if(allocated(salt)) deallocate(salt)

        allocate(sstm(decomp_size(mytask_id+1)))
        allocate(sshm(decomp_size(mytask_id+1)))
        allocate(saltm(decomp_size(mytask_id+1)))
        allocate(saltm_1(decomp_size(mytask_id+1), levlen))
        allocate(vars_3d(decomp_size(mytask_id+1), levlen, num_3d_vars))

        allocate(TAUX(decomp_size(mytask_id+1)))
        allocate(TAUY(decomp_size(mytask_id+1)))
        allocate(PS(decomp_size(mytask_id+1)))
    
    end subroutine scatter_fields

    subroutine deallocate_fields

        implicit none

        if(allocated(TAUX)) deallocate(TAUX)
        if(allocated(TAUY)) deallocate(TAUY)
        if(allocated(PS)) deallocate(PS)
        if(allocated(sst_l)) deallocate(sst_l)
        if(allocated(ssh_l)) deallocate(ssh_l)
        if(allocated(salt_l)) deallocate(salt_l)
        if(allocated(salt_l_1)) deallocate(salt_l_1)
        if(allocated(sstm)) deallocate(sstm)
        if(allocated(sshm)) deallocate(sshm)
        if(allocated(saltm)) deallocate(saltm)
        if(allocated(saltm_1)) deallocate(saltm_1)
        if(allocated(vars_3d)) deallocate(vars_3d)

    end subroutine deallocate_fields

    subroutine scatter_field_2D(global_field, local_field)

        implicit none
        real(kind=4), intent(in)  :: global_field(lonlen, latlen)
        real(kind=4), intent(out) :: local_field(decomp_size(mytask_id+1))
        !----------local variables-----------------------------------
        real(kind=4) gfield(decomp_size(npes), npes)
        real(kind=4) lfield(decomp_size(npes))
        integer :: p, i, j, m
        integer :: displs(1:npes)  !scatter displacements
        integer :: sndcnts(1:npes) !scatter send counts
        integer :: recvcnt(1:npes)  !scatter receive count
        
        gfield = -99999.
        lfield = -99999.
        sndcnts = decomp_size(npes)
        displs(1) = 0
        do p=2, npes
            displs(p) = displs(p-1)+decomp_size(npes)
        end do
        recvcnt = decomp_size(npes)
        if (masterproc) then
            j = 1
            do p=1,npes
                do i=1,decomp_size(p)
                    m = local_grid_cell_index(i,p)
                    gfield(i,p) = global_field(mod(m-1,lonlen)+1, (m-1)/lonlen+1)
                end do
            end do
        end if
        call mpi_scatterv(gfield, sndcnts, displs, mpi_real4, lfield, recvcnt, mpi_real4, 0, mpicomm, ier)
        do i=1,decomp_size(mytask_id+1)
            local_field(i) = lfield(i)
        end do
    end subroutine scatter_field_2D

    subroutine scatter_field_3D(global_field, local_field)

        implicit none
        real(kind=4), intent(in)  :: global_field(lonlen, latlen, levlen)
        real(kind=4), intent(out) :: local_field(decomp_size(mytask_id+1), levlen)
        !----------local variables-----------------------------------
        real(kind=4) gfield(max(decomp_size(mytask_id+1),decomp_size(npes))*levlen, npes)
        real(kind=4) lfield(max(decomp_size(mytask_id+1),decomp_size(npes))*levlen)
        integer :: p, i, j, k, m
        integer :: displs(1:npes)  !scatter displacements
        integer :: sndcnts(1:npes) !scatter send counts
        integer :: recvcnt(1:npes)  !scatter receive count
        
        gfield = -99999.
        lfield = -99999.
        sndcnts = decomp_size*levlen
        displs(1) = 0
        do p=2, npes
            displs(p) = displs(p-1)+max(decomp_size(mytask_id+1),decomp_size(npes))*levlen
        end do
        recvcnt = decomp_size*levlen
        if (masterproc) then
            do p=1,npes
              k = 0  
              do j=1, levlen
                do i=1,decomp_size(p)
                    k = k+1
                    m = local_grid_cell_index(i,p)
                    gfield(k,p) = global_field(mod(m-1,lonlen)+1,(m-1)/lonlen+1, j)
                end do
              end do
            end do
        end if
        call mpi_scatterv(gfield, sndcnts, displs, mpi_real4, lfield, recvcnt, mpi_real4, 0, mpicomm, ier)
        k = 0
        do j=1,levlen
        do i=1,decomp_size(mytask_id+1)
            k = k+1
            local_field(i,j) = lfield(k)
        end do
        end do
    end subroutine scatter_field_3D


end module
