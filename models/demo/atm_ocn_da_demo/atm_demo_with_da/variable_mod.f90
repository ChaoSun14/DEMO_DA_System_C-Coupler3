module variable_mod
    !use mpi
    use spmd_init_mod, only   : masterproc, ier, mpicomm, npes, mytask_id, stdout
    use decomp_init_mod, only : local_grid_cell_index, decomp_size
    use grid_init_mod, only   : latlen, lonlen, levlen

    implicit none
    include 'mpif.h'
    real(kind=4), public, allocatable :: TAUX(:,:), TAUY(:,:), PS(:,:)
    real(kind=4), public, allocatable :: TAUX_l(:), TAUY_l(:), PS_l(:)
    !real(kind=4), public, allocatable :: TAUX_m(:), TAUY_m(:), PS_m(:)
    real(kind=4), public, allocatable :: T(:,:,:), T_l(:,:)
    real(kind=4), public, allocatable :: sst(:), ssh(:), salt(:)
    real(kind=4), public, pointer     :: TAUX_m(:)=>null(), TAUY_m(:)=>null(), PS_m(:)=>null(), T_m(:,:)=>null()

contains
    subroutine variable_init

        call read_input_variables
        call scatter_fields

    end subroutine variable_init

    subroutine read_input_variables

        implicit none
        include "netcdf.inc"

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer :: ncid_input, ret
        integer :: TAUXid, TAUYid, PSid, Tid
        integer :: i, j, k
    
        input_data_dir = ''
        input_file_name = "atm_demo.nc"
        input_file_dir_name = input_data_dir//input_file_name

    
        if (masterproc) then
            ret = nf_open (input_file_name, nf_nowrite, ncid_input)
            ret = nf_inq_varid (ncid_input, "TAUX", TAUXid)
            ret = nf_inq_varid (ncid_input, "TAUY", TAUYid)
            ret = nf_inq_varid (ncid_input, "PS", PSid)
            ret = nf_inq_varid (ncid_input, "T", Tid)

            allocate(TAUX(lonlen, latlen))
            allocate(TAUY(lonlen, latlen))
            allocate(PS(lonlen, latlen))
            allocate(T(lonlen, latlen, levlen))

            ret = nf_get_var_real (ncid_input, TAUXid, TAUX)
            ret = nf_get_var_real (ncid_input, TAUYid, TAUY)
            ret = nf_get_var_real (ncid_input, PSid, PS)
            ret = nf_get_var_real (ncid_input, Tid, T)

        else
            allocate(TAUX(1,1),TAUY(1,1),PS(1,1),T(1,1,1))
        end if
    end subroutine read_input_variables

    subroutine scatter_field_2D(global_field, local_field)

        implicit none
        real(kind=4), intent(in)  :: global_field(lonlen, latlen)
        real(kind=4), intent(out) :: local_field(decomp_size(mytask_id+1))
        !----------local variables-----------------------------------
        real(kind=4) gfield(max(decomp_size(mytask_id+1),decomp_size(npes)), npes)
        real(kind=4) lfield(max(decomp_size(mytask_id+1),decomp_size(npes)))
        integer :: p, i, j, m
        integer :: displs(1:npes)  !scatter displacements
        integer :: sndcnts(1:npes) !scatter send counts
        integer :: recvcnt(1:npes)  !scatter receive count
        
        gfield = -99999.
        lfield = -99999.
        sndcnts = decomp_size
        displs(1) = 0
        do p=2, npes
            displs(p) = displs(p-1)+max(decomp_size(mytask_id+1),decomp_size(npes))
        end do
        recvcnt = decomp_size
        if (masterproc) then
            j = 1
            do p=1,npes
                do i=1,decomp_size(p)
                    m = local_grid_cell_index(i,p)
                    gfield(i,p) = global_field(mod(m-1,lonlen)+1,(m-1)/lonlen+1)
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
    
    subroutine scatter_fields

        implicit none

        allocate(TAUX_l(decomp_size(mytask_id+1)))
        allocate(TAUY_l(decomp_size(mytask_id+1)))
        allocate(PS_l(decomp_size(mytask_id+1)))
        allocate(T_l(decomp_size(mytask_id+1), levlen))

        call scatter_field_2D(TAUX, TAUX_l)
        call scatter_field_2D(TAUY, TAUY_l)
        call scatter_field_2D(PS, PS_l)
        call scatter_field_3D(T, T_l)

        allocate(TAUX_m(decomp_size(mytask_id+1)))
        allocate(TAUY_m(decomp_size(mytask_id+1)))
        allocate(PS_m(decomp_size(mytask_id+1)))
        allocate(T_m(decomp_size(mytask_id+1), levlen))

        allocate(sst(decomp_size(mytask_id+1)))
        allocate(ssh(decomp_size(mytask_id+1)))
        allocate(salt(decomp_size(mytask_id+1)))

    end subroutine scatter_fields

    subroutine deallocate_fields
        implicit none

        if(allocated(TAUX)) deallocate(TAUX)
        if(allocated(TAUY)) deallocate(TAUY)
        if(allocated(PS)) deallocate(PS)
        if(allocated(T)) deallocate(T)
        if(allocated(TAUX_l)) deallocate(TAUX_l)
        if(allocated(TAUY_l)) deallocate(TAUY_l)
        if(allocated(PS_l)) deallocate(PS_l)
        if(allocated(T_l)) deallocate(T_l)
       ! if(allocated(TAUX_m)) deallocate(TAUX_m)
       ! if(allocated(TAUY_m)) deallocate(TAUY_m)
       ! if(allocated(PS_m)) deallocate(PS_m)
       ! if(allocated(T_m)) deallocate(T_m)
        if(allocated(sst)) deallocate(sst)
        if(allocated(ssh)) deallocate(ssh)
        if(allocated(salt)) deallocate(salt)

    end subroutine deallocate_fields

end module

