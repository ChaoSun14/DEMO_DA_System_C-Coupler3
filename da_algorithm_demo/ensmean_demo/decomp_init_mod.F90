module decomp_init_mod_1
    integer, public, allocatable :: local_grid_size(:)
    integer, public, allocatable :: local_grid_index(:,:)
contains
    subroutine decomp_init_da
        use spmd_init_mod_1, only: npes, mytask_id, ier
        use grid_init_mod_1, only: nlat, nlon
        use parse_namelist_mod_1, only : decomp_type_id
        implicit none
include 'mpif.h'
        integer :: i, j, sum
        allocate(local_grid_size(npes))

        local_grid_size = nlat*nlon/npes
        do i = 1, npes-1
            local_grid_size(i) = nlat*nlon/npes
        end do 
        local_grid_size(npes) = nlat*nlon-(nlat*nlon/npes)*(npes-1)

        allocate(local_grid_index(max(local_grid_size(mytask_id+1),local_grid_size(npes)), npes))

        local_grid_index = -999

        if (decomp_type_id == 1) then
            do j = 1, npes
            do i = 1, local_grid_size(j)
                local_grid_index(i,j) = j+(i-1)*npes
            end do
            end do
        else
            do j = 1, npes
            do i = 1, local_grid_size(j)
                if (j.eq.1) then
                    local_grid_index(i,j) = i
                else
                    local_grid_index(i,j) = i+(j-1)*local_grid_size(j-1)
                end if
            end do
            end do
        end if
    end subroutine decomp_init_da

end module decomp_init_mod_1
