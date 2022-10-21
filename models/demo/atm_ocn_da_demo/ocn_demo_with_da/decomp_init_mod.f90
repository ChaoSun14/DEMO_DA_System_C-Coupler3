module decomp_init_mod
    integer, public, allocatable :: decomp_size(:)
    integer, public, allocatable :: local_grid_cell_index(:,:)
contains
    subroutine decomp_init
        use spmd_init_mod, only: npes, mytask_id, ier
        !use mpi
        use grid_init_mod, only: latlen, lonlen
        implicit none

        integer :: i, j, sum
        allocate(decomp_size(npes))

        decomp_size = latlen*lonlen/npes
        do i = 1, npes-1
            decomp_size(i) = latlen*lonlen/npes
        end do 
        decomp_size(npes) = latlen*lonlen-(latlen*lonlen/npes)*(npes-1)

        allocate(local_grid_cell_index(max(decomp_size(mytask_id+1),decomp_size(npes)), npes))

        local_grid_cell_index = -999

       
        do j = 1, npes
        do i = 1, decomp_size(j)
            if (j.eq.1) then
                local_grid_cell_index(i,j) = i
            else
                local_grid_cell_index(i,j) = i+(j-1)*decomp_size(j-1)
            end if
        end do
        end do

    end subroutine decomp_init

end module decomp_init_mod
