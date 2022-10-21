module spmd_init_mod_1
    integer, public :: mytask_id, ier, mpi_comm, npes
    logical, public :: masterpe
contains
    subroutine spmd_init_da(external_mpi_comm)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: external_mpi_comm
        mpi_comm = external_mpi_comm
        if (mpi_comm .eq. MPI_COMM_WORLD) then
            call mpi_init(ier)
        end if
        call mpi_comm_rank(mpi_comm, mytask_id, ier)
        call mpi_comm_size(mpi_comm, npes, ier)
        if (mytask_id == 0) then
            masterpe = .true.
        else
            masterpe = .false.
        end if
    end subroutine spmd_init_da
end module
