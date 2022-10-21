program ocn_demo
    use mpi
    use model_setting_mod, only : ocn_demo_init, ocn_demo_step_on, finalize_ocn_demo
    use spmd_init_mod, only : masterproc, mpicomm

    implicit none
    !include 'mpif.h'
    integer                    :: start_time, end_time, ierr
    real                       :: total_time, max_total_time

    
    call system_clock(start_time)

    call ocn_demo_init
    if (masterproc) then
        print *, "ocn_demo_init finished"
    end if
    call ocn_demo_step_on
    if (masterproc) then
        print *, "ocn_demo finished time integration"
    end if
    call system_clock(end_time)
    total_time = real(end_time-start_time)/10000.
    call mpi_reduce(total_time, max_total_time, 1, MPI_REAL, MPI_MAX, 0, mpicomm, ierr)
    if (masterproc) then
            write(*,*) "[CCPL <OCN_DEMO>] all total time: ", max_total_time
    end if
    call finalize_ocn_demo
    if (masterproc) then
        print *, "ocn_demo has been finalized"
    end if

end program
