program da_algorithm_demo
    use model_setting_mod
    use spmd_init_mod, only : masterpe

    implicit none
include 'mpif.h'
    call da_algorithm_demo_init
    if (masterpe) then
        print *, "da_algorithm_demo_init finished"
    end if
    call da_algorithm_demo_step_on
    if (masterpe) then
        print *, "da_algorithm_demo finished analysis"
    end if
    call finalize_da_algorithm_demo
    if (masterpe) then
        print *, "da_algorithm_demo has been finalized"
    end if
end program
