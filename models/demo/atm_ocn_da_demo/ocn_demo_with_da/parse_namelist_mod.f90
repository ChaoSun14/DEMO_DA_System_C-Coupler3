module parse_namelist_mod
    integer, public :: time_length, time_step, coupling_freq, num_3d_vars
contains
    subroutine parse_namelist
        implicit none
        namelist /ocn_demo_nml/ time_length, time_step, coupling_freq, num_3d_vars
        open(10, file="./ocn_demo.nml")
        read(10, nml=ocn_demo_nml)
        close(10)
    end subroutine parse_namelist
end module
