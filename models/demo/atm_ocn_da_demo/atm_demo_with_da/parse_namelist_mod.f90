module parse_namelist_mod
    integer, public :: time_length, time_step, coupling_freq
contains
    subroutine parse_namelist
        implicit none
        namelist /atm_demo_nml/ time_length, time_step, coupling_freq
        open(10, file="./atm_demo.nml")
        read(10, nml=atm_demo_nml)
        close(10)
    end subroutine parse_namelist
end module
