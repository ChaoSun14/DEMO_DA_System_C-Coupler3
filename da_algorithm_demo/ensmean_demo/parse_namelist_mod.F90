module parse_namelist_mod_1
    character*14      :: ANAL_TIME
    integer, public :: NLATS, NLONS, NLEVS, NANALS, NVARS, decomp_type_id
contains
    subroutine parse_namelist_da
        implicit none
        namelist /da_algorithm_demo_nml/ ANAL_TIME, NLATS, NLONS, NLEVS, NANALS, NVARS, decomp_type_id 
#ifdef ATM_DA
        open(20, file="./atm_da_algorithm_demo.nml")
        read(20, nml=da_algorithm_demo_nml)
        close(20)
#endif
#ifdef OCN_DA
        open(21, file="./ocn_da_algorithm_demo.nml")
        read(21, nml=da_algorithm_demo_nml)
        close(21)
#endif
    end subroutine parse_namelist_da
end module
