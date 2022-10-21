module grid_init_mod_1

    real(kind=4), public, allocatable :: xlat(:), xlon(:), xlev(:)
    integer, public                   :: nlat, nlon, nlev

contains


    subroutine grid_init_all

        use spmd_init_mod_1, only: masterpe, mpi_comm, ier
        implicit none
include 'mpif.h'        
include "netcdf.inc"

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer        :: ncid_input, ret
        integer        :: xlatid, xlonid, xlevid
        integer        :: latdimid, londimid, levdimid 
input_data_dir = ''
#ifdef ATM_DA
        input_file_name = 'atm_demo.nc'
#endif
#ifdef OCN_DA
        input_file_name = 'ocn_demo.nc'
#endif     
        input_file_dir_name = input_data_dir//input_file_name

        if (masterpe) then

            ret = nf_open (input_file_name, nf_nowrite, ncid_input)        
            ret = nf_inq_dimid (ncid_input, "lat", latdimid)
            ret = nf_inq_dimid (ncid_input, "lon", londimid)
            ret = nf_inq_dimid (ncid_input, "lev", levdimid)
            ret = nf_inq_dimlen (ncid_input, latdimid, nlat)
            ret = nf_inq_dimlen (ncid_input, londimid, nlon)
            ret = nf_inq_dimlen (ncid_input, levdimid, nlev)
            
            allocate(xlat(nlat), xlon(nlon), xlev(nlev))

            ret = nf_inq_varid (ncid_input, "lat", xlatid)
            ret = nf_inq_varid (ncid_input, "lon", xlonid)
            ret = nf_inq_varid (ncid_input, "lev", xlevid)
            ret = nf_get_var_real (ncid_input, xlatid, xlat)
            ret = nf_get_var_real (ncid_input, xlonid, xlon)
            ret = nf_get_var_real (ncid_input, xlevid, xlev)
        end if

        call mpi_bcast(nlat, 1, mpi_integer, 0, mpi_comm, ier)
        call mpi_bcast(nlon, 1, mpi_integer, 0, mpi_comm, ier)
        call mpi_bcast(nlev, 1, mpi_integer, 0, mpi_comm, ier)

        if (masterpe == .false.) then
            allocate(xlat(nlat), xlon(nlon),xlev(nlev))
        end if

        call mpi_bcast(xlat, nlat, mpi_real4, 0, mpi_comm, ier)
        call mpi_bcast(xlon, nlon, mpi_real4, 0, mpi_comm, ier)
        call mpi_bcast(xlev, nlev, mpi_real4, 0, mpi_comm, ier)
    end subroutine grid_init_all

end module grid_init_mod_1
