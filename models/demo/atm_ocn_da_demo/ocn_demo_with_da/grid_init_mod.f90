module grid_init_mod

    real(kind=4), public, allocatable :: lat(:), lon(:), lev(:)
    integer, public           :: latlen, lonlen, levlen

contains


    subroutine grid_init

        use spmd_init_mod, only:masterproc, mpicomm, ier
        !use mpi
        implicit none
        include "netcdf.inc"
include 'mpif.h'

        character*1024 :: input_data_dir, input_file_name
        character*1024 :: input_file_dir_name
        integer        :: ncid_input, ret
        integer        :: latid, lonid, levid
        integer        :: latdimid, londimid, levdimid
        integer        :: latndims
        integer, allocatable :: latdimids(:)

        input_data_dir = ''
        input_file_name = 'ocn_demo.nc'
        input_file_dir_name = input_data_dir//input_file_name
        if (masterproc) then
            ret = nf_open (input_file_name, nf_nowrite, ncid_input)
            
            ret = nf_inq_dimid (ncid_input, "lat", latdimid)
            ret = nf_inq_dimid (ncid_input, "lon", londimid)
            ret = nf_inq_dimid (ncid_input, "lev", levdimid)
            ret = nf_inq_dimlen (ncid_input, latdimid, latlen)
            ret = nf_inq_dimlen (ncid_input, londimid, lonlen)
            ret = nf_inq_dimlen (ncid_input, levdimid, levlen)
            
            allocate(lat(latlen), lon(lonlen), lev(levlen))
            ret = nf_inq_varid (ncid_input, "lat", latid)
            ret = nf_inq_varid (ncid_input, "lon", lonid)
            ret = nf_inq_varid (ncid_input, "lev", levid)
            ret = nf_get_var_real (ncid_input, latid, lat)
            ret = nf_get_var_real (ncid_input, lonid, lon)
            ret = nf_get_var_real (ncid_input, levid, lev)
        end if

        call mpi_bcast(latlen, 1, mpi_integer, 0, mpicomm, ier)
        call mpi_bcast(lonlen, 1, mpi_integer, 0, mpicomm, ier)
        call mpi_bcast(levlen, 1, mpi_integer, 0, mpicomm, ier)

        if (masterproc == .false.) then
            allocate(lat(latlen), lon(lonlen), lev(levlen))
        end if

        call mpi_bcast(lat, latlen, mpi_real4, 0, mpicomm, ier)
        call mpi_bcast(lon, lonlen, mpi_real4, 0, mpicomm, ier)
        call mpi_bcast(lev, levlen, mpi_real4, 0, mpicomm, ier)



    end subroutine grid_init

end module grid_init_mod
