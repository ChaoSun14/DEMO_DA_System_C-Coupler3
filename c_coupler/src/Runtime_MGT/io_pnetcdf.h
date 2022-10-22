/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/



#ifndef IO_PNETCDF
#define IO_PNETCDF

#ifdef USE_PARALLEL_IO

#include "io_basis.h"
#include <pnetcdf.h>
#include "cor_global_data.h"
#include "remap_weight_of_strategy_class.h"
#include "remap_common_utils.h"
#include "memory_mgt.h"

#define SCRIP_CENTER_LON_LABEL          "grid_center_lon"
#define SCRIP_CENTER_LAT_LABEL          "grid_center_lat"
#define SCRIP_VERTEX_LON_LABEL          "grid_corner_lon"
#define SCRIP_VERTEX_LAT_LABEL          "grid_corner_lat"
#define SCRIP_MASK_LABEL                "grid_imask"

class IO_pnetcdf: public IO_basis
{
    private:
        int ncfile_id;
        int rcode;
        bool io_with_time_info;
        int time_dim_id;
        int time_count;
        bool is_external_file;
        int pio_proc_num;
        int current_proc_local_id;
        int io_proc_mark;
        MPI_Comm comm;
        
        void write_field_data(Remap_grid_data_class*, Remap_grid_class*, bool, const char*, int, bool, bool);
        void write_field_data(int, Field_mem_info*, Remap_grid_data_class*, Remap_grid_class*, bool, const char*, int, bool, bool, bool);
        void datatype_from_netcdf_to_application(nc_type, char*, const char*);
        void datatype_from_netcdf_to_application(const char*, nc_type*);
        void report_nc_error();

    public:
        IO_pnetcdf(int);
        IO_pnetcdf(int, int, int, MPI_Comm, const char*, const char*, const char*, bool);
        ~IO_pnetcdf();
        
        bool read_data(Remap_data_field*, int, bool);
        void write_grided_data(Remap_grid_data_class*, bool, int, int, bool);
        void write_grid(Remap_grid_class*, bool, bool);
        bool get_io_with_time_info() { return io_with_time_info; }
        bool read_field_data(int, Field_mem_info*, int, bool);

        void write_grided_data(int, Field_mem_info*, bool, int, int, bool);
        long get_dimension_size(const char*, MPI_Comm, bool);
        void put_global_attr(const char*, const void*, const char *, const char *, int);
        void write_grid(int, Field_mem_info*, Remap_grid_class*, bool, bool);
        MPI_Comm get_io_comm() { return comm; }
        bool get_field_datatype(const char*, char*);
        bool get_file_field_attribute(const char*, const char*, char*, char*, bool);

};

#endif

#endif 
