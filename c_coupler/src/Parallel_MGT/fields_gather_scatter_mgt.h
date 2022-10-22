/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef COMP_GATHER_FIELD
#define COMP_GATHER_FIELD


#include "common_utils.h"
#include "memory_mgt.h"
#include "io_netcdf.h"
#include <vector>
#include "datamodel_mgt.h"


class Gather_scatter_rearrange_info
{
    private: 
        char fully_decomp_name[NAME_STR_SIZE];
        char data_type[NAME_STR_SIZE];
        char *io_data_type;
        bool has_IO_field;
        int num_local_procs;
        int num_levels;
        int num_points_in_each_cell;
        int num_total_cells;
        int *counts;
        int *displs;
        void *mpibuf;
        int *rearrange_indexes;

        int host_comp_id;
        int original_decomp_id;
        int new_decomp_id;
        int grid_id;
        int io_grid_id;
        Field_mem_info *IO_field_mem;
        int current_proc_local_id;
        MPI_Comm local_comm;
        Field_mem_info *mirror_local_field_mem;
        Output_handler *output_handler_for_rearrange;
        Input_handler_operator *input_handler_for_rearrange;
        int gather_scatter_type;

    public:
        Gather_scatter_rearrange_info(Field_mem_info*, Field_mem_info*, int);
        ~Gather_scatter_rearrange_info();
        bool match(int, int, int, const char*, int, const char*, int);
        Field_mem_info *gather_field(Field_mem_info*, void*, const char*);
        Field_mem_info *scatter_field(Field_mem_info*, Input_file_time_info*, const char*);
        void scatter_field(Field_mem_info*, bool &);
        Field_mem_info *get_IO_field_mem(Field_mem_info*);
        template <class T> void rearrange_gather_data(T*, T*, int);
        template <class T> void rearrange_scatter_data(T*, T*, int);
        Field_mem_info *get_mirror_local_field_mem() { return mirror_local_field_mem; }
};


class Fields_gather_scatter_mgt
{
    private: 
        std::vector<Gather_scatter_rearrange_info*> gather_scatter_rearrange_infos;
        Gather_scatter_rearrange_info *apply_gather_scatter_rearrange_info(Field_mem_info*, Field_mem_info*, int);

    public:
        Field_mem_info *gather_field(Field_mem_info*, Field_mem_info*, void*, const char*);
        Field_mem_info *scatter_field(Field_mem_info*, Field_mem_info*, Input_file_time_info*, const char*);
        Fields_gather_scatter_mgt() {}
        ~Fields_gather_scatter_mgt();
        void gather_write_field(IO_netcdf*, Field_mem_info*, bool, int, int, bool);
        bool read_scatter_field(IO_netcdf*, Field_mem_info*, const char *, int, bool);
};


#endif