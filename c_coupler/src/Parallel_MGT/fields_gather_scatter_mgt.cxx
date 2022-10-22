/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "fields_gather_scatter_mgt.h"
#include "global_data.h"


template <class T> void Gather_scatter_rearrange_info::rearrange_gather_data(T *buf_src, T *buf_dst, int num_cells_in_fully_decomp)
{
    long i, j, k, m, num_data_in_each_level, num_local_cells, rearrange_indexes_start;
    T *tmp_buf_src, *tmp_buf_dst;


    for (m = 0; m < num_local_procs; m ++) {
        num_data_in_each_level = counts[m] / num_levels;
        num_local_cells = num_data_in_each_level / num_points_in_each_cell;
        rearrange_indexes_start = displs[m] / num_levels / num_points_in_each_cell;
        for (k = 0; k < num_levels; k ++) 
            for (i = 0; i < num_local_cells; i ++) {
                if (rearrange_indexes[rearrange_indexes_start+i] == CCPL_NULL_INT)
                    continue;
                tmp_buf_src = buf_src + displs[m] + k*num_data_in_each_level + i*num_points_in_each_cell;
                tmp_buf_dst = buf_dst + rearrange_indexes[rearrange_indexes_start+i]*num_points_in_each_cell + k*num_cells_in_fully_decomp*num_points_in_each_cell;
                for (j = 0; j < num_points_in_each_cell; j ++)
                    tmp_buf_dst[j] = tmp_buf_src[j];
            }
    }
}


template <class T> void Gather_scatter_rearrange_info::rearrange_scatter_data(T *buf_src, T *buf_dst, int num_cells_in_fully_decomp)
{
    long i, j, k, m, num_data_in_each_level, num_local_cells, rearrange_indexes_start;
    T *tmp_buf_src, *tmp_buf_dst;


    for (m = 0; m < num_local_procs; m ++) {
        num_data_in_each_level = counts[m] / num_levels;
        num_local_cells = num_data_in_each_level / num_points_in_each_cell;
        rearrange_indexes_start = displs[m] / num_levels / num_points_in_each_cell;
        for (k = 0; k < num_levels; k ++)
            for (i = 0; i < num_local_cells; i ++) {
                if (rearrange_indexes[rearrange_indexes_start+i] == CCPL_NULL_INT)
                    continue;
                tmp_buf_src = buf_src + rearrange_indexes[rearrange_indexes_start+i]*num_points_in_each_cell + k*num_cells_in_fully_decomp*num_points_in_each_cell; 
                tmp_buf_dst = buf_dst + displs[m] + k*num_data_in_each_level + i*num_points_in_each_cell;
                for (j = 0; j < num_points_in_each_cell; j ++)
                    tmp_buf_dst[j] = tmp_buf_src[j];
            }
    }  
}


Gather_scatter_rearrange_info::Gather_scatter_rearrange_info(Field_mem_info *local_field, Field_mem_info *io_field, int handler_type)
{
    const char *io_float_datatype;
    const char *io_integer_datatype;
    int io_field_id, io_handler_id;

    host_comp_id = local_field->get_host_comp_id();
    original_decomp_id = local_field->get_decomp_id();
    grid_id = local_field->get_grid_id();
    io_grid_id = io_field != NULL? io_field->get_grid_id() : -1;
    strcpy(data_type, local_field->get_field_data()->get_grid_data_field()->data_type_in_application);
    this->gather_scatter_type = handler_type;

    if (io_field != NULL) {
        io_data_type = strdup(io_field->get_field_data()->get_grid_data_field()->data_type_in_application);
        bool is_io_data_type_float_type = words_are_the_same(io_data_type, DATA_TYPE_FLOAT) || words_are_the_same(io_data_type, DATA_TYPE_DOUBLE);
        io_float_datatype = is_io_data_type_float_type? strdup(io_data_type): strdup("");
        io_integer_datatype = is_io_data_type_float_type? strdup(""): strdup(io_data_type);
    }
    else {
        io_data_type = strdup(data_type);
        io_float_datatype = strdup("");
        io_integer_datatype = strdup("");
    }

    datamodel_mgr->add_handlers_field_mem_buf_mark();
    this->mirror_local_field_mem = memory_manager->alloc_mem(local_field, BUF_MARK_IO_FIELD_MIRROR , datamodel_mgr->get_handlers_field_mem_buf_mark(), local_field->get_data_type(), false, false);
    io_field_id = this->mirror_local_field_mem->get_field_instance_id();

    if (handler_type == OUTPUT_HANDLER) {
        io_handler_id = datamodel_mgr->register_field_instances_output_handler(1, &io_field_id, "file_name", "netcdf", 1, -1, -1, -1, io_grid_id, 0, io_float_datatype, io_integer_datatype, 1, false, 2, "");
        output_handler_for_rearrange = datamodel_mgr->search_output_handler(io_handler_id);
        IO_field_mem = output_handler_for_rearrange->get_unique_IO_field_mem();
    }
    else {
        io_handler_id = datamodel_mgr->register_input_handler_operator(io_field_id, io_grid_id, io_data_type, io_field->get_unit());
        input_handler_for_rearrange = datamodel_mgr->search_input_handler_operator(io_handler_id);
        IO_field_mem = input_handler_for_rearrange->get_unique_IO_field_mem();
    }
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(IO_field_mem->get_data_type(), io_data_type), "Software error in Gather_scatter_rearrange_info::Gather_scatter_rearrange_info: IO_field_mem:%s vs io_data_type:%s, field_name:%s", IO_field_mem->get_data_type(), io_data_type, IO_field_mem->get_field_name());
    has_IO_field = true;


    counts = NULL;
    displs = NULL;
    mpibuf = NULL;
    rearrange_indexes = NULL;
}


bool Gather_scatter_rearrange_info::match(int host_comp_id, int decomp_id, int grid_id, const char *data_type, int io_grid_id, const char *io_field_data_type, int handler_type)
{
    return this->host_comp_id == host_comp_id && this->original_decomp_id == decomp_id && this->grid_id == grid_id && words_are_the_same(this->data_type, data_type) && this->io_grid_id == io_grid_id && words_are_the_same(this->io_data_type, io_field_data_type) && this->gather_scatter_type == handler_type;
}


Field_mem_info *Gather_scatter_rearrange_info::gather_field(Field_mem_info *local_field_mem, void *restart_write_data_file, const char *field_name)
{
	local_field_mem->transformation_between_chunks_array(true);
	mirror_local_field_mem->transformation_between_chunks_array(true);
	memory_manager->copy_field_data_values(mirror_local_field_mem, local_field_mem);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "Successfully memcpy a mirror_local_field_mem");
    output_handler_for_rearrange->execute_handler(true, API_ID_HANDLE_NORMAL_EXPLICIT_OUTPUT, field_name, restart_write_data_file, "");

    return IO_field_mem;
}


Field_mem_info *Gather_scatter_rearrange_info::scatter_field(Field_mem_info *local_field_mem, Input_file_time_info *next_input_file_time, const char *field_name)
{
    char annotation_hint[NAME_STR_SIZE];
    sprintf(annotation_hint, "execute input_handler_operator for scatter field:%s", field_name);
    input_handler_for_rearrange->execute_input_handler_operator(field_name, next_input_file_time, annotation_hint);
    memory_manager->copy_field_data_values(local_field_mem, mirror_local_field_mem);
    mirror_local_field_mem->transformation_between_chunks_array(true);
    local_field_mem->transformation_between_chunks_array(true);
    return IO_field_mem;
}


Field_mem_info *Gather_scatter_rearrange_info::get_IO_field_mem(Field_mem_info *local_field_mem)
{
    if (!has_IO_field)
        return local_field_mem;

    return IO_field_mem;
}


void Gather_scatter_rearrange_info::scatter_field(Field_mem_info *local_field_mem, bool &has_field_in_file)
{
    int has_field_in_file_int = has_field_in_file? 1 : 0;


    MPI_Bcast(&has_field_in_file_int, 1, MPI_INT, 0, local_comm);
    has_field_in_file = (has_field_in_file_int == 1);
    if (!has_field_in_file)
        return;
    
    if (!has_IO_field) {
        MPI_Bcast(local_field_mem->get_data_buf(), local_field_mem->get_size_of_field()*get_data_type_size(data_type), MPI_CHAR, 0, local_comm);
        return;
    }

    if (get_data_type_size(data_type) == 1) {
        if (current_proc_local_id == 0)
            rearrange_scatter_data((char*) IO_field_mem->get_data_buf(), (char*) mpibuf, decomps_info_mgr->get_decomp_info(new_decomp_id)->get_num_local_cells());
        MPI_Scatterv(mpibuf, counts, displs, MPI_CHAR, local_field_mem->get_data_buf(), local_field_mem->get_size_of_field(), MPI_CHAR, 0, local_comm);
    }
    else if (get_data_type_size(data_type) == 2) {
        if (current_proc_local_id == 0)
            rearrange_scatter_data((short*) IO_field_mem->get_data_buf(), (short*) mpibuf, decomps_info_mgr->get_decomp_info(new_decomp_id)->get_num_local_cells());
        MPI_Scatterv(mpibuf, counts, displs, MPI_SHORT, local_field_mem->get_data_buf(), local_field_mem->get_size_of_field(), MPI_SHORT, 0, local_comm);
    }
    else if (get_data_type_size(data_type) == 4) {
        if (current_proc_local_id == 0)
            rearrange_scatter_data((int*) IO_field_mem->get_data_buf(), (int*) mpibuf, decomps_info_mgr->get_decomp_info(new_decomp_id)->get_num_local_cells());
        MPI_Scatterv(mpibuf, counts, displs, MPI_INT, local_field_mem->get_data_buf(), local_field_mem->get_size_of_field(), MPI_INT, 0, local_comm);
    }
    else if (get_data_type_size(data_type) == 8) {
        if (current_proc_local_id == 0)
            rearrange_scatter_data((double*) IO_field_mem->get_data_buf(), (double*) mpibuf, decomps_info_mgr->get_decomp_info(new_decomp_id)->get_num_local_cells());
        MPI_Scatterv(mpibuf, counts, displs, MPI_DOUBLE, local_field_mem->get_data_buf(), local_field_mem->get_size_of_field(), MPI_DOUBLE, 0, local_comm);
    }
    else EXECUTION_REPORT(REPORT_ERROR, false, "C-Coupler error in Gather_scatter_rearrange_info::gather_field\n");
	
	local_field_mem->transformation_between_chunks_array(false);
}


Gather_scatter_rearrange_info::~Gather_scatter_rearrange_info()
{
    if (counts != NULL)
        delete [] counts;
    if (displs != NULL)
        delete [] displs;
    if (mpibuf != NULL)
        delete [] mpibuf;
    if (rearrange_indexes != NULL)
        delete [] rearrange_indexes;
    delete io_data_type;
}


Gather_scatter_rearrange_info *Fields_gather_scatter_mgt::apply_gather_scatter_rearrange_info(Field_mem_info *local_field, Field_mem_info *io_field, int handler_type)
{
    int i;
    Gather_scatter_rearrange_info *rearrange_info;


    for (i = 0; i < gather_scatter_rearrange_infos.size(); i ++)
        if (gather_scatter_rearrange_infos[i]->match(local_field->get_host_comp_id(), local_field->get_decomp_id(), local_field->get_grid_id(), local_field->get_field_data()->get_grid_data_field()->data_type_in_application, io_field != NULL? io_field->get_grid_id() : -1, (io_field != NULL? io_field : local_field)->get_field_data()->get_grid_data_field()->data_type_in_application, handler_type))
            break;

    if (i == gather_scatter_rearrange_infos.size()) {
        rearrange_info = new Gather_scatter_rearrange_info(local_field, io_field, handler_type);
        gather_scatter_rearrange_infos.push_back(rearrange_info);
    }
    else rearrange_info = gather_scatter_rearrange_infos[i];


    if (rearrange_info->get_IO_field_mem(local_field) != local_field && io_field == NULL)
        rearrange_info->get_IO_field_mem(local_field)->copy_in_another_field_attributes(local_field);
    else if (io_field != NULL)
        rearrange_info->get_IO_field_mem(local_field)->copy_in_another_field_attributes(io_field);


    return rearrange_info;
}


Field_mem_info *Fields_gather_scatter_mgt::gather_field(Field_mem_info *local_field, Field_mem_info *io_field, void *restart_write_data_file, const char *field_name)
{
    return apply_gather_scatter_rearrange_info(local_field, io_field, OUTPUT_HANDLER)->gather_field(local_field, restart_write_data_file, field_name);
}


Field_mem_info *Fields_gather_scatter_mgt::scatter_field(Field_mem_info *model_field, Field_mem_info *io_field, Input_file_time_info *next_input_file_time, const char *field_name)
{
    return apply_gather_scatter_rearrange_info(model_field, io_field, INPUT_HANDLER)->scatter_field(model_field, next_input_file_time, field_name);
}


void Fields_gather_scatter_mgt::gather_write_field(IO_netcdf *nc_file, Field_mem_info *local_field, bool write_grid_name, int date, int datesec, bool is_restart_field)
{
    Field_mem_info *global_field = gather_field(local_field, NULL, NULL, NULL);
    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(local_field->get_host_comp_id(), "in gather_write_field") == 0)
        nc_file->write_grided_data(global_field->get_field_data(), write_grid_name, date, datesec, is_restart_field);
}


bool Fields_gather_scatter_mgt::read_scatter_field(IO_netcdf *nc_file, Field_mem_info *local_field, const char *field_IO_name, int time_pos, bool check_existence)
{
    bool has_data_in_file;
    

    Gather_scatter_rearrange_info *rearrage_info = apply_gather_scatter_rearrange_info(local_field, NULL, -1);
    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(local_field->get_host_comp_id(), "in read_scatter_field") == 0) {
        if (field_IO_name != NULL)
            strcpy(rearrage_info->get_IO_field_mem(local_field)->get_field_data()->get_grid_data_field()->field_name_in_IO_file, field_IO_name);
        has_data_in_file = nc_file->read_data(rearrage_info->get_IO_field_mem(local_field)->get_field_data()->get_grid_data_field(), time_pos, check_existence);
        if (field_IO_name != NULL)
            strcpy(rearrage_info->get_IO_field_mem(local_field)->get_field_data()->get_grid_data_field()->field_name_in_IO_file, local_field->get_field_name());
    }
    rearrage_info->scatter_field(local_field, has_data_in_file);
    return has_data_in_file;
}


Fields_gather_scatter_mgt::~Fields_gather_scatter_mgt()
{
    for (int i = 0; i < gather_scatter_rearrange_infos.size(); i ++)
        delete gather_scatter_rearrange_infos[i];
}

