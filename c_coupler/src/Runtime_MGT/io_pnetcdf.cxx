/**********************************************M*****************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Miss Xinzhu Yu and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifdef USE_PARALLEL_IO


#include "io_pnetcdf.h"
#include "execution_report.h"
#include "global_data.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


IO_pnetcdf::IO_pnetcdf(int ncfile_id)
{
    this->io_with_time_info = false;
    strcpy(this->object_name, "NULL");
    strcpy(this->file_type, FILE_TYPE_NETCDF);
    strcpy(this->file_name, "NULL");
    strcpy(this->open_format, "NULL");
    this->is_external_file = true;
    this->ncfile_id = ncfile_id;
}


IO_pnetcdf::IO_pnetcdf(int host_comp_id, int pio_proc_num, int io_proc_mark, MPI_Comm comm, const char *object_name, const char *file_name, const char *format, bool io_with_time_info)
{
    this->io_with_time_info = io_with_time_info;
    this->pio_proc_num = pio_proc_num;
    strcpy(this->object_name, object_name);
    strcpy(this->file_type, FILE_TYPE_NETCDF);
    strcpy(this->file_name, file_name);
    strcpy(this->open_format, format);
    this->comm = comm;
    this->is_external_file = false;
    this->current_proc_local_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "in parallel IO");
    this->io_proc_mark = io_proc_mark;

    if (words_are_the_same(format, "r"))
        rcode = ncmpi_open(comm, file_name, NC_NOWRITE, MPI_INFO_NULL, &ncfile_id);
    else if (words_are_the_same(format, "w")) {
        rcode = ncmpi_create(comm, file_name, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile_id);
        report_nc_error();
        ncmpi_enddef(ncfile_id);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "the format of openning netcdf file must be read or write (\"r\" or \"w\")\n");
    report_nc_error();

    rcode = ncmpi_close(ncfile_id);
    report_nc_error();
}


IO_pnetcdf::~IO_pnetcdf()
{
}


void IO_pnetcdf::report_nc_error()
{
    EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", ncmpi_strerror(rcode), file_name);
}

void IO_pnetcdf::write_grid(Remap_grid_class *associated_grid, bool write_grid_name, bool use_script_format) {

}

void IO_pnetcdf::write_field_data(Remap_grid_data_class *field_data, Remap_grid_class *interchange_grid, bool is_grid_data, const char *grid_field_type, int dim_ncid_num_vertex, bool write_grid_name, bool use_script_format) {

}

bool IO_pnetcdf::read_data(Remap_data_field *read_data_field, int time_pos, bool check_existence) {

}

void IO_pnetcdf::write_grided_data(Remap_grid_data_class *grided_data, bool write_grid_name, int date, int datesec, bool is_restart_field) {

}

void IO_pnetcdf::datatype_from_netcdf_to_application(nc_type nc_data_type, char *application_data_type, const char *field_name)
{
    switch (nc_data_type) {
        case NC_BYTE:
            strcpy(application_data_type, DATA_TYPE_BOOL);
            break;
        case NC_CHAR:
            strcpy(application_data_type, DATA_TYPE_CHAR);
            break;
        case NC_INT:
            strcpy(application_data_type, DATA_TYPE_INT);
            break;
        case NC_SHORT:
            strcpy(application_data_type, DATA_TYPE_SHORT);
            break;
        case NC_FLOAT:
            strcpy(application_data_type, DATA_TYPE_FLOAT);
            break;
        case NC_DOUBLE:
            strcpy(application_data_type, DATA_TYPE_DOUBLE);
            break;
        default:
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "nc data type is not supported when reading field \"%s\" from netcdf file\n", field_name);
            break;
    }
}


void IO_pnetcdf::datatype_from_netcdf_to_application(const char *application_data_type, nc_type *nc_data_type)
{
    if (words_are_the_same(application_data_type, DATA_TYPE_BOOL))
        *nc_data_type = NC_INT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_CHAR))
        *nc_data_type = NC_CHAR;
    else if (words_are_the_same(application_data_type, DATA_TYPE_INT))
        *nc_data_type = NC_INT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_SHORT))
        *nc_data_type = NC_SHORT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_FLOAT))
        *nc_data_type = NC_FLOAT;
    else if (words_are_the_same(application_data_type, DATA_TYPE_DOUBLE))
        *nc_data_type = NC_DOUBLE;
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in datatype_from_netcdf_to_application \"%s\"\n", application_data_type);
}


bool IO_pnetcdf::read_field_data(int host_comp_id, Field_mem_info *io_field_instance, int time_pos, bool check_existence)
{
    int i, num_attributes, num_dimensions, variable_id, dimension_ids[256];
    int num_sized_sub_grids, total_required_data_size, req, st;
    Remap_grid_class *sized_sub_grids[256];
    char variable_name[256];
    nc_type nc_data_type;
    unsigned long data_size;
    long long dimension_size;
    Remap_field_attribute field_attribute;
    MPI_Offset starts[256], counts[256], field_attribute_size;

    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(host_comp_id);
    int local_proc_id = comp_node->get_current_proc_local_id();
    Remap_grid_data_class *grided_data = io_field_instance->get_field_data();
    Remap_data_field *read_data_field = grided_data->get_grid_data_field();
    grided_data->get_coord_value_grid()->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);

    rcode = ncmpi_open(comm, file_name, NC_NOWRITE, MPI_INFO_NULL, &ncfile_id);
    report_nc_error();

    rcode = ncmpi_inq_varid(ncfile_id, read_data_field->field_name_in_IO_file, &variable_id);
    if (!check_existence && rcode == NC_ENOTVAR) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Does not find the field \"%s\" in the data file \"%s\"", read_data_field->field_name_in_IO_file, file_name);
        rcode = ncmpi_close(ncfile_id);
        report_nc_error();
        return false;
    }
    report_nc_error();
    rcode = ncmpi_inq_var(ncfile_id, variable_id, variable_name, &nc_data_type, &num_dimensions, dimension_ids, &num_attributes);
    report_nc_error();

    EXECUTION_REPORT(REPORT_ERROR, -1, num_dimensions == num_sized_sub_grids, "the number of dimensions in netcdf file (which is %d) is different from the dimensions defined by the variable (\"%s\") in model (which is %d).", num_dimensions, io_field_instance->get_field_name(), num_sized_sub_grids);//1. need modified

    int local_required_data_size = read_data_field->required_data_size;//2. total size
    MPI_Allreduce(&local_required_data_size, &total_required_data_size, 1, MPI_INT, MPI_SUM, comm);

    for (i = 0, data_size = 1; i < num_dimensions; i ++) {
        rcode = ncmpi_inq_dimlen(ncfile_id, dimension_ids[i], &dimension_size);
        report_nc_error();

        if (time_pos != -1 && i == 0) {
            EXECUTION_REPORT(REPORT_ERROR, -1, time_pos >= 0 && time_pos < dimension_size, "C-Coupler error in IO_pnetcdf::read_data");
            starts[0] = time_pos;
            counts[0] = 1;
        }
        else {
            if (sized_sub_grids[num_sized_sub_grids-1-i]->get_is_sphere_grid()) {
                int local_h2d_size = sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size();
                int total_h2d_size;
                MPI_Allreduce(&local_h2d_size, &total_h2d_size, 1, MPI_INT, MPI_SUM, comm);
                EXECUTION_REPORT(REPORT_ERROR, -1, dimension_size == total_h2d_size, "the dimension_size in netcdf file (which is %d) does not match the gird size defined by the variable (\"%s\") in model (which is %d)", dimension_size, io_field_instance->get_field_name(), sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size());
                const int *grid_local_cell_global_index = decomps_info_mgr->get_decomp_info(io_field_instance->get_decomp_id())->get_local_cell_global_indx();
                starts[i] = grid_local_cell_global_index[0];
                counts[i] = sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size();
            }
            else {
                EXECUTION_REPORT(REPORT_ERROR, -1, dimension_size == sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size(), "the dimension_size in netcdf file (which is %d) does not match the gird size defined by the variable (\"%s\") in model (which is %d)", dimension_size, io_field_instance->get_field_name(), sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size());
                if (local_proc_id == 0) {
                    starts[i] = 0;
                    counts[i] = sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size();
                }
                else {
                    starts[i] = 0;
                    counts[i] = 0;
                }
            }
            data_size *= sized_sub_grids[num_sized_sub_grids-1-i]->get_grid_size();
        }
    }

    read_data_field->read_data_size = data_size;
    if (read_data_field->data_buf != NULL)
        EXECUTION_REPORT(REPORT_ERROR, -1, read_data_field->required_data_size == read_data_field->read_data_size, 
                     "the data size of field \"%s\" in netcdf file (which is %d) is different from of the data size of field \"%s\" determined by grid (which is %d)\n",
                     read_data_field->field_name_in_IO_file, read_data_field->required_data_size, 
                     read_data_field->field_name_in_application, read_data_field->read_data_size);
    else {
        read_data_field->required_data_size = data_size;
        read_data_field->data_buf = new char [data_size*get_data_type_size(read_data_field->data_type_in_application)];
    }

    datatype_from_netcdf_to_application(nc_data_type, read_data_field->data_type_in_IO_file, read_data_field->field_name_in_IO_file);

    
    for (i = 0; i < num_attributes; i ++) {
        rcode = ncmpi_inq_attname(ncfile_id, variable_id, i, field_attribute.attribute_name);
        report_nc_error();
        rcode = ncmpi_inq_att(ncfile_id, variable_id, field_attribute.attribute_name, &nc_data_type, &field_attribute_size);
        report_nc_error();
        field_attribute.attribute_size = (unsigned long) field_attribute_size;
        datatype_from_netcdf_to_application(nc_data_type, field_attribute.attribute_type, field_attribute.attribute_name);
        EXECUTION_REPORT(REPORT_ERROR, -1, get_data_type_size(field_attribute.attribute_type)*field_attribute.attribute_size <= 8192, 
                     "value of attribute \"%s\" is out-of-bound (8192 bytes) when reading info of field \"%s\" from netcdf file\n", 
                     field_attribute.attribute_name, read_data_field->field_name_in_IO_file);
        switch (nc_data_type) {
            case NC_BYTE:
            case NC_CHAR:
                rcode = ncmpi_get_att_text(ncfile_id, variable_id, field_attribute.attribute_name, field_attribute.attribute_value);
                field_attribute.attribute_value[field_attribute.attribute_size] = '\0';
                break;
            case NC_SHORT:
                rcode = ncmpi_get_att_short(ncfile_id, variable_id, field_attribute.attribute_name, (short*)field_attribute.attribute_value);
                break;
            case NC_INT:
                rcode = ncmpi_get_att_int(ncfile_id, variable_id, field_attribute.attribute_name, (int*)field_attribute.attribute_value);
                break;
            case NC_FLOAT:
                rcode = ncmpi_get_att_float(ncfile_id, variable_id, field_attribute.attribute_name, (float*)field_attribute.attribute_value);
                break;
            case NC_DOUBLE:
                rcode = ncmpi_get_att_double(ncfile_id, variable_id, field_attribute.attribute_name, (double*)field_attribute.attribute_value);
                break;
            default:
                EXECUTION_REPORT(REPORT_ERROR, -1, false, "Netcdf file may be corrupted: netcdf data type is unknown\n");
        }
        report_nc_error();
        read_data_field->field_attributes.push_back(field_attribute);
    }

    
    if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_CHAR)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_CHAR), 
                     "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) must be the same (char)\n",
                     read_data_field->field_name_in_application, read_data_field->field_name_in_IO_file);
        rcode = ncmpi_iget_vara_uchar(ncfile_id, variable_id, starts, counts, (unsigned char *) read_data_field->data_buf, &req);
    }
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_BOOL)) {
        char *temp_buffer = new char[read_data_field->required_data_size*8];
        if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE)) {
            rcode = ncmpi_iget_vara_double(ncfile_id, variable_id, starts, counts, (double*) temp_buffer, &req);            
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((double *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_FLOAT)) {
            rcode = ncmpi_iget_vara_float(ncfile_id, variable_id, starts, counts, (float*) temp_buffer, &req);
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((float *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_INT)) {
            rcode = ncmpi_iget_vara_int(ncfile_id, variable_id, starts, counts, (int*) temp_buffer, &req);
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((int *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_SHORT)) {
            rcode = ncmpi_iget_vara_short(ncfile_id, variable_id, starts, counts, (short*) temp_buffer, &req);
            for (i = 0; i < read_data_field->required_data_size; i ++)
                ((bool *) read_data_field->data_buf)[i] = ((short *) temp_buffer)[i] != 0;
        }
        else if (words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_CHAR))
            rcode = ncmpi_iget_vara_uchar(ncfile_id, variable_id, starts, counts, (unsigned char *) read_data_field->data_buf, &req);
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, 
                          "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) in netcdf is \"%s\", which is not support\n",
                          read_data_field->field_name_in_application, 
                          read_data_field->field_name_in_IO_file,
                          read_data_field->data_type_in_IO_file);

        delete [] temp_buffer;
    }
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_FLOAT)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_FLOAT) || 
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_LONG) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_INT) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_SHORT), 
                     "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) in netcdf file must be float, double, long, int or short\n",
                     read_data_field->field_name_in_application, 
                     read_data_field->field_name_in_IO_file);
        rcode = ncmpi_iget_vara_float(ncfile_id, variable_id, starts, counts, (float*) read_data_field->data_buf, &req);
    }
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_INT)) 
        rcode = ncmpi_iget_vara_int(ncfile_id, variable_id, starts, counts, (int *) read_data_field->data_buf, &req);
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_SHORT)) 
        rcode = ncmpi_iget_vara_short(ncfile_id, variable_id, starts, counts, (short *) read_data_field->data_buf, &req);
    else if (words_are_the_same(read_data_field->data_type_in_application, DATA_TYPE_DOUBLE)) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_FLOAT) || 
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_DOUBLE) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_LONG) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_INT) ||
                     words_are_the_same(read_data_field->data_type_in_IO_file, DATA_TYPE_SHORT), 
                     "the data type of field (\"%s\" in program and  \"%s\" in netcdf file) in netcdf file must be float, double, long, int or short\n",
                     read_data_field->field_name_in_application, 
                     read_data_field->field_name_in_IO_file);
        rcode = ncmpi_iget_vara_double(ncfile_id, variable_id, starts, counts, (double*) read_data_field->data_buf, &req);
    }
    report_nc_error();

    ncmpi_wait_all(ncfile_id, 1, &req, &st);
    report_nc_error();
    rcode = ncmpi_close(ncfile_id);
    report_nc_error();

    return true;
}


void IO_pnetcdf::write_grid(int host_comp_id, Field_mem_info *output_field_instance, Remap_grid_class *associated_grid, bool write_grid_name, bool use_script_format)
{
    int num_sized_sub_grids, num_leaf_grids, num_masked_sub_grids, num_sphere_leaf_grids, i, dim_ncid, num_h2d_sized_sub_grids;
    Remap_grid_class *sized_sub_grids[256], *leaf_grids[256], *masked_sub_grids[256], *h2d_sized_sub_grids[256];
    Remap_grid_data_class *grid_data_field;
    char tmp_string[256];


    if (associated_grid == NULL)
        return;

    rcode = ncmpi_open(comm, file_name, NC_WRITE, MPI_INFO_NULL, &ncfile_id);
    report_nc_error();

    rcode = ncmpi_redef(ncfile_id);
    report_nc_error();
    associated_grid->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
    for (i = 0; i < num_sized_sub_grids; i ++)
        if (sized_grids_map.find(sized_sub_grids[i]) == sized_grids_map.end()) {
            if (write_grid_name) {
                if (sized_sub_grids[i]->get_num_dimensions() == 1)
                    sprintf(tmp_string, "dim_%s_%d", sized_sub_grids[i]->get_coord_label(), recorded_grids.size());
                else if (sized_sub_grids[i]->get_is_sphere_grid())
                    sprintf(tmp_string, "dim_H2D_%d", recorded_grids.size());
                else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error in IO_pnetcdf::write_grid");
            }
            else if (sized_sub_grids[i]->get_num_dimensions() == 1)
                sprintf(tmp_string, "%s", sized_sub_grids[i]->get_coord_label());
            else sprintf(tmp_string, "grid_size_%s", sized_sub_grids[i]->get_grid_name());
            if (sized_sub_grids[i]->get_is_sphere_grid())
                rcode = ncmpi_def_dim(ncfile_id, tmp_string, decomps_info_mgr->get_decomp_info(output_field_instance->get_decomp_id())->get_num_global_cells(), &dim_ncid);
            else rcode = ncmpi_def_dim(ncfile_id, tmp_string, sized_sub_grids[i]->get_grid_size(), &dim_ncid);
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "define dim %s for grid \"%s\" (%lx) in ncfile %s", tmp_string, sized_sub_grids[i]->get_grid_name(), sized_sub_grids[i], file_name);
            report_nc_error();
            sized_grids_map[sized_sub_grids[i]] = dim_ncid;
            recorded_grids.push_back(sized_sub_grids[i]);
        }
    if (use_script_format) {
        if (sized_grids_map.find(associated_grid) == sized_grids_map.end()) {
            rcode = ncmpi_def_dim(ncfile_id, "grid_size", associated_grid->get_grid_size(), &dim_ncid);
            report_nc_error();
            sized_grids_map[associated_grid] = dim_ncid;
            recorded_grids.push_back(associated_grid);
        }
        rcode = ncmpi_def_dim(ncfile_id, "grid_rank", num_sized_sub_grids, &dim_ncid);
        report_nc_error();
        int grid_dim_size[256], dims_ncid;
        for (int i = 0; i < num_sized_sub_grids; i ++)
            grid_dim_size[i] = sized_sub_grids[i]->get_grid_size();
        rcode = ncmpi_def_var(ncfile_id, "grid_dims", NC_INT, 1, &dim_ncid, &dims_ncid);
        report_nc_error();
        rcode = ncmpi_enddef(ncfile_id);
        report_nc_error();
        ncmpi_put_var_int(ncfile_id, dims_ncid, grid_dim_size);
        report_nc_error();
        rcode = ncmpi_redef(ncfile_id);
        report_nc_error();
    }
    ncmpi_enddef(ncfile_id);
    report_nc_error();

    associated_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, associated_grid);
    for (i = 0; i < num_leaf_grids; i ++) {
        EXECUTION_REPORT(REPORT_LOG,host_comp_id, true, "leaf_grids size:%d", leaf_grids[i]->get_grid_size());
        grid_data_field = leaf_grids[i]->get_grid_center_field();
        if (grid_data_field != NULL) 
            write_field_data(host_comp_id, output_field_instance, grid_data_field, associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format, !words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LEV));
        /*{
            if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LEV))
                write_field_data(host_comp_id, output_field_instance, grid_data_field, associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format, false);
            else write_field_data(host_comp_id, output_field_instance, grid_data_field, associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format, true);
        }*/
        if (leaf_grids[i]->get_sigma_grid_sigma_value_field() != NULL) {
            EXECUTION_REPORT(REPORT_ERROR, -1, grid_data_field == NULL, "Software error in IO_pnetcdf::write_grid");
            write_field_data(host_comp_id, output_field_instance, leaf_grids[i]->get_sigma_grid_sigma_value_field(), associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format, false);
            if (leaf_grids[i]->get_hybrid_grid_coefficient_field() != NULL)
                write_field_data(host_comp_id, output_field_instance, leaf_grids[i]->get_hybrid_grid_coefficient_field(), associated_grid, true, GRID_CENTER_LABEL, -1, write_grid_name, use_script_format, false);
            sprintf(tmp_string, "P0");
            if (write_grid_name)
                sprintf(tmp_string, "grid_%d_P0", get_recorded_grid_num(leaf_grids[i]));
            double P0 = leaf_grids[i]->get_sigma_grid_top_value();
            rcode = ncmpi_redef(ncfile_id);
            report_nc_error();
            rcode = ncmpi_put_att_double(ncfile_id, NC_GLOBAL, tmp_string, NC_DOUBLE, 1, &P0);
            report_nc_error();
            ncmpi_enddef(ncfile_id);
            report_nc_error();
        }
        grid_data_field = leaf_grids[i]->get_grid_vertex_field();
        if (grid_data_field != NULL && !grid_data_field->get_coord_value_grid()->get_are_vertex_values_set_in_default()) {
            if (grid_data_field->get_coord_value_grid()->get_num_dimensions() == 1)
               sprintf(tmp_string, "num_vertexes_%s", grid_data_field->get_coord_value_grid()->get_coord_label());
            else if (use_script_format)
                sprintf(tmp_string, "grid_corners");
            else sprintf(tmp_string, "num_vertexes_H2D");
            rcode = ncmpi_inq_dimid(ncfile_id, tmp_string, &dim_ncid);
            if (rcode == NC_EBADDIM) {
                rcode = ncmpi_redef(ncfile_id);
                report_nc_error();
                rcode = ncmpi_def_dim(ncfile_id, tmp_string, grid_data_field->get_coord_value_grid()->get_num_vertexes(), &dim_ncid);
                ncmpi_enddef(ncfile_id);
                report_nc_error();
            }
            write_field_data(host_comp_id, output_field_instance, grid_data_field, associated_grid, true, GRID_VERTEX_LABEL, dim_ncid, write_grid_name, use_script_format, false);
        }
    }

    associated_grid->get_masked_sub_grids(&num_masked_sub_grids, masked_sub_grids);
    for (i = 0; i < num_masked_sub_grids; i ++)
        write_field_data(host_comp_id, output_field_instance, masked_sub_grids[i]->get_grid_mask_field(), associated_grid, true, GRID_MASK_LABEL, -1, write_grid_name, use_script_format, true);
    if (associated_grid->get_grid_imported_area() != NULL)
        write_field_data(host_comp_id, output_field_instance, associated_grid->get_grid_imported_area(), associated_grid, true, "area", -1, write_grid_name, use_script_format, true);

    rcode = ncmpi_close(ncfile_id);
    report_nc_error();
}


void IO_pnetcdf::write_field_data(int host_comp_id, Field_mem_info *output_field_instance, Remap_grid_data_class *field_data, 
                                Remap_grid_class *interchange_grid,
                                bool is_grid_data, 
                                const char *grid_field_type, 
                                int dim_ncid_num_vertex,
                                bool write_grid_name,
                                bool use_script_format,
                                bool h2d_related)
{
    int num_sized_sub_grids, num_dims, i, num_total_procs, local_required_data_size, total_required_data_size=-1;
    unsigned long io_data_size;
    MPI_Offset dimension_size[256];
    MPI_Offset starts[256], counts[256];
    Remap_grid_class *sized_sub_grids[256];
    char tmp_string[256];
    int var_ncid, dim_ncids[256], req;
    nc_type nc_data_type;
    int *temp_buffer, st;
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(host_comp_id);
    int local_proc_id = comp_node->get_current_proc_local_id();

    if (!is_grid_data)
        field_data->set_masked_cell_to_missing_value();

    tmp_string[0] = '\0';
    if (is_grid_data && write_grid_name) {
        if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application, "mask"))
            sprintf(tmp_string, "grid_%d_%s", get_recorded_grid_num(field_data->get_coord_value_grid()), field_data->get_grid_data_field()->field_name_in_application);
        else {
            if (words_are_the_same(grid_field_type, GRID_VERTEX_LABEL))
                sprintf(tmp_string, "grid_%d_%s_%s", get_recorded_grid_num(field_data->get_coord_value_grid()), grid_field_type, field_data->get_grid_data_field()->field_name_in_application);
            else sprintf(tmp_string, "grid_%d_%s", get_recorded_grid_num(field_data->get_coord_value_grid()), field_data->get_grid_data_field()->field_name_in_application);
        }
    }
    else {
        if (!words_are_the_same(field_data->get_grid_data_field()->field_name_in_IO_file, "\0") && !is_grid_data)
            sprintf(tmp_string, "%s", field_data->get_grid_data_field()->field_name_in_IO_file);
        else {
            if (words_are_the_same(grid_field_type, GRID_VERTEX_LABEL))
                if (!use_script_format)
                    sprintf(tmp_string , "%s_%s", grid_field_type, field_data->get_grid_data_field()->field_name_in_application);
                else if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application,COORD_LABEL_LON))
                    sprintf(tmp_string , SCRIP_VERTEX_LON_LABEL);
                else sprintf(tmp_string , SCRIP_VERTEX_LAT_LABEL);
            else {
                sprintf(tmp_string , "%s", field_data->get_grid_data_field()->field_name_in_application);
                if (use_script_format) {
                    if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application, COORD_LABEL_LON))
                        sprintf(tmp_string , SCRIP_CENTER_LON_LABEL);
                    else if (words_are_the_same(field_data->get_grid_data_field()->field_name_in_application, COORD_LABEL_LAT))
                        sprintf(tmp_string , SCRIP_CENTER_LAT_LABEL);
                    else if (words_are_the_same(grid_field_type, GRID_MASK_LABEL))
                        sprintf(tmp_string , SCRIP_MASK_LABEL);
                }
            }
        }
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "IO field name is %s", tmp_string);
    }

    rcode = ncmpi_inq_varid(ncfile_id, tmp_string, &var_ncid);
    if (rcode != NC_ENOTVAR) {
        if (is_grid_data)
            return;
        else EXECUTION_REPORT(REPORT_WARNING, -1, io_with_time_info,
                            "field data \"%s\" has been written to netcdf file \"%s\" before. The old data will be overwritten\n",
                            field_data->get_grid_data_field()->field_name_in_application, file_name);
    }

    if (interchange_grid != NULL) {
        field_data->interchange_grid_data(interchange_grid);
        if (sized_grids_map.find(field_data->get_coord_value_grid()) != sized_grids_map.end()) {
            num_sized_sub_grids = 1;
            sized_sub_grids[0] = field_data->get_coord_value_grid();
        }
        else field_data->get_coord_value_grid()->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
        for (i = 0; i < num_sized_sub_grids; i ++)
            dim_ncids[i] = sized_grids_map[sized_sub_grids[num_sized_sub_grids-1-i]];

		if (report_error_enabled) {
	        for (i = 0, io_data_size = 1; i < num_sized_sub_grids; i ++) {
	            rcode = ncmpi_inq_dimlen(ncfile_id, dim_ncids[i], &dimension_size[i]);
	            report_nc_error();
	            io_data_size *= dimension_size[i];
	        }
	        if (io_data_size == 1) io_data_size = 0;
	        local_required_data_size = field_data->get_grid_data_field()->required_data_size;
	        MPI_Allreduce(&local_required_data_size, &total_required_data_size, 1, MPI_INT, MPI_SUM, comm);
	        if (h2d_related)
	            EXECUTION_REPORT(REPORT_ERROR, -1, total_required_data_size == io_data_size, "C-Coupler error: the data size in field %s for writing and IO file must be the same: %ld : %ld", field_data->get_grid_data_field()->field_name_in_application, total_required_data_size, io_data_size);
	        else
	            EXECUTION_REPORT(REPORT_ERROR, -1, field_data->get_grid_data_field()->required_data_size == io_data_size, "C-Coupler error: the data size in field for writing and IO file must be the same: %ld : %ld", field_data->get_grid_data_field()->required_data_size, io_data_size);
		}

        for (i = 0; i < num_sized_sub_grids; i ++) {
            EXECUTION_REPORT(REPORT_ERROR, -1, sized_grids_map.find(sized_sub_grids[i]) != sized_grids_map.end(), "remap software error1 in write_field_data\n");
            if (h2d_related) {
                if (sized_sub_grids[i]->get_is_sphere_grid()) {
					const int *grid_local_cell_global_index = decomps_info_mgr->get_decomp_info(output_field_instance->get_decomp_id())->get_local_cell_global_indx();
					int start = grid_local_cell_global_index[0];
                    starts[num_sized_sub_grids-1-i] = start;
                }
                else starts[num_sized_sub_grids-1-i] = 0;
				counts[num_sized_sub_grids-1-i] = sized_sub_grids[i]->get_grid_size();
            }
            else {
                if (local_proc_id == 0) {
                    starts[num_sized_sub_grids-1-i] = 0;
                    counts[num_sized_sub_grids-1-i] = sized_sub_grids[i]->get_grid_size();
                }
                else {
                    starts[num_sized_sub_grids-1-i] = 0;
                    counts[num_sized_sub_grids-1-i] = 0;
                }
            }
        }
        num_dims = num_sized_sub_grids;
    }
    else num_dims = 0;

    if (is_grid_data) {
        if (dim_ncid_num_vertex != -1) {
            starts[num_dims] = 0;
            counts[num_dims] = field_data->get_coord_value_grid()->get_num_vertexes();
            dim_ncids[num_dims++] = dim_ncid_num_vertex;
        }
    }

    if (!is_grid_data && io_with_time_info) {
        for (i = num_dims; i > 0; i --) {
            dim_ncids[i] = dim_ncids[i-1];
            starts[i] = starts[i-1];
            counts[i] = counts[i-1];
        }
        num_dims ++;
        dim_ncids[0] = time_dim_id;
        starts[0] = time_count - 1;
        counts[0] = 1;
    }

    rcode = ncmpi_inq_varid(ncfile_id, tmp_string, &var_ncid);
    if (rcode == NC_ENOTVAR) {
        rcode = ncmpi_redef(ncfile_id);
        report_nc_error();
        datatype_from_netcdf_to_application(field_data->get_grid_data_field()->data_type_in_IO_file, &nc_data_type);
        rcode = ncmpi_def_var(ncfile_id, tmp_string, nc_data_type, num_dims, dim_ncids, &var_ncid);
        report_nc_error();
        for (i = 0; i < field_data->get_grid_data_field()->field_attributes.size(); i ++) {
            datatype_from_netcdf_to_application(field_data->get_grid_data_field()->field_attributes[i].attribute_type, &nc_data_type);
            switch (nc_data_type) {
                case NC_BYTE:
                case NC_CHAR:
                    rcode = ncmpi_put_att_text(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name,
                                            field_data->get_grid_data_field()->field_attributes[i].attribute_size, field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_SHORT:
                    rcode = ncmpi_put_att_short(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (short*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_INT:
                    rcode = ncmpi_put_att_int(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (int*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_FLOAT:
                    rcode = ncmpi_put_att_float(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (float*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                case NC_DOUBLE:
                    rcode = ncmpi_put_att_double(ncfile_id, var_ncid, field_data->get_grid_data_field()->field_attributes[i].attribute_name, nc_data_type,
                                             field_data->get_grid_data_field()->field_attributes[i].attribute_size, (double*)field_data->get_grid_data_field()->field_attributes[i].attribute_value);
                    break;
                default:
                    EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error2 in write_field_data\n");
            }
            report_nc_error();
        }
        rcode = ncmpi_enddef(ncfile_id);
        report_nc_error();
    }

    if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL)) {
        temp_buffer = new int [field_data->get_grid_data_field()->required_data_size];
        for (long i = 0; i < field_data->get_grid_data_field()->required_data_size; i ++)
            if (((bool*)field_data->get_grid_data_field()->data_buf)[i])
                temp_buffer[i] = 1;
            else temp_buffer[i] = 0;
        rcode = ncmpi_iput_vara_int(ncfile_id, var_ncid, starts, counts, temp_buffer, &req);
    }
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_CHAR))
        rcode = ncmpi_iput_vara_schar(ncfile_id, var_ncid, starts, counts, (signed char *) field_data->get_grid_data_field()->data_buf, &req);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT))
        rcode = ncmpi_iput_vara_float(ncfile_id, var_ncid, starts, counts, (float*) field_data->get_grid_data_field()->data_buf, &req);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_INT))
        rcode = ncmpi_iput_vara_int(ncfile_id, var_ncid, starts, counts, (int *) field_data->get_grid_data_field()->data_buf, &req);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_SHORT))
        rcode = ncmpi_iput_vara_short(ncfile_id, var_ncid, starts, counts, (short *) field_data->get_grid_data_field()->data_buf, &req);
    else if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE)) {
        rcode = ncmpi_iput_vara_double(ncfile_id, var_ncid, starts, counts, (double*) field_data->get_grid_data_field()->data_buf, &req);
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "remap software error3 in write_field_data\n");
    report_nc_error();
    ncmpi_wait_all(ncfile_id, 1, &req, &st);
    report_nc_error();
    if (words_are_the_same(field_data->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL))
        delete [] temp_buffer;
}


void IO_pnetcdf::write_grided_data(int host_comp_id, Field_mem_info *output_field_instance, bool write_grid_name, int date, int datesec, bool is_restart_field)
{
    MPI_Offset starts, counts, dim_len;
    int current_date, current_datesec;
    int time_var_id, date_var_id, datesec_var_id, tmp_var_id;;
    Remap_grid_data_class *tmp_field_data_for_io;
    Remap_grid_data_class *grided_data = output_field_instance->get_field_data();

    if (execution_phase_number == 0)
        return;
    
    rcode = ncmpi_open(comm, file_name, NC_WRITE, MPI_INFO_NULL, &ncfile_id);
    report_nc_error();

    if (!io_with_time_info)
        EXECUTION_REPORT(REPORT_ERROR, -1, date == -1 && datesec == -1, "remap software error in write_grided_data \n");
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, date > 0 && datesec >= 0, "remap software error in write_grided_data \n");
        rcode = ncmpi_inq_dimid(ncfile_id, "time", &time_dim_id); 
        if (rcode == NC_EBADDIM) {
            rcode = ncmpi_redef(ncfile_id);
            report_nc_error();
            rcode = ncmpi_def_dim(ncfile_id, "time", NC_UNLIMITED, &time_dim_id);
            report_nc_error();
            rcode = ncmpi_def_var(ncfile_id, "time", NC_INT, 1, &time_dim_id, &time_var_id);
            report_nc_error();
            rcode = ncmpi_def_var(ncfile_id, "date", NC_INT, 1, &time_dim_id, &date_var_id);
            report_nc_error();
            rcode = ncmpi_def_var(ncfile_id, "datesec", NC_INT, 1, &time_dim_id, &datesec_var_id);
            report_nc_error();
            ncmpi_enddef(ncfile_id);
            report_nc_error();
            time_count = 0;
            current_date = -1;
            current_datesec = -1;
        }
        else {
            rcode = ncmpi_inq_dimlen(ncfile_id, time_dim_id, &dim_len);
            time_count = dim_len;
            report_nc_error();
            starts = time_count - 1;
            counts = 1;
            rcode = ncmpi_inq_varid(ncfile_id, "date", &date_var_id);
            report_nc_error();
            rcode = ncmpi_get_vara_int_all(ncfile_id, date_var_id, &starts, &counts, &current_date);
            report_nc_error();
            rcode = ncmpi_inq_varid(ncfile_id, "datesec", &datesec_var_id);
            report_nc_error();
            rcode = ncmpi_get_vara_int_all(ncfile_id, datesec_var_id, &starts, &counts, &current_datesec);
            report_nc_error();
        }
        if (current_date != date || current_datesec != datesec) {
            time_count ++;
            starts = time_count - 1;
            counts = 1;
            rcode = ncmpi_inq_varid(ncfile_id, "time", &time_var_id);  
            report_nc_error();
            EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "starts:%d, time_count:%d,counts:%d, time_var_id:%x,nc_file_id:%x", starts, time_count, counts, time_var_id,ncfile_id);
            rcode = ncmpi_put_vara_int_all(ncfile_id, time_var_id, &starts, &counts, &time_count);
            report_nc_error();
            EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "finished writing time2");
            rcode = ncmpi_inq_varid(ncfile_id, "date", &date_var_id);  
            report_nc_error();
            rcode = ncmpi_put_vara_int_all(ncfile_id, date_var_id, &starts, &counts, &date);   
            report_nc_error();
            rcode = ncmpi_inq_varid(ncfile_id, "datesec", &datesec_var_id);  
            report_nc_error();
            rcode = ncmpi_put_vara_int_all(ncfile_id, datesec_var_id, &starts, &counts, &datesec);
            report_nc_error();
        }
    }

    rcode = ncmpi_close(ncfile_id);
    report_nc_error();

    write_grid(host_comp_id, output_field_instance, grided_data->get_coord_value_grid(), write_grid_name, false);

    rcode = ncmpi_open(comm, file_name, NC_WRITE, MPI_INFO_NULL, &ncfile_id);
    report_nc_error();

    if (strlen(grided_data->get_grid_data_field()->data_type_in_IO_file) == 0)
        strcpy(grided_data->get_grid_data_field()->data_type_in_IO_file, grided_data->get_grid_data_field()->data_type_in_application);
    tmp_field_data_for_io = generate_field_data_for_IO(grided_data, is_restart_field);
    write_field_data(host_comp_id, output_field_instance, tmp_field_data_for_io, grided_data->get_coord_value_grid(), false, "", -1, write_grid_name, false, true);
    if (tmp_field_data_for_io != grided_data)
        delete tmp_field_data_for_io;

    rcode = ncmpi_close(ncfile_id);
    report_nc_error();
}


long IO_pnetcdf::get_dimension_size(const char *dim_name, MPI_Comm comm, bool is_root_proc)
{
    int dimension_id;
    long dimension_size = -1;


    if (is_root_proc) {
        rcode = ncmpi_open(comm, file_name, NC_NOWRITE, MPI_INFO_NULL, &ncfile_id);
        report_nc_error();
        rcode = ncmpi_inq_dimid(ncfile_id, dim_name, &dimension_id);
        if (rcode == NC_NOERR) {
            rcode = ncmpi_inq_dimlen(ncfile_id, dimension_id, (MPI_Offset *)(&dimension_size));
            report_nc_error();   
        }
        rcode = ncmpi_close(ncfile_id);
        report_nc_error();
    }

    if (comm != MPI_COMM_NULL)
        MPI_Bcast(&dimension_size, 1, MPI_LONG, 0, comm);

    return dimension_size;
}

void IO_pnetcdf::put_global_attr(const char *text_title, const void *attr_value, const char *local_data_type, const char *nc_data_type, int size)
{
    int nc_datatype;

    
    if (!is_external_file) {
        rcode = ncmpi_open(comm, file_name, NC_WRITE, MPI_INFO_NULL, &ncfile_id);
        report_nc_error();
    }
    rcode = ncmpi_redef(ncfile_id);
    report_nc_error();

    if (words_are_the_same(nc_data_type, DATA_TYPE_STRING))
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(local_data_type, DATA_TYPE_STRING), "software error in IO_pnetcdf::put_global_attr: miss match of data type");
    else if (words_are_the_same(nc_data_type, DATA_TYPE_FLOAT))
        nc_datatype = NC_FLOAT;
    else if (words_are_the_same(nc_data_type, DATA_TYPE_DOUBLE))
        nc_datatype = NC_DOUBLE;
    else if (words_are_the_same(nc_data_type, DATA_TYPE_INT))
        nc_datatype = NC_INT;
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "software error in IO_pnetcdf::put_global_attr: wrong nc data type %s", nc_datatype);
    
    if (words_are_the_same(local_data_type, DATA_TYPE_STRING))
        rcode = ncmpi_put_att_text(ncfile_id, NC_GLOBAL, text_title, strlen((const char*)attr_value), (const char*)attr_value);
    else if (words_are_the_same(local_data_type, DATA_TYPE_FLOAT))
        rcode = ncmpi_put_att_float(ncfile_id, NC_GLOBAL, text_title, NC_FLOAT, size, (const float*)attr_value);
    else if (words_are_the_same(local_data_type, DATA_TYPE_DOUBLE))
        rcode = ncmpi_put_att_double(ncfile_id, NC_GLOBAL, text_title, NC_DOUBLE, size, (const double*)attr_value);
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "software error in IO_pnetcdf::put_global_attr: wrong local data type %s", local_data_type);
    report_nc_error();
    ncmpi_enddef(ncfile_id);
    report_nc_error();
    if (!is_external_file) {
        rcode = ncmpi_close(ncfile_id);
        report_nc_error();
    }
}


bool IO_pnetcdf::get_field_datatype(const char *field_name, char *var_data_type)
{
    int variable_id = NC_GLOBAL;
    nc_type nc_data_type;

    rcode = ncmpi_open(comm, file_name, NC_NOWRITE, MPI_INFO_NULL, &ncfile_id);
    report_nc_error();

    if (field_name != NULL) {
        rcode = ncmpi_inq_varid(ncfile_id, field_name, &variable_id);
        if (rcode != NC_NOERR) {
            rcode = nc_close(ncfile_id);
            report_nc_error();
            return false;
        }
    }
    rcode = ncmpi_inq_vartype(ncfile_id, variable_id, &nc_data_type);
    if (rcode != NC_NOERR) {
        rcode = ncmpi_close(ncfile_id);
        report_nc_error();
        return false;
    }
    datatype_from_netcdf_to_application(nc_data_type, var_data_type, field_name);
    rcode = ncmpi_close(ncfile_id);
    report_nc_error();
}


bool IO_pnetcdf::get_file_field_attribute(const char *field_name, const char *attribute_name, char *attribute_value, char *data_type, bool neccessity)
{
    int variable_id = NC_GLOBAL;
    nc_type nc_data_type;
    MPI_Offset attribute_size;

    rcode = ncmpi_open(comm, file_name, NC_NOWRITE, MPI_INFO_NULL, &ncfile_id);
    report_nc_error();

    if (field_name != NULL) {
        rcode = ncmpi_inq_varid(ncfile_id, field_name, &variable_id);
        if (rcode != NC_NOERR) {
            rcode = ncmpi_close(ncfile_id);
            report_nc_error();
            return false;
        }
    }
    rcode = ncmpi_inq_att(ncfile_id, variable_id, attribute_name, &nc_data_type, &attribute_size);
    if (rcode != NC_NOERR) {
        rcode = ncmpi_close(ncfile_id);
        report_nc_error();
        return false;
    }
    switch (nc_data_type) {
        case NC_BYTE:
        case NC_CHAR:
            strcpy(data_type, DATA_TYPE_STRING);
            rcode = ncmpi_get_att_text(ncfile_id, variable_id, attribute_name, attribute_value);
            attribute_value[attribute_size] = '\0';
            break;
        case NC_SHORT:
            rcode = ncmpi_get_att_short(ncfile_id, variable_id, attribute_name, (short*)attribute_value);
            strcpy(data_type, DATA_TYPE_SHORT);
            break;
        case NC_INT:
            strcpy(data_type, DATA_TYPE_INT);
            rcode = ncmpi_get_att_int(ncfile_id, variable_id, attribute_name, (int*)attribute_value);
            break;
        case NC_FLOAT:
            strcpy(data_type, DATA_TYPE_FLOAT);
            rcode = ncmpi_get_att_float(ncfile_id, variable_id, attribute_name, (float*)attribute_value);
            break;
        case NC_DOUBLE:
            strcpy(data_type, DATA_TYPE_DOUBLE);
            rcode = ncmpi_get_att_double(ncfile_id, variable_id, attribute_name, (double*)attribute_value);
            break;
        default:
            rcode = ncmpi_close(ncfile_id);
            report_nc_error();
            return false;
    }
    if (rcode != NC_NOERR) {
        rcode = ncmpi_close(ncfile_id);
        if (neccessity) report_nc_error();
        return false;
    }

    rcode = ncmpi_close(ncfile_id);
    report_nc_error();

    return true;
}

#endif

