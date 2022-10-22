/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu.
  *  If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "original_grid_mgt.h"
#include <netcdf.h>
#include <unistd.h>
#include <math.h>



Original_grid_info::Original_grid_info(int set_comp_id, Original_grid_info *member_original_grid, int grid_id)
{
    this->ensemble_member_grid = member_original_grid;
    this->comp_id = set_comp_id;
    this->grid_id = grid_id;
    this->grid_name = strdup(member_original_grid->get_grid_name());
    this->bottom_field_name = NULL;
    if (member_original_grid->bottom_field_name != NULL)
        this->bottom_field_name = strdup(member_original_grid->bottom_field_name);
    this->comp_full_name = strdup(comp_comm_group_mgt_mgr->search_global_node(set_comp_id)->get_comp_full_name());
    this->original_CoR_grid = member_original_grid->original_CoR_grid;
    this->V3D_sub_CoR_grid = member_original_grid->V3D_sub_CoR_grid;
    this->H2D_sub_CoR_grid = member_original_grid->H2D_sub_CoR_grid;
    this->V1D_sub_CoR_grid = member_original_grid->V1D_sub_CoR_grid;
    this->Time1D_sub_CoR_grid = member_original_grid->Time1D_sub_CoR_grid;
    this->Tracer1D_sub_CoR_grid = member_original_grid->Tracer1D_sub_CoR_grid;
    this->max_sub_grid_under_V3D = NULL;
    original_grid_mgr->push_back_grid(this);
    if (member_original_grid->max_sub_grid_under_V3D != NULL)
        if (member_original_grid->max_sub_grid_under_V3D == member_original_grid)
            this->max_sub_grid_under_V3D = this;
        else this->max_sub_grid_under_V3D = original_grid_mgr->promote_ensemble_member_grid_to_set(set_comp_id, member_original_grid->max_sub_grid_under_V3D, false);
    this->Time1D_sub_grid = NULL;
    if (member_original_grid->Time1D_sub_grid != NULL)
        if (member_original_grid->Time1D_sub_grid == member_original_grid)
            this->Time1D_sub_grid == this;
        else this->Time1D_sub_grid = original_grid_mgr->promote_ensemble_member_grid_to_set(set_comp_id, member_original_grid->Time1D_sub_grid, false);
    this->Tracer1D_sub_grid = NULL;
    if (member_original_grid->Tracer1D_sub_grid != NULL)
        if (member_original_grid->Tracer1D_sub_grid == member_original_grid)
            this->Tracer1D_sub_grid = this;
        else this->Tracer1D_sub_grid = original_grid_mgr->promote_ensemble_member_grid_to_set(set_comp_id, member_original_grid->Tracer1D_sub_grid, false);
    this->H2D_sub_grid_order = member_original_grid->H2D_sub_grid_order;
    this->V1D_sub_grid_order = member_original_grid->V1D_sub_grid_order;
    this->T1D_sub_grid_order = member_original_grid->T1D_sub_grid_order;
    this->tracer_sub_grid_order = member_original_grid->tracer_sub_grid_order;
    this->bottom_field_id = member_original_grid->bottom_field_id;
    this->bottom_field_variation_type = member_original_grid->bottom_field_id;   // 0: static; 1: dynamic; 2: external
    this->mid_point_grid = NULL;
    if (member_original_grid->mid_point_grid != NULL)
        this->mid_point_grid = original_grid_mgr->promote_ensemble_member_grid_to_set(set_comp_id, member_original_grid->mid_point_grid, false);
    this->interface_level_grid = NULL;
    if (member_original_grid->interface_level_grid != NULL)
        this->interface_level_grid = original_grid_mgr->promote_ensemble_member_grid_to_set(set_comp_id, member_original_grid->interface_level_grid, false);
    this->checksum_H2D_mask = member_original_grid->checksum_H2D_mask;
    this->checksum_H2D_center_lon = member_original_grid->checksum_H2D_center_lon;
    this->checksum_H2D_center_lat = member_original_grid->checksum_H2D_center_lat;
    this->used_in_md_grid = member_original_grid->used_in_md_grid;
    if (original_CoR_grid->get_num_dimensions() > 2)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, member_original_grid->sub_grids_id.size() > 0, "Software error in Original_grid_info::Original_grid_info");
    for (int i = 0; i < member_original_grid->sub_grids_id.size(); i ++)
        this->sub_grids_id.push_back(original_grid_mgr->promote_ensemble_member_grid_to_set(set_comp_id, original_grid_mgr->search_grid_info(member_original_grid->sub_grids_id[i]), false)->get_grid_id());
}


Original_grid_info::Original_grid_info(int comp_id, int grid_id, const char *grid_name, const char *annotation, Remap_grid_class *original_CoR_grid, bool model_registration, bool recursively, bool calculate_checksum_H2D_mask)
{
    this->comp_id = comp_id;
    this->used_in_md_grid = false;
    this->original_CoR_grid = original_CoR_grid;
    this->bottom_field_variation_type = BOTTOM_FIELD_VARIATION_UNSET;
    this->V3D_lev_field_variation_type = BOTTOM_FIELD_VARIATION_UNSET;
    this->bottom_field_name = NULL;
    this->V3D_lev_field_id = -1;
    this->bottom_field_id = -1;
    this->checksum_H2D_mask = -1;
    this->checksum_H2D_center_lon = -1;
    this->checksum_H2D_center_lat = -1;
    this->mid_point_grid = NULL;
    this->interface_level_grid = NULL;
    this->grid_name = strdup(grid_name);
    this->ensemble_member_grid = NULL;
    comp_full_name = strdup(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "Original_grid_info")->get_full_name());
    generate_remapping_grids();

    max_sub_grid_under_V3D = NULL;
    Time1D_sub_grid = NULL;
    Tracer1D_sub_grid = NULL;
    if (recursively) {
        if (Time1D_sub_CoR_grid != NULL)
            Time1D_sub_grid = original_grid_mgr->search_or_add_grid_info(this, comp_id, Time1D_sub_CoR_grid);
        if (Tracer1D_sub_CoR_grid != NULL)
            Tracer1D_sub_grid = original_grid_mgr->search_or_add_grid_info(this, comp_id, Tracer1D_sub_CoR_grid);
        if (V3D_sub_CoR_grid != NULL)
            max_sub_grid_under_V3D = original_grid_mgr->search_or_add_grid_info(this, comp_id, V3D_sub_CoR_grid);
        else if (H2D_sub_CoR_grid != NULL)
            max_sub_grid_under_V3D = original_grid_mgr->search_or_add_grid_info(this, comp_id, H2D_sub_CoR_grid);
        else if (V1D_sub_CoR_grid != NULL)
            max_sub_grid_under_V3D = original_grid_mgr->search_or_add_grid_info(this, comp_id, V1D_sub_CoR_grid);
    }

    this->grid_id = original_grid_mgr->get_num_original_grids() | TYPE_GRID_LOCAL_ID_PREFIX;
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish allocating the original grid \"%s\" with ID %x for comp %x with CoR grid %s", grid_name, grid_id, comp_id, original_CoR_grid->get_grid_name());
    annotation_mgr->add_annotation(this->grid_id, "grid_registration", annotation);
}


void Original_grid_info::calculate_H2D_grid_checksums()
{
    if (H2D_sub_CoR_grid == NULL || checksum_H2D_mask != -1)
        return;

    MPI_Comm comm = comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comm_group();
    Distributed_H2D_grid_engine *distributed_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(comp_full_name, H2D_sub_CoR_grid);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, distributed_grid != NULL, "Software error in Original_grid_info::calculate_H2D_grid_checksums: %s %s", comp_full_name, grid_name);
    void *data_array = distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_grid_mask_field()->get_grid_data_field()->data_buf;
    checksum_H2D_mask = calculate_checksum_of_array(data_array, distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_grid_size(), sizeof(bool), distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_local_cell_global_indexes(), comm);
    data_array = distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf;
    checksum_H2D_center_lon = calculate_checksum_of_array(data_array, distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_grid_size(), sizeof(double), distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_local_cell_global_indexes(), comm);
    data_array = distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf;
    checksum_H2D_center_lat = calculate_checksum_of_array(data_array, distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_grid_size(), sizeof(double), distributed_grid->get_basic_decomp_grid()->get_decomp_grid()->get_local_cell_global_indexes(), comm);
}


void Original_grid_info::write_H2D_grid()
{
    if (!output_H2D_grid_enabled)
        return;
    if (H2D_sub_CoR_grid != NULL && V1D_sub_CoR_grid == NULL && Time1D_sub_CoR_grid == NULL && Tracer1D_sub_CoR_grid == NULL && !H2D_sub_CoR_grid->get_is_empty_grid()) {
        char nc_file_name[NAME_STR_SIZE];
        sprintf(nc_file_name, "%s/%s@%s.nc", comp_comm_group_mgt_mgr->get_internal_H2D_grids_dir(), grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "Original_grid_info")->get_full_name());
        if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in register_h2d_grid_with_data") == 0) {
            IO_netcdf *netcdf_file_object = new IO_netcdf("H2D_grid_data", nc_file_name, "w", false);
            netcdf_file_object->write_grid(H2D_sub_CoR_grid, false, true);
            netcdf_file_object->put_global_attr("edge_type", "LON_LAT", DATA_TYPE_STRING, DATA_TYPE_STRING, -1);   // to be modified
            Remap_grid_class *leaf_grids[256];
            int num_leaf_grids;
            H2D_sub_CoR_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, H2D_sub_CoR_grid);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(leaf_grids[0]->get_coord_label(), COORD_LABEL_LON), "software error in Original_grid_info::Original_grid_info");
            double temp_double = H2D_sub_CoR_grid->get_boundary_min_lon();
            netcdf_file_object->put_global_attr("min_lon", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            temp_double = H2D_sub_CoR_grid->get_boundary_max_lon();
            netcdf_file_object->put_global_attr("max_lon", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            temp_double = H2D_sub_CoR_grid->get_boundary_min_lat();
            netcdf_file_object->put_global_attr("min_lat", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            temp_double = H2D_sub_CoR_grid->get_boundary_max_lat();
            netcdf_file_object->put_global_attr("max_lat", &temp_double, DATA_TYPE_DOUBLE, DATA_TYPE_DOUBLE, 1);
            netcdf_file_object->put_global_attr("mask_checksum", &checksum_H2D_mask, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
            netcdf_file_object->put_global_attr("center_lon_checksum", &checksum_H2D_center_lon, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
            netcdf_file_object->put_global_attr("center_lat_checksum", &checksum_H2D_center_lat, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
            if (leaf_grids[0]->get_grid_cyclic())
                netcdf_file_object->put_global_attr("cyclic_or_acyclic", "cyclic", DATA_TYPE_STRING, DATA_TYPE_STRING, -1);
            else netcdf_file_object->put_global_attr("cyclic_or_acyclic", "acyclic", DATA_TYPE_STRING, DATA_TYPE_STRING, -1);
            netcdf_file_object->put_global_attr("title", grid_name, DATA_TYPE_STRING, DATA_TYPE_STRING, -1);
            delete netcdf_file_object;
            char status_file_name[NAME_STR_SIZE];
            sprintf(status_file_name, "%s/%s@%s.nc.end", comp_comm_group_mgt_mgr->get_comps_ending_config_status_dir(), grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "Original_grid_info")->get_full_name());
            FILE *status_file = fopen(status_file_name, "w+");
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, status_file != NULL, "Software error in Comp_comm_group_mgt_node::merge_comp_comm_info: configuration ending status file cannot be created: %s", status_file_name);
            fclose(status_file);
        }
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Original_grid_info::Original_grid_info"));
    }
}

Original_grid_info::~Original_grid_info()
{
    delete [] grid_name;
    delete [] comp_full_name;
    if (bottom_field_name != NULL)
        delete [] bottom_field_name;
}


void Original_grid_info::reset_grid_data()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    if (original_CoR_grid != NULL) {
        delete original_CoR_grid;
        original_CoR_grid = NULL;
    }
}



const char *Original_grid_info::get_annotation()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    return annotation_mgr->get_annotation(this->grid_id, "grid_registration");
}


void Original_grid_info::link_sub_grids()
{
    if (sub_grids_id.size() > 0)
        return;

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli1 link sub grids %s %lx %lx", grid_name, this, get_max_sub_grid_under_V3D());

    if (get_Time1D_sub_CoR_grid() != NULL && get_Time1D_sub_CoR_grid() != original_CoR_grid)
        sub_grids_id.push_back(original_grid_mgr->search_or_add_grid_info(this, comp_id, get_Time1D_sub_CoR_grid())->get_grid_id());

    if (get_V1D_sub_CoR_grid() != NULL && get_V1D_sub_CoR_grid() != original_CoR_grid)
        sub_grids_id.push_back(original_grid_mgr->search_or_add_grid_info(this, comp_id, get_V1D_sub_CoR_grid())->get_grid_id());

    if (get_H2D_sub_CoR_grid() != NULL && get_H2D_sub_CoR_grid() != original_CoR_grid)
        sub_grids_id.push_back(original_grid_mgr->search_or_add_grid_info(this, comp_id, get_H2D_sub_CoR_grid())->get_grid_id());

    if (Tracer1D_sub_CoR_grid != NULL && Tracer1D_sub_CoR_grid != original_CoR_grid)
        sub_grids_id.push_back(original_grid_mgr->search_or_add_grid_info(this, comp_id, Tracer1D_sub_CoR_grid)->get_grid_id());
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli2 link sub grids %s", grid_name);
}


void Original_grid_info::generate_remapping_grids()
{
    Remap_grid_class *lon_sub_grid = NULL, *lat_sub_grid = NULL, *leaf_grids[256];
    int num_leaf_grids;


    V3D_sub_CoR_grid = NULL;
    H2D_sub_CoR_grid = NULL;
    V1D_sub_CoR_grid = NULL;
    Time1D_sub_CoR_grid = NULL;
    Tracer1D_sub_CoR_grid = NULL;
    H2D_sub_grid_order = -1;
    V1D_sub_grid_order = -1;
    T1D_sub_grid_order = -1;
    tracer_sub_grid_order = -1;

    if (original_CoR_grid->get_num_dimensions() == 2 && (original_CoR_grid->has_grid_coord_label(COORD_LABEL_LON) || original_CoR_grid->has_grid_coord_label(COORD_LABEL_LAT)))
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_CoR_grid->get_is_H2D_grid(), "Software error in Original_grid_info::generate_remapping_grids: not a horizontal grid");
    if (original_CoR_grid->get_num_dimensions() == 1)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(original_CoR_grid->get_coord_label(), COORD_LABEL_LEV) || words_are_the_same(original_CoR_grid->get_coord_label(), COORD_LABEL_TIME) || words_are_the_same(original_CoR_grid->get_coord_label(), COORD_LABEL_TRACER), "Software error in Original_grid_info::generate_remapping_grids: not a vertical grid, not a time grid, or not a tracer grid");
    H2D_sub_CoR_grid = original_CoR_grid->get_sphere_sub_grid();
    V3D_sub_CoR_grid = original_CoR_grid->get_V3D_sub_grid();
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, V3D_sub_CoR_grid != NULL == original_CoR_grid->has_V3D_sub_grid(), "Software error in Original_grid_info::generate_remapping_grids: not a horizontal grid");

    original_CoR_grid->get_leaf_grids(&num_leaf_grids, leaf_grids, original_CoR_grid);
    for (int i = 0; i < num_leaf_grids; i ++) {
        if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LON)) {
            lon_sub_grid = leaf_grids[i];
            H2D_sub_grid_order = i;
        }
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LAT))
            lat_sub_grid = leaf_grids[i];
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_LEV)) {
            V1D_sub_CoR_grid = leaf_grids[i];
            V1D_sub_grid_order = i;
        }
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_TIME)) {
            Time1D_sub_CoR_grid = leaf_grids[i];
            T1D_sub_grid_order = i;
        }
        else if (words_are_the_same(leaf_grids[i]->get_coord_label(), COORD_LABEL_TRACER)) {
            Tracer1D_sub_CoR_grid = leaf_grids[i];
            tracer_sub_grid_order = i;
        }
    }

    if (H2D_sub_CoR_grid == NULL && (lon_sub_grid != NULL || lat_sub_grid != NULL)) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lon_sub_grid != NULL && lat_sub_grid != NULL, "Software error in Original_grid_info::generate_remapping_grids: wrong sub grids for a horizontal grid");
        leaf_grids[0] = lon_sub_grid;
        leaf_grids[1] = lat_sub_grid;
        char H2D_grid_name[NAME_STR_SIZE];
        sprintf(H2D_grid_name, "H2D_grid_for_%s", grid_name);
        H2D_sub_CoR_grid = new Remap_grid_class(H2D_grid_name, 2, leaf_grids, 0);
    }

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL || V1D_sub_CoR_grid != NULL || Time1D_sub_CoR_grid != NULL || Tracer1D_sub_CoR_grid, "Software error in Original_grid_info::generate_remapping_grids: empty grid");
}


bool Original_grid_info::is_V1D_sub_grid_after_H2D_sub_grid()
{
    if (H2D_sub_CoR_grid == NULL || V1D_sub_CoR_grid == NULL)
        return true;

    return H2D_sub_grid_order < V1D_sub_grid_order;
}


void Original_grid_info::set_unique_3D_lev_field(int field_id, const char *type, const char *annotation)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, bottom_field_variation_type == BOTTOM_FIELD_VARIATION_UNSET && V3D_lev_field_variation_type == BOTTOM_FIELD_VARIATION_UNSET, "Software error in Original_grid_info::set_unique_3D_lev_field");
    annotation_mgr->add_annotation(grid_id, "set 3-D level field", annotation);
    V3D_lev_field_variation_type = words_are_the_same(type, "constant") ? 0 : 1;
    V3D_lev_field_id = field_id;
}


void Original_grid_info::set_unique_bottom_field(int field_id, int type, const char *annotation)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, bottom_field_variation_type == BOTTOM_FIELD_VARIATION_UNSET, "Software error in Original_grid_info::set_unique_bottom_field");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !used_in_md_grid, "Error happens when setting the unique surface field of the 3D grid \"%s\": this grid has already been used to generate a multi-dimensional grid. Please check the model code with the annotation \"%s\".", grid_name, annotation);
    annotation_mgr->add_annotation(grid_id, "set surface field", annotation);
    bottom_field_variation_type = type;
    if (type != BOTTOM_FIELD_VARIATION_EXTERNAL) {
        Field_mem_info *field_inst = memory_manager->get_field_instance(field_id);
        bottom_field_name = strdup(field_inst->get_field_name());
        bottom_field_id = field_id;
        Remap_grid_class *decomp_grid = decomp_grids_mgr->search_decomp_grid_info(field_inst->get_decomp_id(), get_original_CoR_grid(), true)->get_decomp_grid();
        if (decomp_grid != NULL)
            decomp_grid->set_level_V3D_coord_dynamic_trigger_field(field_inst->get_field_data());
    }
}


void Original_grid_info::write_grid_into_array(char **temp_array_buffer, long &buffer_max_size, long &buffer_content_size)
{
    get_original_CoR_grid()->write_grid_into_array(temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&V3D_lev_field_variation_type, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&bottom_field_variation_type, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&checksum_H2D_mask, sizeof(long), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&checksum_H2D_center_lon, sizeof(long), temp_array_buffer, buffer_max_size, buffer_content_size);
    write_data_into_array_buffer(&checksum_H2D_center_lat, sizeof(long), temp_array_buffer, buffer_max_size, buffer_content_size);
}


void Original_grid_info::get_grid_data(const char *host_comp_full_name, int decomp_id, int chunk_index, const char *label, const char *data_type, int array_size, char *grid_data, const char *annotation, const char *API_label)
{
    Remap_grid_data_class *grid_field;
    Remap_grid_class *field_CoR_grid;
    Decomp_info *decomp_info = NULL;
    int offset = 0, grid_field_size = array_size;
    Comp_comm_group_mgt_node *host_comp_node = comp_comm_group_mgt_mgr->search_global_node(host_comp_full_name != NULL ? host_comp_full_name : comp_full_name);


    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), words_are_the_same(label, COORD_LABEL_LON) || words_are_the_same(label, COORD_LABEL_LAT) || words_are_the_same(label, GRID_MASK_LABEL), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the label (currently is \"%s\") is wrong (must be \"lon\", \"lat\" or \"mask\"). Please verify the model code with the annotation \"%s\".", API_label, grid_name, label, annotation);
    if (words_are_the_same(label, GRID_MASK_LABEL)) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), words_are_the_same(data_type, DATA_TYPE_INT), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the data type (currently is %s) corresponding to the parameter of \"grid_data\" is wrong when the label is \"mask\" (the right data type must be integer). Please verify the model code with the annotation \"%s\".", API_label, grid_name, data_type, annotation);
    }
    else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), words_are_the_same(data_type, DATA_TYPE_FLOAT) || words_are_the_same(data_type, DATA_TYPE_DOUBLE), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the data type (currently is %s) corresponding to the parameter of \"grid_data\" is wrong when the label (currently is \"%s\") is not \"mask\" (the right data type must be floating-point). Please verify the model code with the annotation \"%s\".", API_label, grid_name, data_type, label, annotation);

    EXECUTION_REPORT(REPORT_ERROR, host_comp_node->get_comp_id(), H2D_sub_CoR_grid != NULL && !H2D_sub_CoR_grid->get_is_empty_grid(), "Software error in Original_grid_info::get_grid_data");

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), array_size > 0, "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array of \"grid_data\" has not been allocated. Please verify the model code with the annotation \"%s\".", API_label, grid_name, annotation);
    if (decomp_id != -1) {
        field_CoR_grid = decomp_grids_mgr->search_decomp_grid_info(decomp_id, original_CoR_grid, false)->get_decomp_grid();
        decomp_info = decomps_info_mgr->get_decomp_info(decomp_id);
        grid_field_size = decomp_info->get_num_local_cells();
        if (decomp_info->get_num_chunks() == 0)
            chunk_index = -1;
        else if (chunk_index != -1)
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, decomp_info->get_comp_id(), chunk_index >= 1 && chunk_index <= decomp_info->get_num_chunks(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the given parallel decomposition (\"%s\") has %d chunks while the given \"chunk_index\" (%d) is out of the range [1,%d]. Please verify the model code with the annotation \"%s\".", API_label, grid_name, decomp_info->get_decomp_name(), decomp_info->get_num_chunks(), chunk_index, decomp_info->get_num_chunks(), annotation);
        if (chunk_index == -1) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), array_size == decomp_info->get_num_local_cells(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array size (currently is %d) of the parameter of \"grid_data\" is different from the size corresponding to the parallel decomposition (currently is %d). Please verify the model code with the annotation \"%s\".", API_label, grid_name, array_size, decomp_info->get_num_local_cells(), annotation);
        }
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), array_size == decomp_info->get_chunk_size(chunk_index - 1), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array size (currently is %d) of the parameter of \"grid_data\" is different from the size corresponding to the parallel decomposition (currently is %d). Please verify the model code with the annotation \"%s\".", API_label, grid_name, array_size, decomp_info->get_chunk_size(chunk_index - 1), annotation);
        for (int i = 0; i < chunk_index - 1; i ++)
            offset += decomp_info->get_chunk_size(i);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), original_grid_mgr->get_original_CoR_grid(decomp_info->get_grid_id()) == original_CoR_grid, "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the grid_id and decomp_id do not correspond to the same H2D grid. Please verify the model code with the annotation \"%s\".", API_label, grid_name, annotation);
    }
    else {
        chunk_index = -1;
        Distributed_H2D_grid_engine *distributed_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(host_comp_full_name != NULL ? host_comp_full_name : comp_full_name, H2D_sub_CoR_grid);
        if (distributed_grid != NULL)
            field_CoR_grid = distributed_grid->get_basic_decomp_grid()->get_decomp_grid();
        else field_CoR_grid = H2D_sub_CoR_grid;
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_node->get_comp_id(), array_size == H2D_sub_CoR_grid->get_grid_size(), "Error happens when calling the API \"%s\" to get the grid data of an H2D grid \"%s\": the array size (currently is %d) of the parameter of \"grid_data\" is different from the size corresponding to the parallel decomposition (currently is %d). Please verify the model code with the annotation \"%s\".", API_label, grid_name, array_size, field_CoR_grid->get_grid_size(), annotation);
    }

    if (!words_are_the_same(label, GRID_MASK_LABEL)) {
        grid_field = field_CoR_grid->get_grid_center_field(label);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_field != NULL, "Software error in Original_grid_info::get_grid_data: NULL grid center field");
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(grid_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "Software error in Original_grid_info::get_grid_data: wrong center field data type");
        if (!grid_field->get_coord_value_grid()->get_is_H2D_grid()) {
            grid_field = field_CoR_grid->expand_to_generate_full_coord_value(grid_field);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid == field_CoR_grid, "Software error in Original_grid_info::get_grid_data: NULL grid center field");
        }
        double *grid_field_data_buf = (double*)grid_field->get_grid_data_field()->data_buf;
        if (field_CoR_grid != H2D_sub_CoR_grid && decomp_id == -1) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, decomp_id == -1, "Software error in Original_grid_info::get_grid_data: NULL grid center field");
            int aggregated_grid_data_size;
            decomp_info = decomps_info_mgr->get_decomp_info(distributed_H2D_grid_mgr->search_distributed_H2D_grid(host_comp_full_name != NULL ? host_comp_full_name : comp_full_name, H2D_sub_CoR_grid)->get_basic_decomp_grid()->get_decomp_id());
            for (int i = 0; i < field_CoR_grid->get_grid_size(); i ++)
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_CoR_grid->get_local_cell_global_indexes()[i] == decomp_info->get_local_cell_global_indx()[i], "Software error in Original_grid_info::get_grid_data: %d: %d vs %d", i, field_CoR_grid->get_local_cell_global_indexes()[i], decomp_info->get_local_cell_global_indx()[i]);
            grid_field_data_buf = (double*)check_and_aggregate_local_grid_data(host_comp_node->get_comp_id(), API_ID_GRID_MGT_GET_H2D_GRID_DATA, host_comp_node->get_comm_group(), "Original_grid_info::get_grid_data", H2D_sub_CoR_grid->get_grid_size(), field_CoR_grid->get_grid_size(), sizeof(double), (char*) grid_field->get_grid_data_field()->data_buf, label, field_CoR_grid->get_grid_size(), field_CoR_grid->get_local_cell_global_indexes(), aggregated_grid_data_size, "in Original_grid_info::get_grid_data");
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, aggregated_grid_data_size == H2D_sub_CoR_grid->get_grid_size(), "Software error in Original_grid_info::get_grid_data: NULL grid center field");
        }
        if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
            transform_datatype_of_arrays(grid_field_data_buf, (float*)grid_data, array_size);
        else transform_datatype_of_arrays(grid_field_data_buf, (double*)grid_data, array_size);
        if (grid_field != field_CoR_grid->get_grid_center_field(label))
            delete grid_field;
        if (grid_field_data_buf != (double*)grid_field->get_grid_data_field()->data_buf)
            delete [] grid_field_data_buf;
    }
    else {
        grid_field = field_CoR_grid->get_grid_mask_field();
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_field != NULL, "Software error in Original_grid_info::get_grid_data: NULL grid mask field");
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(grid_field->get_grid_data_field()->data_type_in_application, DATA_TYPE_BOOL), "Software error in Original_grid_info::get_grid_data: wrong mask field data type");
        bool *grid_field_data_buf = (bool*)grid_field->get_grid_data_field()->data_buf;
        if (field_CoR_grid != H2D_sub_CoR_grid && decomp_id == -1) {
            int aggregated_grid_data_size;
            grid_field_data_buf = (bool*)check_and_aggregate_local_grid_data(host_comp_node->get_comp_id(), API_ID_GRID_MGT_GET_H2D_GRID_DATA, host_comp_node->get_comm_group(), "Original_grid_info::get_grid_data", H2D_sub_CoR_grid->get_grid_size(), field_CoR_grid->get_grid_size(), sizeof(bool), (char*) grid_field->get_grid_data_field()->data_buf, label, field_CoR_grid->get_grid_size(), field_CoR_grid->get_local_cell_global_indexes(), aggregated_grid_data_size, "in Original_grid_info::get_grid_data");
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, aggregated_grid_data_size == H2D_sub_CoR_grid->get_grid_size(), "Software error in Original_grid_info::get_grid_data: NULL grid center field");
        }
        transform_datatype_of_arrays(grid_field_data_buf, (int*)grid_data, array_size);
        if (grid_field_data_buf != (bool*)grid_field->get_grid_data_field()->data_buf)
            delete [] grid_field_data_buf;
    }
}


void Original_grid_info::set_mid_point_grid(Original_grid_info *new_grid)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    this->mid_point_grid = new_grid;
    new_grid->interface_level_grid = this;
    if (this->bottom_field_variation_type != BOTTOM_FIELD_VARIATION_UNSET)
        new_grid->set_unique_bottom_field(this->bottom_field_id, this->bottom_field_variation_type, annotation_mgr->get_annotation(this->grid_id, "set surface field"));
}


void Original_grid_info::set_grid_checksum(long checksum_mask, long checksum_lon, long checksum_lat)
{
    this->checksum_H2D_mask = checksum_mask;
    this->checksum_H2D_center_lon = checksum_lon;
    this->checksum_H2D_center_lat = checksum_lat;
}


double *Original_grid_info::get_center_lon_values(const char *host_comp_full_name)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL, "Software error in Original_grid_info::get_center_lon_values");
    double *center_lon_values = new double [H2D_sub_CoR_grid->get_grid_size()];
    get_grid_data(host_comp_full_name, -1, -1, "lon", DATA_TYPE_DOUBLE, H2D_sub_CoR_grid->get_grid_size(), (char*) center_lon_values, "internal", "internal");
    return center_lon_values;
}


double *Original_grid_info::get_center_lat_values(const char *host_comp_full_name)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL, "Software error in Original_grid_info::get_center_lat_values");
    double *center_lat_values = new double [H2D_sub_CoR_grid->get_grid_size()];
    get_grid_data(host_comp_full_name, -1, -1, "lat", DATA_TYPE_DOUBLE, H2D_sub_CoR_grid->get_grid_size(), (char*) center_lat_values, "internal", "internal");
    return center_lat_values;
}


bool Original_grid_info::is_H2D_grid_and_the_same_as_another_grid(const char *host_comp_full_name, Original_grid_info *another_grid)
{
    bool check_result, distributed_check_result;
    double *this_center_lon, *this_center_lat, *another_center_lon, *another_center_lat;
    Distributed_H2D_grid_engine *this_distributed_H2D_grid, *another_distributed_H2D_grid;


    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->H2D_sub_CoR_grid == NULL && another_grid->H2D_sub_CoR_grid == NULL || this->H2D_sub_CoR_grid != NULL && another_grid->H2D_sub_CoR_grid != NULL, "Software error in Original_grid_info::is_V1D_sub_grid_the_same_as_another_grid");

    if (this->H2D_sub_CoR_grid == NULL)
        return true;

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->get_H2D_sub_grid()->get_checksum_H2D_mask() != -1 && another_grid->get_H2D_sub_grid()->get_checksum_H2D_mask() != -1 && this->get_H2D_sub_grid()->get_checksum_H2D_center_lon() != -1 && another_grid->get_H2D_sub_grid()->get_checksum_H2D_center_lon() != -1 && this->get_H2D_sub_grid()->get_checksum_H2D_center_lat() != -1 && another_grid->get_H2D_sub_grid()->get_checksum_H2D_center_lat() != -1, "Software error in Original_grid_info::is_H2D_grid_and_the_same_as_another_grid");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->get_H2D_sub_grid() != NULL && another_grid->get_H2D_sub_grid() != NULL, "Software error in Original_grid_info::is_H2D_grid_and_the_same_as_another_grid");

    return this->get_H2D_sub_grid()->get_checksum_H2D_mask() == another_grid->get_H2D_sub_grid()->get_checksum_H2D_mask() && this->get_H2D_sub_grid()->get_checksum_H2D_center_lon() == another_grid->get_H2D_sub_grid()->get_checksum_H2D_center_lon() && this->get_H2D_sub_grid()->get_checksum_H2D_center_lat() == another_grid->get_H2D_sub_grid()->get_checksum_H2D_center_lat();

    this_distributed_H2D_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(host_comp_full_name, this->H2D_sub_CoR_grid);
    another_distributed_H2D_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(host_comp_full_name, another_grid->H2D_sub_CoR_grid);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, another_distributed_H2D_grid != NULL, "Software error in Original_grid_info::is_H2D_grid_and_the_same_as_another_grid: \"%s\"  \"%s\"", another_grid->comp_full_name, another_grid->H2D_sub_CoR_grid->get_grid_name());
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this_distributed_H2D_grid != NULL, "Software error in Original_grid_info::is_H2D_grid_and_the_same_as_another_grid: \"%s\" \"%s\"", this->comp_full_name, this->H2D_sub_CoR_grid->get_grid_name());
    distributed_check_result = this_distributed_H2D_grid->is_the_same_as_another_distributed_H2D_grid(another_distributed_H2D_grid);

    if (report_error_enabled) {
        this_center_lon = this->get_center_lon_values(host_comp_full_name);
        another_center_lon = another_grid->get_center_lon_values(host_comp_full_name);
        check_result = are_two_coord_arrays_same(this_center_lon, another_center_lon, this->get_H2D_sub_CoR_grid()->get_grid_size(), another_grid->get_H2D_sub_CoR_grid()->get_grid_size(), false);
        delete [] this_center_lon;
        delete [] another_center_lon;
        if (!check_result) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !distributed_check_result, "Software error in Original_grid_info::is_H2D_grid_and_the_same_as_another_grid");
            return false;
        }
        this_center_lat = this->get_center_lat_values(host_comp_full_name);
        another_center_lat = another_grid->get_center_lat_values(host_comp_full_name);
        check_result = are_two_coord_arrays_same(this_center_lat, another_center_lat, this->get_H2D_sub_CoR_grid()->get_grid_size(), another_grid->get_H2D_sub_CoR_grid()->get_grid_size(), false);
        delete [] this_center_lat;
        delete [] another_center_lat;
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, check_result && distributed_check_result || !check_result && !distributed_check_result, "Software error in Original_grid_info::is_H2D_grid_and_the_same_as_another_grid");
    }

    return false;
}


bool Original_grid_info::is_V1D_sub_grid_the_same_as_another_grid(Original_grid_info *another_grid)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->V1D_sub_CoR_grid == NULL && another_grid->V1D_sub_CoR_grid == NULL || this->V1D_sub_CoR_grid != NULL && another_grid->V1D_sub_CoR_grid != NULL, "Software error in Original_grid_info::is_V1D_sub_grid_the_same_as_another_grid");

    if (this->V1D_sub_CoR_grid == NULL)
        return true;

    //if (this->V1D_sub_CoR_grid->is_sigma_grid() || another_grid->V1D_sub_CoR_grid->is_sigma_grid())
    if (this->get_V3D_sub_CoR_grid() != NULL || another_grid->get_V3D_sub_CoR_grid() != NULL)
        if (this->get_V3D_sub_CoR_grid()->is_sigma_grid() || another_grid->get_V3D_sub_CoR_grid()->is_sigma_grid())
            return false;

    if (!words_are_the_same(this->V1D_sub_CoR_grid->get_coord_unit(), another_grid->V1D_sub_CoR_grid->get_coord_unit()))
        return false;

    if (this->V1D_sub_CoR_grid->get_grid_center_field() == NULL || another_grid->V1D_sub_CoR_grid->get_grid_center_field() == NULL)
        return false;

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(this->V1D_sub_CoR_grid->get_grid_center_field()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE) && words_are_the_same(another_grid->V1D_sub_CoR_grid->get_grid_center_field()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "Software error in Original_grid_info::is_V1D_sub_grid_the_same_as_another_grid");

    return are_two_coord_arrays_same((double*)(this->V1D_sub_CoR_grid->get_grid_center_field()->get_grid_data_field()->data_buf), (double*)(another_grid->V1D_sub_CoR_grid->get_grid_center_field()->get_grid_data_field()->data_buf), this->V1D_sub_CoR_grid->get_grid_size(), another_grid->V1D_sub_CoR_grid->get_grid_size(), true);
}


bool Original_grid_info::is_Time_sub_grid_the_same_as_another_grid(Original_grid_info *another_grid)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, Time1D_sub_CoR_grid == NULL, "Time grid is not supported yet");
    return true;
}


bool Original_grid_info::is_Tracer_sub_grid_the_same_as_another_grid(Original_grid_info *another_grid)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->Tracer1D_sub_CoR_grid == NULL && another_grid->Tracer1D_sub_CoR_grid == NULL || this->Tracer1D_sub_CoR_grid != NULL && another_grid->Tracer1D_sub_CoR_grid != NULL, "Software error in Original_grid_info::is_Tracer_sub_grid_the_same_as_another_grid");

    if (this->Tracer1D_sub_CoR_grid == NULL)
        return true;

    return this->Tracer1D_sub_CoR_grid->get_grid_size() == another_grid->Tracer1D_sub_CoR_grid->get_grid_size();
}


bool Original_grid_info::is_the_same_as_another_grid(const char *host_comp_full_name, Original_grid_info *another_grid)
{
    if (this->H2D_sub_grid_order != another_grid->H2D_sub_grid_order || this->V1D_sub_grid_order != another_grid->V1D_sub_grid_order ||
            this->T1D_sub_grid_order != another_grid->T1D_sub_grid_order || this->tracer_sub_grid_order != another_grid->tracer_sub_grid_order)
        return false;

    return is_H2D_grid_and_the_same_as_another_grid(host_comp_full_name, another_grid) && is_V1D_sub_grid_the_same_as_another_grid(another_grid) && is_Time_sub_grid_the_same_as_another_grid(another_grid) && is_Tracer_sub_grid_the_same_as_another_grid(another_grid);
}


void Original_grid_info::copy_bottom_field_variation_type(Original_grid_info *dst_grid)
{
    if (bottom_field_variation_type != BOTTOM_FIELD_VARIATION_UNSET) {
        dst_grid->bottom_field_variation_type = bottom_field_variation_type;
        dst_grid->bottom_field_id = bottom_field_id;
    }
}


bool Original_grid_info::has_sigma_sub_grid()
{
    return V3D_sub_CoR_grid != NULL && V3D_sub_CoR_grid->is_sigma_grid();
}


bool Original_grid_info::does_use_V3D_level_coord()
{
    return V3D_sub_CoR_grid != NULL && V3D_sub_CoR_grid->does_use_V3D_level_coord();
}


int Original_grid_info::get_bottom_field_id()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    if (bottom_field_id != -1)
        return bottom_field_id;

    return original_grid_mgr->get_bottom_field_id_for_grid(this);
}


int Original_grid_info::get_V3D_lev_field_id()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    if (V3D_lev_field_id != -1)
        return V3D_lev_field_id;

    return original_grid_mgr->get_V3D_lev_field_id_for_grid(this);
}


Original_grid_info *Original_grid_info::get_interface_level_grid()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    return interface_level_grid;
}


Original_grid_info *Original_grid_info::get_mid_point_grid()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ensemble_member_grid == NULL);
    return mid_point_grid;
}


Original_grid_info *Original_grid_info::get_H2D_sub_grid()
{
    if (is_H2D_grid())
        return this;

    for (int i = 0; i < sub_grids_id.size(); i ++) {
        Original_grid_info *H2D_sub_grid = original_grid_mgr->get_original_grid(sub_grids_id[i])->get_H2D_sub_grid();
        if (H2D_sub_grid != NULL)
            return H2D_sub_grid;
    }

    if (interface_level_grid != NULL)
        return interface_level_grid->get_H2D_sub_grid();

    return NULL;
}


bool Original_grid_info::get_H2D_sub_grid_full_name(char *full_name)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL, "Software error in Original_grid_info::get_H2D_sub_grid_full_name");

    if (is_H2D_grid()) {
        sprintf(full_name, "%s@%s", grid_name, comp_full_name);
        return true;
    }

    sprintf(full_name, "H2D_of_%s@%s", grid_name, comp_full_name);
    return true;
}


Original_grid_info *Original_grid_info::get_sub_original_grid_corresponding_to_another(Original_grid_info *another_grid)
{
    if (another_grid->get_V1D_sub_CoR_grid() != NULL && this->get_V1D_sub_CoR_grid() == NULL ||
            another_grid->get_H2D_sub_CoR_grid() != NULL && this->get_H2D_sub_CoR_grid() == NULL ||
            another_grid->get_Tracer1D_sub_grid() != NULL && this->get_Tracer1D_sub_grid() == NULL ||
            another_grid->get_Time1D_sub_grid() != NULL && this->get_V1D_sub_CoR_grid() == NULL)
        return NULL;

    if ((another_grid->get_V1D_sub_CoR_grid() == NULL && this->get_V1D_sub_CoR_grid() == NULL || another_grid->get_V1D_sub_CoR_grid() != NULL && this->get_V1D_sub_CoR_grid() != NULL) &&
            (another_grid->get_H2D_sub_CoR_grid() == NULL && this->get_H2D_sub_CoR_grid() == NULL || another_grid->get_H2D_sub_CoR_grid() != NULL && this->get_H2D_sub_CoR_grid() != NULL) &&
            (another_grid->get_Tracer1D_sub_grid() == NULL && this->get_Tracer1D_sub_grid() == NULL || another_grid->get_Tracer1D_sub_grid() != NULL && this->get_Tracer1D_sub_grid() != NULL) &&
            (another_grid->get_Time1D_sub_grid() == NULL && this->get_Time1D_sub_grid() == NULL || another_grid->get_Time1D_sub_grid() != NULL && this->get_Time1D_sub_grid() != NULL))
        return this;

    for (int i = 0; i < sub_grids_id.size(); i ++) {
        Original_grid_info *sub_original_grid = original_grid_mgr->search_grid_info(sub_grids_id[i])->get_sub_original_grid_corresponding_to_another(another_grid);
        if (sub_original_grid != NULL)
            return sub_original_grid;
    }

    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Encounter a special case in Original_grid_info::get_sub_original_grid_corresponding_to_another. Please ask Dr. Liu (liuli-cess@tsinghua.edu.cn) for specific help.");

    return NULL;
}

void Original_grid_info::rename_grid_name(char *new_grid_name) {
    delete [] grid_name;
    grid_name = const_cast<char*>(strdup(new_grid_name));
};


long Original_grid_info::get_checksum_H2D_mask()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL && checksum_H2D_mask != -1, "Software error in Original_grid_info::get_checksum_H2D_mask: %s", grid_name);
    return checksum_H2D_mask;
}


long Original_grid_info::get_checksum_H2D_center_lon()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL && checksum_H2D_center_lon != -1, "Software error in Original_grid_info::checksum_H2D_center_lon");
    return checksum_H2D_center_lon;
}


long Original_grid_info::get_checksum_H2D_center_lat()
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, H2D_sub_CoR_grid != NULL && checksum_H2D_center_lat != -1, "Software error in Original_grid_info::checksum_H2D_center_lat");
    return checksum_H2D_center_lat;
}


Original_grid_mgt::Original_grid_mgt()
{
    original_grids.clear();
    CoR_script_name[0] = '\0';
    CoR_grids = NULL;
}


void Original_grid_mgt::initialize_CoR_grids()
{
    if (CoR_grids != NULL)
        return;

    sprintf(CoR_script_name, "%s/CCPL_grid.cor", comp_comm_group_mgt_mgr->get_root_comp_config_dir());
    FILE *fp = fopen(CoR_script_name, "r");
    if (fp == NULL)
        CoR_script_name[0] = '\0';
    else fclose(fp);
    if (strlen(CoR_script_name) != 0) {
        char current_dir[NAME_STR_SIZE], grids_dir[NAME_STR_SIZE];
        EXECUTION_REPORT(REPORT_ERROR, -1, getcwd(current_dir, NAME_STR_SIZE) != NULL, "Cannot get the current working directory for running the model");
        sprintf(grids_dir, "%s/grids_weights", comp_comm_group_mgt_mgr->get_root_comp_config_dir());
        EXECUTION_REPORT(REPORT_ERROR, -1, chdir(grids_dir) == 0, "Fail to access the directory of the CoR grid data files: \"%s\". Please verify.", grids_dir);
        CoR_grids = new Remap_mgt(CoR_script_name);
        chdir(current_dir);
    }
    else CoR_grids = new Remap_mgt(NULL);
}


Original_grid_mgt::~Original_grid_mgt()
{
    for (int i = 0; i < original_grids.size(); i ++)
        if (original_grids[i] != NULL)
            delete original_grids[i];

    delete CoR_grids;
}


void Original_grid_mgt::delete_external_original_grids()
{
    for (int i = 0; i < original_grids.size(); i ++) {
        if (original_grids[i] == NULL)
            continue;
        Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grids[i]->get_comp_full_name());
        if (comp_node == NULL || comp_node->get_current_proc_local_id() == -1) {
            original_grids[i]->reset_grid_data();
            delete original_grids[i];
            original_grids[i] = NULL;
        }
    }
}


Original_grid_info *Original_grid_mgt::search_grid_info(const char *grid_name, int comp_id)
{
    for (int i = 0; i < original_grids.size(); i ++)
        if (original_grids[i] != NULL)
            if (words_are_the_same(original_grids[i]->get_grid_name(), grid_name) && original_grids[i]->get_comp_id() == comp_id)
                return original_grids[i];

    return NULL;
}


Original_grid_info *Original_grid_mgt::search_grid_info(int comp_id, Remap_grid_class *CoR_grid)
{
    for (int i = 0; i < original_grids.size(); i ++)
        if (original_grids[i] != NULL)
            if (original_grids[i]->get_original_CoR_grid() == CoR_grid && original_grids[i]->get_comp_id() == comp_id)
                return original_grids[i];

    return NULL;
}



Original_grid_info *Original_grid_mgt::search_or_add_grid_info(Original_grid_info *root_original_grid, int comp_id, Remap_grid_class *CoR_grid)
{
    char grid_name[NAME_STR_SIZE];


    if (root_original_grid->get_original_CoR_grid() == CoR_grid)
        return root_original_grid;

    for (int i = 0; i < original_grids.size(); i ++)
        if (original_grids[i] != NULL)
            if (original_grids[i]->get_comp_id() == comp_id && original_grids[i]->get_original_CoR_grid() == CoR_grid)
                return original_grids[i];

    sprintf(grid_name, "original_sub_grid_%s", CoR_grid->get_grid_name());
    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, "registering temp original grid", CoR_grid, false, true, false));
    return original_grids[original_grids.size() - 1];
}


Original_grid_info *Original_grid_mgt::search_grid_info(int grid_id)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id), "software error in Original_grid_mgt::search_grid_info based on grid_id");
    return original_grids[grid_id & TYPE_ID_SUFFIX_MASK];
}


int Original_grid_mgt::register_H2D_grid_via_comp(int comp_id, const char *grid_name, const char *annotation)
{
    char XML_file_name[NAME_STR_SIZE], nc_file_name[NAME_STR_SIZE], status_file_name[NAME_STR_SIZE];
    const char *another_comp_full_name = NULL, *another_comp_grid_name = NULL;
    int line_number;


    sprintf(XML_file_name, "%s/all/coupling_connections/%s.coupling_connections.xml", comp_comm_group_mgt_mgr->get_config_root_dir(), comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in register_H2D_grid_via_comp")->get_full_name());
    TiXmlDocument *XML_file = open_XML_file_to_read(comp_id, XML_file_name, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_H2D_grid_via_comp"), false);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, XML_file != NULL, "Error happens when calling the API \"CCPL_register_H2D_grid_from_another_component\" to register an H2D grid \"%s\": the coupling connection configuration file (\"%s\") does not exist. The API call is at the model code with the annotation \"%s\". ", grid_name, XML_file_name, annotation);

    TiXmlElement *root_XML_element;
    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(comp_id, XML_file_name, XML_file);
    for (; root_XML_element_node != NULL; root_XML_element_node = root_XML_element_node->NextSibling()) {
        if (root_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
            continue;
        root_XML_element = root_XML_element_node->ToElement();
        if (words_are_the_same(root_XML_element->Value(), "local_grids"))
            break;
    }
    if (root_XML_element_node != NULL) {
        for (TiXmlNode *grid_XML_element_node = root_XML_element->FirstChild(); grid_XML_element_node != NULL; grid_XML_element_node = grid_XML_element_node->NextSibling()) {
            if (grid_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                continue;
            TiXmlElement *grid_XML_element = grid_XML_element_node->ToElement();
            const char *xml_grid_name = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, grid_XML_element, "local_grid_name", XML_file_name, line_number, "grid name of the current component", "the coupling connection configuration file", true);
            if (words_are_the_same(xml_grid_name, grid_name)) {
                another_comp_full_name = get_XML_attribute(comp_id, 512, grid_XML_element, "another_comp_full_name", XML_file_name, line_number, "the full name of the another component", "the coupling connection configuration file", true);
                another_comp_grid_name = get_XML_attribute(comp_id, CCPL_NAME_STR_LEN, grid_XML_element, "another_comp_grid_name", XML_file_name, line_number, "the grid name of the another component", "the coupling connection configuration file", true);
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, strlen(another_comp_grid_name) > 0, "Error happens when calling the API \"CCPL_register_H2D_grid_from_another_component\" to register an H2D grid \"%s\": the coupling connection configuration file (\"%s\") specifies an empty name of the remote grid. Please check the XML file around line number %d", grid_name, XML_file_name, line_number);
                break;
            }
        }
    }

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, another_comp_full_name != NULL && another_comp_grid_name != NULL, "Error happens when calling the API \"CCPL_register_H2D_grid_from_another_component\" to register an H2D grid \"%s\": the coupling connection configuration file (\"%s\") does not contain the information for this grid. The API call is at the model code with the annotation \"%s\". ", grid_name, XML_file_name, annotation);
    sprintf(nc_file_name, "%s/%s@%s.nc", comp_comm_group_mgt_mgr->get_internal_H2D_grids_dir(), another_comp_grid_name, another_comp_full_name);
    sprintf(status_file_name, "%s/%s@%s.nc.end", comp_comm_group_mgt_mgr->get_comps_ending_config_status_dir(), another_comp_grid_name, another_comp_full_name);
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Wait to read NetCDF file \"%s\" to register H2D grid \"%s\" based on the grid \"%s\" of remote component \"%s\". Dead wait will be encounted if the full name of the remote component is wrong. So please make sure the full name of the remote component is correct in the the coupling connection configuration file (\"%s\")", nc_file_name, grid_name, another_comp_grid_name, another_comp_full_name, XML_file_name);
    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in register_H2D_grid_via_comp") == 0) {
        while (true) {
            FILE *status_file = fopen(status_file_name, "r");
            if (status_file == NULL) {
                if (comp_comm_group_mgt_mgr->has_comp_ended_configuration(another_comp_full_name))
                    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Fail to read NetCDF file \"%s\" to register H2D grid \"%s\": the remote component model \"%s\" has ended its coupling configuration stage without registering the required grid (\"%s\") before. Please check the configuration file (\"%s\") or the corresponding model code.", nc_file_name, grid_name, another_comp_full_name, another_comp_grid_name, XML_file_name);
                continue;
            }
            fclose(status_file);
            break;
        }
    }

    delete XML_file;
    synchronize_comp_processes_for_API(comp_id, API_ID_GRID_MGT_REG_H2D_GRID_VIA_COMP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_H2D_grid_via_comp"), "register_H2D_grid_via_comp", annotation);
    if (!report_error_enabled)
        MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_H2D_grid_via_comp"));

    return register_H2D_grid_via_file(comp_id, grid_name, nc_file_name, annotation);
}


void Original_grid_mgt::calculate_min_max_H2D_coord_value(int comp_id, char *center_coord, char *vertex_coord, int size_center_coord, int size_vertex_coord, const char *data_type, double &min_coord_value, double &max_coord_value)
{
    if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
        float float_min_value = (float)NULL_COORD_VALUE, float_max_value = (float)NULL_COORD_VALUE;
        get_min_value_in_array((float*)center_coord, size_center_coord, true, (float)NULL_COORD_VALUE, float_min_value);
        get_min_value_in_array((float*)vertex_coord, size_vertex_coord, true, (float)NULL_COORD_VALUE, float_min_value);
        get_max_value_in_array((float*)center_coord, size_center_coord, true, (float)NULL_COORD_VALUE, float_max_value);
        get_max_value_in_array((float*)vertex_coord, size_vertex_coord, true, (float)NULL_COORD_VALUE, float_max_value);
        min_coord_value = float_min_value;
        max_coord_value = float_max_value;
    }
    else {
        min_coord_value = NULL_COORD_VALUE, max_coord_value = NULL_COORD_VALUE;
        get_min_value_in_array((double*)center_coord, size_center_coord, true, NULL_COORD_VALUE, min_coord_value);
        get_min_value_in_array((double*)vertex_coord, size_vertex_coord, true, NULL_COORD_VALUE, min_coord_value);
        get_max_value_in_array((double*)center_coord, size_center_coord, true, NULL_COORD_VALUE, max_coord_value);
        get_max_value_in_array((double*)vertex_coord, size_vertex_coord, true, NULL_COORD_VALUE, max_coord_value);
    }
    double *all_min_coord_values = new double [comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_num_procs()];
    double *all_max_coord_values = new double [comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_num_procs()];
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The min and max values calculated by C-Coupler is %lf and %lf", min_coord_value, max_coord_value);
    MPI_Gather(&min_coord_value, 1, MPI_DOUBLE, all_min_coord_values, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, ""));
    MPI_Gather(&max_coord_value, 1, MPI_DOUBLE, all_max_coord_values, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, ""));
    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "") == 0)
        for (int i = 0; i < comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_num_procs(); i ++) {
            if (are_floating_values_equal(NULL_COORD_VALUE, min_coord_value) || !are_floating_values_equal(NULL_COORD_VALUE, all_min_coord_values[i]) && min_coord_value > all_min_coord_values[i])
                min_coord_value = all_min_coord_values[i];
            if (are_floating_values_equal(NULL_COORD_VALUE, max_coord_value) || !are_floating_values_equal(NULL_COORD_VALUE, all_max_coord_values[i]) && max_coord_value < all_max_coord_values[i])
                max_coord_value = all_max_coord_values[i];
        }
    MPI_Bcast(&min_coord_value, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, ""));
    MPI_Bcast(&max_coord_value, 1, MPI_DOUBLE, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, ""));
    delete [] all_min_coord_values;
    delete [] all_max_coord_values;
}


void Original_grid_mgt::common_checking_for_H2D_registration_via_data(int comp_id, const char *grid_name, const char *edge_type, const char *coord_unit, char *cyclic_or_acyclic, const char *data_type, int size_mask, int size_center_lon,
        int size_center_lat, int size_vertex_lon, int size_vertex_lat, int *mask, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, char *vertex_lon, char *vertex_lat, const char *annotation, int API_id)
{
    char API_label[NAME_STR_SIZE], hint[NAME_STR_SIZE];
    double eps = 1.0000001;
    double min_lon_value, max_lon_value, min_lat_value, max_lat_value, real_min_lon_value, real_max_lon_value, real_min_lat_value, real_max_lat_value;


    get_API_hint(comp_id, API_id, API_label);
    sprintf(hint, "registering an H2D grid \"%s\"", grid_name);

    transform_datatype_of_arrays(min_lon, (char*)(&min_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(max_lon, (char*)(&max_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(min_lat, (char*)(&min_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
    transform_datatype_of_arrays(max_lat, (char*)(&max_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, edge_type, "edge_type", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, cyclic_or_acyclic, "cyclic_or_acyclic", annotation);
    check_API_parameter_string(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, data_type, "implicit data type", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, min_lon_value, "\"min_lon\"", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, max_lon_value, "\"max_lon\"", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, min_lat_value, "\"min_lat\"", annotation);
    check_API_parameter_double(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, max_lat_value, "\"max_lat\"", annotation);

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(edge_type, "LON_LAT") || words_are_the_same(edge_type, "GREAT_ARC") || words_are_the_same(edge_type, "XY") || words_are_the_same(edge_type, "TriPolar"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is \"%s\") of parameter \"edge_type\" is not \"LON_LAT\", \"GREAT_ARC\", \"TriPolar\" or \"XY\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, edge_type, annotation);
    if (words_are_the_same(edge_type, "LON_LAT") || words_are_the_same(edge_type, "GREAT_ARC") || words_are_the_same(edge_type, "TriPolar"))
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(coord_unit, COORD_UNIT_DEGREES) || words_are_the_same(coord_unit, COORD_UNIT_RADIANS), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the parameter \"coord_unit\" (currently is \"%s\") is not \"degrees\" or \"radians\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, coord_unit, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(cyclic_or_acyclic, "cyclic") || words_are_the_same(cyclic_or_acyclic, "acyclic"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is \"%s\") of parameter \"cyclic_or_acyclic\" is not \"cyclic\" or \"acyclic\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, cyclic_or_acyclic, annotation);
    if (words_are_the_same(edge_type, "XY")) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(cyclic_or_acyclic, "acyclic"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is \"%s\") of parameter \"cyclic_or_acyclic\" is not \"acyclic\" when the parameter \"edge_type\" has been set to \"XY\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, cyclic_or_acyclic, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !words_are_the_same(coord_unit, COORD_UNIT_DEGREES) && !words_are_the_same(coord_unit, COORD_UNIT_RADIANS), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the parameter \"coord_unit\" (currently is \"%s\") can not be \"degrees\" or \"radians\" when the \"edge_type\" is \"XY\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, coord_unit, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the current \"edge_type\" \"XY\" has not been fully supported in C-Coupler. Please contact Dr. Li Liu for help.", grid_name, API_label);
    }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type, DATA_TYPE_FLOAT) || words_are_the_same(data_type, DATA_TYPE_DOUBLE), "software error in register_h2d_grid_with_data: wrong implicit data type");

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries("integer", mask, size_mask, 0, 1, 0, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"mask\" are wrong (not 0 and 1). Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_floating_values_equal(NULL_COORD_VALUE, min_lon_value) && are_floating_values_equal(NULL_COORD_VALUE, max_lon_value) || !are_floating_values_equal(NULL_COORD_VALUE, min_lon_value) && !are_floating_values_equal(NULL_COORD_VALUE, max_lon_value), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": \"min_lon\" (currently is %lf) and \"max_lon\" (currently is %lf) must be CCPL_NULL_COORD_VALUE (%lf) or not at the same time. Please check the model code related to the annotation \"%s\"", grid_name, API_label, min_lon_value, max_lon_value, NULL_COORD_VALUE, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_floating_values_equal(NULL_COORD_VALUE, min_lat_value) && are_floating_values_equal(NULL_COORD_VALUE, max_lat_value) || !are_floating_values_equal(NULL_COORD_VALUE, min_lat_value) && !are_floating_values_equal(NULL_COORD_VALUE, max_lat_value), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": \"min_lat\" (currently is %lf) and \"max_lat\" (currently is %lf) must be CCPL_NULL_COORD_VALUE (%lf) or not at the same time. Please check the model code related to the annotation \"%s\"", grid_name, API_label, min_lat_value, max_lat_value, NULL_COORD_VALUE, annotation);
    if (size_center_lon != 0) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lon < 0 || size_vertex_lon % size_center_lon == 0, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size of the parameter \"vertex_lon\" must be an integer multiple of the array size of the parameter \"center_lon\" (the current two array sizes are %d and %d respectively). Please check the model code related to the annotation \"%s\"", grid_name, API_label, size_vertex_lon, size_center_lon, annotation);
    }
    else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lon <= 0, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size of the parameter \"vertex_lon\" must be an integer multiple of the array size of the parameter \"center_lon\" (the current two array sizes are %d and %d respectively). Please check the model code related to the annotation \"%s\"", grid_name, API_label, size_vertex_lon, size_center_lon, annotation);
    calculate_min_max_H2D_coord_value(comp_id, center_lat, vertex_lat, size_center_lat, size_vertex_lat, data_type, real_min_lat_value, real_max_lat_value);
    calculate_min_max_H2D_coord_value(comp_id, center_lon, vertex_lon, size_center_lon, size_vertex_lon, data_type, real_min_lon_value, real_max_lon_value);
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The min and max latitude values of the grid \"%s\" calculated by C-Coupler is %lf and %lf", grid_name, real_min_lat_value, real_max_lat_value);
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The min and max longitude values of the grid \"%s\" calculated by C-Coupler is %lf and %lf", grid_name, real_min_lon_value, real_max_lon_value);

    if (!are_floating_values_equal(NULL_COORD_VALUE, min_lat_value)) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, max_lat_value > min_lat_value, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": \"min_lat\" (currently is %lf) is not smaller than \"max_lat\" (currently is %lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, max_lat_value, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) center_lat, size_center_lat, min_lat_value, max_lat_value, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"center_lat\" are not between \"min_lat\" (%lf) and \"max_lat\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, max_lat_value, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) vertex_lat, size_vertex_lat, min_lat_value, max_lat_value, (double) NULL_COORD_VALUE, true), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"vertex_lat\" are not between \"min_lat\" (%lf) and \"max_lat\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, max_lat_value, annotation);
    }
    else {
        min_lat_value = real_min_lat_value;
        max_lat_value = real_max_lat_value;
    }
    if (!(words_are_the_same(coord_unit, COORD_UNIT_DEGREES) && are_floating_values_equal((double) - 90.0, min_lat_value) ||
            words_are_the_same(coord_unit, COORD_UNIT_RADIANS) && are_floating_values_equal((double) - PI / 2, min_lat_value)))
        min_lat_value = real_min_lat_value;
    if (!(words_are_the_same(coord_unit, COORD_UNIT_DEGREES) && are_floating_values_equal((double)90.0, max_lat_value) ||
            words_are_the_same(coord_unit, COORD_UNIT_RADIANS) && are_floating_values_equal((double)PI / 2, max_lat_value)))
        max_lat_value = real_max_lat_value;

    if (words_are_the_same(coord_unit, COORD_UNIT_DEGREES)) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lat_value, 1, (double) - 90.0 * eps, (double) 90.0 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lat\" is wrong (not between -90 and 90). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lat_value, 1, (double) - 90.0 * eps, (double) 90.0 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lat\" is wrong (not between -90 and 90). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lat_value, annotation);
        if (!are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lon_value, 1, (double) - 360.0 * eps, (double) 360.0 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lon\" is wrong (not between -360 and 360). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, annotation);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lon_value, 1, (double) - 360.0 * eps, (double) 360.0 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lon\" is wrong (not between -360 and 360). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lon_value, annotation);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, fabs(max_lon_value - min_lon_value) <= ((double)360.0)*eps, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the difference between \"min_lon\" and \"max_lon\" (%lf) is wrong (not between -360 and 360 degrees). Please check the model code related to the annotation \"%s\".", grid_name, API_label, fabs(max_lon_value - min_lon_value), annotation);
        }
        if (are_floating_values_equal((double) - 90.0, min_lat_value) || are_floating_values_equal((double)90.0, max_lat_value))
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(cyclic_or_acyclic, "cyclic"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the H2D grid must be but currently is not a cyclic grid when it covers north or south pole. Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
        if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
            min_lon_value = -360;
            max_lon_value = 360;
        }
    }
    else if (words_are_the_same(coord_unit, COORD_UNIT_RADIANS)) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lat_value, 1, -((double)PI) / 2 * eps, ((double)PI) / 2 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lat\" are wrong (not between -PI/2 and PI/2). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lat_value, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lat_value, 1, -((double)PI) / 2 * eps, ((double)PI) / 2 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lat\" are wrong (not between -PI/2 and PI/2). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lat_value, annotation);
        if (!are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &min_lon_value, 1, -((double)PI) * 2 * eps, ((double)PI) * 2 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"min_lon\" are wrong (not between -2PI and 2PI). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, annotation);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(DATA_TYPE_DOUBLE, &max_lon_value, 1, -((double)PI) * 2 * eps, ((double)PI) * 2 * eps, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the specified value (%lf) of the parameter \"max_lon\" are wrong (not between -2PI and 2PI). Please check the model code related to the annotation \"%s\".", grid_name, API_label, max_lon_value, annotation);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, fabs(max_lon_value - min_lon_value) <= ((double)PI) * 2 * eps, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the difference between \"min_lon\" and \"max_lon\" is wrong (not between -2PI and 2PI). Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
        }
        if (are_floating_values_equal((double) - PI / 2, min_lat_value) || are_floating_values_equal((double)PI / 2, max_lat_value))
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(cyclic_or_acyclic, "cyclic"), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the H2D grid must be but currently is not a cyclic grid when it covers north or south pole. Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
        if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
            min_lon_value = -((double)PI) * 2;
            max_lon_value = ((double)PI) * 2;
        }
    }
    if (are_floating_values_equal(NULL_COORD_VALUE, min_lon_value)) {
        min_lon_value = real_min_lon_value;
        max_lon_value = real_max_lon_value;
    }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) center_lon, size_center_lon, min_lon_value, max_lon_value, (double) NULL_COORD_VALUE, false), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"center_lon\" are not between \"min_lon\" (%lf) and \"max_lon\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, max_lon_value, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries(data_type, (double*) vertex_lon, size_vertex_lon, min_lon_value, max_lon_value, (double) NULL_COORD_VALUE, true), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": some values of the parameter \"vertex_lon\" are not between \"min_lon\" (%lf) and \"max_lon\" (%lf). Please check the model code related to the annotation \"%s\".", grid_name, API_label, min_lon_value, max_lon_value, annotation);

    if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
        if (words_are_the_same(coord_unit, COORD_UNIT_DEGREES)) {
            min_lon_value = 0;
            max_lon_value = 360;
        }
        else if (words_are_the_same(coord_unit, COORD_UNIT_RADIANS)) {
            min_lon_value = 0;
            max_lon_value = ((double)PI) * 2;
        }
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Software error in Original_grid_mgt::common_checking_for_H2D_registration_via_data");
    }
    else if (min_lon_value < max_lon_value) {
        min_lon_value = real_min_lon_value;
        max_lon_value = real_max_lon_value;
    }
    transform_datatype_of_arrays((char*)(&min_lat_value), min_lat, DATA_TYPE_DOUBLE, data_type, 1);
    transform_datatype_of_arrays((char*)(&max_lat_value), max_lat, DATA_TYPE_DOUBLE, data_type, 1);
    transform_datatype_of_arrays((char*)(&min_lon_value), min_lon, DATA_TYPE_DOUBLE, data_type, 1);
    transform_datatype_of_arrays((char*)(&max_lon_value), max_lon, DATA_TYPE_DOUBLE, data_type, 1);
}


int Original_grid_mgt::register_H2D_grid_empty(int comp_id, const char *grid_name, int grid_size, const char *annotation)
{
    check_API_parameter_int(comp_id, API_ID_GRID_MGT_REG_H2D_GRID_WITHOUT_DATA, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, grid_size, "\"grid_size\"", annotation);
    int grid_id = create_H2D_grid_from_global_data(comp_id, grid_name, "degrees", "acyclic", NULL, grid_size, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, false, annotation);
}


int Original_grid_mgt::register_H2D_grid_via_local_data(int comp_id, const char *grid_name, const char *edge_type, const char *coord_unit, char *cyclic_or_acyclic, const char *data_type, int grid_size, int num_local_cells, int size_local_cells_global_index, int size_center_lon, int size_center_lat,
        int size_mask, int size_area, int size_vertex_lon, int size_vertex_lat, int *local_cells_global_index, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *decomp_name, int *decomp_id, const char *annotation, int API_id)
{
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data");
    char API_label[NAME_STR_SIZE], hint[NAME_STR_SIZE], local_grid_name[NAME_STR_SIZE], local_decomp_name[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
    int global_area_size = -1, global_mask_size = -1;
    int num_local_vertexes = 0, num_global_vertexes, local_mask_size = 0, local_area_size = 0, num_total_cells;
    int local_temp_int_array[4], global_temp_int_array[4];


    get_API_hint(comp_id, API_id, API_label);
    sprintf(hint, "registering an H2D grid \"%s\"", grid_name);

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_cells >= 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the specified number (%d) of local grid cells is wrong (smaller than 0). Please check the model code related to \"%s\"", API_label, grid_name, num_local_cells, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_cells <= size_local_cells_global_index, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size (currently is %d) of \"local_cells_global_index\" is smaller than \"num_local_cells\" (currently is %d). Please check the model code related to \"%s\"", API_label, grid_name, num_local_cells, size_local_cells_global_index, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, *decomp_id == 0 && strlen(decomp_name) > 0 || *decomp_id == -1 && strlen(decomp_name) == 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": parameters \"decomp_name\" and \"decomp_id\" are not speicified/unspecified at the same time. Please check the model code related to the annotation \"%s\".", API_label, grid_name, annotation);
    if (size_mask >= 0) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_cells == size_mask, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size of the parameter \"mask\" must be the same as the value of the parameter \"num_local_cells\" (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, num_local_cells, size_mask, annotation);
        local_mask_size = size_mask;
    }
    if (size_area >= 0) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_cells == size_area, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size of the parameter \"area\" must be the same as the value of the parameter \"num_local_cells\" (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, num_local_cells, size_area, annotation);
        local_area_size = size_area;
    }
    local_temp_int_array[0] = num_local_cells;
    local_temp_int_array[1] = local_mask_size;
    local_temp_int_array[2] = local_area_size;
    MPI_Allreduce(local_temp_int_array, global_temp_int_array, 3, MPI_INT, MPI_SUM, comp_node->get_comm_group());
    num_total_cells = global_temp_int_array[0];
    global_mask_size = global_temp_int_array[1];
    global_area_size = global_temp_int_array[2];

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_total_cells == grid_size, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the total sum of the parameter \"num_local_cells\" among all processes of the corresponding component model is not the same with the parameter \"grid_size\" (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, num_total_cells, grid_size, annotation);

    common_checking_for_H2D_registration_via_data(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, data_type, size_mask, size_center_lon, size_center_lat,
            size_vertex_lon, size_vertex_lat, mask, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, vertex_lon, vertex_lat, annotation, API_id);

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, grid_size > 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the \"grid_size\" (currently is %d) is too small (smaller than 1). Please check the model code related to \"%s\"", API_label, grid_name, grid_size, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_cells == size_center_lon, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size of the parameter \"center_lon\" must be the same as the value of the parameter \"num_local_cells\" (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, num_local_cells, size_center_lon, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_cells == size_center_lat, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size of the parameter \"center_lat\" must be the same as the value of the parameter \"num_local_cells\" (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, num_local_cells, size_center_lat, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, global_mask_size == 0 || global_mask_size == grid_size, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the parameter \"mask\" must be specified (or not specified) at the same time consistently among all processes of the corresponding component model. Please check the model code related to the annotation \"%s\".", API_label, grid_name, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, global_area_size == 0 || global_area_size == grid_size, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the parameter \"area\" must be specified (or not specified) at the same time consistently among all processes of the corresponding component model. Please check the model code related to the annotation \"%s\".", API_label, grid_name, annotation);
    if (global_mask_size == 0)
        global_mask_size = -1;
    if (global_area_size == 0)
        global_area_size = -1;

    if (num_local_cells == 0) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lon <= 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the parameter \"vertex_lon\" should not be specified or its array size should be 0 (currently is %d) when the parameter \"num_local_cells\" is 0. Please check the model code related to the annotation \"%s\".", API_label, grid_name, size_vertex_lon, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lat <= 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the parameter \"vertex_lat\" should not be specified or its array size should be 0 (currently is %d) when the parameter \"num_local_cells\" is 0. Please check the model code related to the annotation \"%s\".", API_label, grid_name, size_vertex_lat, annotation);
    }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lon == size_vertex_lat, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array sizes of the parameters \"vertex_lon\" and \"vertex_lat\" must be the same (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, size_vertex_lon, size_vertex_lat, annotation);
    if (size_vertex_lon >= 0)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lon % num_local_cells == 0, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the array size of the parameter \"vertex_lon\" must be an integer multiple of the parameter \"num_local_cells\" (currently are different: %d and %d respectively). Please check the model code related to the annotation \"%s\".", API_label, grid_name, size_vertex_lon, num_local_cells, annotation);
    if (num_local_cells != 0 && size_vertex_lon >= 0)
        num_local_vertexes = size_vertex_lon / num_local_cells;
    MPI_Allreduce(&num_local_vertexes, &num_global_vertexes, 1, MPI_INT, MPI_MAX, comp_node->get_comm_group());
    if (num_global_vertexes > 0 && num_local_cells > 0)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_local_vertexes == num_global_vertexes, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": the parameters \"vertex_lon\" and \"vertex_lat\" must be specified (or not specified) at the same time consistently among all processes of the corresponding component model. Please check the model code related to the annotation \"%s\".", API_label, grid_name, annotation);

    check_API_parameter_int(comp_id, API_id, comm, NULL, grid_size, "\"grid_size\"", annotation);
    for (int i = 0; i < num_local_cells; i ++)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, local_cells_global_index[i] >= 1 && local_cells_global_index[i] <= grid_size, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": some values in \"local_cells_global_index\" (for example, %d) are out of bound (smaller than 1 or larger than \"grid_size\" (%d)). Please check the model code related to the annotation \"%s\".", API_label, grid_name, local_cells_global_index[i], grid_size, annotation);

    sprintf(local_grid_name, "local_grid_of_%s", grid_name);
    Remap_grid_class *local_grid = generate_H2D_CoR_grid(comp_id, local_grid_name, coord_unit, cyclic_or_acyclic, data_type, num_local_cells, 0, num_global_vertexes, size_center_lon, size_center_lat, size_mask,
                                   size_area, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, false, annotation);
    Remap_grid_class *empty_global_grid = generate_H2D_CoR_grid(comp_id, grid_name, coord_unit, cyclic_or_acyclic, data_type, grid_size, 0, num_global_vertexes, grid_size, grid_size,
                                          global_mask_size, global_area_size, min_lon, max_lon, min_lat, max_lat, NULL, NULL, NULL, NULL, NULL, NULL, true, annotation);
    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, empty_global_grid, true, true, false));
    int grid_id = original_grids[original_grids.size() - 1]->get_grid_id();

    if (strlen(decomp_name) != 0)
        strcpy(local_decomp_name, decomp_name);
    else sprintf(local_decomp_name, "register_decomp_of_%s", grid_name);
    *decomp_id = decomps_info_mgr->register_H2D_parallel_decomposition(local_decomp_name, grid_id, -1, num_local_cells, local_cells_global_index, 0, NULL, annotation);
    Decomp_grid_info *registered_decomp_grid = decomp_grids_mgr->add_decomp_grid_info(*decomp_id, get_original_grid(grid_id)->get_original_CoR_grid(), local_grid);
    Distributed_H2D_grid_engine *distributed_H2D_grid = distributed_H2D_grid_mgr->generate_distributed_H2D_grid_engine(comp_node->get_full_name(), grid_name, empty_global_grid, grid_size, 0, 0, get_original_grid(grid_id), false);

    if (report_error_enabled) {
        Distribute_merge_sort<int> *local_cell_global_index_sort = new Distribute_merge_sort<int>(-1, -1, NULL, NULL);
        common_sort_struct<int> *global_index_sort_structs = NULL;
        int aver_decomp_size = grid_size / comp_node->get_num_procs();
        int num_local_sorted_global_indexes = num_local_cells;
        if (grid_size % comp_node->get_num_procs() > 0)
            aver_decomp_size ++;
        if (num_local_cells > 0)
            global_index_sort_structs = new common_sort_struct<int> [num_local_cells];
        for (int i = 0; i < num_local_cells; i ++) {
            global_index_sort_structs[i].content = local_cells_global_index[i] - 1;
            global_index_sort_structs[i].target_proc_id = global_index_sort_structs[i].content / aver_decomp_size;
        }
        local_cell_global_index_sort->do_data_sorting_with_target_process_id(&global_index_sort_structs, decomps_info_mgr->get_decomp_info(*decomp_id), comp_node, &num_local_sorted_global_indexes);
        for (int i = 0; i < num_local_sorted_global_indexes; i ++)
            global_index_sort_structs[i].key = global_index_sort_structs[i].content;
        do_quick_sort(global_index_sort_structs, (int*)NULL, 0, num_local_sorted_global_indexes - 1);
        for (int i = 0; i < num_local_sorted_global_indexes - 1; i ++)
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, global_index_sort_structs[i].content != global_index_sort_structs[i + 1].content, "Error happens when calling the API \"%s\" to register an H2D grid \"%s\": there are multiple copies of the No. %d cell in the parameter \"local_cells_global_index\" among processes (each grid cell can only have one copy). Please check the model code related to the annotation \"%s\".", API_label, grid_name, global_index_sort_structs[i].content + 1, annotation);
        if (global_index_sort_structs != NULL)
            delete [] global_index_sort_structs;
        delete local_cell_global_index_sort;
    }

    Remap_grid_data_class *basic_decomp_cor_grid_lon_field = distributed_H2D_grid->generate_basic_decomp_grid_from_decomp_grid(comp_node->get_full_name(), registered_decomp_grid, comp_node->get_full_name())->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON);
    distributed_H2D_grid->get_basic_decomp_grid()->get_decomp_grid()->check_coord_values_range();
    for (int i = 0; i < basic_decomp_cor_grid_lon_field->get_grid_data_field()->required_data_size; i ++)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, ((double*)basic_decomp_cor_grid_lon_field->get_grid_data_field()->data_buf)[i] != NULL_COORD_VALUE);
//  distributed_H2D_grid_mgr->generate_vertexes_via_distributed_triangulation(distributed_H2D_grid_mgr->search_distributed_H2D_grid(comp_node->get_full_name(), get_original_grid(grid_id)));
    if (output_H2D_grid_enabled)
        distributed_H2D_grid_mgr->generate_full_grid_data_from_distributed_grid(comp_id, empty_global_grid, distributed_H2D_grid->get_basic_decomp_grid()->get_decomp_grid(), false);
    get_original_grid(grid_id)->calculate_H2D_grid_checksums();
    get_original_grid(grid_id)->write_H2D_grid();

    return grid_id;
}


Remap_grid_class *Original_grid_mgt::generate_H2D_CoR_grid(int comp_id, const char *grid_name, const char *coord_unit, const char *cyclic_or_acyclic, const char *data_type, int dim_size1, int dim_size2, int num_vertex, int size_center_lon,
        int size_center_lat, int size_mask, int size_area, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, bool is_global_grid, const char *annotation)
{
    char true_H2D_grid_name[NAME_STR_SIZE], true_lon_grid_name[NAME_STR_SIZE], true_lat_grid_name[NAME_STR_SIZE];
    Remap_grid_class *CoR_H2D_grid, *CoR_lon_grid, *CoR_lat_grid, *sub_grids[256];
    long grid_size = dim_size2 == 0 ? dim_size1 : dim_size1 * dim_size2;
    double min_lon_value, max_lon_value, min_lat_value, max_lat_value;

    sprintf(true_H2D_grid_name, "%s@%s", grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, annotation)->get_full_name());
    sprintf(true_lon_grid_name, "lon_%s", true_H2D_grid_name);
    sprintf(true_lat_grid_name, "lat_%s", true_H2D_grid_name);
    if (dim_size2 == 0) {
        CoR_lon_grid = new Remap_grid_class(true_lon_grid_name, "lon", coord_unit, cyclic_or_acyclic, 0);
        CoR_lat_grid = new Remap_grid_class(true_lat_grid_name, "lat", coord_unit, cyclic_or_acyclic, 0);
    }
    else {
        CoR_lon_grid = new Remap_grid_class(true_lon_grid_name, "lon", coord_unit, cyclic_or_acyclic, dim_size1);
        CoR_lat_grid = new Remap_grid_class(true_lat_grid_name, "lat", coord_unit, cyclic_or_acyclic, dim_size2);
    }
    sub_grids[0] = CoR_lon_grid;
    sub_grids[1] = CoR_lat_grid;
    CoR_H2D_grid = new Remap_grid_class(true_H2D_grid_name, 2, sub_grids, grid_size);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, size_center_lat >= 0 && size_center_lon >= 0, "Software error in Original_grid_mgt::generate_H2D_CoR_grid");
    if (size_center_lon == 0 && false) {
        CoR_lat_grid->set_empty_grid();
        CoR_lon_grid->set_empty_grid();
        CoR_H2D_grid->set_empty_grid();
    }
    else {
        if (size_mask >= 0 || center_lon == NULL)
            CoR_H2D_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)mask, 0);
        else {
            int *tmp_mask = NULL;
            if (CoR_H2D_grid->get_grid_size() > 0) {
                tmp_mask = new int [CoR_H2D_grid->get_grid_size()];
                for (int i = 0; i < CoR_H2D_grid->get_grid_size(); i ++)
                    tmp_mask[i] = 1;
            }
            CoR_H2D_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)tmp_mask, 0);
            if (tmp_mask != NULL)
                delete [] tmp_mask;
        }
        if (size_area >= 0)
            CoR_H2D_grid->read_grid_data_from_array("area", "area", data_type, (const char*)area, 0);
        if (size_center_lon == grid_size) {
            CoR_H2D_grid->read_grid_data_from_array("center", COORD_LABEL_LON, data_type, (const char*)center_lon, 0);
            CoR_H2D_grid->read_grid_data_from_array("center", COORD_LABEL_LAT, data_type, (const char*)center_lat, 0);
            if (num_vertex > 0) {
                CoR_H2D_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, data_type, (const char*)vertex_lon, num_vertex);
                CoR_H2D_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, data_type, (const char*)vertex_lat, num_vertex);
            }
        }
        else {
            CoR_lon_grid->read_grid_data_from_array("center", COORD_LABEL_LON, data_type, (const char*)center_lon, 0);
            CoR_lat_grid->read_grid_data_from_array("center", COORD_LABEL_LAT, data_type, (const char*)center_lat, 0);
            if (num_vertex > 0) {
                CoR_lon_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, data_type, (const char*)vertex_lon, num_vertex);
                CoR_lat_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, data_type, (const char*)vertex_lat, num_vertex);
            }
        }

        transform_datatype_of_arrays(min_lon, (char*)(&min_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
        transform_datatype_of_arrays(max_lon, (char*)(&max_lon_value), data_type, DATA_TYPE_DOUBLE, 1);
        transform_datatype_of_arrays(min_lat, (char*)(&min_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
        transform_datatype_of_arrays(max_lat, (char*)(&max_lat_value), data_type, DATA_TYPE_DOUBLE, 1);
        if (words_are_the_same(coord_unit, COORD_UNIT_RADIANS)) {
            min_lon_value = min_lon_value * 180 / PI;
            max_lon_value = max_lon_value * 180 / PI;
            min_lat_value = min_lat_value * 180 / PI;
            max_lat_value = max_lat_value * 180 / PI;
        }
        if (words_are_the_same(cyclic_or_acyclic, "cyclic")) {
            min_lon_value = 0.0;
            max_lon_value = 360.0;
        }
        CoR_H2D_grid->set_grid_boundary(min_lon_value, max_lon_value, min_lat_value, max_lat_value);
        CoR_H2D_grid->end_grid_definition_stage(NULL);
    }

    if (is_global_grid) {
        remap_grid_manager->add_remap_grid(CoR_lon_grid);
        remap_grid_manager->add_remap_grid(CoR_lat_grid);
        remap_grid_manager->add_remap_grid(CoR_H2D_grid);
        if (center_lon == NULL) {
            CoR_lon_grid->transform_coord_unit_from_radian_to_degrees();
            CoR_lat_grid->transform_coord_unit_from_radian_to_degrees();
            CoR_H2D_grid->transform_coord_unit_from_radian_to_degrees();
        }
    }

    return CoR_H2D_grid;
}


int Original_grid_mgt::create_H2D_grid_from_global_data(int comp_id, const char *grid_name, const char *coord_unit, const char *cyclic_or_acyclic, const char *data_type, int dim_size1, int dim_size2, int num_vertex, int size_center_lon,
        int size_center_lat, int size_mask, int size_area, int size_vertex_lon, int size_vertex_lat, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, bool generate_basic_decomp_grid, const char *annotation)
{
    Remap_grid_class *CoR_H2D_grid = generate_H2D_CoR_grid(comp_id, grid_name, coord_unit, cyclic_or_acyclic, data_type, dim_size1, dim_size2, num_vertex, size_center_lon, size_center_lat, size_mask,
                                     size_area, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, true, annotation);
    Original_grid_info *original_grid = new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, CoR_H2D_grid, true, true, true);
    original_grids.push_back(original_grid);
    int grid_id = original_grids[original_grids.size() - 1]->get_grid_id();

    if (dim_size2 == 0)
        dim_size1 = 0;

    if (center_lon != NULL) {
        original_grid->write_H2D_grid();
        Distributed_H2D_grid_engine *distributed_H2D_grid = distributed_H2D_grid_mgr->generate_distributed_H2D_grid_engine(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, annotation)->get_full_name(), grid_name, CoR_H2D_grid, CoR_H2D_grid->get_grid_size(), dim_size1, dim_size2, original_grids[original_grids.size() - 1], generate_basic_decomp_grid);
//      if (generate_basic_decomp_grid)
//          distributed_H2D_grid_mgr->generate_vertexes_via_distributed_triangulation(distributed_H2D_grid);
        original_grid->calculate_H2D_grid_checksums();
        CoR_H2D_grid->empty_grid_data();
    }

    return grid_id;
}


void Original_grid_mgt::generate_distributed_H2D_grid_between_comps(const char *src_comp_full_name, Original_grid_info *src_original_grid, const char *dst_comp_full_name, Original_grid_info *dst_original_grid)
{
    Distributed_H2D_grid_engine *src_distributed_H2D_grid = NULL, *dst_distributed_H2D_grid = NULL;
    Comp_comm_group_mgt_node *src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_comp_full_name);
    Comp_comm_group_mgt_node *dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_comp_full_name);
    Decomp_grid_info *src_decomp_grid = NULL, *dst_decomp_grid;
    int has_distributed_H2D_grid = 0, *pointer_has_distributed_H2D_grid = new int [1];
    long array_size = sizeof(int);


    if (src_original_grid != NULL) {
        src_distributed_H2D_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(src_comp_full_name, src_original_grid);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_distributed_H2D_grid != NULL, "Software error in Original_grid_mgt::generate_distributed_H2D_grid_between_comps: %s", src_original_grid->get_grid_name());
        src_decomp_grid = src_distributed_H2D_grid->get_basic_decomp_grid();
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "generate_distributed_H2D_grid_between_comps from \"%s\" (\"%s\") to \"%s\"", src_comp_full_name, src_original_grid->get_original_CoR_grid()->get_grid_name(), dst_comp_full_name);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_decomp_grid != NULL, "Software error in Original_grid_mgt::generate_distributed_H2D_grid_between_comps");
    }
    if (dst_original_grid != NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "generate_distributed_H2D_grid_between_comps from \"%s\" to \"%s\" (\"%s\")", src_comp_full_name, dst_comp_full_name, dst_original_grid->get_original_CoR_grid()->get_grid_name());
        if (distributed_H2D_grid_mgr->search_distributed_H2D_grid(dst_comp_full_name, dst_original_grid) != NULL)
            has_distributed_H2D_grid = 1;
        if (report_error_enabled) {
            int bcasted_has_distributed_H2D_grid = has_distributed_H2D_grid;
            EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(&bcasted_has_distributed_H2D_grid, 1, MPI_INT, 0, dst_comp_node->get_comm_group()) == MPI_SUCCESS);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, bcasted_has_distributed_H2D_grid == has_distributed_H2D_grid, "Software error in Original_grid_mgt::generate_distributed_H2D_grid_between_comps");
        }
    }

    pointer_has_distributed_H2D_grid[0] = has_distributed_H2D_grid;
    transfer_array_from_one_comp_to_another(dst_comp_node->get_current_proc_local_id(), dst_comp_node->get_local_proc_global_id(0), src_comp_node->get_current_proc_local_id(), src_comp_node->get_local_proc_global_id(0), src_comp_node->get_comm_group(), (char**)(&pointer_has_distributed_H2D_grid), array_size);
    has_distributed_H2D_grid = pointer_has_distributed_H2D_grid[0];
    delete [] pointer_has_distributed_H2D_grid;

    if (has_distributed_H2D_grid == 1)
        return;

    if (dst_original_grid != NULL)
        /* num_lons and num_lats are not used yet*/
        dst_distributed_H2D_grid = distributed_H2D_grid_mgr->generate_distributed_H2D_grid_engine(dst_comp_full_name, dst_original_grid->get_grid_name(), dst_original_grid->get_original_CoR_grid(), dst_original_grid->get_original_CoR_grid()->get_grid_size(), 0, 0, dst_original_grid, false);

    dst_decomp_grid = (dst_distributed_H2D_grid != NULL ? dst_distributed_H2D_grid : src_distributed_H2D_grid)->generate_basic_decomp_grid_from_decomp_grid(src_comp_full_name, src_decomp_grid, dst_comp_full_name);
}


int Original_grid_mgt::register_H2D_grid_via_global_data(int comp_id, const char *grid_name, const char *edge_type, const char *coord_unit, char *cyclic_or_acyclic, const char *data_type, int dim_size1, int dim_size2, int size_center_lon, int size_center_lat,
        int size_mask, int size_area, int size_vertex_lon, int size_vertex_lat, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *annotation, int API_id)
{
    int data_type_size, grid_size, num_vertex = 0;
    char true_H2D_grid_name[NAME_STR_SIZE], true_lon_grid_name[NAME_STR_SIZE], true_lat_grid_name[NAME_STR_SIZE];
    char coord_label[NAME_STR_SIZE], coord_name[NAME_STR_SIZE], API_label[NAME_STR_SIZE], hint[NAME_STR_SIZE];
    Remap_grid_class *CoR_H2D_grid, *CoR_lon_grid, *CoR_lat_grid, *sub_grids[256];


    get_API_hint(comp_id, API_id, API_label);
    sprintf(hint, "registering an H2D grid \"%s\"", grid_name);

    common_checking_for_H2D_registration_via_data(comp_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, data_type, size_mask, size_center_lon, size_center_lat,
            size_vertex_lon, size_vertex_lat, mask, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, vertex_lon, vertex_lat, annotation, API_id);

    data_type_size = get_data_type_size(data_type);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, size_center_lon, data_type_size, (const char*)center_lon, "center_lon", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, size_center_lat, data_type_size, (const char*)center_lat, "center_lat", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, size_mask, sizeof(int), (const char*)mask, "mask", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, size_area, data_type_size, (const char*)area, "area", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, size_vertex_lon, data_type_size, (const char*)vertex_lon, "vertex_lon", annotation);
    check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), hint, size_vertex_lat, data_type_size, (const char*)vertex_lat, "vertex_lat", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, dim_size1, "dim_size1", annotation);
    check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), NULL, dim_size2, "dim_size2", annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, dim_size1 > 3, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is %d) of the parameter \"dim_size1\" is wrong. It must be larger than 3. Please check the model code related to the annotation \"%s\".", grid_name, API_label, dim_size1, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, dim_size2 == 0 || dim_size2 > 3, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the value (currently is %d) of the parameter \"dim_size2\" is wrong. It must be 0 or a positive value larger than 3. Please check the model code related to the annotation \"%s\".", grid_name, API_label, dim_size2, annotation);
    if (dim_size2 == 0) {
        grid_size = dim_size1;
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_center_lon == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lon\" is different from the grid size that is determined by \"dim_size1\" (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lon, grid_size, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_center_lat == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lat\" is different from the grid size that is determined by \"dim_size1\" (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lat, grid_size, annotation);
    }
    else {
        grid_size = (dim_size1) * (dim_size2);
        if (size_center_lon == dim_size1) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_center_lat == dim_size2, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lat\" is different from \"dim_size2\" (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lat, dim_size2, annotation);
        }
        else {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_center_lon == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lon\" is different from \"dim_size1\" (currently is %d) and the grid size (currently is %d) that is determined by the product of \"dim_size1\" and \"dim_size2\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lon, dim_size1, grid_size, annotation);
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_center_lat == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"center_lat\" is different from the grid size (currently is %d) that is determined by the product of \"dim_size1\" and \"dim_size2\" (the array size of \"center_lon\" is the same as the grid size). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_center_lat, grid_size, annotation);
        }
    }
    if (size_mask > 0)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_mask == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"mask\" is different from the grid size (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_mask, grid_size, annotation);
    if (size_area > 0)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_area == grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size (currently is %d) of \"area\" is different from the grid size (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_area, grid_size, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, (size_vertex_lon == -1 && size_vertex_lat == -1) || (size_vertex_lon > 0 && size_vertex_lat > 0), "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the optional parameters \"vertex_lon\" and \"vertex_lat\" are not set/unset at the same time. Please check the model code related to the annotation \"%s\".", grid_name, API_label, annotation);
    if (size_vertex_lon > 0) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, (size_vertex_lon) % (size_center_lon) == 0, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size of \"vertex_lon\" (currently is %d) is not an integral multiple of the array size (currently is %d) of \"center_lon\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_vertex_lon, size_center_lon, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, (size_vertex_lat) % (size_center_lat) == 0, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array size of \"vertex_lat\" (currently is %d) is not an integral multiple of the array size (currently is %d) of \"center_lat\". Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_vertex_lat, size_center_lat, annotation);
        int num_vertex_lon = (size_vertex_lon) / (size_center_lon);
        int num_vertex_lat = (size_vertex_lat) / (size_center_lat);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_vertex_lon == num_vertex_lat, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the numbers of vertexes determined by \"vertex_lon\" (currently is %d) and \"vertex_lat\" (currently is %d) respectively are not the same. Please check the model code related to the annotation \"%s\".", grid_name, API_label, num_vertex_lon, num_vertex_lat, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_vertex_lon > grid_size && size_vertex_lat > grid_size || size_vertex_lon < grid_size && size_vertex_lat < grid_size, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the array sizes of \"vertex_lon\" and \"vertex_lat\" (currently are %d and %d respectively) are not bigger/smaller than the grid size at the same time. Please check the model code related to the annotation \"%s\".", grid_name, API_label, size_vertex_lon, size_vertex_lat, annotation);
        if (size_center_lon == grid_size) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_vertex_lon >= 3, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the number of vertexes (currently is %d) is wrong as it is smaller than 3. Please check the model code related to the annotation \"%s\".", grid_name, API_label, num_vertex_lon, annotation);
        }
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id , num_vertex_lon == 2, "Error happens when registering an H2D grid \"%s\" through the API \"%s\": the number of vertexes (currently is %d) is wrong as it is not 2. Please check the model code related to the annotation \"%s\".", grid_name, API_label, num_vertex_lon, annotation);
        num_vertex = num_vertex_lon;
    }

    return create_H2D_grid_from_global_data(comp_id, grid_name, coord_unit, cyclic_or_acyclic, data_type, dim_size1, dim_size2, num_vertex, size_center_lon, size_center_lat,
                                            size_mask, size_area, size_vertex_lon, size_vertex_lat, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, true, annotation);
}


int Original_grid_mgt::register_H2D_grid_via_file(int comp_id, const char *grid_name, const char *data_file_name, const char *annotation)
{
    int rcode, ncfile_id, grid_id;
    int size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat;
    long dim_lon_size, dim_lat_size, dim_H2D_size, dim_size1, dim_size2;
    char *center_lon, *center_lat, *vertex_lon, *vertex_lat, *area;
    char min_lon[NAME_STR_SIZE], max_lon[NAME_STR_SIZE], min_lat[NAME_STR_SIZE], max_lat[NAME_STR_SIZE];
    int *mask;
    char data_type_for_center_lat[NAME_STR_SIZE], data_type_for_center_lon[NAME_STR_SIZE], data_type_for_vertex_lon[NAME_STR_SIZE], data_type_for_vertex_lat[NAME_STR_SIZE], data_type_for_mask[NAME_STR_SIZE], data_type_for_area[NAME_STR_SIZE];
    char data_type_temp[NAME_STR_SIZE];
    char edge_type[NAME_STR_SIZE], cyclic_or_acyclic[NAME_STR_SIZE], unit_center_lon[NAME_STR_SIZE], unit_center_lat[NAME_STR_SIZE], unit_vertex_lon[NAME_STR_SIZE], unit_vertex_lat[NAME_STR_SIZE];
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_H2D_grid_via_file");
    bool is_root_proc = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "in register_H2D_grid_via_file") == 0;


    check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE, comm, "registering an H2D grid", data_file_name, "data_file_name", annotation);

    IO_netcdf *netcdf_file_object = new IO_netcdf("H2D_grid_data", data_file_name, "r", false);
    dim_lon_size = netcdf_file_object->get_dimension_size(COORD_LABEL_LON, comm, is_root_proc);
    dim_lat_size = netcdf_file_object->get_dimension_size(COORD_LABEL_LAT, comm, is_root_proc);
    dim_H2D_size = netcdf_file_object->get_dimension_size("grid_size", comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_CENTER_LON_LABEL, (void**)(&center_lon), &size_center_lon, data_type_for_center_lon, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_CENTER_LAT_LABEL, (void**)(&center_lat), &size_center_lat, data_type_for_center_lat, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_VERTEX_LON_LABEL, (void**)(&vertex_lon), &size_vertex_lon, data_type_for_vertex_lon, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_VERTEX_LAT_LABEL, (void**)(&vertex_lat), &size_vertex_lat, data_type_for_vertex_lat, comm, is_root_proc);
    netcdf_file_object->read_file_field("area", (void**)(&area), &size_area, data_type_for_area, comm, is_root_proc);
    netcdf_file_object->read_file_field(SCRIP_MASK_LABEL, (void**)(&mask), &size_mask, data_type_for_mask, comm, is_root_proc);
    if (dim_lon_size > 0 && dim_lat_size > 0 && dim_H2D_size > 0)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, dim_H2D_size == dim_lon_size * dim_lat_size, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the size (currently is %d) of dimension \"grid_size\" is different from the multiple (currently is %d) of sizes of dimensions \"lon\" and \"lat\"", grid_name, annotation, dim_H2D_size, dim_lon_size * dim_lat_size, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, dim_H2D_size > 0 || (dim_lon_size > 0 && dim_lat_size > 0), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: the dimension size (dimensions \"lon\" and \"lat\" in the file) or the grid size (dimension \"grid_size\" in the file) is not correctly specified in the file \"%s\". Please verify.", grid_name, annotation, data_file_name);
    if (dim_lon_size > 0 && dim_lat_size > 0) {
        dim_size1 = dim_lon_size;
        dim_size2 = dim_lat_size;
    }
    else {
        dim_size1 = dim_H2D_size;
        dim_size2 = 0;
    }

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, center_lon != NULL, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: the longitude value for the center of each grid point (variable \"lon\" in the file) is not specified in the data file \"%s\". ",
                                      grid_name, annotation, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, center_lat != NULL, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: the latitude value for the center of each grid point (variable \"lat\" in the file) is not specified in the data file \"%s\". ",
                                      grid_name, annotation, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, vertex_lon != NULL && vertex_lat != NULL || vertex_lon == NULL && vertex_lat == NULL, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the longitude and latitude values for each vertex (variables \"vertex_lon\" and \"vertex_lat\" in the file) of each grid point must be specified/unspecified at the same time",
                                      grid_name, annotation, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_center_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variables \"lon\" and \"lat\" are not the same", grid_name, annotation, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon, DATA_TYPE_FLOAT) || words_are_the_same(data_type_for_center_lon, DATA_TYPE_DOUBLE), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variables \"lon\" is not floating-point", grid_name, annotation, data_file_name);
    if (vertex_lon != NULL) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_vertex_lon), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"vertex_lon\" is different from the data type of variable \"lon\".", grid_name, annotation, data_file_name);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_vertex_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"vertex_lat\" is different from the data type of variable \"lat\".", grid_name, annotation, data_file_name);
    }
    if (area != NULL)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_area), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"area\" is different from the data type of variable \"lon\".", grid_name, annotation, data_file_name);
    if (mask != NULL)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(data_type_for_mask, DATA_TYPE_INT), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the data type of variable \"mask\" is not \"integer\".", grid_name, annotation, data_file_name);

    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "edge_type", edge_type, data_type_temp, comm, is_root_proc) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the global attribute \"edge_type\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "cyclic_or_acyclic", cyclic_or_acyclic, data_type_temp, comm, is_root_proc) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the global attribute \"cyclic_or_acyclic\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, strlen(edge_type) > 0, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", \"edge_type\" is not specified as a global attribute", grid_name, annotation, data_file_name);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, strlen(cyclic_or_acyclic) > 0, "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", \"cyclic_or_acyclic\" is not specified as a global attribute", grid_name, annotation, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, (netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LON_LABEL, "unit", unit_center_lon, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LON_LABEL, "units", unit_center_lon, data_type_temp, comm, is_root_proc)) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"%s\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, SCRIP_CENTER_LON_LABEL, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, (netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LAT_LABEL, "unit", unit_center_lat, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_CENTER_LAT_LABEL, "units", unit_center_lat, data_type_temp, comm, is_root_proc)) && words_are_the_same(data_type_temp, DATA_TYPE_STRING), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"%s\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, SCRIP_CENTER_LAT_LABEL, data_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "min_lon", min_lon, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"min_lon\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(min_lon, min_lon, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "max_lon", max_lon, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"max_lon\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(max_lon, max_lon, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "min_lat", min_lat, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"min_lat\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(min_lat, min_lat, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(NULL, "max_lat", max_lat, data_type_temp, comm, is_root_proc) && (words_are_the_same(data_type_temp, DATA_TYPE_FLOAT) || words_are_the_same(data_type_temp, DATA_TYPE_DOUBLE)), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the value of the global attribute \"max_lat\" from the data file \"%s\": it does not exist or its type is not float or double", grid_name, annotation, data_file_name);
    transform_datatype_of_arrays(max_lat, max_lat, data_type_temp, data_type_for_center_lon, 1);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(unit_center_lon, unit_center_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the units of variables \"lon\" and \"lat\" are different", grid_name, annotation, data_file_name);
    if (vertex_lon != NULL) {
        EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LON_LABEL, "unit", unit_vertex_lon, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LON_LABEL, "units", unit_vertex_lon, data_type_temp, comm, is_root_proc), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"vertex_lon\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LAT_LABEL, "unit", unit_vertex_lat, data_type_temp, comm, is_root_proc) || netcdf_file_object->get_file_field_string_attribute(SCRIP_VERTEX_LAT_LABEL, "units", unit_vertex_lat, data_type_temp, comm, is_root_proc), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: fail to get the unit of variable \"vertex_lat\" from the data file \"%s\": it does not exist or its type is not string", grid_name, annotation, data_file_name);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(unit_center_lon, unit_vertex_lon), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the units of variables \"lon\" and \"vertex_lon\" are different", grid_name, annotation, data_file_name);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(unit_center_lat, unit_vertex_lat), "Error happens when registering an H2D grid \"%s\" (the corresponding model code annotation is \"%s\") through the API CCPL_register_H2D_grid_via_data_file: in the data file \"%s\", the units of variables \"lat\" and \"vertex_lat\" are different", grid_name, annotation, data_file_name);
    }

    grid_id = register_H2D_grid_via_global_data(comp_id, grid_name, edge_type, unit_center_lon, cyclic_or_acyclic, data_type_for_center_lon, dim_size1, dim_size2, size_center_lon, size_center_lat,
              size_mask, size_area, size_vertex_lon, size_vertex_lat, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation, API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE);

    delete [] center_lon;
    delete [] center_lat;
    if (vertex_lon != NULL) {
        delete [] vertex_lon;
        delete [] vertex_lat;
    }
    if (mask != NULL)
        delete [] mask;
    if (area != NULL)
        delete [] area;

    return grid_id;
}


int Original_grid_mgt::register_V1D_grid_via_data(int API_id, int comp_id, const char *grid_name, int grid_type, const char *coord_unit, int grid_size,
        double value1, const double *value2, const double *value3, bool enabling_check, const char *annotation)
{
    char full_grid_name[NAME_STR_SIZE];
    Remap_grid_class *CoR_V1D_grid;


    if (enabling_check) {
        check_API_parameter_int(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_V1D_grid_via_data"), NULL, grid_size, "implicit grid size", annotation);
        check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_V1D_grid_via_data"), "registering a V1D grid", 1, sizeof(double), (const char*)(&value1), "floating-point parameters", annotation);
        if (value2 != NULL)
            check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_V1D_grid_via_data"), "registering a V1D grid", grid_size, sizeof(double), (const char*)(value2), "floating-point parameters", annotation);
        if (value3 != NULL)
            check_API_parameter_data_array(comp_id, API_id, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_V1D_grid_via_data"), "registering a V1D grid", grid_size, sizeof(double), (const char*)(value3), "floating-point parameters", annotation);
    }

    sprintf(full_grid_name, "%s@%s", grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "in register_V1D_grid_via_data")->get_full_name());
    if (grid_type == 0)
        CoR_V1D_grid = new Remap_grid_class(full_grid_name, COORD_LABEL_TRACER, "", NULL, grid_size);
    else if (grid_type <= 4)
        CoR_V1D_grid = new Remap_grid_class(full_grid_name, COORD_LABEL_LEV, coord_unit, NULL, grid_size);
    switch (grid_type) {
    case 1:
        CoR_V1D_grid->read_grid_data_from_array("center", COORD_LABEL_LEV, DATA_TYPE_DOUBLE, (const char*)value2, 0);
        break;
    case 2:
        CoR_V1D_grid->set_lev_grid_sigma_info(value1, value2, NULL, 1.0);
        break;
    case 3:
        CoR_V1D_grid->set_lev_grid_sigma_info(value1, value2, value3, 1.0);
        break;
    default:
        break;
    }

    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, CoR_V1D_grid, true, true, false));

    remap_grid_manager->add_remap_grid(CoR_V1D_grid);

    return original_grids[original_grids.size() - 1]->get_grid_id();
}


int Original_grid_mgt::register_md_grid_via_multi_grids(int comp_id, const char *grid_name, int sub_grid1_id, int sub_grid2_id, int sub_grid3_id, int size_mask, int *mask, bool enabling_check, const char *annotation)
{
    int num_sub_grids = 2, num_H2D_grid = 0, num_V1D_grid = 0, num_T1D_grid = 0, num_tracer_grid = 0, i;
    char full_grid_name[NAME_STR_SIZE];
    Remap_grid_class *CoR_MD_grid, *sub_grids[3];


    MPI_Comm local_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Original_grid_mgt::register_md_grid_via_multi_grids");
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(sub_grid1_id), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"sub_grid1_id\" (currently 0x%x) is wrong. Please check the model code related to the annotation \"%s\".", grid_name, sub_grid1_id, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(sub_grid1_id) == comp_id, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the component model corresponding to the grid with id of \"sub_grid1_id\" is different from the current component with the id of comp_id. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(sub_grid2_id), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"sub_grid2_id\" is wrong. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(sub_grid2_id) == comp_id, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the component corresponding to the grid with id of \"sub_grid2_id\" is different from the current component with the id of comp_id. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    sub_grids[0] = original_grid_mgr->get_original_CoR_grid(sub_grid1_id);
    sub_grids[1] = original_grid_mgr->get_original_CoR_grid(sub_grid2_id);
    original_grid_mgr->get_original_grid(sub_grid1_id)->set_used_in_md_grid();
    original_grid_mgr->get_original_grid(sub_grid2_id)->set_used_in_md_grid();
    int have_sub_grid3 = sub_grid3_id != -1 ? 1 : 0;
    if (enabling_check) {
        check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, "for registering a multi-dimension grid", original_grid_mgr->get_name_of_grid(sub_grid1_id), "\"sub_grid1_id\" (the corresponding grid name)", annotation);
        check_API_parameter_data_array(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in register_h2d_grid_with_data"), "registering a multi-dimension grid", size_mask, sizeof(int), (const char*)mask, "mask", annotation);
        check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, "for registering a multi-dimension grid", original_grid_mgr->get_name_of_grid(sub_grid2_id), "\"sub_grid2_id\" (the corresponding grid name)", annotation);
        check_API_parameter_int(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, NULL, have_sub_grid3, "\"sub_grid3_id\" (specified or not)", annotation);
    }

    if (sub_grid3_id != -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(sub_grid3_id), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"sub_grid3_id\" is wrong. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
        if (enabling_check)
            check_API_parameter_string(comp_id, API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS, local_comm, "for registering a multi-dimension grid", original_grid_mgr->get_name_of_grid(sub_grid3_id), "\"sub_grid3_id\" (the corresponding grid name)", annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->get_comp_id_of_grid(sub_grid3_id) == comp_id, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the component corresponding to the grid with id of \"sub_grid3_id\" is different from the current component with the id of comp_id. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
        sub_grids[num_sub_grids++] = original_grid_mgr->get_original_CoR_grid(sub_grid3_id);
        original_grid_mgr->get_original_grid(sub_grid3_id)->set_used_in_md_grid();
    }

    for (i = 0; i < num_sub_grids; i ++)
        if (sub_grids[i]->get_is_H2D_grid())
            num_H2D_grid ++;
        else if (sub_grids[i]->has_grid_coord_label(COORD_LABEL_LEV))
            num_V1D_grid ++;
        else if (sub_grids[i]->has_grid_coord_label(COORD_LABEL_TIME))
            num_T1D_grid ++;
        else if (sub_grids[i]->has_grid_coord_label(COORD_LABEL_TRACER))
            num_tracer_grid ++;
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": grid \"%s\" cannot be used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, sub_grids[i]->get_grid_name(), annotation);

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_H2D_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one H2D grid are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_V1D_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one V1D grid are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_T1D_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one T1D grid are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, num_tracer_grid <= 1, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": more than one normal 1D grid (not time grid or vertical grid) are used to generate a multi-dimension grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);

    sprintf(full_grid_name, "%s@%s", grid_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "register_md_grid_via_multi_grids")->get_full_name());
    CoR_MD_grid = new Remap_grid_class(full_grid_name, num_sub_grids, sub_grids, 0);
    if (CoR_MD_grid->has_V3D_sub_grid() && CoR_MD_grid->get_num_dimensions() > 3) {
        for (i = 0; i < num_sub_grids; i ++)
            if (sub_grids[i]->has_V3D_sub_grid())
                break;
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, i < num_sub_grids, "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": This multi-dimensional grid has a V3D sub grid, while no input sub grid is a V3D grid or has a V3D sub grid. In such a case, a V3D sub grid should be registered in advance. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
    }
    if (size_mask > 0)    {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !(CoR_MD_grid->has_grid_coord_label(COORD_LABEL_TIME) || CoR_MD_grid->has_grid_coord_label(COORD_LABEL_TRACER)), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": \"mask\" should not be given a time grid or a normal 1D grid is used for generating a multi-dimensional grid. Please check the model code related to the annotation \"%s\".", grid_name, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, size_mask == CoR_MD_grid->get_grid_size(), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": the array size (currently is %d) of \"mask\" is different from the grid size (currently is %d). Please check the model code related to the annotation \"%s\".", grid_name, size_mask, CoR_MD_grid->get_grid_size(), annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, are_array_values_between_boundaries("integer", mask, size_mask, 0, 1, 0, false), "Error happens when calling the API \"CCPL_register_MD_grid_via_multi_grids\" to register a multi-dimension grid \"%s\": some values of the parameter \"mask\" are wrong (not 0 and 1). Please check the model code related to the annotation \"%s\".", grid_name, annotation);
        CoR_MD_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)mask, 0);
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "the mask field for multi-dimensional grid is not fully supported in C-Coupler. Please contact Dr. Li Liu for help. ");
    }

    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, CoR_MD_grid, true, true, false));
    original_grids[original_grids.size() - 1]->set_md_grid();
    original_grids[original_grids.size() - 1]->write_H2D_grid();
    remap_grid_manager->add_remap_grid(CoR_MD_grid);
    original_grid_mgr->get_original_grid(sub_grid1_id)->copy_bottom_field_variation_type(original_grids[original_grids.size() - 1]);
    original_grid_mgr->get_original_grid(sub_grid2_id)->copy_bottom_field_variation_type(original_grids[original_grids.size() - 1]);
    if (have_sub_grid3)
        original_grid_mgr->get_original_grid(sub_grid3_id)->copy_bottom_field_variation_type(original_grids[original_grids.size() - 1]);
    original_grids[original_grids.size() - 1]->add_sub_grid_id(sub_grid1_id);
    original_grids[original_grids.size() - 1]->add_sub_grid_id(sub_grid2_id);
    if (have_sub_grid3)
        original_grids[original_grids.size() - 1]->add_sub_grid_id(sub_grid3_id);
    if (enabling_check && original_grids[original_grids.size() - 1]->get_H2D_sub_grid() != NULL)
        original_grids[original_grids.size() - 1]->set_grid_checksum(original_grids[original_grids.size() - 1]->get_H2D_sub_grid()->get_checksum_H2D_mask(), original_grids[original_grids.size() - 1]->get_H2D_sub_grid()->get_checksum_H2D_center_lon(), original_grids[original_grids.size() - 1]->get_H2D_sub_grid()->get_checksum_H2D_center_lat());
    return original_grids[original_grids.size() - 1]->get_grid_id();
}


void Original_grid_mgt::set_3D_grid_3D_vertical_coord_field_inst(int grid_id, int field_id, const char *static_or_dynamic, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    Original_grid_info *original_grid;
    Field_mem_info *field_inst;
    int comp_id;
    MPI_Comm local_comm;


    get_API_hint(-1, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, API_label);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id), "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of a 3-D grid: the parameter of \"grid_id\" is wrong. Please verify the model code with the annotation \"%s\".", annotation);
    original_grid = original_grid_mgr->get_original_grid(grid_id);
    comp_id = original_grid->get_comp_id();
    local_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Original_grid_mgt::set_3D_grid_3D_vertical_coord_field_inst");
    check_API_parameter_string(comp_id, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, local_comm, "setting the 3-D vertical coordinate field of a 3-D grid", original_grid->get_grid_name(), "\"grid_id\" (the name of the 3-D grid)", annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid->is_V3D_grid(), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": this grid is not a 3-D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !original_grid->get_original_CoR_grid()->is_sigma_grid(), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": cannot set a 3-D vertical field to this grid because its V1D sub grid is a SIGMA or HYBRID grid. Please verify the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, !original_grid->get_original_CoR_grid()->does_use_V3D_level_coord(), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": cannot set a 3-D vertical field to this grid because its vertical coordinate values have already been set. Please verify the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, memory_manager->check_is_legal_field_instance_id(field_id), "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the grid \"%s\": the parameter of \"field_id\" is wrong. Please verify the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    check_API_parameter_field_instance(comp_id, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, local_comm, "setting the 3-D vertical coordinate field of a 3-D grid", field_id, "\"field_id\" (the corresponding 3-D field)", annotation);
    field_inst = memory_manager->get_field_instance(field_id);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, field_inst->get_grid_id() == grid_id, "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": the given field instance is not on the given 3-D grid. Please verify the model code with the annotation \"%s\".", original_grid->get_grid_name(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(original_grid->get_original_CoR_grid()->get_a_leaf_grid(COORD_LABEL_LEV)->get_coord_unit(), field_inst->get_unit()), "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": The unit (\"%s\") of the grid coordinate is different from the unit (\"%s\") of the given field instance (\"%s\"). Please verify the model code with the annotation \"%s\".", original_grid->get_grid_name(), original_grid->get_original_CoR_grid()->get_a_leaf_grid(COORD_LABEL_LEV)->get_coord_unit(), field_inst->get_unit(), field_inst->get_field_name(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(static_or_dynamic, "constant") || words_are_the_same(static_or_dynamic, "variable"), "Error happens when calling the API \"CCPL_set_3D_grid_3D_vertical_coord_field\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": the input parameter of \"label\" (the current value is \"%s\") must be \"constant\" or \"variable\". Please verify the model code with the annotation \"%s\".", original_grid->get_grid_name(), static_or_dynamic, annotation);
    check_API_parameter_string(comp_id, API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD, local_comm, "setting the 3-D vertical coordinate field of a 3-D grid", original_grid->get_grid_name(), "\"label\" (specification for constant or variable of 3-D field)", annotation);
    if (original_grid->get_V3D_lev_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"%s\" to set the 3-D vertical coordinate field of the 3-D grid \"%s\": the 3-D vertical coordinate field has been set before at the model code with the annotation \"%s\" and cannot be set again at the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation_mgr->get_annotation(grid_id, "set 3-D level field"), annotation);
    original_grid->set_unique_3D_lev_field(field_id, static_or_dynamic, annotation);
    original_grid->get_original_CoR_grid()->set_using_V3D_level_coord();
    decomp_grids_mgr->set_decomp_grids_using_3D_level_coord(original_grid->get_original_CoR_grid());
}


void Original_grid_mgt::set_3d_grid_bottom_field(int comp_id, int grid_id, int field_id, int static_or_dynamic_or_external, int API_id, const char *API_label, const char *annotation)
{
    Original_grid_info *original_grid = original_grid_mgr->get_original_grid(grid_id);
    Original_grid_info *field_original_grid;
    MPI_Comm local_comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Original_grid_mgt::set_3d_grid_bottom_field");


    synchronize_comp_processes_for_API(comp_id, API_id, local_comm, "setting the surface field of a 3-D grid", annotation);
    check_API_parameter_string(comp_id, API_id, local_comm, "setting the surface field of a 3-D grid", original_grid->get_grid_name(), "\"grid_id\" (the name of the 3-D grid)", annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid->is_V3D_grid(), "Error happens when calling the API \"%s\" to set the surface field of a 3-D grid \"%s\": this grid is not a 3-D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid->get_original_CoR_grid()->is_sigma_grid(), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": cannot set the surface field to this grid because its V1D sub grid is not a SIGMA or HYBRID grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation);
    if (original_grid->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the surface field has been set before at the model code with the annotation \"%s\" and cannot be set again at the model code with the annotation \"%s\".", API_label, original_grid->get_grid_name(), annotation_mgr->get_annotation(grid_id, "set surface field"), annotation);
    if (static_or_dynamic_or_external != BOTTOM_FIELD_VARIATION_EXTERNAL) {
        check_API_parameter_field_instance(comp_id, API_id, local_comm, "setting the surface field of a 3-D grid", field_id, "\"field_id\" (the corresponding surface field)", annotation);
        Field_mem_info *field_inst = memory_manager->get_field_instance(field_id);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, field_inst->get_grid_id() != -1, "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the surface field \"%s\" is not on a grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
        field_original_grid = original_grid_mgr->get_original_grid(field_inst->get_grid_id());
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, field_original_grid->is_H2D_grid(), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the grid \"%s\" corresponding to the surface field \"%s\" is not a H2D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, field_original_grid->get_original_CoR_grid()->is_subset_of_grid(original_grid->get_original_CoR_grid()), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the grid \"%s\" corresponding to the surface field \"%s\" is not a sub grid of the 3-D grid. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, words_are_the_same(field_inst->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_FLOAT) || words_are_the_same(field_inst->get_field_data()->get_grid_data_field()->data_type_in_application, DATA_TYPE_DOUBLE), "Error happens when calling the API \"%s\" to set the surface field of the 3-D grid \"%s\": the data type of the surface field \"%s\" is not floating-point. Please check the model code related to the annotation \"%s\".", API_label, original_grid->get_grid_name(), field_inst->get_field_name(), annotation);
    }
    original_grid->set_unique_bottom_field(field_id, static_or_dynamic_or_external, annotation);
    if (original_grid->get_mid_point_grid() != NULL)
        original_grid->get_mid_point_grid()->set_unique_bottom_field(field_id, static_or_dynamic_or_external, annotation);
    if (original_grid->get_interface_level_grid() != NULL)
        original_grid->get_interface_level_grid()->set_unique_bottom_field(field_id, static_or_dynamic_or_external, annotation);
}


int Original_grid_mgt::get_CoR_defined_grid(int comp_id, const char *grid_name, const char *CoR_grid_name, const char *annotation)
{
    Original_grid_info *original_grid, *original_H2D_sub_grid = NULL;
    Remap_grid_class *original_CoR_grid, *script_CoR_grid;
    char CoR_script_file_name[NAME_STR_SIZE];


    sprintf(CoR_script_file_name, "%s/CCPL_grid.cor", comp_comm_group_mgt_mgr->get_root_comp_config_dir());
    script_CoR_grid = remap_grid_manager->search_remap_grid_with_grid_name(CoR_grid_name);
    if (script_CoR_grid == NULL)
        if (strlen(CoR_script_name) > 0) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not define the corresponding CoR grid. Please check the CoR script or the model code corresponding to the annotation \"%s\"",
                                              grid_name, CoR_grid_name, CoR_script_name, annotation);
        }
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, false, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not exist. Please verify",
                                                   grid_name, CoR_grid_name, CoR_script_file_name);
    EXECUTION_REPORT(REPORT_ERROR, comp_id, script_CoR_grid->format_sub_grids(script_CoR_grid), "Please modify the definition of grid \"%s\" in the CoR script \"%s\". We propose to order the dimensions of the grid into the order such as lon, lat, level and time");
    if (script_CoR_grid->get_is_H2D_grid() && are_floating_values_equal(NULL_COORD_VALUE, script_CoR_grid->get_boundary_min_lon()))
        script_CoR_grid->set_grid_boundary(0, 360.0, -90.0, 90.0);

    script_CoR_grid->end_grid_definition_stage(NULL);
    original_CoR_grid = script_CoR_grid;
//  original_CoR_grid = script_CoR_grid->duplicate_grid(script_CoR_grid);
    original_CoR_grid->end_grid_definition_stage(NULL);
    original_grid = new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, annotation, original_CoR_grid, true, true, true);
    original_grids.push_back(original_grid);

    if (original_grid->get_H2D_sub_CoR_grid() != NULL) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_CoR_grid->get_grid_center_field(COORD_LABEL_LON) != NULL, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not specify the center values of longitude (X) of the CoR grid. Please verify.", grid_name, CoR_grid_name, CoR_script_file_name);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_CoR_grid->get_grid_center_field(COORD_LABEL_LAT) != NULL, "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not specify the center values of latitude (Y) of the CoR grid. Please verify.", grid_name, CoR_grid_name, CoR_script_file_name);
    }
    if (original_grid->get_V1D_sub_CoR_grid() != NULL) {
        Remap_grid_data_class *level_field = original_CoR_grid->get_grid_center_field(COORD_LABEL_LEV);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, level_field != NULL && level_field->get_coord_value_grid()->get_grid_size() == original_grid->get_V1D_sub_CoR_grid()->get_grid_size(),
                                          "Error happens when calling the API \"CCPL_register_CoR_defined_grid\" to register grid \"%s\" based on the CoR grid \"%s\": the CoR script \"%s\" does not specify the center values of the vertical coordinate of the CoR grid or the vertical coordinate is not a Z grid. Please verify.", grid_name, CoR_grid_name, CoR_script_file_name);
    }
    if (original_grid->get_H2D_sub_CoR_grid() != NULL) {
        if (original_CoR_grid->get_is_H2D_grid())
            original_H2D_sub_grid = original_grid;
        else original_H2D_sub_grid = search_grid_info(add_original_grid(comp_id, original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), original_grid->get_H2D_sub_CoR_grid()));
        if (distributed_H2D_grid_mgr->search_distributed_H2D_grid(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, annotation)->get_full_name(), original_H2D_sub_grid) == NULL) {
            original_H2D_sub_grid->get_original_CoR_grid()->remove_vertexes_info();
            Distributed_H2D_grid_engine *distributed_H2D_grid = distributed_H2D_grid_mgr->generate_distributed_H2D_grid_engine(comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, annotation)->get_full_name(), original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), original_grid->get_H2D_sub_CoR_grid(), original_grid->get_H2D_sub_CoR_grid()->get_grid_size(), 0, 0, original_H2D_sub_grid, true);
//          distributed_H2D_grid_mgr->generate_vertexes_via_distributed_triangulation(distributed_H2D_grid);
        }
        original_H2D_sub_grid->calculate_H2D_grid_checksums();
        original_grid->set_grid_checksum(original_H2D_sub_grid->get_checksum_H2D_mask(), original_H2D_sub_grid->get_checksum_H2D_center_lon(), original_H2D_sub_grid->get_checksum_H2D_center_lat());
    }
    original_grid->link_sub_grids();
    original_grid->write_H2D_grid();

//  if (original_grid->get_H2D_sub_CoR_grid() != NULL)
//      original_grid->get_H2D_sub_CoR_grid()->empty_grid_data();

    return original_grid->get_grid_id();
}


bool Original_grid_mgt::is_grid_id_legal(int grid_id) const
{
    int true_grid_id = grid_id & TYPE_ID_SUFFIX_MASK;


    if ((grid_id & TYPE_ID_PREFIX_MASK) != TYPE_GRID_LOCAL_ID_PREFIX)
        return false;

    if (true_grid_id < 0 || true_grid_id >= original_grids.size())
        return false;

    return true;
}


int Original_grid_mgt::get_comp_id_of_grid(int grid_id) const
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id & TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_comp_id_of_grid");
    return original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_comp_id();
}


const char *Original_grid_mgt::get_name_of_grid(int grid_id) const
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id & TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_name_of_grid");
    return original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_grid_name();
}


Remap_grid_class *Original_grid_mgt::get_original_CoR_grid(int grid_id) const
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id & TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_original_CoR_grid");
    return original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_original_CoR_grid();
}


Original_grid_info *Original_grid_mgt::get_original_grid(int grid_id) const
{
    if (!is_grid_id_legal(grid_id))
        return NULL;

    return original_grids[grid_id & TYPE_ID_SUFFIX_MASK];
}

int Original_grid_mgt::get_grid_size(int grid_id, const char *annotation) const
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id), "Error happens when getting the size of a grid: the parameter of grid ID is wrong. Please verify the model code corresponding to the annotation \"%s\"", annotation);
    return get_original_CoR_grid(grid_id)->get_grid_size();
}


int Original_grid_mgt::get_grid_id(int comp_id, const char *grid_name, const char *annotation)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->is_legal_local_comp_id(comp_id, false), "Error happens when getting the ID of a grid: the parameter of component ID is wrong. Please verify the model code corresponding to the annotation \"%s\"", annotation);
    Original_grid_info *original_grid = search_grid_info(grid_name, comp_id);
    if (original_grid == NULL)
        return -1;
    else return original_grid->get_grid_id();
}


int Original_grid_mgt::add_original_grid(int comp_id, const char *grid_name, Remap_grid_class *original_CoR_grid)
{
    Original_grid_info *existing_grid = search_grid_info(grid_name, comp_id);
    if (existing_grid != NULL) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_CoR_grid == existing_grid->get_original_CoR_grid(), "Software error in Original_grid_mgt::add_original_grid");
        return existing_grid->get_grid_id();
    }

    original_grids.push_back(new Original_grid_info(comp_id, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, grid_name, "Original_grid_mgt::add_original_grid", original_CoR_grid, false, true, false));
    return original_grids[original_grids.size() - 1]->get_grid_id();
}


int Original_grid_mgt::get_total_grid_size_beyond_H2D(int grid_id)
{
    int total_grid_size_beyond_H2D = 1;

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id & TYPE_ID_SUFFIX_MASK] != NULL, "Software error in Original_grid_mgt::get_total_grid_size_beyond_H2D: wrong grid id");
    if (original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_V1D_sub_CoR_grid() != NULL)
        total_grid_size_beyond_H2D *= original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_V1D_sub_CoR_grid()->get_grid_size();
    if (original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_Time1D_sub_CoR_grid() != NULL)
        total_grid_size_beyond_H2D *= original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_Time1D_sub_CoR_grid()->get_grid_size();
    if (original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_Tracer1D_sub_grid() != NULL)
        total_grid_size_beyond_H2D *= original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->get_Tracer1D_sub_grid()->get_original_CoR_grid()->get_grid_size();

    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, total_grid_size_beyond_H2D > 0, "Software error in Original_grid_mgt::get_total_grid_size_beyond_H2D: wrong size of vertical grid");
    return total_grid_size_beyond_H2D;
}


bool Original_grid_mgt::is_V1D_sub_grid_after_H2D_sub_grid(int grid_id)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_grid_id_legal(grid_id) && original_grids[grid_id & TYPE_ID_SUFFIX_MASK] != NULL, "Original_grid_info::is_V1D_sub_grid_after_H2D_sub_grid: wrong grid id");
    return original_grids[grid_id & TYPE_ID_SUFFIX_MASK]->is_V1D_sub_grid_after_H2D_sub_grid();
}


void Original_grid_mgt::register_mid_point_grid(int level_3D_grid_id, int *mid_3D_grid_id, int *mid_1D_grid_id, int size_mask, const int *mask, const char *annotation, const char *API_label)
{
    Original_grid_info *level_3D_grid, *mid_3D_grid, *mid_1D_grid;
    Remap_grid_class *level_1D_CoR_grid, *level_H2D_CoR_grid, *mid_1D_CoR_grid, *mid_3D_CoR_grid, *sub_grids[256], *sized_sub_grids[256], *existing_level_1D_CoR_grid;
    char grid_name[NAME_STR_SIZE];
    int num_sized_sub_grids;


    level_3D_grid = get_original_grid(level_3D_grid_id);
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(level_3D_grid->get_comp_id(), "in Original_grid_mgt::register_mid_point_grid");
    synchronize_comp_processes_for_API(level_3D_grid->get_comp_id(), API_ID_GRID_MGT_REG_MID_POINT_GRID, comm, "registering a mid-point grid", annotation);
    check_API_parameter_string(level_3D_grid->get_comp_id(), API_ID_GRID_MGT_REG_MID_POINT_GRID, comm, "registering a mid-point grid", level_3D_grid->get_grid_name(), "\"level_3D_grid_id\" (the name of the interface-level grid)", annotation);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, level_3D_grid->get_comp_id(), level_3D_grid->is_V3D_grid(), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the grid specified by the parameter \"level_3D_grid_id\" is not a 3-D grid. Please verify the model code with the annotation \"%s\".", API_label, annotation);
    if (level_3D_grid->get_mid_point_grid() != NULL)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, level_3D_grid->get_comp_id(), false, "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the mid-point grid of the grid \"%s\" has been registered before (at the model code with the annotation \"%s\"). Please verify the model code with the annotation \"%s\".", API_label, level_3D_grid->get_grid_name(), level_3D_grid->get_mid_point_grid()->get_annotation(), annotation);
    if (level_3D_grid->get_interface_level_grid() != NULL)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, level_3D_grid->get_comp_id(), false, "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the specified interface-level grid itself is a mid-point grid (registered at the model code with the annotation \"%s\") of the grid \"%s\". It cannot be used to register another mid-point grid. Please verify the model code with the annotation \"%s\".", API_label, level_3D_grid->get_annotation(), level_3D_grid->get_interface_level_grid()->get_grid_name(), annotation);
    check_API_parameter_data_array(level_3D_grid->get_comp_id(), API_ID_GRID_MGT_REG_MID_POINT_GRID, comm, "registering a mid-point grid", size_mask, sizeof(int), (const char*)mask, "\"mask\"", annotation);
    level_1D_CoR_grid = level_3D_grid->get_V1D_sub_CoR_grid();
    level_H2D_CoR_grid = level_3D_grid->get_H2D_sub_CoR_grid();
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, level_1D_CoR_grid != NULL && level_H2D_CoR_grid != NULL, "Software error in Original_grid_mgt::register_mid_point_grid: NULL level_1D_CoR_grid or NULL level_H2D_CoR_grid");
    existing_level_1D_CoR_grid = level_1D_CoR_grid->get_mid_point_grid();
    mid_1D_CoR_grid = level_1D_CoR_grid->generate_mid_point_grid();
    level_3D_grid->get_original_CoR_grid()->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
    if (sized_sub_grids[0]->has_grid_coord_label(COORD_LABEL_LON) || sized_sub_grids[0]->has_grid_coord_label(COORD_LABEL_LAT)) {
        sub_grids[0] = level_H2D_CoR_grid;
        sub_grids[1] = mid_1D_CoR_grid;
    }
    else {
        sub_grids[0] = mid_1D_CoR_grid;
        sub_grids[1] = level_H2D_CoR_grid;
    }
    sprintf(grid_name, "mid_grid_for_%s", level_3D_grid->get_original_CoR_grid()->get_grid_name());
    mid_3D_CoR_grid = level_3D_grid->get_original_CoR_grid()->get_mid_point_grid();
    if (mid_3D_CoR_grid == NULL) {
        mid_3D_CoR_grid = new Remap_grid_class(grid_name, 2, sub_grids, 0);
        remap_grid_manager->add_remap_grid(mid_3D_CoR_grid);
    }
    mid_1D_CoR_grid->end_grid_definition_stage(NULL);
    mid_3D_CoR_grid->end_grid_definition_stage(NULL);
    mid_1D_grid = search_grid_info(mid_1D_CoR_grid->get_grid_name(), level_3D_grid->get_comp_id());
    if (mid_1D_grid == NULL) {
        mid_1D_grid = new Original_grid_info(level_3D_grid->get_comp_id(), original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, mid_1D_CoR_grid->get_grid_name(), annotation, mid_1D_CoR_grid, true, true, false);
        original_grids.push_back(mid_1D_grid);
        if (existing_level_1D_CoR_grid == NULL)
            remap_grid_manager->add_remap_grid(mid_1D_grid->get_original_CoR_grid());
    }
    *mid_1D_grid_id = mid_1D_grid->get_grid_id();
    mid_3D_grid = search_grid_info(mid_3D_CoR_grid->get_grid_name(), level_3D_grid->get_comp_id());
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, mid_3D_grid == NULL, "Software error in Original_grid_mgt::register_mid_point_grid: mid_3D_grid exists");
    mid_3D_grid = new Original_grid_info(level_3D_grid->get_comp_id(), original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX, mid_3D_CoR_grid->get_grid_name(), annotation, mid_3D_CoR_grid, true, true, false);
    original_grids.push_back(mid_3D_grid);
    *mid_3D_grid_id = mid_3D_grid->get_grid_id();
    level_3D_grid->set_mid_point_grid(mid_3D_grid);
    if (size_mask > 0) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, level_3D_grid->get_comp_id(), size_mask == mid_3D_CoR_grid->get_grid_size(), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: the array size of \"mask\" is different from the size of the mid-point grid. Please verify the model code with the annotation \"%s\".", API_label, annotation);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, level_3D_grid->get_comp_id(), are_array_values_between_boundaries("integer", mask, size_mask, 0, 1, 0, false), "Error happens when calling the API \"%s\" to register the mid-point grid of a grid: some values of the parameter \"mask\" are wrong (not 0 and 1). Please check the model code related to the annotation \"%s\".", API_label, annotation);
        mid_3D_CoR_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_INT, (const char*)mask, 0);
    }
    if (level_1D_CoR_grid->get_super_grid_of_setting_coord_values() == level_1D_CoR_grid)
        mid_1D_CoR_grid->set_super_grid_of_setting_coord_values(mid_1D_CoR_grid);
    else if (level_1D_CoR_grid->get_super_grid_of_setting_coord_values() != NULL && level_1D_CoR_grid->get_super_grid_of_setting_coord_values()->get_num_dimensions() == 3)
        mid_1D_CoR_grid->set_super_grid_of_setting_coord_values(mid_3D_CoR_grid);
    level_3D_grid->get_original_CoR_grid()->set_mid_point_grid(mid_3D_CoR_grid);
    for (int i = 0; i < level_3D_grid->sub_grids_id.size(); i ++) {
        if (search_grid_info(level_3D_grid->sub_grids_id[i])->get_original_CoR_grid()->has_grid_coord_label(COORD_LABEL_LEV)) {

            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, search_grid_info(level_3D_grid->sub_grids_id[i])->get_original_CoR_grid()->get_num_dimensions() == 1, "Software error in Original_grid_mgt::register_mid_point_grid");
            mid_3D_grid->add_sub_grid_id(*mid_1D_grid_id);
        }
        else mid_3D_grid->add_sub_grid_id(level_3D_grid->sub_grids_id[i]);
    }
    mid_3D_grid->set_grid_checksum(level_3D_grid->get_checksum_H2D_mask(), level_3D_grid->get_checksum_H2D_center_lon(), level_3D_grid->get_checksum_H2D_center_lat());
}


Original_grid_info *Original_grid_mgt::promote_ensemble_member_grid_to_set(int set_comp_id, Original_grid_info *member_grid, bool at_top_level)
{
    Original_grid_info *set_original_grid = original_grid_mgr->search_grid_info(member_grid->get_grid_name(), set_comp_id);
    Comp_comm_group_mgt_node *member_comp_node = comp_comm_group_mgt_mgr->search_global_node(member_grid->get_comp_id());
    Comp_comm_group_mgt_node *set_comp_node = comp_comm_group_mgt_mgr->search_global_node(set_comp_id);
    int i;
    int *member_root_proc_ids = NULL;
    int member_root_proc_set_id = -1;
    int num_ensemble_members;
    char *ensemble_members_full_names;
    MPI_Status status;


    if (at_top_level)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->search_global_node(set_comp_id)->is_given_comp_a_child(comp_comm_group_mgt_mgr->search_global_node(member_grid->get_comp_id())), "Original_grid_mgt::promote_ensemble_member_grid_to_set");
    if (set_original_grid != NULL) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, set_original_grid->get_comp_id() == set_comp_id, "Original_grid_mgt::promote_ensemble_member_grid_to_set");
        return set_original_grid;
    }
    set_original_grid = new Original_grid_info(set_comp_id, member_grid, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX);

    if (at_top_level) {
        if (member_comp_node->get_current_proc_local_id() == 0)
            member_root_proc_set_id = set_comp_node->get_current_proc_local_id();
        if (set_comp_node->get_current_proc_local_id() == 0)
            member_root_proc_ids = new int [set_comp_node->get_num_procs()];
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Gather(&member_root_proc_set_id, 1, MPI_INT, member_root_proc_ids, 1, MPI_INT, 0, set_comp_node->get_comm_group()) == MPI_SUCCESS);
        if (set_comp_node->get_current_proc_local_id() == 0) {
            num_ensemble_members = 0;
            for (int i = 0; i < set_comp_node->get_num_procs(); i ++)
                if (member_root_proc_ids[i] > -1)
                    num_ensemble_members ++;
        }
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(&num_ensemble_members, 1, MPI_INT, 0, set_comp_node->get_comm_group()) == MPI_SUCCESS);
        ensemble_members_full_names = new char [NAME_STR_SIZE * num_ensemble_members];
        if (set_comp_node->get_current_proc_local_id() == 0) {
            i = 0;
            num_ensemble_members = 0;
            if (member_comp_node->get_current_proc_local_id() == 0) {
                strcpy(ensemble_members_full_names, member_comp_node->get_full_name());
                num_ensemble_members ++;
                i ++;
            }
            for (; i < set_comp_node->get_num_procs(); i ++)
                if (member_root_proc_ids[i] > -1) {
                    EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Recv(ensemble_members_full_names + NAME_STR_SIZE * num_ensemble_members, NAME_STR_SIZE, MPI_CHAR, i, i, set_comp_node->get_comm_group(), &status) == MPI_SUCCESS);
                    num_ensemble_members ++;
                }
        }
        if (member_comp_node->get_current_proc_local_id() == 0 && set_comp_node->get_current_proc_local_id() != 0)
            EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Send((char*)member_comp_node->get_full_name(), NAME_STR_SIZE, MPI_CHAR, 0, set_comp_node->get_current_proc_local_id(), set_comp_node->get_comm_group()) == MPI_SUCCESS);
        EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(ensemble_members_full_names, NAME_STR_SIZE * num_ensemble_members, MPI_CHAR, 0, set_comp_node->get_comm_group()) == MPI_SUCCESS);
        for (i = 0; i < num_ensemble_members; i ++) {
            Comp_comm_group_mgt_node *current_member_comp_node = comp_comm_group_mgt_mgr->search_global_node(ensemble_members_full_names + i * NAME_STR_SIZE);
            if (report_error_enabled) {
                int has_current_member_comp_node = current_member_comp_node == NULL ? 0 : 1;
                EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(&has_current_member_comp_node, 1, MPI_INT, 0, member_comp_node->get_comm_group()) == MPI_SUCCESS);
                EXECUTION_REPORT(REPORT_ERROR, -1, has_current_member_comp_node == 0 && current_member_comp_node == NULL || has_current_member_comp_node != 0 && current_member_comp_node != NULL, "Software error in Original_grid_mgt::promote_ensemble_member_grid_to_set");
            }
            if (current_member_comp_node == NULL)
                current_member_comp_node = comp_comm_group_mgt_mgr->load_comp_info_from_XML(member_comp_node->get_comp_id(), ensemble_members_full_names + i * NAME_STR_SIZE, member_comp_node->get_comm_group());
            if (original_grid_mgr->search_grid_info(member_grid->get_grid_name(), current_member_comp_node->get_comp_id()) == NULL) {
                EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_member_comp_node != member_comp_node, "Software error in Original_grid_mgt::promote_ensemble_member_grid_to_set");
                new Original_grid_info(current_member_comp_node->get_comp_id(), member_grid, original_grids.size() | TYPE_GRID_LOCAL_ID_PREFIX);
            }
        }

        if (set_original_grid->get_H2D_sub_grid() != NULL) {
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, member_grid->get_H2D_sub_grid() != NULL, "Software error in Original_grid_mgt::promote_ensemble_member_grid_to_set");
            Comp_comm_group_mgt_node *member_comp_node = comp_comm_group_mgt_mgr->search_global_node(member_grid->get_comp_id());
            Comp_comm_group_mgt_node *set_comp_node = comp_comm_group_mgt_mgr->search_global_node(set_comp_id);
            char first_member_full_name[NAME_STR_SIZE];
            strcpy(first_member_full_name, member_comp_node->get_full_name());
            EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Bcast(first_member_full_name, NAME_STR_SIZE, MPI_CHAR, 0, set_comp_node->get_comm_group()) == MPI_SUCCESS);
            bool is_first_member = words_are_the_same(first_member_full_name, member_comp_node->get_full_name());
            Decomp_grid_info *member_decomp_grid = NULL;
            Distributed_H2D_grid_engine *dst_distributed_H2D_grid = distributed_H2D_grid_mgr->generate_distributed_H2D_grid_engine(set_comp_node->get_full_name(), set_original_grid->get_H2D_sub_grid()->get_grid_name(), set_original_grid->get_H2D_sub_grid()->get_original_CoR_grid(), set_original_grid->get_H2D_sub_grid()->get_original_CoR_grid()->get_grid_size(), 0, 0, set_original_grid->get_H2D_sub_grid(), false);
            if (is_first_member)
                member_decomp_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(member_comp_node->get_full_name(), member_grid->get_H2D_sub_CoR_grid())->get_basic_decomp_grid();
            dst_distributed_H2D_grid->generate_basic_decomp_grid_from_decomp_grid(first_member_full_name, member_decomp_grid, set_comp_node->get_full_name());
            EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Barrier(set_comp_node->get_comm_group()) == MPI_SUCCESS);
        }

        if (member_root_proc_ids != NULL)
            delete [] member_root_proc_ids;
        delete [] ensemble_members_full_names;
    }
    return set_original_grid;
}


int Original_grid_mgt::get_V3D_lev_field_id_for_grid(Original_grid_info *original_grid)
{
    for (int i = 0; i < original_grids.size(); i ++) {
        if (original_grids[i] == original_grid)
            continue;
        if (original_grids[i]->V3D_lev_field_id != -1 && original_grids[i]->get_original_CoR_grid()->is_similar_grid_with(original_grid->get_original_CoR_grid())) {
            original_grid->V3D_lev_field_id = original_grids[i]->V3D_lev_field_id;
            original_grid->V3D_lev_field_variation_type = original_grids[i]->V3D_lev_field_variation_type;
            return original_grids[i]->V3D_lev_field_id;
        }
    }

    return -1;
}


int Original_grid_mgt::get_bottom_field_id_for_grid(Original_grid_info *original_grid)
{
    for (int i = 0; i < original_grids.size(); i ++) {
        if (original_grids[i] == original_grid)
            continue;
        if (original_grids[i]->bottom_field_id != -1 && original_grids[i]->get_original_CoR_grid()->is_similar_grid_with(original_grid->get_original_CoR_grid())) {
            original_grid->bottom_field_id = original_grids[i]->bottom_field_id;
            original_grid->bottom_field_variation_type = original_grids[i]->bottom_field_variation_type;
            original_grid->bottom_field_name = strdup(original_grids[i]->bottom_field_name);
            return original_grids[i]->bottom_field_id;
        }
    }

    return -1;
}

