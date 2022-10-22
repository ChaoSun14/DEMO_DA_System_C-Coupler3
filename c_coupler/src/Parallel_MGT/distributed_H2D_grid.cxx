/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "distributed_H2D_grid.h"
#include "triangulation.h"
#include "remap_utils_nearest_points.h"
#include <math.h>


void Distributed_H2D_grid_engine::initialize_data(const char *comp_full_name, const char *grid_name, long global_grid_size, int num_lons, int num_lats, Remap_grid_class *CoR_grid_initialization, Original_grid_info *original_grid)
{
	this->comp_full_name = strdup(comp_full_name);
	this->grid_name = strdup(grid_name);
	this->num_lons = num_lons;
	this->num_lats = num_lats;
	this->basic_decomp_grid = NULL;
	this->global_grid_size = global_grid_size;
	this->CoR_grid_initialization = CoR_grid_initialization;
	this->original_grid = original_grid;
	if (num_lons > 0)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, global_grid_size == num_lats * num_lons, "Software error in Distributed_H2D_grid_engine::initialize_data");	
}


/* Default distributed grid from global CoR grid */
Distributed_H2D_grid_engine::Distributed_H2D_grid_engine(const char *comp_full_name, const char *grid_name, Remap_grid_class *CoR_grid_initialization, long global_grid_size, int num_lons, int num_lats, Original_grid_info *original_grid, bool generate_basic_decomp_grid)
{
	initialize_data(comp_full_name, grid_name, global_grid_size, num_lons, num_lats, CoR_grid_initialization, original_grid);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, CoR_grid_initialization->get_is_H2D_grid(), "Software error in Distributed_H2D_grid_engine::Distributed_H2D_grid_engine");
	if (generate_basic_decomp_grid) {
		Decomp_info *default_decomp = decomps_info_mgr->generate_default_parallel_decomp(original_grid, comp_comm_group_mgt_mgr->search_global_node(comp_full_name)->get_comp_id());
		basic_decomp_grid = decomp_grids_mgr->search_decomp_grid_info(default_decomp->get_decomp_id(), CoR_grid_initialization, true);
	}
}


Distributed_H2D_grid_engine::~Distributed_H2D_grid_engine()
{
	delete [] comp_full_name;
	delete [] grid_name;
}


Decomp_grid_info *Distributed_H2D_grid_engine::get_basic_decomp_grid()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, basic_decomp_grid != NULL, "Software error in Distributed_H2D_grid_engine::get_basic_decomp_grid: %s %lx", grid_name, CoR_grid_initialization);
	return basic_decomp_grid;
}


void Distributed_H2D_grid_engine::allocate_grid_field_inst(Remap_grid_data_class *grid_data_field, int grid_id, int decomp_id, bool change_grid_data_buf, Field_mem_info **grid_field_insts, int &num_grid_field_insts, int connection_id)
{
	if (grid_data_field == NULL)
		return;
	
	Field_mem_info *field_inst = memory_manager->alloc_mem(grid_data_field->get_grid_data_field()->field_name_in_application, decomp_id, grid_id, BUF_MARK_GRID_FIELD-connection_id, grid_data_field->get_grid_data_field()->data_type_in_application, "unit_neglected", "C-Coupler internal", false, true);

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "allocate field inst of grid %s: %s  %s", grid_name, grid_data_field->get_grid_data_field()->field_name_in_application, grid_data_field->get_grid_data_field()->data_type_in_application);
	if (change_grid_data_buf) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field->get_grid_data_field()->required_data_size == field_inst->get_size_of_field(), "Software error in Distributed_H2D_grid_engine::allocate_grid_field_inst: %d  %d", grid_data_field->get_grid_data_field()->required_data_size, field_inst->get_size_of_field());
		field_inst->reset_field_data(grid_data_field);
//		field_inst->reset_mem_buf(grid_data_field->get_grid_data_field()->data_buf, true, REG_FIELD_TAG_NONE);
	}
	grid_field_insts[num_grid_field_insts++] = field_inst;
	grid_fields_insts.push_back(field_inst);
}


Decomp_grid_info *Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer(const char *src_host_comp_full_name, const char *src_comp_full_name, Decomp_grid_info *src_decomp_grid, const char *dst_host_comp_full_name, const char *dst_comp_full_name, int dst_decomp_id, Remap_grid_class *dst_CoR_grid)
{
	Comp_comm_group_mgt_node *src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_host_comp_full_name);
	Comp_comm_group_mgt_node *dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_host_comp_full_name);
	Remap_grid_class *src_decomp_CoR_grid = NULL, *src_original_CoR_grid = NULL;
	int src_original_grid_id = -1, src_comp_id = -1;
	int src_num_grid_field_insts=0, dst_num_grid_field_insts=0;
	int grid_id_for_num_vertex_src=-1, grid_id_for_vertexes_src=-1, grid_id_for_num_vertex_dst=-1, grid_id_for_vertexes_dst=-1;
	char grid_name_for_num_vertex[NAME_STR_SIZE], grid_name_for_vertexes[NAME_STR_SIZE];
	Remap_grid_data_class *grid_data_field;
	Runtime_trans_algorithm *runtime_send_algorithm, *runtime_recv_algorithm;
	Field_mem_info *src_grid_field_insts[1024], *dst_grid_field_insts[1024];
	Remap_grid_data_class *dst_center_lons=NULL, *dst_center_lats=NULL, *dst_vertex_lons=NULL, *dst_vertex_lats=NULL, *dst_grid_mask=NULL, *dst_area = NULL;
	Decomp_grid_info *decomp_grid = NULL;
	int connection_id = coupling_generator->apply_connection_id();
	double time1, time2, time3;


	wtime(&time1);
	if (src_comp_node != NULL && src_comp_node->get_current_proc_local_id() != -1) {
		if (src_decomp_grid == NULL)
			src_decomp_grid = basic_decomp_grid; 
		src_decomp_CoR_grid = src_decomp_grid->get_decomp_grid();
		src_original_CoR_grid = src_decomp_grid->get_original_grid();
		src_original_grid_id = decomps_info_mgr->get_decomp_info(src_decomp_grid->get_decomp_id())->get_grid_id();
		src_comp_id = original_grid_mgr->get_original_grid(src_original_grid_id)->get_comp_id();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_decomp_CoR_grid != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
	}

	if (src_comp_node != NULL && src_comp_node->get_current_proc_local_id() != -1) {
		EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Start to generate decomp_grid via data transfer from component \"%s\" (\"%s\") (current component) to \"%s\" (\"%s\") : %s.", src_comp_full_name, src_host_comp_full_name, dst_comp_full_name, dst_host_comp_full_name, src_original_CoR_grid->get_grid_name());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_CoR_grid != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
	}
	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1)
		EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Start to generate decomp_grid via data transfer from component \"%s\" (\"%s\") to \"%s\" (\"%s\") (current component): %s", src_comp_full_name, src_host_comp_full_name, dst_comp_full_name, dst_host_comp_full_name, dst_CoR_grid->get_grid_name());
	if (src_comp_node == dst_comp_node) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_comp_node->get_current_proc_local_id() != -1, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_CoR_grid != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_CoR_grid != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_CoR_grid == dst_CoR_grid, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer: %s (%lx)   %s (%lx)", src_original_CoR_grid->get_grid_name(), src_original_CoR_grid, dst_CoR_grid->get_grid_name(), dst_CoR_grid);
	}
	
	sprintf(grid_name_for_num_vertex, "temp_grid_for_num_vertex_of_%s", grid_name);
	sprintf(grid_name_for_vertexes, "temp_grid_for_vertexes_of_%s", grid_name);
	if (src_comp_node != NULL && src_comp_node->get_current_proc_local_id() != -1 && src_original_CoR_grid->get_num_vertexes() > 0) {
		grid_id_for_num_vertex_src = original_grid_mgr->register_V1D_grid_via_data(API_ID_GRID_MGT_REG_NORMAL_1D_GRID_NO_DATA, src_comp_id, grid_name_for_num_vertex, 4, "none", src_original_CoR_grid->get_num_vertexes(), 0.0, NULL, NULL, false, "C-Coupler internal");
		grid_id_for_vertexes_src = original_grid_mgr->register_md_grid_via_multi_grids(src_comp_id, grid_name_for_vertexes, grid_id_for_num_vertex_src, src_original_grid_id, -1, -1, NULL, false, "C-Coupler internal");
		grid_data_field = src_decomp_CoR_grid->get_grid_vertex_field(COORD_LABEL_LON);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		if (src_decomp_CoR_grid->get_grid_size() > 0)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field->get_grid_data_field()->data_buf!= NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer: %d");
		allocate_grid_field_inst(grid_data_field, grid_id_for_vertexes_src, src_decomp_grid->get_decomp_id(), true, src_grid_field_insts, src_num_grid_field_insts, connection_id);
		grid_data_field = src_decomp_CoR_grid->get_grid_vertex_field(COORD_LABEL_LAT);
		allocate_grid_field_inst(grid_data_field, grid_id_for_vertexes_src, src_decomp_grid->get_decomp_id(), true, src_grid_field_insts, src_num_grid_field_insts, connection_id);
	}
	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1 && dst_CoR_grid->get_num_vertexes() > 0) {
		grid_id_for_num_vertex_dst = original_grid_mgr->register_V1D_grid_via_data(API_ID_GRID_MGT_REG_NORMAL_1D_GRID_NO_DATA, original_grid->get_comp_id(), grid_name_for_num_vertex, 4, "none", dst_CoR_grid->get_num_vertexes(), 0.0, NULL, NULL, false, "C-Coupler internal");
		grid_id_for_vertexes_dst = original_grid_mgr->register_md_grid_via_multi_grids(original_grid->get_comp_id(), grid_name_for_vertexes, grid_id_for_num_vertex_dst, original_grid->get_grid_id(), -1, -1, NULL, false, "C-Coupler internal");
		grid_data_field = dst_CoR_grid->get_grid_vertex_field(COORD_LABEL_LON);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		allocate_grid_field_inst(grid_data_field, grid_id_for_vertexes_dst, dst_decomp_id, false, dst_grid_field_insts, dst_num_grid_field_insts, connection_id);
		dst_vertex_lons = dst_grid_field_insts[dst_num_grid_field_insts-1]->get_field_data();
		grid_data_field = dst_CoR_grid->get_grid_vertex_field(COORD_LABEL_LAT);
		allocate_grid_field_inst(grid_data_field, grid_id_for_vertexes_dst, dst_decomp_id, false, dst_grid_field_insts, dst_num_grid_field_insts, connection_id);
		dst_vertex_lats = dst_grid_field_insts[dst_num_grid_field_insts-1]->get_field_data();
	}
	if (src_comp_node != NULL && src_comp_node->get_current_proc_local_id() != -1) {
		allocate_grid_field_inst(src_decomp_CoR_grid->get_grid_center_field(COORD_LABEL_LON), src_original_grid_id, src_decomp_grid->get_decomp_id(), true, src_grid_field_insts, src_num_grid_field_insts, connection_id);
		allocate_grid_field_inst(src_decomp_CoR_grid->get_grid_center_field(COORD_LABEL_LAT), src_original_grid_id, src_decomp_grid->get_decomp_id(), true, src_grid_field_insts, src_num_grid_field_insts, connection_id);
		allocate_grid_field_inst(src_decomp_CoR_grid->get_grid_mask_field(), src_original_grid_id, src_decomp_grid->get_decomp_id(), true, src_grid_field_insts, src_num_grid_field_insts, connection_id);
		allocate_grid_field_inst(src_decomp_CoR_grid->get_grid_imported_area(), src_original_grid_id, src_decomp_grid->get_decomp_id(), true, src_grid_field_insts, src_num_grid_field_insts, connection_id);
	}
	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1) {
		grid_data_field = dst_CoR_grid->get_grid_center_field(COORD_LABEL_LON);
		allocate_grid_field_inst(grid_data_field, original_grid->get_grid_id(), dst_decomp_id, false, dst_grid_field_insts, dst_num_grid_field_insts, connection_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		dst_center_lons = dst_grid_field_insts[dst_num_grid_field_insts-1]->get_field_data();
		grid_data_field = dst_CoR_grid->get_grid_center_field(COORD_LABEL_LAT);
		allocate_grid_field_inst(grid_data_field, original_grid->get_grid_id(), dst_decomp_id, false, dst_grid_field_insts, dst_num_grid_field_insts, connection_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		dst_center_lats = dst_grid_field_insts[dst_num_grid_field_insts-1]->get_field_data();
		grid_data_field = dst_CoR_grid->get_grid_mask_field();
		allocate_grid_field_inst(grid_data_field, original_grid->get_grid_id(), dst_decomp_id, false, dst_grid_field_insts, dst_num_grid_field_insts, connection_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field != NULL, "Software error in Distributed_H2D_grid_engine::generate_decomp_grid_via_data_transfer");
		dst_grid_mask = dst_grid_field_insts[dst_num_grid_field_insts-1]->get_field_data();
		allocate_grid_field_inst(dst_CoR_grid->get_grid_imported_area(), original_grid->get_grid_id(), dst_decomp_id, false, dst_grid_field_insts, dst_num_grid_field_insts, connection_id);
		if (dst_CoR_grid->get_grid_imported_area() != NULL)
			dst_area = dst_grid_field_insts[dst_num_grid_field_insts-1]->get_field_data();
	}
	wtime(&time3);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "start to generate_simple_data_transfer_connection");

	coupling_generator->generate_simple_data_transfer_connection(src_host_comp_full_name, src_comp_full_name, src_num_grid_field_insts, src_grid_field_insts, dst_host_comp_full_name, dst_comp_full_name, dst_num_grid_field_insts, dst_grid_field_insts, &runtime_send_algorithm, &runtime_recv_algorithm);
	wtime(&time2);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "prepare coupling connection for grid distributed %lf vs %lf", time2 - time1, time3 - time1);

	if (runtime_send_algorithm != NULL) {
		runtime_send_algorithm->pass_transfer_parameters(-1, -1);
		runtime_send_algorithm->run(true);
	}
	if (runtime_recv_algorithm != NULL) {
		runtime_recv_algorithm->pass_transfer_parameters(-1, -1);
		runtime_recv_algorithm->run(true);
	}
	if (runtime_send_algorithm != NULL)
		runtime_send_algorithm->wait_sending_data();
	wtime(&time3);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "use %lf seconds for grid distributed transfer", time3 - time2);

	if (runtime_send_algorithm != NULL)
		delete runtime_send_algorithm;
	if (runtime_recv_algorithm != NULL)
		delete runtime_recv_algorithm;

	for (int i = 0; i < dst_num_grid_field_insts; i ++)
		dst_grid_field_insts[i]->transformation_between_chunks_array(true);

	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1)
		decomp_grid = new Decomp_grid_info(dst_decomp_id, dst_CoR_grid, dst_center_lons, dst_center_lats, dst_vertex_lons, dst_vertex_lats, dst_grid_mask, dst_area);

	if (src_comp_node != NULL && src_comp_node->get_current_proc_local_id() != -1)
		EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Finish generating decomp_grid via data transfer from component \"%s\" (current component) to \"%s\".", src_comp_full_name, dst_comp_full_name);
	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1)
		EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Finish generating decomp_grid via data transfer from component \"%s\" (current component) to \"%s\".", src_comp_full_name, dst_comp_full_name);

	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1) {
		Decomp_info *dst_decomp_info = decomps_info_mgr->get_decomp_info(dst_decomp_id);
		for (int i = 0; i < dst_decomp_info->get_num_local_cells(); i ++) {
			if ((dst_decomp_info->get_local_cell_global_indx())[i] == CCPL_NULL_INT) {
				((double*)dst_center_lons->get_grid_data_field()->data_buf)[i] = NULL_COORD_VALUE;
				((double*)dst_center_lats->get_grid_data_field()->data_buf)[i] = NULL_COORD_VALUE;
				((bool*)dst_grid_mask->get_grid_data_field()->data_buf)[i] = false;
				if (dst_vertex_lats != NULL)
					for (int j = 0; j < decomp_grid->get_decomp_grid()->get_num_vertexes(); j ++) {
						((double*)dst_vertex_lons->get_grid_data_field()->data_buf)[i*decomp_grid->get_decomp_grid()->get_num_vertexes()+j] = NULL_COORD_VALUE;
						((double*)dst_vertex_lats->get_grid_data_field()->data_buf)[i*decomp_grid->get_decomp_grid()->get_num_vertexes()+j] = NULL_COORD_VALUE;
					}
			}
		}
	}

	if (src_comp_node != NULL && src_comp_node->get_current_proc_local_id() != -1)
		EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_node->get_comp_id(), true, "Finish generating decomp_grid via data transfer from component \"%s\" (current component) to \"%s\": %s %lx.", src_comp_full_name, dst_comp_full_name, src_original_CoR_grid->get_grid_name(), src_original_CoR_grid);
	if (dst_comp_node != NULL && dst_comp_node->get_current_proc_local_id() != -1)
		EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_node->get_comp_id(), true, "Finish generating decomp_grid via data transfer from component \"%s\" to \"%s\" (current component): %s %lx", src_comp_full_name, dst_comp_full_name, dst_CoR_grid->get_grid_name(), dst_CoR_grid);

	for (int i = 0; i < grid_fields_insts.size(); i ++) {
		grid_fields_insts[i]->empty_grided_field_data();
		memory_manager->delete_field_inst(grid_fields_insts[i]);
	}
	grid_fields_insts.clear();

	return decomp_grid;
}


Decomp_grid_info *Distributed_H2D_grid_engine::generate_basic_decomp_grid_from_decomp_grid(const char *src_comp_full_name, Decomp_grid_info *src_decomp_grid, const char *dst_comp_full_name)
{
	int dst_decomp_id = -1;
	Remap_grid_class *dst_CoR_grid = NULL;
	char *dst_decomp_comp_full_name = NULL;

	
	if (words_are_the_same(comp_full_name, dst_comp_full_name)) {
		dst_decomp_id = decomps_info_mgr->generate_default_parallel_decomp(original_grid, comp_comm_group_mgt_mgr->search_global_node(dst_comp_full_name)->get_comp_id())->get_decomp_id();
		dst_CoR_grid = original_grid->get_original_CoR_grid();
		dst_decomp_comp_full_name = (char*) (comp_comm_group_mgt_mgr->search_global_node(decomps_info_mgr->get_decomp_info(dst_decomp_id)->get_comp_id())->get_full_name());
	}

	Decomp_grid_info *new_decomp_grid = generate_decomp_grid_via_data_transfer(src_comp_full_name, src_comp_full_name, src_decomp_grid, dst_comp_full_name, dst_decomp_comp_full_name, dst_decomp_id, dst_CoR_grid);

	if (words_are_the_same(comp_full_name, dst_comp_full_name) && basic_decomp_grid != NULL) {
		Remap_grid_data_class *old_center_lon = basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON);
		Remap_grid_data_class *old_center_lat = basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LAT);
		Remap_grid_data_class *new_center_lon = new_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON);
		Remap_grid_data_class *new_center_lat = new_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LAT);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, are_two_coord_arrays_same((double*)old_center_lon->get_grid_data_field()->data_buf, (double*)new_center_lon->get_grid_data_field()->data_buf, basic_decomp_grid->get_decomp_grid()->get_grid_size(), new_decomp_grid->get_decomp_grid()->get_grid_size(), true), "Software error in Distributed_H2D_grid_engine::generate_basic_decomp_grid_from_decomp_grid");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, are_two_coord_arrays_same((double*)old_center_lat->get_grid_data_field()->data_buf, (double*)new_center_lat->get_grid_data_field()->data_buf, basic_decomp_grid->get_decomp_grid()->get_grid_size(), new_decomp_grid->get_decomp_grid()->get_grid_size(), true), "Software error in Distributed_H2D_grid_engine::generate_basic_decomp_grid_from_decomp_grid");
	}
	else if (words_are_the_same(comp_full_name, dst_comp_full_name)) {
		basic_decomp_grid = new_decomp_grid;
		decomp_grids_mgr->add_decomp_grid_info(basic_decomp_grid);
	}

	return new_decomp_grid;
}


void Distributed_H2D_grid_engine::update_basic_decomp_grid(Decomp_grid_info * src_decomp_grid)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, basic_decomp_grid != NULL, "Software error in Distributed_H2D_grid_engine::update_basic_decomp_grid");
	decomp_grids_mgr->delete_decomp_grid(basic_decomp_grid);
	basic_decomp_grid = NULL;
	generate_basic_decomp_grid_from_decomp_grid(comp_full_name, src_decomp_grid, comp_full_name);
	if (report_error_enabled) {
		basic_decomp_grid->get_decomp_grid()->check_center_vertex_values_consistency_2D(true);
		double *local_center_lons = (double*) (basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
		double *local_center_lats = (double*) (basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);
		const int *local_cell_global_indexes = basic_decomp_grid->get_decomp_grid()->get_local_cell_global_indexes();
		double *global_center_lons = (double*) (CoR_grid_initialization->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
		double *global_center_lats = (double*) (CoR_grid_initialization->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);	
		if (CoR_grid_initialization == CoR_grid_initialization->get_grid_center_field(COORD_LABEL_LON)->get_coord_value_grid()) {
			if (global_center_lons != NULL)
				for (int i = 0; i < basic_decomp_grid->get_decomp_grid()->get_grid_size(); i ++) 
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, local_center_lons[i] == global_center_lons[local_cell_global_indexes[i]] && local_center_lats[i] == global_center_lats[local_cell_global_indexes[i]], "Software error in Distributed_H2D_grid_engine::update_basic_decomp_grid");
		}
		else {
			int dim1 = CoR_grid_initialization->get_a_leaf_grid(COORD_LABEL_LON)->get_grid_size();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dim1 != 0, "Software error in Distributed_H2D_grid_engine::update_basic_decomp_grid");
			for (int i = 0; i < basic_decomp_grid->get_decomp_grid()->get_grid_size(); i ++) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, local_center_lons[i] == global_center_lons[local_cell_global_indexes[i]%dim1] && local_center_lats[i] == global_center_lats[local_cell_global_indexes[i]/dim1], "Software error in Distributed_H2D_grid_engine::update_basic_decomp_grid");
			}			
		}
	}
}


bool Distributed_H2D_grid_engine::is_the_same_as_another_distributed_H2D_grid(Distributed_H2D_grid_engine *another)
{
	int local_check_result, overall_check_result;
	double *this_center_lon, *another_center_lon, *this_center_lat, *another_center_lat; 
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(this->comp_full_name);


	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(this->comp_full_name, another->comp_full_name), "Software error in Distributed_H2D_grid_engine::is_the_same_as_another_distributed_H2D_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->basic_decomp_grid != NULL && this->basic_decomp_grid->get_decomp_grid()->get_grid_size() > 0 && another->basic_decomp_grid != NULL && another->basic_decomp_grid->get_decomp_grid()->get_grid_size() > 0, "Software error in Distributed_H2D_grid_engine::is_the_same_as_another_distributed_H2D_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->CoR_grid_initialization->get_grid_size() == this->global_grid_size && another->CoR_grid_initialization->get_grid_size() == another->global_grid_size, "Software error in Distributed_H2D_grid_engine::is_the_same_as_another_distributed_H2D_grid");

	if (this->global_grid_size != another->global_grid_size)
		return false;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this->basic_decomp_grid->get_decomp_grid()->get_grid_size() == another->basic_decomp_grid->get_decomp_grid()->get_grid_size(), "Software error in Distributed_H2D_grid_engine::is_the_same_as_another_distributed_H2D_grid");	

	this_center_lon = (double*)this->basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf;
	another_center_lon = (double*)another->basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf;
	local_check_result = are_two_coord_arrays_same(this_center_lon, another_center_lon, this->basic_decomp_grid->get_decomp_grid()->get_grid_size(), another->basic_decomp_grid->get_decomp_grid()->get_grid_size(), false)? 1 : 0;
	if (local_check_result == 1) {
		this_center_lat = (double*)this->basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf;
		another_center_lat = (double*)another->basic_decomp_grid->get_decomp_grid()->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf;
		local_check_result = are_two_coord_arrays_same(this_center_lat, another_center_lat, this->basic_decomp_grid->get_decomp_grid()->get_grid_size(), another->basic_decomp_grid->get_decomp_grid()->get_grid_size(), false)? 1 : 0;
	}

	EXECUTION_REPORT(REPORT_ERROR, -1, MPI_Allreduce(&local_check_result, &overall_check_result, 1, MPI_INT, MPI_SUM, comp_node->get_comm_group()) == MPI_SUCCESS);
	return overall_check_result == comp_node->get_num_procs();
}


Remapping_grid_domain_decomp_engine::Remapping_grid_domain_decomp_engine(int src_original_grid_id, int comp_id, Remap_grid_class *src_original_cor_grid, Remap_grid_class *dst_original_cor_grid, double whole_domain_min_lon, double whole_domain_max_lon, double whole_domain_min_lat, double whole_domain_max_lat)
{
	double whole_domain_lon_diff, whole_domain_lat_diff, lon_lat_diff_ratio;
	double pole_region_min_width_fator = 3, pole_region_max_width_fator = 10.0, nonpole_region_max_width_fator = 45.0, workload_redundant_copy_factor = 2.0, pole_region_width, pole_region_lat_value, whole_area_pre_process, pole_workload_shrink_fator = 1.0;
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	double time1, time2, time3;

	wtime(&time1);
	this->comp_id = comp_id;
	this->num_total_procs = comp_node->get_num_procs();
	this->src_original_grid_id = src_original_grid_id;
	this->src_original_cor_grid = src_original_cor_grid;
	this->dst_original_cor_grid = dst_original_cor_grid;
	this->src_grid_original_num_vertex = src_original_cor_grid->get_num_vertexes();
	if (dst_original_cor_grid != NULL)
		this->dst_grid_original_num_vertex = dst_original_cor_grid->get_num_vertexes();
	this->whole_domain_min_lon = whole_domain_min_lon;
	this->whole_domain_max_lon = whole_domain_max_lon;
	this->whole_domain_min_lat = whole_domain_min_lat;
	this->whole_domain_max_lat = whole_domain_max_lat;
	this->nonPole_region_min_lat = whole_domain_min_lat;
	this->nonPole_region_max_lat = whole_domain_max_lat;
	this->num_northPole_subdomain = 0;
	this->num_southPole_subdomain = 0;
	this->src_grid_one_sided_comm_data_buf = NULL;
	this->dst_grid_one_sided_comm_data_buf = NULL;

	whole_domain_lat_diff = whole_domain_max_lat - whole_domain_min_lat;
	if (whole_domain_min_lon < whole_domain_max_lon)
		whole_domain_lon_diff = whole_domain_max_lon - whole_domain_min_lon;
	else whole_domain_lon_diff = whole_domain_max_lon - whole_domain_min_lon + 360.0;	

	if ((whole_domain_max_lat == 90.0 || whole_domain_min_lat == -90.0) && words_are_the_same(src_original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES)) {
		whole_area_pre_process = fabs(sin(DEGREE_TO_RADIAN(whole_domain_max_lat))-sin(DEGREE_TO_RADIAN(whole_domain_min_lat))) / comp_node->get_num_procs() / workload_redundant_copy_factor;
		pole_region_lat_value = 90. - RADIAN_TO_DEGREE(asin(1 - whole_area_pre_process/pole_workload_shrink_fator));
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "pole_region_lat_value is %lf", pole_region_lat_value);
		pole_region_width = MAX(MIN(pole_region_max_width_fator, pole_region_lat_value), pole_region_min_width_fator);
		if (whole_domain_max_lat == 90.0) {
			num_northPole_subdomain = 1;
			nonPole_region_max_lat = whole_domain_max_lat - pole_region_width;
			if (src_original_cor_grid->get_boundary_max_lat() != 90.0)
				nonPole_region_max_lat = src_original_cor_grid->get_boundary_max_lat() + (90.0 - src_original_cor_grid->get_boundary_max_lat()) * 0.1;
			if (dst_original_cor_grid != NULL && dst_original_cor_grid->get_boundary_max_lat() != 90.0)
				nonPole_region_max_lat = dst_original_cor_grid->get_boundary_max_lat() + (90.0 - dst_original_cor_grid->get_boundary_max_lat()) * 0.1;
		}
		if (whole_domain_min_lat == -90.0) {
			num_southPole_subdomain = 1;
			nonPole_region_min_lat = whole_domain_min_lat + pole_region_width;
			if (src_original_cor_grid->get_boundary_min_lat() != -90.0)
				nonPole_region_min_lat = src_original_cor_grid->get_boundary_min_lat() + (-90.0 - src_original_cor_grid->get_boundary_min_lat()) * 0.1;
			if (dst_original_cor_grid != NULL && dst_original_cor_grid->get_boundary_min_lat() != -90.0)
				nonPole_region_min_lat = dst_original_cor_grid->get_boundary_min_lat() + (-90.0 - dst_original_cor_grid->get_boundary_min_lat()) * 0.1;
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, nonPole_region_min_lat < whole_domain_max_lat && nonPole_region_max_lat > whole_domain_min_lat, "Software error in remapping_grid_domain_decom_engine::remapping_grid_domain_decom_engine");
	}

	int num_remained_procs = MAX(1, comp_node->get_num_procs()-num_northPole_subdomain-num_southPole_subdomain);
	whole_domain_lat_diff = nonPole_region_max_lat - nonPole_region_min_lat;
	lon_lat_diff_ratio = whole_domain_lon_diff / whole_domain_lat_diff;
	double num_nonPole_region_y_subdomains_double = sqrt((double)(num_remained_procs) / lon_lat_diff_ratio);
	double num_nonPole_region_x_subdomains_double = sqrt(lon_lat_diff_ratio)*num_nonPole_region_y_subdomains_double;
	num_nonPole_region_y_subdomains = MIN(MAX(1, (int)num_nonPole_region_y_subdomains_double), num_remained_procs);
	num_nonPole_region_x_subdomains = num_remained_procs/num_nonPole_region_y_subdomains;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_nonPole_region_x_subdomains >= 1 && num_nonPole_region_y_subdomains >= 1 && num_nonPole_region_x_subdomains*num_nonPole_region_y_subdomains <= num_remained_procs && num_nonPole_region_x_subdomains*num_nonPole_region_y_subdomains*2 > num_remained_procs, "Software error in remapping_grid_domain_decom_engine::remapping_grid_domain_decom_engine");
	if (words_are_the_same(src_original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES)) {
		if (whole_domain_lat_diff/nonpole_region_max_width_fator > num_nonPole_region_y_subdomains)
			num_nonPole_region_y_subdomains = (int) (whole_domain_lat_diff/nonpole_region_max_width_fator) + 1;
		if (whole_domain_lon_diff/nonpole_region_max_width_fator > num_nonPole_region_x_subdomains)
			num_nonPole_region_x_subdomains = (int) (whole_domain_lon_diff/nonpole_region_max_width_fator) + 1;
	}

	nonPole_subdomain_x_side_length = whole_domain_lon_diff / num_nonPole_region_x_subdomains;
	nonPole_subdomain_y_side_length = whole_domain_lat_diff / num_nonPole_region_y_subdomains;
	num_total_subdomains = num_northPole_subdomain + num_southPole_subdomain + num_nonPole_region_x_subdomains*num_nonPole_region_y_subdomains;

	EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Decompose the common grid domain into %d (%d %d %dx%d) subdomains", num_total_subdomains, num_northPole_subdomain,num_southPole_subdomain,num_nonPole_region_x_subdomains, num_nonPole_region_y_subdomains);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "start to generate_subdomains_decomp_grid src");
	src_subdomain_decomp_grid = generate_subdomains_decomp_grid(src_original_cor_grid, src_subdomains_CoR_grids);
	dst_subdomain_decomp_grid = NULL;

	wtime(&time2);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "The time for generate_subdomains_decomp_grid src is %lf", time2 - time1);
	time1 = time2;

	if (dst_original_cor_grid != NULL) {
		dst_subdomain_decomp_grid = generate_subdomains_decomp_grid(dst_original_cor_grid, dst_subdomains_CoR_grids);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_subdomains_CoR_grids.size() == dst_subdomains_CoR_grids.size(), "Remapping_grid_domain_decomp_engine::Remapping_grid_domain_decomp_engine");
	}

	wtime(&time2);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "The time for generate_subdomains_decomp_grid dst is %lf", time2 - time1);

	confirm_or_generate_grid_vertexes(src_original_cor_grid, src_subdomains_CoR_grids, comp_node, true);
	if (dst_original_cor_grid != NULL)
		confirm_or_generate_grid_vertexes(dst_original_cor_grid, dst_subdomains_CoR_grids, comp_node, false);

	if (dst_original_cor_grid == NULL)
		return;

	wtime(&time3);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "The time for generating vertexes is %lf", time3 - time2);

	generate_grid_one_sided_comm_info(src_subdomains_CoR_grids, &src_grid_one_sided_comm_data_buf, &src_grid_MPI_Win);
	initialize_subdomains_halo_info(src_original_cor_grid, src_subdomains_CoR_grids, comp_node);
	wtime(&time1);
	src_grid_original_num_vertex = 4;
	for (int i = 0; i < num_total_subdomains; i ++) {
		if (i % num_total_procs == comp_node->get_current_proc_local_id()) {
			expand_subdomain_halo_grid(i, src_original_cor_grid, src_subdomains_CoR_grids, true);
			if (src_subdomains_expanded_CoR_grids[i/num_total_procs] != NULL) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_subdomains_expanded_CoR_grids[i/num_total_procs]->get_num_vertexes() > 0, "Software error in Remapping_grid_domain_decomp_engine::Remapping_grid_domain_decomp_engine");
				if (report_error_enabled)
					src_subdomains_expanded_CoR_grids[i/num_total_procs]->check_center_vertex_values_consistency_2D(true);
			}
		}
	}	

	wtime(&time2);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "The time for expanding subdomains is %lf vs %lf", time2 - time3, time1-time3);

	//destroy_grid_one_sided_comm_data_buf(true);
}


void Remapping_grid_domain_decomp_engine::clear_subdomains_halo_info()
{
	for (int i = 0; i < src_subdomains_expanded_CoR_grids.size(); i ++)
		if (src_subdomains_expanded_CoR_grids[i] != NULL)
			delete src_subdomains_expanded_CoR_grids[i];
	for (int i = 0; i < subdomains_final_triangles.size(); i ++) 
		if (subdomains_final_triangles[i] != NULL)
			delete [] subdomains_final_triangles[i];

	destroy_remote_subdomain_CoR_grids();
	subdomains_halo_num_levels.clear();
	subdomains_halo_subdomains_CoR_grids.clear();
	subdomains_halo_subdomains_IDs.clear();
	subdomains_halo_lon_bounds.clear();
	subdomains_halo_lat_bounds.clear();
	src_subdomains_expanded_CoR_grids.clear();
	subdomains_final_triangles.clear();
	local_subdomains_index.clear();
	num_subdomains_final_triangles.clear();
	remote_subdomain_CoR_grids.clear();
}


void Remapping_grid_domain_decomp_engine::initialize_subdomains_halo_info(Remap_grid_class *original_cor_grid, std::vector<std::pair<int, Remap_grid_class*> >& subdomains_CoR_grids, Comp_comm_group_mgt_node *comp_node)
{
	clear_subdomains_halo_info();

	for (int i = 0; i < subdomains_CoR_grids.size(); i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_CoR_grids[i].first%num_total_procs == comp_node->get_current_proc_local_id() && i == subdomains_CoR_grids[i].first/num_total_procs, "Software error in Remapping_grid_domain_decomp_engine::Remapping_grid_domain_decomp_engine");
		subdomains_halo_num_levels.push_back(0);
		std::vector<Remap_grid_class*> CoR_grids;
		std::vector<int> subdomains_IDs;
		std::pair<double, double> lon_lat_bounds;
		subdomains_halo_subdomains_CoR_grids.push_back(CoR_grids);
		subdomains_halo_subdomains_IDs.push_back(subdomains_IDs);
		subdomains_halo_lon_bounds.push_back(lon_lat_bounds);
		subdomains_halo_lat_bounds.push_back(lon_lat_bounds);
		src_subdomains_expanded_CoR_grids.push_back(NULL);
		calculate_subdomain_halo_bounds(original_cor_grid, subdomains_CoR_grids[i].first, 0, subdomains_halo_lon_bounds[subdomains_CoR_grids[i].first/num_total_procs], subdomains_halo_lat_bounds[subdomains_CoR_grids[i].first/num_total_procs]); 
	}
}


void Remapping_grid_domain_decomp_engine::confirm_or_generate_grid_vertexes(Remap_grid_class *original_cor_grid, std::vector<std::pair<int, Remap_grid_class*> >& subdomains_CoR_grids, Comp_comm_group_mgt_node *comp_node, bool is_src)
{
	int grid_original_num_vertex = is_src? src_grid_original_num_vertex : dst_grid_original_num_vertex;
	double time1, time2, time3, time4;


	if (grid_original_num_vertex > 0) {
		for (int i = 0; i < subdomains_CoR_grids.size(); i ++)
			if (subdomains_CoR_grids[i].second != NULL && subdomains_CoR_grids[i].second->get_grid_size() > 0)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_CoR_grids[i].second->get_num_vertexes() > 0, "Software error in Remapping_grid_domain_decomp_engine::confirm_or_generate_grid_vertexes");
		return;
	}

	wtime(&time1);

	for (int i = 0; i < subdomains_CoR_grids.size(); i ++)
		if (subdomains_CoR_grids[i].second != NULL)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_CoR_grids[i].second->get_num_vertexes() == 0, "Software error in Remapping_grid_domain_decomp_engine::confirm_or_generate_grid_vertexes");

	if (is_src)
		generate_grid_one_sided_comm_info(src_subdomains_CoR_grids, &src_grid_one_sided_comm_data_buf, &src_grid_MPI_Win);
	else generate_grid_one_sided_comm_info(dst_subdomains_CoR_grids, &dst_grid_one_sided_comm_data_buf, &dst_grid_MPI_Win);

	initialize_subdomains_halo_info(original_cor_grid, subdomains_CoR_grids, comp_node);

	triangulation_comm_data_buf = NULL;
	subdomain_triangulation_info = new int [subdomains_CoR_grids.size()*3];
	triangulation_comm_data_buf_max_size = 0;
	triangulation_comm_data_buf_content_size = 0;
	write_data_into_array_buffer(subdomain_triangulation_info, subdomains_CoR_grids.size()*3*sizeof(int), &triangulation_comm_data_buf, triangulation_comm_data_buf_max_size, triangulation_comm_data_buf_content_size);
	
	for (int i = 0; i < num_total_subdomains; i ++) {
		if (i % num_total_procs == comp_node->get_current_proc_local_id()) {
			int subdomain_index_intra_proc = i / comp_node->get_num_procs();
			subdomain_triangulation_info[subdomain_index_intra_proc*3+0] = i;
			subdomain_triangulation_info[subdomain_index_intra_proc*3+1] = 0;
			subdomain_triangulation_info[subdomain_index_intra_proc*3+2] = 0;
			subdomains_final_triangles.push_back(NULL);
			local_subdomains_index.push_back(i);
			num_subdomains_final_triangles.push_back(0);
			expand_subdomain_halo_grid(i, original_cor_grid, subdomains_CoR_grids, is_src);
		}
	}

	long temp_long = 0;
	MPI_Win_create(triangulation_comm_data_buf, triangulation_comm_data_buf_content_size, sizeof(char), MPI_INFO_NULL, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comm_group(), &triangulation_comm_win);
	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, comp_node->get_current_proc_local_id(), 0, triangulation_comm_win);
	write_data_into_array_buffer(subdomain_triangulation_info, subdomains_CoR_grids.size()*3*sizeof(int), &triangulation_comm_data_buf, triangulation_comm_data_buf_max_size, temp_long);
	MPI_Win_unlock(comp_node->get_current_proc_local_id(), triangulation_comm_win);

	temp_long = 0;

	MPI_Barrier(comp_node->get_comm_group());

	bool is_sphere_grid = original_cor_grid->get_is_sphere_grid() && words_are_the_same(original_cor_grid->get_sphere_grid_coord_unit(), COORD_UNIT_DEGREES);
	for (int current_subdomain_index = 0; current_subdomain_index < local_subdomains_index.size(); current_subdomain_index ++)
		check_halo_subdomain_triangulation_checksum(original_cor_grid, current_subdomain_index, is_sphere_grid, comp_node);
	
	MPI_Win_free(&triangulation_comm_win);

	if (triangulation_comm_data_buf != NULL)
		delete [] triangulation_comm_data_buf;
	if (subdomain_triangulation_info != NULL)
		delete [] subdomain_triangulation_info;

	destroy_grid_one_sided_comm_data_buf(is_src);
	destroy_remote_subdomain_CoR_grids();
	clear_subdomains_halo_info();

	for (int i = 0; i < subdomains_CoR_grids.size(); i ++)
		if (subdomains_CoR_grids[i].second != NULL)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_CoR_grids[i].second->get_num_vertexes() > 0, "Software error in Remapping_grid_domain_decomp_engine::confirm_or_generate_grid_vertexes");

	wtime(&time4);
	EXECUTION_REPORT_LOG(REPORT_LOG, comp_node->get_comp_id(), true, "Spend %lf seconds for distributed generating vertexes for grid \"%s\" with %d cells", time4-time1, original_cor_grid->get_grid_name(), original_cor_grid->get_grid_size());
}


void Remapping_grid_domain_decomp_engine::check_halo_subdomain_triangulation_checksum(Remap_grid_class *original_cor_grid, int current_subdomain_index, bool is_sphere_grid, Comp_comm_group_mgt_node *comp_node)
{
	std::vector<int> subdomain_halo_subdomain_IDs;
	int remote_subdomain_owner_id, subdomain_index_intra_remote_process;
	int halo_subdomain_triangulation_info[3];
	int temp_int, point_index, triangle_index;
	Triangle_inline* remote_result_triangles;

	get_subdomain_halo_subdomains_IDs(original_cor_grid, local_subdomains_index[current_subdomain_index], 1, subdomain_halo_subdomain_IDs);
	for(int halo_index = 0; halo_index < subdomain_halo_subdomain_IDs.size(); halo_index ++ ) {
		remote_subdomain_owner_id = subdomain_halo_subdomain_IDs[halo_index] % comp_node->get_num_procs();
		subdomain_index_intra_remote_process = subdomain_halo_subdomain_IDs[halo_index] / comp_node->get_num_procs();
		MPI_Win_lock(MPI_LOCK_SHARED, remote_subdomain_owner_id, 0, triangulation_comm_win);
		MPI_Get(halo_subdomain_triangulation_info, sizeof(int)*3, MPI_CHAR, remote_subdomain_owner_id, subdomain_index_intra_remote_process*sizeof(int)*3, sizeof(int)*3, MPI_CHAR, triangulation_comm_win);
		MPI_Win_unlock(remote_subdomain_owner_id, triangulation_comm_win);

		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_halo_subdomain_IDs[halo_index] == halo_subdomain_triangulation_info[0], "Software error1 in Remapping_grid_domain_decomp_engine::Remapping_grid_domain_decomp_engine, %d and %d", subdomain_halo_subdomain_IDs[halo_index], halo_subdomain_triangulation_info[0]);

		if(halo_subdomain_triangulation_info[1] > 0) {
			long len_remote_result_triangle_points_buffer = halo_subdomain_triangulation_info[1]*(6*sizeof(double)+3*sizeof(int));
			char *remote_result_triangle_points_buffer = new char [len_remote_result_triangle_points_buffer];
			MPI_Win_lock(MPI_LOCK_SHARED, remote_subdomain_owner_id, 0, triangulation_comm_win);
			MPI_Get(remote_result_triangle_points_buffer, len_remote_result_triangle_points_buffer, MPI_CHAR, remote_subdomain_owner_id, halo_subdomain_triangulation_info[2], len_remote_result_triangle_points_buffer, MPI_CHAR, triangulation_comm_win);
			MPI_Win_unlock(remote_subdomain_owner_id, triangulation_comm_win);

			remote_result_triangles = new Triangle_inline [halo_subdomain_triangulation_info[1]];
			for(triangle_index = halo_subdomain_triangulation_info[1] - 1; triangle_index >= 0; triangle_index --) {
				for (point_index = 2; point_index >= 0; point_index --) {
					read_data_from_array_buffer(&(remote_result_triangles[triangle_index].v[point_index].y), sizeof(double), remote_result_triangle_points_buffer, len_remote_result_triangle_points_buffer, true);
					read_data_from_array_buffer(&(remote_result_triangles[triangle_index].v[point_index].x), sizeof(double), remote_result_triangle_points_buffer, len_remote_result_triangle_points_buffer, true);
					read_data_from_array_buffer(&temp_int, sizeof(int), remote_result_triangle_points_buffer, len_remote_result_triangle_points_buffer, true);
					remote_result_triangles[triangle_index].v[point_index].id = temp_int;
				}
				if (is_sphere_grid) {
					for (point_index = 0; point_index < 3; point_index ++) {
						if ( fabs(remote_result_triangles[triangle_index].v[point_index].x - remote_result_triangles[triangle_index].v[(point_index+1)%3].x) > PAT_CYCLIC_EDGE_THRESHOLD) {
							remote_result_triangles[triangle_index].is_cyclic = true;
							break;
						}
					}
					if (point_index == 3)
						remote_result_triangles[triangle_index].is_cyclic = false;
				}
			}

			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, len_remote_result_triangle_points_buffer == 0, "the len_remote_result_triangle_points_buffer is wrong, %d.", len_remote_result_triangle_points_buffer);

			check_boundary_triangles_checksum(original_cor_grid, current_subdomain_index, subdomain_halo_subdomain_IDs[halo_index], remote_result_triangles, halo_subdomain_triangulation_info[1]);

			delete [] remote_result_triangles;
			delete [] remote_result_triangle_points_buffer;
		}
	}
}


void Remapping_grid_domain_decomp_engine::check_boundary_triangles_checksum(Remap_grid_class *original_cor_grid, int current_subdomain_index, int halo_subdomain_id, Triangle_inline* remote_result_triangles, int num_remote_result_triangles)
{
	int subdomain_id = local_subdomains_index[current_subdomain_index];
	PatCC_Point *head, *tail;
	unsigned long checksum_subdomain, checksum_remote_subdomain;
	double subdomain_min_lon, subdomain_max_lon, subdomain_min_lat, subdomain_max_lat;
	double halo_subdomain_min_lon, halo_subdomain_max_lon, halo_subdomain_min_lat, halo_subdomain_max_lat;
	std::pair<double, double> subdomain_bound_lon, subdomain_bound_lat;
	bool is_sphere_grid = original_cor_grid->get_is_sphere_grid() && words_are_the_same(original_cor_grid->get_sphere_grid_coord_unit(), COORD_UNIT_DEGREES);
	
	calculate_subdomain_halo_bounds(original_cor_grid, halo_subdomain_id, 0, subdomain_bound_lon, subdomain_bound_lat); 
	halo_subdomain_min_lon = subdomain_bound_lon.first;
	halo_subdomain_max_lon = subdomain_bound_lon.second;
	halo_subdomain_min_lat = subdomain_bound_lat.first;
	halo_subdomain_max_lat = subdomain_bound_lat.second;
	calculate_subdomain_halo_bounds(original_cor_grid, subdomain_id, 0, subdomain_bound_lon, subdomain_bound_lat); 
	subdomain_min_lon = subdomain_bound_lon.first;
	subdomain_max_lon = subdomain_bound_lon.second;
	subdomain_min_lat = subdomain_bound_lat.first;
	subdomain_max_lat = subdomain_bound_lat.second;

	if(is_north_pole_subdomain(subdomain_id)) {
		head = new PatCC_Point(halo_subdomain_min_lon, halo_subdomain_max_lat, false);
		tail = new PatCC_Point(halo_subdomain_max_lon, halo_subdomain_max_lat, false);
	}
	else if(is_south_pole_subdomain(subdomain_id)) {
		head = new PatCC_Point(halo_subdomain_min_lon, halo_subdomain_min_lat, false);
		tail = new PatCC_Point(halo_subdomain_max_lon, halo_subdomain_min_lat, false);
	}
	else {
		int subdomain_x_index, subdomain_y_index, halo_x_index, halo_y_index;
		get_subdomain_XY_index(subdomain_id, subdomain_x_index, subdomain_y_index);
		get_subdomain_XY_index(halo_subdomain_id, halo_x_index, halo_y_index);

		if(subdomain_x_index == halo_x_index && subdomain_y_index == halo_y_index-1) {// UP
			head = new PatCC_Point(subdomain_min_lon, subdomain_max_lat, false);
			tail = new PatCC_Point(subdomain_max_lon, subdomain_max_lat, false);
		}
		else if(subdomain_x_index == halo_x_index && subdomain_y_index == halo_y_index+1) {// DOWN
			head = new PatCC_Point(subdomain_min_lon, subdomain_min_lat, false);
			tail = new PatCC_Point(subdomain_max_lon, subdomain_min_lat, false);
		}
		else if(subdomain_x_index == halo_x_index+1 && subdomain_y_index == halo_y_index) { // LEFT
			head = new PatCC_Point(subdomain_min_lon, subdomain_min_lat, false);
			tail = new PatCC_Point(subdomain_min_lon, subdomain_max_lat, false);
		}
		else if (subdomain_x_index == halo_x_index-1 && subdomain_y_index == halo_y_index) { // RIGHT
			head = new PatCC_Point(subdomain_max_lon, subdomain_min_lat, false);
			tail = new PatCC_Point(subdomain_max_lon, subdomain_max_lat, false);
		}
		else {
			head = new PatCC_Point(0.0, 0.0, false);
			tail = new PatCC_Point(0.0, 0.0, false);
		}
	}

	checksum_subdomain = cal_triangles_checksum_on_boundary(subdomains_final_triangles[current_subdomain_index], num_subdomains_final_triangles[current_subdomain_index], *head, *tail, subdomain_id);
	checksum_remote_subdomain = cal_triangles_checksum_on_boundary(remote_result_triangles, num_remote_result_triangles, *head, *tail, halo_subdomain_id);

	EXECUTION_REPORT(REPORT_ERROR, -1, checksum_subdomain == checksum_remote_subdomain, "Inconsitent distributed triangulation detected. Please try through decreasing the number of processes of component models or contact Dr. Li Liu (liuli-cess@tsinghua.edu.cn).\n Wrong checksum is %lu(%d) %lu(%d), (%lf %lf %lf %lf)(%lf %lf %lf %lf) (%lf %lf) (%lf %lf)", checksum_subdomain, local_subdomains_index[current_subdomain_index], checksum_remote_subdomain, halo_subdomain_id, subdomain_min_lon, subdomain_max_lon, subdomain_min_lat, subdomain_max_lat, halo_subdomain_min_lon, halo_subdomain_max_lon, halo_subdomain_min_lat, halo_subdomain_max_lat, head->x, head->y, tail->x, tail->y);

	delete head;
	delete tail;
}


Remapping_grid_domain_decomp_engine::~Remapping_grid_domain_decomp_engine()
{
	delete src_subdomain_decomp_grid;
	if (dst_original_cor_grid != NULL)
		delete dst_subdomain_decomp_grid;

	for (int i = 0; i < src_subdomains_CoR_grids.size(); i ++)
		if (src_subdomains_CoR_grids[i].second != NULL)
			delete src_subdomains_CoR_grids[i].second;

	for (int i = 0; i < dst_subdomains_CoR_grids.size(); i ++)
		if (dst_subdomains_CoR_grids[i].second != NULL)
			delete dst_subdomains_CoR_grids[i].second;

	destroy_grid_one_sided_comm_data_buf(true);
	destroy_grid_one_sided_comm_data_buf(false);
	destroy_remote_subdomain_CoR_grids();
	clear_subdomains_halo_info();
}


void Remapping_grid_domain_decomp_engine::destroy_grid_one_sided_comm_data_buf(bool is_src)
{
	if (is_src && src_grid_one_sided_comm_data_buf != NULL) {
		MPI_Win_free(&src_grid_MPI_Win);
		delete [] src_grid_one_sided_comm_data_buf;
		src_grid_one_sided_comm_data_buf = NULL;
	}
	if (!is_src && dst_grid_one_sided_comm_data_buf != NULL) {
		MPI_Win_free(&dst_grid_MPI_Win);
		delete [] dst_grid_one_sided_comm_data_buf;
		dst_grid_one_sided_comm_data_buf = NULL;
	}
}


void Remapping_grid_domain_decomp_engine::destroy_remote_subdomain_CoR_grids()
{
	for (int i = 0; i < remote_subdomain_CoR_grids.size(); i ++)
		if (remote_subdomain_CoR_grids[i] != NULL)
			delete remote_subdomain_CoR_grids[i];
	remote_subdomain_CoR_grids.clear();
}

int Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID(double lon_value, double lat_value)
{
	int x_index, y_index;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lat_value != NULL_COORD_VALUE && lon_value != NULL_COORD_VALUE, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID1");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lat_value >= whole_domain_min_lat && lat_value <= whole_domain_max_lat, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID2: %lf vs (%lf, %lf) ", lat_value, whole_domain_min_lat, whole_domain_max_lat);
	if (whole_domain_min_lon < whole_domain_max_lon) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lon_value >= whole_domain_min_lon && lon_value <= whole_domain_max_lon, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID3");
	}
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lon_value >= whole_domain_min_lon && lon_value <= 360.0 || lon_value >= 0.0 && lon_value <= whole_domain_max_lon, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID4");
		
	if (lat_value > nonPole_region_max_lat) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_northPole_subdomain > 0, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID5");
		return 0;
	}
	if (lat_value < nonPole_region_min_lat) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_southPole_subdomain > 0, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID6");
		return num_northPole_subdomain;
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_nonPole_region_x_subdomains > 0 && num_nonPole_region_y_subdomains > 0, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID7");
	y_index = (int) ((lat_value - nonPole_region_min_lat) / nonPole_subdomain_y_side_length);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, y_index >= 0 && y_index <= num_nonPole_region_y_subdomains, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID8");
	y_index = MIN(y_index, num_nonPole_region_y_subdomains - 1);

	if (whole_domain_min_lon > whole_domain_max_lon && lon_value < whole_domain_min_lon)
		x_index = (int) ((lon_value + 360.0 - whole_domain_min_lon) / nonPole_subdomain_x_side_length);
	else x_index = (int) ((lon_value - whole_domain_min_lon) / nonPole_subdomain_x_side_length);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, x_index >= 0 && x_index <= num_nonPole_region_x_subdomains, "Software error in Remapping_grid_domain_decomp_engine::calculate_grid_point_subdomain_ID9");
	x_index = MIN(x_index, num_nonPole_region_x_subdomains - 1);
	
	return num_southPole_subdomain + num_northPole_subdomain + y_index*num_nonPole_region_x_subdomains+x_index;
}


Decomp_grid_info *Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid(Remap_grid_class *original_cor_grid, std::vector<std::pair<int, Remap_grid_class*> > &subdomains_CoR_grids)
{
	Comp_comm_group_mgt_node *host_comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	Distributed_H2D_grid_engine *distributed_H2D_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(host_comp_node->get_comp_full_name(), original_cor_grid);
	if (distributed_H2D_grid == NULL)
		distributed_H2D_grid = distributed_H2D_grid_mgr->search_distributed_H2D_grid(NULL, original_cor_grid);
	Decomp_grid_info *basic_decomp_grid = distributed_H2D_grid->get_basic_decomp_grid();
	Decomp_grid_info *subdomains_decomp_grid;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, basic_decomp_grid != NULL && decomps_info_mgr->get_decomp_info(basic_decomp_grid->get_decomp_id())->get_host_comp_id() == comp_id, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
	Remap_grid_class *decomp_cor_grid = basic_decomp_grid->get_decomp_grid();
	const int *basic_local_cell_global_indx = decomps_info_mgr->get_decomp_info(basic_decomp_grid->get_decomp_id())->get_local_cell_global_indx();
	double *decomp_lon_values = (double*) decomp_cor_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf;
	double *decomp_lat_values = (double*) decomp_cor_grid->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, decomp_cor_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->required_data_size == decomp_cor_grid->get_grid_size(), "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, decomp_cor_grid->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->required_data_size == decomp_cor_grid->get_grid_size(), "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
	Grid_point_subdomain_mapping_info *grid_points_subdomain_map = NULL;
	Grid_point_subdomain_mapping_info *grid_points_whole_map = NULL;
	int subdomains_decomp_id, total_num_cells_in_local_subdomains = 0;
	int *cell_indexes_in_local_subdomains = NULL, *cell_indexes_in_local_subdomains_tmp = NULL;
	long global_size;
	char decomp_name[4*NAME_STR_SIZE], subdomain_decomp_name[4*NAME_STR_SIZE], src_grid_name[NAME_STR_SIZE];
	Distribute_merge_sort<Grid_point_subdomain_mapping_info> *Grid_point_subdomain_mapping_info_sort = new Distribute_merge_sort<Grid_point_subdomain_mapping_info>(-1, -1, NULL, NULL);
	common_sort_struct<Grid_point_subdomain_mapping_info> *Grid_point_subdomain_mapping_info_data = NULL;
	int num_procs_in_sorting = Grid_point_subdomain_mapping_info_sort->calculate_max_power2(host_comp_node->get_num_procs());	
	int num_local_sorted_data = decomp_cor_grid->get_grid_size();
	int send_recv_mark, proc_id_send_to, proc_id_recv_from;
	int num_local_subdomains, *num_cells_in_local_subdomains, *displ_cells_in_local_subdomains;
	std::pair<double, double> lon_bounds, lat_bounds;
	double time1, time2, time3;


	if (decomp_cor_grid->get_grid_size() > 0)
		Grid_point_subdomain_mapping_info_data = new common_sort_struct<Grid_point_subdomain_mapping_info> [decomp_cor_grid->get_grid_size()];
	for (int i = 0; i < decomp_cor_grid->get_grid_size(); i ++) {
		Grid_point_subdomain_mapping_info_data[i].content.global_grid_cell_index = basic_local_cell_global_indx[i];
		Grid_point_subdomain_mapping_info_data[i].content.subdomain_id = calculate_grid_point_subdomain_ID(decomp_lon_values[i], decomp_lat_values[i]);
		Grid_point_subdomain_mapping_info_data[i].content.owner_proc_id = Grid_point_subdomain_mapping_info_data[i].content.subdomain_id % num_total_procs;
		Grid_point_subdomain_mapping_info_data[i].target_proc_id = Grid_point_subdomain_mapping_info_data[i].content.owner_proc_id;
	}
	Grid_point_subdomain_mapping_info_sort->do_data_sorting_with_target_process_id(&Grid_point_subdomain_mapping_info_data, decomps_info_mgr->get_decomp_info(basic_decomp_grid->get_decomp_id()), host_comp_node, &num_local_sorted_data);
	if (num_local_sorted_data > 0)
		cell_indexes_in_local_subdomains = new int [num_local_sorted_data];
	for (int i = 0; i < num_local_sorted_data; i ++) 
		Grid_point_subdomain_mapping_info_data[i].key = Grid_point_subdomain_mapping_info_data[i].content.global_grid_cell_index+1;
	do_quick_sort(Grid_point_subdomain_mapping_info_data, (int*)NULL, 0, num_local_sorted_data-1);		
	for (int i = 0; i < num_local_sorted_data; i ++) 
		cell_indexes_in_local_subdomains[i] = Grid_point_subdomain_mapping_info_data[i].content.global_grid_cell_index+1;	

	if (report_error_enabled) {
		if (decomp_cor_grid->get_grid_size() > 0)
			grid_points_subdomain_map = new Grid_point_subdomain_mapping_info [decomp_cor_grid->get_grid_size()];
		for (int i = 0; i < decomp_cor_grid->get_grid_size(); i ++) {
			grid_points_subdomain_map[i].global_grid_cell_index = basic_local_cell_global_indx[i];
			grid_points_subdomain_map[i].subdomain_id = calculate_grid_point_subdomain_ID(decomp_lon_values[i], decomp_lat_values[i]);
			grid_points_subdomain_map[i].owner_proc_id = grid_points_subdomain_map[i].subdomain_id % num_total_procs;
		}
		gather_array_in_one_comp(host_comp_node->get_num_procs(), host_comp_node->get_current_proc_local_id(), grid_points_subdomain_map, decomp_cor_grid->get_grid_size(), sizeof(Grid_point_subdomain_mapping_info), NULL, (void**)(&grid_points_whole_map), global_size, host_comp_node->get_comm_group());
		if (host_comp_node->get_current_proc_local_id() == 0)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_cor_grid->get_grid_size() == global_size, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
		global_size *= sizeof(Grid_point_subdomain_mapping_info);
		bcast_array_in_one_comp(host_comp_node->get_current_proc_local_id(), (char**)(&grid_points_whole_map), global_size, host_comp_node->get_comm_group());
		for (int i = 0; i < original_cor_grid->get_grid_size(); i ++) {
			if (grid_points_whole_map[i].owner_proc_id == host_comp_node->get_current_proc_local_id())
				total_num_cells_in_local_subdomains ++;
		}
		if (total_num_cells_in_local_subdomains > 0)
			cell_indexes_in_local_subdomains_tmp = new int [total_num_cells_in_local_subdomains];
		total_num_cells_in_local_subdomains = 0;
		for (int i = 0; i < original_cor_grid->get_grid_size(); i ++)
			if (grid_points_whole_map[i].owner_proc_id == host_comp_node->get_current_proc_local_id())
				cell_indexes_in_local_subdomains_tmp[total_num_cells_in_local_subdomains++] = grid_points_whole_map[i].global_grid_cell_index + 1;
		EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The grid \"%s\" has %d local grid cells after subdomain decomposition for remapping", original_cor_grid->get_grid_name(), total_num_cells_in_local_subdomains);
		do_quick_sort(cell_indexes_in_local_subdomains_tmp, (int*)NULL, 0, total_num_cells_in_local_subdomains-1);	
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, total_num_cells_in_local_subdomains == num_local_sorted_data, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid: %d vs %d", total_num_cells_in_local_subdomains, num_local_sorted_data);
		for (int i = 0; i < num_local_sorted_data; i ++) 
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cell_indexes_in_local_subdomains_tmp[i] == cell_indexes_in_local_subdomains[i], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid: %d: %d vs %d", i, cell_indexes_in_local_subdomains_tmp[i], cell_indexes_in_local_subdomains[i]);
		if (grid_points_subdomain_map != NULL)
			delete [] grid_points_subdomain_map;
		delete [] grid_points_whole_map;
		if (cell_indexes_in_local_subdomains_tmp != NULL)
			delete [] cell_indexes_in_local_subdomains_tmp;
		cell_indexes_in_local_subdomains_tmp = NULL;
	}

	if (num_local_sorted_data > 0)
		cell_indexes_in_local_subdomains_tmp = new int [num_local_sorted_data];
	for (int i = 0; i < num_local_sorted_data; i ++) 
		Grid_point_subdomain_mapping_info_data[i].key = Grid_point_subdomain_mapping_info_data[i].content.subdomain_id;
	do_quick_sort(Grid_point_subdomain_mapping_info_data, (int*)NULL, 0, num_local_sorted_data-1);		
	for (int i = 0; i < num_local_sorted_data; i ++) {
		cell_indexes_in_local_subdomains[i] = Grid_point_subdomain_mapping_info_data[i].content.global_grid_cell_index+1;	
		cell_indexes_in_local_subdomains_tmp[i] = cell_indexes_in_local_subdomains[i];
	}

	delete Grid_point_subdomain_mapping_info_sort;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->search_grid_info(src_original_grid_id) != NULL, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->search_grid_info(src_original_grid_id)->get_H2D_sub_grid() != NULL, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid: %s  %x", original_grid_mgr->search_grid_info(src_original_grid_id)->get_grid_name(), original_grid_mgr->search_grid_info(src_original_grid_id)->get_H2D_sub_CoR_grid());
	sprintf(src_grid_name, "%s@%s", original_grid_mgr->search_grid_info(src_original_grid_id)->get_H2D_sub_grid()->get_grid_name(), original_grid_mgr->search_grid_info(src_original_grid_id)->get_H2D_sub_grid()->get_comp_full_name());
	if (dst_original_cor_grid != NULL) {
		if (original_cor_grid == src_original_cor_grid)
			sprintf(decomp_name, "decomp_of_%s_remapping_to_%s_at_%s", src_grid_name, dst_original_cor_grid->get_grid_name(), host_comp_node->get_full_name());
		else sprintf(decomp_name, "decomp_of_%s_remapped_from_%s_at_%s", dst_original_cor_grid->get_grid_name(), src_grid_name, host_comp_node->get_full_name());
	}
	else sprintf(decomp_name, "internal_decomp_of_%s_at_%s", src_grid_name, host_comp_node->get_full_name());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, distributed_H2D_grid->get_original_grid()->get_original_CoR_grid() == original_cor_grid, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
	if (decomps_info_mgr->search_decomp_info(decomp_name, distributed_H2D_grid->get_original_grid()->get_comp_id()) != NULL)
		subdomains_decomp_id = decomps_info_mgr->search_decomp_info(decomp_name, distributed_H2D_grid->get_original_grid()->get_comp_id())->get_decomp_id();
	else subdomains_decomp_id = decomps_info_mgr->register_H2D_parallel_decomposition(decomp_name, distributed_H2D_grid->get_original_grid()->get_grid_id(), decomps_info_mgr->get_decomp_info(basic_decomp_grid->get_decomp_id())->get_host_comp_id(), num_local_sorted_data, cell_indexes_in_local_subdomains, 0, NULL, "");
		
	EXECUTION_REPORT(REPORT_LOG, -1, true, "start to generate_decomp_grid_via_data_transfer");
	const char *comp_full_name = distributed_H2D_grid->get_original_grid()->get_comp_full_name();
	wtime(&time1);
	subdomains_decomp_grid = distributed_H2D_grid->generate_decomp_grid_via_data_transfer(host_comp_node->get_comp_full_name(), comp_full_name, NULL, host_comp_node->get_full_name(), comp_full_name, subdomains_decomp_id, original_cor_grid);
	wtime(&time2);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "the time for generate_decomp_grid_via_data_transfer is %lf", time2 - time1);

	num_local_subdomains = num_total_subdomains/num_total_procs;
	if (num_total_subdomains%num_total_procs > host_comp_node->get_current_proc_local_id())
		num_local_subdomains ++;
	num_cells_in_local_subdomains = new int [num_local_subdomains];
	displ_cells_in_local_subdomains = new int [num_local_subdomains];
	for (int i = 0; i < num_local_subdomains; i ++)
		num_cells_in_local_subdomains[i] = 0;
	for (int i = 0; i < num_local_sorted_data; i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, Grid_point_subdomain_mapping_info_data[i].content.subdomain_id%num_total_procs == host_comp_node->get_current_proc_local_id() && Grid_point_subdomain_mapping_info_data[i].content.subdomain_id/num_total_procs < num_local_subdomains, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
		num_cells_in_local_subdomains[Grid_point_subdomain_mapping_info_data[i].content.subdomain_id/num_total_procs] ++;
	}
	for (int i=0, j=0; i < num_local_subdomains; i ++) {
		displ_cells_in_local_subdomains[i] = j;
		j += num_cells_in_local_subdomains[i];
		num_cells_in_local_subdomains[i] = 0;
	}
	for (int i = 0; i < num_local_sorted_data; i ++) {
		int local_subdomain_id = Grid_point_subdomain_mapping_info_data[i].content.subdomain_id/num_total_procs;
		cell_indexes_in_local_subdomains[displ_cells_in_local_subdomains[local_subdomain_id]+num_cells_in_local_subdomains[local_subdomain_id]] = i;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cell_indexes_in_local_subdomains_tmp[displ_cells_in_local_subdomains[local_subdomain_id]+num_cells_in_local_subdomains[local_subdomain_id]] == Grid_point_subdomain_mapping_info_data[i].content.global_grid_cell_index+1, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
		num_cells_in_local_subdomains[local_subdomain_id] ++;
	}

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Will push %d subdomain grids", num_local_subdomains);

	for (int i = 0; i < num_local_subdomains; i ++) {
		std::pair<int, Remap_grid_class*> subdomain_CoR_grid;
		if (num_cells_in_local_subdomains[i] == 0) {
			subdomain_CoR_grid.first = host_comp_node->get_current_proc_local_id() + i*num_total_procs;
			subdomain_CoR_grid.second = NULL;
		} else {
			subdomain_CoR_grid.first = Grid_point_subdomain_mapping_info_data[displ_cells_in_local_subdomains[i]].content.subdomain_id;
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_CoR_grid.first%num_total_procs == host_comp_node->get_current_proc_local_id(), "Software error in Remapping_grid_domain_decomp_engine::generate_subdomains_decomp_grid");
			sprintf(subdomain_decomp_name, "%s_at_subdomain_%d", decomp_name, Grid_point_subdomain_mapping_info_data[displ_cells_in_local_subdomains[i]].content.subdomain_id);
			subdomain_CoR_grid.second = subdomains_decomp_grid->get_decomp_grid()->generate_decomp_grid(cell_indexes_in_local_subdomains+displ_cells_in_local_subdomains[i], num_cells_in_local_subdomains[i], subdomain_decomp_name);
			subdomain_CoR_grid.second->set_local_cell_global_indexes(cell_indexes_in_local_subdomains_tmp+displ_cells_in_local_subdomains[i], true);
			calculate_subdomain_halo_bounds(original_cor_grid, subdomain_CoR_grid.first, 0, lon_bounds, lat_bounds);
			subdomain_CoR_grid.second->set_grid_boundary(lon_bounds.first, lon_bounds.second, lat_bounds.first, lat_bounds.second);
		}
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Push subdomain grid %d (%d) with displacement %d with size %d", subdomain_CoR_grid.first, i, displ_cells_in_local_subdomains[i], num_cells_in_local_subdomains[i]);
		subdomains_CoR_grids.push_back(subdomain_CoR_grid);
	}

	if (cell_indexes_in_local_subdomains != NULL)
		delete [] cell_indexes_in_local_subdomains;
	if (cell_indexes_in_local_subdomains_tmp != NULL)
		delete [] cell_indexes_in_local_subdomains_tmp;
	if (Grid_point_subdomain_mapping_info_data != NULL)
		delete [] Grid_point_subdomain_mapping_info_data;
	delete [] num_cells_in_local_subdomains;
	delete [] displ_cells_in_local_subdomains;

	return subdomains_decomp_grid;
}


void Remapping_grid_domain_decomp_engine::generate_grid_one_sided_comm_info(std::vector<std::pair<int, Remap_grid_class*> >& subdomains_CoR_grids, char **grid_one_sided_comm_data_buf, MPI_Win *grid_MPI_Win)
{
	int last_displ, current_displ = 0, *subdomains_grids_array_info = new int [subdomains_CoR_grids.size()*3];
	char *temp_array_buffer = NULL;
	long buffer_max_size, buffer_content_size = 0;
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);


	if ((*grid_one_sided_comm_data_buf) != NULL)
		return;

	for (int i = 0; i < subdomains_CoR_grids.size(); i ++) {
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "the subdomain id is %d , index %d, num vertexes: %d.", subdomains_CoR_grids[i].first, i, subdomains_CoR_grids[i].second == NULL? 0 : subdomains_CoR_grids[i].second->get_num_vertexes());
		subdomains_grids_array_info[i*3+0] = subdomains_CoR_grids[i].first;
		subdomains_grids_array_info[i*3+2] = current_displ + subdomains_CoR_grids.size()*3*sizeof(int);

		if (subdomains_CoR_grids[i].second == NULL) {
			subdomains_grids_array_info[i*3+1] = 0;
		}
		else {
			last_displ = buffer_content_size;
			subdomains_CoR_grids[i].second->write_grid_into_array(&temp_array_buffer, buffer_max_size, buffer_content_size);
			current_displ = buffer_content_size;
			subdomains_grids_array_info[i*3+1] = current_displ - last_displ;
		}
	}
	// grid information compression
	*grid_one_sided_comm_data_buf = new char [subdomains_CoR_grids.size()*3*sizeof(int)+buffer_content_size];
	MPI_Win_create(*grid_one_sided_comm_data_buf, subdomains_CoR_grids.size()*3*sizeof(int)+buffer_content_size, sizeof(char), MPI_INFO_NULL, comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_comm_group(), grid_MPI_Win);
	MPI_Win_lock(MPI_LOCK_EXCLUSIVE, comp_node->get_current_proc_local_id(), 0, *grid_MPI_Win);
	memcpy(*grid_one_sided_comm_data_buf, subdomains_grids_array_info, subdomains_CoR_grids.size()*3*sizeof(int));
	memcpy(*grid_one_sided_comm_data_buf+subdomains_CoR_grids.size()*3*sizeof(int), temp_array_buffer, buffer_content_size);
	MPI_Win_unlock(comp_node->get_current_proc_local_id(), *grid_MPI_Win);
	MPI_Barrier(comp_node->get_comm_group());
	delete [] subdomains_grids_array_info;
	if (temp_array_buffer != NULL)
		delete [] temp_array_buffer;
}


Remap_grid_class *Remapping_grid_domain_decomp_engine::get_subdomain_grid(int subdomain_id, std::vector<std::pair<int, Remap_grid_class*> >& subdomains_CoR_grids, bool is_src)
{
	MPI_Win grid_MPI_Win = is_src ? src_grid_MPI_Win : dst_grid_MPI_Win;
	Remap_grid_class *subdomain_grid = NULL;
	int owner_process_id, subdomain_index_intra_owner_process;
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	double time1, time2, time3, time4, time5, time6, time7;

	wtime(&time1);
	owner_process_id = subdomain_id % comp_node->get_num_procs();
	subdomain_index_intra_owner_process = subdomain_id / comp_node->get_num_procs();
	if (owner_process_id == comp_node->get_current_proc_local_id()) {
		subdomain_grid = subdomains_CoR_grids[subdomain_index_intra_owner_process].second;
		if (subdomain_grid != NULL)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_CoR_grids[subdomain_index_intra_owner_process].first == subdomain_id, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_grid");
 	}
	else {
		int subdomain_transfer_info[3];
		MPI_Win_lock(MPI_LOCK_SHARED, owner_process_id, 0, grid_MPI_Win);
		wtime(&time6);
		MPI_Get(subdomain_transfer_info, sizeof(int)*3, MPI_CHAR, owner_process_id, subdomain_index_intra_owner_process*sizeof(int)*3, sizeof(int)*3, MPI_CHAR, grid_MPI_Win);
		wtime(&time7);
		MPI_Win_unlock(owner_process_id, grid_MPI_Win);
		wtime(&time4);
		if (subdomain_transfer_info[1] > 0) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_transfer_info[0] == subdomain_id, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_grid: %d: %d vs %d", owner_process_id, subdomain_transfer_info[0], subdomain_id);
			char *remote_subdomain_grid_data_buf = new char [subdomain_transfer_info[1]];
			long buf_size = subdomain_transfer_info[1];
			wtime(&time5);
			MPI_Win_lock(MPI_LOCK_SHARED, owner_process_id, 0, grid_MPI_Win);
			MPI_Get(remote_subdomain_grid_data_buf, subdomain_transfer_info[1], MPI_CHAR, owner_process_id, subdomain_transfer_info[2], subdomain_transfer_info[1], MPI_CHAR, grid_MPI_Win);
			MPI_Win_unlock(owner_process_id, grid_MPI_Win);
			wtime(&time2);
			subdomain_grid = new Remap_grid_class(NULL, NULL, remote_subdomain_grid_data_buf, buf_size, false);
			wtime(&time3);
			EXECUTION_REPORT(REPORT_LOG, -1, true, "MPI get for %d bytes use %lf vs %lf vs %lf vs %lf vs %lf vs %lf vs %lf seconds from %d process", subdomain_transfer_info[1], time4-time1, time5-time4, time2-time1, time3-time1, time6-time1, time7-time6, time4-time7, owner_process_id);
			delete [] remote_subdomain_grid_data_buf;
		}
		remote_subdomain_CoR_grids.push_back(subdomain_grid);
	}
	
	if (subdomain_grid != NULL && subdomain_grid->get_grid_size() > 0)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_grid->get_local_cell_global_indexes() != NULL, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_grid");
	return subdomain_grid;
}


void Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains_IDs(Remap_grid_class *original_cor_grid, int subdomain_ID, int level, std::vector<int> &halo_subdomains_IDs)
{
	int subdomain_x_index, subdomain_y_index, subdomain_halo_right_bound, subdomain_halo_left_bound, subdomain_halo_top_bound, subdomain_halo_bottom_bound;
	int x_loop_begin, x_loop_end, y_loop_begin, y_loop_end, current_subdomain_id, last_subdomain_id = -CCPL_NULL_INT, current_i_index, current_j_index;

	halo_subdomains_IDs.clear();

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, level > 0, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains_IDs");

	if (is_north_pole_subdomain(subdomain_ID)) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, level <= num_nonPole_region_y_subdomains+num_southPole_subdomain, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains");
		if (level <= num_nonPole_region_y_subdomains) {
			for (int i = 0; i < num_nonPole_region_x_subdomains; i ++)
				halo_subdomains_IDs.push_back((num_nonPole_region_y_subdomains-level)*num_nonPole_region_x_subdomains+i+num_northPole_subdomain+num_southPole_subdomain);
		}
		else halo_subdomains_IDs.push_back(num_northPole_subdomain);
		return;
	}

	if (is_south_pole_subdomain(subdomain_ID)) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, level <= num_nonPole_region_y_subdomains+num_northPole_subdomain, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains");
		if (level <= num_nonPole_region_y_subdomains) {
			for (int i = 0; i < num_nonPole_region_x_subdomains; i ++)
				halo_subdomains_IDs.push_back((level-1)*num_nonPole_region_x_subdomains+i+num_northPole_subdomain+num_southPole_subdomain);
		}
		else halo_subdomains_IDs.push_back(0);
		return;
	}
	
	get_subdomain_XY_index(subdomain_ID, subdomain_x_index, subdomain_y_index);
	subdomain_halo_top_bound = subdomain_y_index - level;
	subdomain_halo_bottom_bound = subdomain_y_index + level;
	subdomain_halo_left_bound = subdomain_x_index - level;
	subdomain_halo_right_bound = subdomain_x_index + level;
	if (whole_domain_min_lon == 0.0 && whole_domain_max_lon	== 360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES)) {
		if (2*level < num_nonPole_region_x_subdomains) {
			subdomain_halo_left_bound = (subdomain_halo_left_bound+num_nonPole_region_x_subdomains)%num_nonPole_region_x_subdomains;
			subdomain_halo_right_bound = subdomain_halo_right_bound % num_nonPole_region_x_subdomains;
		}
		else if (2*level == num_nonPole_region_x_subdomains) {
			subdomain_halo_left_bound = (subdomain_halo_left_bound+num_nonPole_region_x_subdomains)%num_nonPole_region_x_subdomains;
			subdomain_halo_right_bound = -CCPL_NULL_INT;				
		}
		else {
			subdomain_halo_left_bound = -CCPL_NULL_INT;
			subdomain_halo_right_bound = -CCPL_NULL_INT;
		}
	}
	x_loop_begin = subdomain_halo_left_bound;
	x_loop_end = subdomain_halo_right_bound;
	if (x_loop_begin == -CCPL_NULL_INT)
		x_loop_begin = 0;
	if (x_loop_end == -CCPL_NULL_INT)
		x_loop_end = num_nonPole_region_x_subdomains - 1;
	if (x_loop_begin > x_loop_end)
		x_loop_end += num_nonPole_region_x_subdomains;
	y_loop_begin = subdomain_halo_top_bound+1 >= 0? subdomain_halo_top_bound+1 : 0;
	y_loop_end = subdomain_halo_bottom_bound < num_nonPole_region_y_subdomains? subdomain_halo_bottom_bound : num_nonPole_region_y_subdomains;	

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "subdomain_halo_left_bound  is %d %d, top bound is %d %d", subdomain_halo_left_bound, subdomain_halo_right_bound, subdomain_halo_top_bound, subdomain_halo_bottom_bound);
	for (int j = 0; j < 2; j ++) {
		current_j_index = (j == 0)? subdomain_halo_top_bound : subdomain_halo_bottom_bound;
		if (current_j_index != -CCPL_NULL_INT)
			for (int i = x_loop_begin; i <= x_loop_end; i ++) {
				int temp_i = i;				
				if (whole_domain_min_lon == 0.0 && whole_domain_max_lon == 360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
					temp_i = i%num_nonPole_region_x_subdomains;
				current_subdomain_id = get_subdomain_ID_from_XY_index(temp_i, current_j_index);
				if (last_subdomain_id != current_subdomain_id && current_subdomain_id != -CCPL_NULL_INT)
					halo_subdomains_IDs.push_back(current_subdomain_id);
				last_subdomain_id = current_subdomain_id;
			}
	}
	for (int i = 0; i < 2; i ++) {
		current_i_index = (i == 0)? subdomain_halo_left_bound : subdomain_halo_right_bound;
		if (current_i_index != -CCPL_NULL_INT)
			for (int j = y_loop_begin; j < y_loop_end; j ++) {
				current_subdomain_id = get_subdomain_ID_from_XY_index(current_i_index, j);
				if (last_subdomain_id != current_subdomain_id && current_subdomain_id != -CCPL_NULL_INT) {
					halo_subdomains_IDs.push_back(current_subdomain_id);
				}
				last_subdomain_id = current_subdomain_id;
			}
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, halo_subdomains_IDs.size() > 0, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains: %d %d %d: %d %d %d %d: %d %d %d %d", subdomain_ID, num_total_subdomains, level, subdomain_halo_top_bound, subdomain_halo_bottom_bound, subdomain_halo_right_bound, subdomain_halo_left_bound, num_northPole_subdomain, num_southPole_subdomain, num_nonPole_region_x_subdomains, num_nonPole_region_y_subdomains);
	for (int i = 0; i < halo_subdomains_IDs.size(); i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, halo_subdomains_IDs[i] >= 0 && halo_subdomains_IDs[i] < num_total_subdomains, "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains");
		for (int j = i+1; j < halo_subdomains_IDs.size(); j ++)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, halo_subdomains_IDs[i] != halo_subdomains_IDs[j], "Software error in Remapping_grid_domain_decomp_engine::get_subdomain_halo_subdomains: %d %d", j, halo_subdomains_IDs[i]);
	}
}


int Remapping_grid_domain_decomp_engine::get_subdomain_ID_from_XY_index(int x_index, int y_index)
{
	if (y_index < -num_southPole_subdomain || y_index >= num_nonPole_region_y_subdomains+num_northPole_subdomain || x_index < 0 || x_index >= num_nonPole_region_x_subdomains)
		return -CCPL_NULL_INT;
	if (y_index == num_nonPole_region_y_subdomains)
		return 0;
	if (y_index < 0)
		return num_northPole_subdomain;
	
	return num_northPole_subdomain + num_southPole_subdomain + y_index*num_nonPole_region_x_subdomains + x_index;
}


void Remapping_grid_domain_decomp_engine::get_subdomain_XY_index(int subdomain_ID, int &x_index, int &y_index)
{
	if (is_north_pole_subdomain(subdomain_ID)) {
		y_index = num_nonPole_region_y_subdomains;
		x_index = 0;
	}
	else if (is_south_pole_subdomain(subdomain_ID)) {
		y_index = -1;
		x_index = 0;
	}
	else {
		x_index = (subdomain_ID-num_northPole_subdomain-num_southPole_subdomain)%num_nonPole_region_x_subdomains;
		y_index = (subdomain_ID-num_northPole_subdomain-num_southPole_subdomain)/num_nonPole_region_x_subdomains;
	}
}


int Remapping_grid_domain_decomp_engine::calculate_two_subdomains_distance(Remap_grid_class *original_cor_grid, int subdomain1_ID, int subdomain2_ID)
{
	int subdomain1_x_index, subdomain1_y_index, subdomain2_x_index, subdomain2_y_index;
	int x_dist, y_dist;

	
	get_subdomain_XY_index(subdomain1_ID, subdomain1_x_index, subdomain1_y_index);
	get_subdomain_XY_index(subdomain2_ID, subdomain2_x_index, subdomain2_y_index);
	if (subdomain1_y_index == -1 || subdomain1_y_index == num_nonPole_region_y_subdomains || subdomain2_y_index == -1 || subdomain2_y_index == num_nonPole_region_y_subdomains)
		return abs(subdomain1_y_index-subdomain2_y_index);
	
	y_dist = abs(subdomain1_y_index-subdomain2_y_index);
	x_dist = abs(subdomain1_x_index-subdomain2_x_index);
	if (whole_domain_min_lon == 0.0 && whole_domain_max_lon	== 360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
		x_dist = MIN(x_dist, num_nonPole_region_x_subdomains-x_dist);
	return MAX(y_dist, x_dist);
}


Remap_grid_class *Remapping_grid_domain_decomp_engine::expand_src_subdomain_halo_grid(int subdomain_index)
{
	expand_subdomain_halo_grid(src_subdomains_CoR_grids[subdomain_index].first, src_original_cor_grid, src_subdomains_CoR_grids, true);
	return src_subdomains_expanded_CoR_grids[subdomain_index];
}


void Remapping_grid_domain_decomp_engine::expand_subdomain_halo_grid(int subdomain_ID, Remap_grid_class *original_cor_grid, std::vector<std::pair<int, Remap_grid_class*> >& subdomains_CoR_grids, bool is_src)
{
	int subdomain_index = subdomain_ID / num_total_procs;
	int halo_level = (++subdomains_halo_num_levels[subdomain_index]);
	int k, j;
	std::vector<int> halo_subdomains_IDs;
	double time1, time2, time3;
	
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_CoR_grids[subdomain_index].first == subdomain_ID, "software error in Remapping_grid_domain_decomp_engine::expand_subdomain_halo_grid");
	get_subdomain_halo_subdomains_IDs(original_cor_grid, subdomain_ID, halo_level, halo_subdomains_IDs);

	EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "expand %d subdomain of grid \%s\"", subdomain_ID, original_cor_grid->get_grid_name());
	if(report_error_enabled) {
		for (j = 0; j < halo_subdomains_IDs.size(); j ++)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, calculate_two_subdomains_distance(original_cor_grid, halo_subdomains_IDs[j], subdomain_ID) == halo_level, "software error in Remapping_grid_domain_decomp_engine::expand_subdomain_halo_grid");
		for (k = 0; k < num_total_subdomains; k ++) {
			if (calculate_two_subdomains_distance(original_cor_grid, k, subdomain_ID) == halo_level) {
				for (j = 0; j < halo_subdomains_IDs.size(); j ++)
					if (k == halo_subdomains_IDs[j])
						break;
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, j < halo_subdomains_IDs.size(), "software error in Remapping_grid_domain_decomp_engine::expand_subdomain_halo_grid");
			}
		}
		for (j = 0; j < halo_subdomains_IDs.size(); j ++)
			for (k = 0; k < subdomains_halo_subdomains_IDs[subdomain_index].size(); k ++)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, halo_subdomains_IDs[j] != subdomains_halo_subdomains_IDs[subdomain_index][k], "software error in Remapping_grid_domain_decomp_engine::expand_subdomain_halo_grid");
	}

	for (j = 0; j < halo_subdomains_IDs.size(); j ++)
		subdomains_halo_subdomains_IDs[subdomain_index].push_back(halo_subdomains_IDs[j]);

	wtime(&time1);
	for (j = 0; j < halo_subdomains_IDs.size(); j ++) {
		Remap_grid_class *halo_subdomain_CoR_grid = get_subdomain_grid(halo_subdomains_IDs[j], subdomains_CoR_grids, is_src);
		if (halo_subdomain_CoR_grid != NULL)
			subdomains_halo_subdomains_CoR_grids[subdomain_index].push_back(halo_subdomain_CoR_grid);
	}
	wtime(&time2);
	calculate_subdomain_halo_bounds(original_cor_grid, subdomain_ID, halo_level, subdomains_halo_lon_bounds[subdomain_ID/num_total_procs], subdomains_halo_lat_bounds[subdomain_ID/num_total_procs]);
	generate_subdomain_expanded_CoR_grid_with_triangulation(subdomain_ID, subdomains_CoR_grids, original_cor_grid, is_src);
	wtime(&time3);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "time to expand is %lf vs %lf", time2-time1, time3-time2);
}


void Remapping_grid_domain_decomp_engine::get_src_current_expanded_subdomain_boundaries(int subdomain_index, double &min_lon, double &max_lon, double &min_lat, double &max_lat)
{
	std::pair<double, double> subdomains_halo_lon_bound, subdomains_halo_lat_bound;

	calculate_subdomain_halo_bounds(src_original_cor_grid, src_subdomains_CoR_grids[subdomain_index].first, subdomains_halo_num_levels[subdomain_index], subdomains_halo_lon_bound, subdomains_halo_lat_bound);
	min_lon = subdomains_halo_lon_bound.first;
	max_lon = subdomains_halo_lon_bound.second;
	min_lat = subdomains_halo_lat_bound.first;
	max_lat = subdomains_halo_lat_bound.second;
}


void Remapping_grid_domain_decomp_engine::calculate_subdomain_halo_bounds(Remap_grid_class *original_cor_grid, int subdomain_ID, int level, std::pair<double, double>& subdomains_halo_lon_bound, std::pair<double, double>& subdomains_halo_lat_bound)
{
	double subdomain_halo_min_lon, subdomain_halo_max_lon, subdomain_halo_min_lat, subdomain_halo_max_lat;
	int subdomain_x_index, subdomain_y_index;


	if (is_north_pole_subdomain(subdomain_ID)) {
		subdomain_halo_min_lon = 0.0;
		subdomain_halo_max_lon = 360.0;
		subdomain_halo_max_lat = 90.0;
		subdomain_halo_min_lat = MAX(whole_domain_min_lat, nonPole_region_max_lat-level*nonPole_subdomain_y_side_length);
	}
	else if (is_south_pole_subdomain(subdomain_ID)) {
		subdomain_halo_min_lon = 0.0;
		subdomain_halo_max_lon = 360.0;
		subdomain_halo_min_lat = -90.0;
		subdomain_halo_max_lat = MIN(whole_domain_max_lat, nonPole_region_min_lat+level*nonPole_subdomain_y_side_length);
	}
	else {
		get_subdomain_XY_index(subdomain_ID, subdomain_x_index, subdomain_y_index);
		subdomain_halo_min_lat = MAX(whole_domain_min_lat, nonPole_region_min_lat+(subdomain_y_index-level)*nonPole_subdomain_y_side_length);
		subdomain_halo_max_lat = MIN(whole_domain_max_lat, nonPole_region_min_lat+(subdomain_y_index+1+level)*nonPole_subdomain_y_side_length);
		if (subdomain_halo_min_lat <= -90.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
			subdomain_halo_min_lat = -89.5;
		if (subdomain_halo_max_lat >= 90.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
			subdomain_halo_max_lat = 89.5;
		if (whole_domain_min_lon == 0.0 && whole_domain_max_lon == 360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES)) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, level < num_nonPole_region_x_subdomains/2, "Software error in Remapping_grid_domain_decomp_engine::calculate_subdomain_halo_bounds: %d %d", level, num_nonPole_region_x_subdomains);
			subdomain_halo_min_lon = whole_domain_min_lon + (subdomain_x_index-level)*nonPole_subdomain_x_side_length;
			subdomain_halo_max_lon = whole_domain_min_lon + (subdomain_x_index+1+level)*nonPole_subdomain_x_side_length;
			if (subdomain_halo_min_lon < 0.0)
				subdomain_halo_min_lon += 360.0;
			if (subdomain_halo_max_lon > 360.0)
				subdomain_halo_max_lon -= 360.0;
		}
		else {
			double temp_whole_domain_max_lon = whole_domain_max_lon;
			if (whole_domain_min_lon > whole_domain_max_lon)
				temp_whole_domain_max_lon += 360.0;
			subdomain_halo_min_lon = MAX(whole_domain_min_lon, whole_domain_min_lon+(subdomain_x_index-level)*nonPole_subdomain_x_side_length);
			subdomain_halo_max_lon = MIN(temp_whole_domain_max_lon, whole_domain_min_lon+(subdomain_x_index+1+level)*nonPole_subdomain_x_side_length);
			if (subdomain_halo_max_lon > 360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
				subdomain_halo_max_lon -= 360.0;
			if (subdomain_halo_min_lon > 360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
				subdomain_halo_min_lon -= 360.0;
		}
	}

	subdomains_halo_lon_bound.first = subdomain_halo_min_lon;
	subdomains_halo_lon_bound.second = subdomain_halo_max_lon;
	subdomains_halo_lat_bound.first = subdomain_halo_min_lat;
	subdomains_halo_lat_bound.second = subdomain_halo_max_lat;
}


bool Remapping_grid_domain_decomp_engine::is_point_in_a_domain(double point_lon, double point_lat, double domain_bound_min_lon, double domain_bound_max_lon, double domain_bound_min_lat, double domain_bound_max_lat)
{
	double temp_point_lon = point_lon, temp_domain_bound_max_lon = domain_bound_max_lon;

	if (domain_bound_min_lon > domain_bound_max_lon) {
		temp_domain_bound_max_lon += 360.0;
		if (point_lon <= domain_bound_max_lon)
			temp_point_lon += 360.0;
	}

	return temp_point_lon >= domain_bound_min_lon && temp_point_lon <= temp_domain_bound_max_lon && point_lat >= domain_bound_min_lat && point_lat <= domain_bound_max_lat;
}


void *Remapping_grid_domain_decomp_engine::generate_expanded_grid_field_data(int active_size, std::vector<Remap_grid_class*> &subdomain_CoR_grids, std::vector<const void*> &subdomain_grid_data_bufs, bool *expanded_domain_mask, int data_type_size, std::vector<int> &num_values_per_cell)
{
	void *expanded_grid_field_data;
	int all_cells_iter = 0, active_cells_iter = 0, i, max_num_values_per_cell = -1;


	for (i = 0; i < subdomain_CoR_grids.size(); i ++)
		if (subdomain_grid_data_bufs[i] != NULL)
			break;
	if (i == subdomain_CoR_grids.size() || num_values_per_cell.size() == 0)
		return NULL;

	for (i = 0; i < num_values_per_cell.size(); i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_values_per_cell[i] > 0, "Software error in Remapping_grid_domain_decomp_engine::generate_expanded_grid_field_data");
		if (max_num_values_per_cell < num_values_per_cell[i])
			max_num_values_per_cell = num_values_per_cell[i];
	}
	if (max_num_values_per_cell > 1)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, data_type_size == 8, "Software error in Remapping_grid_domain_decomp_engine::generate_expanded_grid_field_data");

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, active_size > 0 && data_type_size > 0, "Software error in Remapping_grid_domain_decomp_engine::generate_expanded_grid_field_data %d vs %d vs %d", active_size, data_type_size);
	expanded_grid_field_data = new char [active_size*max_num_values_per_cell*data_type_size];
	for (int i = 0; i < subdomain_grid_data_bufs.size(); i ++) {
		if (subdomain_grid_data_bufs[i] == NULL) 
			continue;
		char *src_data_buf = (char*) subdomain_grid_data_bufs[i];
		char *dst_data_buf = (char*) expanded_grid_field_data;		
		for (int j = 0; j < subdomain_CoR_grids[i]->get_grid_size(); j ++) {
			if (expanded_domain_mask[all_cells_iter]) {
				for (int t = 0; t < num_values_per_cell[i]*data_type_size; t ++)
					dst_data_buf[active_cells_iter*max_num_values_per_cell*data_type_size+t] = src_data_buf[j*num_values_per_cell[i]*data_type_size+t];
				for (int t = num_values_per_cell[i]; t < max_num_values_per_cell; t ++)
					((double*)dst_data_buf)[active_cells_iter*max_num_values_per_cell+t] = NULL_COORD_VALUE;
				active_cells_iter ++;
			}
			all_cells_iter ++;
		}
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, active_size == active_cells_iter, "Software error in Remapping_grid_domain_decomp_engine::generate_expanded_grid_field_data");

	return expanded_grid_field_data;
}


unsigned long Remapping_grid_domain_decomp_engine::cal_triangles_checksum_on_boundary(Triangle_inline* result_triangles, int num_result_triangles, PatCC_Point head, PatCC_Point tail, int subdomain_ID, double threshold)
{

//	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "begin cal_triangles_checksum_on_boundary, triangle is %d", num_result_triangles);
	if (float_eq(head.x, tail.x) && float_eq(head.y, tail.y))
		return 0;

	PatCC_Point sibling_head, sibling_tail;
	bool match_sibling = head.x > tail.x;

	if (match_sibling) {
		sibling_head.x = head.x - 360.0;
		sibling_tail.x = tail.x;
		sibling_head.y = head.y;
		sibling_tail.y = tail.y;

		tail.x += 360.0;
	}

	unsigned long checksum = 0;
	unsigned count = 0;

	for(unsigned i = 0; i < num_result_triangles; i ++) {
		if (result_triangles[i].is_cyclic) {
			double tmp_x[3];

			for(int t = 0; t < 3; t ++)
				tmp_x[t] = result_triangles[i].v[t].x >= 180 ? result_triangles[i].v[t].x - 360 : result_triangles[i].v[t].x;

			Triangle_inline tmp_triangle1(PatCC_Point(tmp_x[0], result_triangles[i].v[0].y, result_triangles[i].v[0].id, false), PatCC_Point(tmp_x[1], result_triangles[i].v[1].y, result_triangles[i].v[1].id, false), PatCC_Point(tmp_x[2], result_triangles[i].v[2].y, result_triangles[i].v[2].id, false), true);

			for(int t = 0; t < 3; t ++)
				tmp_x[t] = result_triangles[i].v[t].x < 180 ? result_triangles[i].v[t].x + 360 : result_triangles[i].v[t].x;

			Triangle_inline tmp_triangle2(PatCC_Point(tmp_x[0], result_triangles[i].v[0].y, result_triangles[i].v[0].id, false), PatCC_Point(tmp_x[1], result_triangles[i].v[1].y, result_triangles[i].v[1].id, false), PatCC_Point(tmp_x[2], result_triangles[i].v[2].y, result_triangles[i].v[2].id, false), true);

			if (is_triangle_intersecting_with_segment(&tmp_triangle1, head, tail, threshold) || 
				is_triangle_intersecting_with_segment(&tmp_triangle2, head, tail, threshold) || (match_sibling && (is_triangle_intersecting_with_segment(&tmp_triangle1, sibling_head, sibling_tail, threshold) || is_triangle_intersecting_with_segment(&tmp_triangle2, sibling_head, sibling_tail, threshold)))){
//				EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "the triangle is %d(%lf, %lf) %d(%lf, %lf) %d(%lf, %lf)", result_triangles[i].v[0].id, result_triangles[i].v[0].x, result_triangles[i].v[0].y, result_triangles[i].v[1].id, result_triangles[i].v[1].x, result_triangles[i].v[1].y, result_triangles[i].v[2].id, result_triangles[i].v[2].x, result_triangles[i].v[2].y);	
				checksum += hash_triangle_by_id(result_triangles[i]);
					count++;
			}
		} else {
			if (is_triangle_intersecting_with_segment(&result_triangles[i], head, tail, threshold) || (match_sibling && is_triangle_intersecting_with_segment(&result_triangles[i], sibling_head, sibling_tail, threshold))) {
//				EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "the triangle is %d(%lf, %lf) %d(%lf, %lf) %d(%lf, %lf)", result_triangles[i].v[0].id, result_triangles[i].v[0].x, result_triangles[i].v[0].y, result_triangles[i].v[1].id, result_triangles[i].v[1].x, result_triangles[i].v[1].y, result_triangles[i].v[2].id, result_triangles[i].v[2].x, result_triangles[i].v[2].y);	
				checksum += hash_triangle_by_id(result_triangles[i]);
				count++;
			}
		}
	}
	checksum *= count;
	checksum &= 0x0FFFFFFFFFFFFFFF;

//	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "end cal_triangles_checksum_on_boundary");
	return checksum;
}


void Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation(int subdomain_ID, std::vector<std::pair<int, Remap_grid_class*> >& origin_subdomains_CoR_grids, Remap_grid_class *original_cor_grid, bool is_src)
{
	int total_size = 0, active_size, k;
	std::vector<Remap_grid_class *> *subdomain_halo_subdomains_CoR_grids = &(subdomains_halo_subdomains_CoR_grids[subdomain_ID/num_total_procs]);
	std::vector<Remap_grid_class *> expanded_subdomains_CoR_grids;
	Remap_grid_class *subdomain_CoR_grid = origin_subdomains_CoR_grids[subdomain_ID/num_total_procs].second, *sub_grids[2], *expanded_CoR_grid;
	std::vector<const void *> subdomains_center_lon_fields, subdomains_center_lat_fields, subdomains_vertex_lon_fields, subdomains_vertex_lat_fields, subdomains_mask_fields, subdomains_global_cell_indexes;
	bool *expanded_domain_mask;
	int subdomain_type; // 0 for common subdomain, 1 for north pole, -1 for south pole.
	double expanded_bound_min_lon, expanded_bound_max_lon, expanded_bound_min_lat, expanded_bound_max_lat;
	char expanded_H2D_grid_name[NAME_STR_SIZE], expanded_lon_grid_name[NAME_STR_SIZE], expanded_lat_grid_name[NAME_STR_SIZE];
	Triangle_inline *result_triangles;
	int subdomain_triangulation_info_index;
	std::pair<double, double> lon_bounds, lat_bounds;
	std::vector<int> num_vertexes_in_grid, num_values_per_cell;
	int grid_original_num_vertex = is_src? src_grid_original_num_vertex : dst_grid_original_num_vertex;
	int max_num_vertexes_in_halos = -1;
	double time1, time2, time3, time4, time5, time6, time7;
	double *new_vertex_lon_values, *new_vertex_lat_values;
	int new_max_num_voronoi_diagram_vertex;

	wtime(&time1);
	if (src_subdomains_expanded_CoR_grids[subdomain_ID/num_total_procs] != NULL) {
		delete src_subdomains_expanded_CoR_grids[subdomain_ID/num_total_procs];
		src_subdomains_expanded_CoR_grids[subdomain_ID/num_total_procs] = NULL;
	}

	expanded_bound_min_lon = subdomains_halo_lon_bounds[subdomain_ID/num_total_procs].first;
	expanded_bound_max_lon = subdomains_halo_lon_bounds[subdomain_ID/num_total_procs].second;
	expanded_bound_min_lat = subdomains_halo_lat_bounds[subdomain_ID/num_total_procs].first;
	expanded_bound_max_lat = subdomains_halo_lat_bounds[subdomain_ID/num_total_procs].second;
	if (subdomain_CoR_grid != NULL && subdomain_CoR_grid->get_grid_size() > 0) {
		expanded_subdomains_CoR_grids.push_back(subdomain_CoR_grid);
		num_vertexes_in_grid.push_back(subdomain_CoR_grid->get_num_vertexes());
		max_num_vertexes_in_halos = max_num_vertexes_in_halos > subdomain_CoR_grid->get_num_vertexes()? max_num_vertexes_in_halos : subdomain_CoR_grid->get_num_vertexes();
		num_values_per_cell.push_back(1);
		subdomains_center_lon_fields.push_back(subdomain_CoR_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
		subdomains_center_lat_fields.push_back(subdomain_CoR_grid->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);
		if (grid_original_num_vertex > 0 && subdomain_CoR_grid->get_grid_vertex_field(COORD_LABEL_LON) != NULL)
			subdomains_vertex_lon_fields.push_back(subdomain_CoR_grid->get_grid_vertex_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
		else subdomains_vertex_lon_fields.push_back(NULL);
		if (grid_original_num_vertex > 0 && subdomain_CoR_grid->get_grid_vertex_field(COORD_LABEL_LAT) != NULL)
			subdomains_vertex_lat_fields.push_back(subdomain_CoR_grid->get_grid_vertex_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);
		else subdomains_vertex_lat_fields.push_back(NULL);
		subdomains_mask_fields.push_back(subdomain_CoR_grid->get_grid_mask_field()->get_grid_data_field()->data_buf);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_CoR_grid->get_grid_mask_field() != NULL, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
		subdomains_global_cell_indexes.push_back(subdomain_CoR_grid->get_local_cell_global_indexes());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_CoR_grid->get_local_cell_global_indexes() != NULL, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
		total_size = subdomain_CoR_grid->get_grid_size();
	}
	wtime(&time2);
	for (int i = 0; i < (*subdomain_halo_subdomains_CoR_grids).size(); i ++) {
		if ((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_size() == 0) {			
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The bypass NO. %d halo sub grid out of %d for the subdomain %d of grid \"%s\"", i, (*subdomain_halo_subdomains_CoR_grids).size(), subdomain_ID, (*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_name());
			continue;
		}
		expanded_subdomains_CoR_grids.push_back((*subdomain_halo_subdomains_CoR_grids)[i]);
		num_vertexes_in_grid.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_num_vertexes());
		max_num_vertexes_in_halos = max_num_vertexes_in_halos > (*subdomain_halo_subdomains_CoR_grids)[i]->get_num_vertexes()? max_num_vertexes_in_halos : (*subdomain_halo_subdomains_CoR_grids)[i]->get_num_vertexes();
		num_values_per_cell.push_back(1);
		total_size += (*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_size();
		subdomains_center_lon_fields.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
		subdomains_center_lat_fields.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);
		if (grid_original_num_vertex > 0 && (*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_vertex_field(COORD_LABEL_LON) != NULL)
			subdomains_vertex_lon_fields.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_vertex_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
		else subdomains_vertex_lon_fields.push_back(NULL);
		if (grid_original_num_vertex > 0 && (*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_vertex_field(COORD_LABEL_LAT) != NULL)
			subdomains_vertex_lat_fields.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_vertex_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);
		else subdomains_vertex_lat_fields.push_back(NULL);
		subdomains_mask_fields.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_mask_field()->get_grid_data_field()->data_buf);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_mask_field() != NULL, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
		if ((*subdomain_halo_subdomains_CoR_grids)[i]->get_grid_size() > 0)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (*subdomain_halo_subdomains_CoR_grids)[i]->get_local_cell_global_indexes() != NULL, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
		subdomains_global_cell_indexes.push_back((*subdomain_halo_subdomains_CoR_grids)[i]->get_local_cell_global_indexes());
	}
	wtime(&time3);

	if (total_size == 0)
		return;

	expanded_domain_mask = new bool [total_size];
	active_size = 0;
	k = 0;
	for (int i = 0; i < subdomains_center_lon_fields.size(); i ++) {
		double *lon_data_buf = (double*) subdomains_center_lon_fields[i];
		double *lat_data_buf = (double*) subdomains_center_lat_fields[i];
		for (int j = 0; j < expanded_subdomains_CoR_grids[i]->get_grid_size(); j ++) {
			expanded_domain_mask[k] = is_point_in_a_domain(lon_data_buf[j], lat_data_buf[j], expanded_bound_min_lon, expanded_bound_max_lon, expanded_bound_min_lat, expanded_bound_max_lat);
			if (subdomain_CoR_grid != NULL && i == 0)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, expanded_domain_mask[k], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
			if (expanded_domain_mask[k])
				active_size ++;
			k ++;
		}
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, k == total_size, "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");

	if (active_size == 0) {
		delete [] expanded_domain_mask;
		return;
	}
	wtime(&time4);
	double *expanded_CoR_grid_center_lon_data = (double*) generate_expanded_grid_field_data(active_size, expanded_subdomains_CoR_grids, subdomains_center_lon_fields, expanded_domain_mask, sizeof(double), num_values_per_cell);
	double *expanded_CoR_grid_center_lat_data = (double*) generate_expanded_grid_field_data(active_size, expanded_subdomains_CoR_grids, subdomains_center_lat_fields, expanded_domain_mask, sizeof(double), num_values_per_cell);
	double *expanded_CoR_grid_vertex_lon_data = (double*) generate_expanded_grid_field_data(active_size, expanded_subdomains_CoR_grids, subdomains_vertex_lon_fields, expanded_domain_mask, sizeof(double), num_vertexes_in_grid);
	double *expanded_CoR_grid_vertex_lat_data = (double*) generate_expanded_grid_field_data(active_size, expanded_subdomains_CoR_grids, subdomains_vertex_lat_fields, expanded_domain_mask, sizeof(double), num_vertexes_in_grid);
	bool *expanded_CoR_grid_mask_data = (bool*) generate_expanded_grid_field_data(active_size, expanded_subdomains_CoR_grids, subdomains_mask_fields, expanded_domain_mask, sizeof(bool), num_values_per_cell);
	int *expanded_CoR_grid_global_cell_indexes = (int*) generate_expanded_grid_field_data(active_size, expanded_subdomains_CoR_grids, subdomains_global_cell_indexes, expanded_domain_mask, sizeof(int), num_values_per_cell);
	sprintf(expanded_H2D_grid_name, "expanded_subdomain_%d_of_%s", subdomain_ID, original_cor_grid->get_grid_name());
	sprintf(expanded_lon_grid_name, "lon_of_%s", expanded_H2D_grid_name);
	sprintf(expanded_lat_grid_name, "lat_of_%s", expanded_H2D_grid_name);
	wtime(&time5);
	if (expanded_bound_min_lon==0.0 && expanded_bound_max_lon==360.0 && words_are_the_same(original_cor_grid->get_coord_unit(), COORD_UNIT_DEGREES))
		sub_grids[0] = new Remap_grid_class(expanded_lon_grid_name, "lon", original_cor_grid->get_coord_unit(), "cyclic", 0);
	else sub_grids[0] = new Remap_grid_class(expanded_lon_grid_name, "lon", original_cor_grid->get_coord_unit(), "acyclic", 0);
	sub_grids[1] = new Remap_grid_class(expanded_lat_grid_name, "lat", original_cor_grid->get_coord_unit(), "acyclic", 0);
	expanded_CoR_grid = new Remap_grid_class(expanded_H2D_grid_name, 2, sub_grids, active_size);
	expanded_CoR_grid->set_generated_from_duplication();
	expanded_CoR_grid->set_local_cell_global_indexes(expanded_CoR_grid_global_cell_indexes, false);
	expanded_CoR_grid->set_grid_boundary(expanded_bound_min_lon, expanded_bound_max_lon, expanded_bound_min_lat, expanded_bound_max_lat);
	expanded_CoR_grid->read_grid_data_from_array("center", COORD_LABEL_LON, DATA_TYPE_DOUBLE, (const char*)expanded_CoR_grid_center_lon_data, 0);
	expanded_CoR_grid->read_grid_data_from_array("center", COORD_LABEL_LAT, DATA_TYPE_DOUBLE, (const char*)expanded_CoR_grid_center_lat_data, 0);
	expanded_CoR_grid->read_grid_data_from_array("mask", "mask", DATA_TYPE_BOOL, (const char*)expanded_CoR_grid_mask_data, 0);
	int count = 0;
	for (int i = 0; i < expanded_CoR_grid->get_grid_size(); i ++)
		if (expanded_CoR_grid_mask_data[i])
			count ++;

	if (expanded_CoR_grid_vertex_lon_data != NULL) {
		if (is_src)
			expanded_CoR_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, DATA_TYPE_DOUBLE, (const char*)expanded_CoR_grid_vertex_lon_data, max_num_vertexes_in_halos);
		else expanded_CoR_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, DATA_TYPE_DOUBLE, (const char*)expanded_CoR_grid_vertex_lon_data, max_num_vertexes_in_halos);
	}
	if (expanded_CoR_grid_vertex_lat_data != NULL) {
		if (is_src)
			expanded_CoR_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, DATA_TYPE_DOUBLE, (const char*)expanded_CoR_grid_vertex_lat_data, max_num_vertexes_in_halos);
		else
			expanded_CoR_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, DATA_TYPE_DOUBLE, (const char*)expanded_CoR_grid_vertex_lat_data, max_num_vertexes_in_halos);
		if (report_error_enabled)
			expanded_CoR_grid->check_center_vertex_values_consistency_2D(true);
	}

	if (report_error_enabled && original_cor_grid->get_grid_center_field(COORD_LABEL_LON) != NULL && original_cor_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->required_data_size > 0) {
		double *lon_values_global = (double*) original_cor_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf;
		double *lat_values_global = (double*) original_cor_grid->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf;
		bool *masks_global = NULL;
		if (original_cor_grid->get_grid_mask_field() != NULL)
			masks_global = (bool*) original_cor_grid->get_grid_mask_field()->get_grid_data_field()->data_buf;
		int dim1 = original_cor_grid->get_a_leaf_grid(COORD_LABEL_LON)->get_grid_size();
		if (original_cor_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->required_data_size == original_cor_grid->get_grid_size()) {
			for (int i = 0; i < expanded_CoR_grid->get_grid_size(); i ++)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, expanded_CoR_grid_center_lon_data[i] == lon_values_global[expanded_CoR_grid_global_cell_indexes[i]] && expanded_CoR_grid_center_lat_data[i] == lat_values_global[expanded_CoR_grid_global_cell_indexes[i]], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation: (%lf, %lf) vs (%lf, %lf)", expanded_CoR_grid_center_lon_data[i], lon_values_global[expanded_CoR_grid_global_cell_indexes[i]], expanded_CoR_grid_center_lat_data[i], lat_values_global[expanded_CoR_grid_global_cell_indexes[i]]);
			if (masks_global != NULL) {
				for (int i = 0; i < expanded_CoR_grid->get_grid_size(); i ++)
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, expanded_CoR_grid_mask_data[i] == masks_global[expanded_CoR_grid_global_cell_indexes[i]], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
			}
			else {
				for (int i = 0; i < expanded_CoR_grid->get_grid_size(); i ++)
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, expanded_CoR_grid_mask_data[i], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation");
			}
		}
		else {
			for (int i = 0; i < expanded_CoR_grid->get_grid_size(); i ++)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, expanded_CoR_grid_center_lon_data[i] == lon_values_global[expanded_CoR_grid_global_cell_indexes[i]%dim1] && expanded_CoR_grid_center_lat_data[i] == lat_values_global[expanded_CoR_grid_global_cell_indexes[i]/dim1], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation: (%lf, %lf) vs (%lf, %lf)", expanded_CoR_grid_center_lon_data[i], lon_values_global[expanded_CoR_grid_global_cell_indexes[i]%dim1], expanded_CoR_grid_center_lat_data[i], lat_values_global[expanded_CoR_grid_global_cell_indexes[i]/dim1]);
			for (int i = 0; i < expanded_CoR_grid->get_grid_size(); i ++)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, expanded_CoR_grid_mask_data[i] == masks_global[expanded_CoR_grid_global_cell_indexes[i]], "Software error in Remapping_grid_domain_decomp_engine::generate_subdomain_expanded_CoR_grid_with_triangulation: (%lf, %lf) vs (%lf, %lf)", expanded_CoR_grid_center_lon_data[i], lon_values_global[expanded_CoR_grid_global_cell_indexes[i]%dim1], expanded_CoR_grid_center_lat_data[i], lat_values_global[expanded_CoR_grid_global_cell_indexes[i]/dim1]);
		}
	}

	delete [] expanded_CoR_grid_center_lon_data;
	delete [] expanded_CoR_grid_center_lat_data;
	delete [] expanded_CoR_grid_mask_data;
	if (expanded_CoR_grid_vertex_lon_data != NULL)
		delete [] expanded_CoR_grid_vertex_lon_data;
	if (expanded_CoR_grid_vertex_lat_data != NULL)
		delete [] expanded_CoR_grid_vertex_lat_data;
	delete [] expanded_domain_mask;

	wtime(&time6);

	if (expanded_CoR_grid->get_num_vertexes() == 0 && subdomain_CoR_grid != NULL && subdomain_CoR_grid->get_grid_size() >= 0) {  // enabled at any time for more testing
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Trianglulation will be conducted for src expanded grid %s %ld: %0.18f", original_cor_grid->get_grid_name(), expanded_CoR_grid->get_grid_size(), expanded_bound_max_lat);

		if (is_north_pole_subdomain(subdomain_ID))
			subdomain_type = 1;
		else if (is_south_pole_subdomain(subdomain_ID))
			subdomain_type = -1;
		else 
			subdomain_type = 0;

		calculate_subdomain_halo_bounds(original_cor_grid, subdomain_ID, 0, lon_bounds, lat_bounds);
		PatCC_Delaunay_Voronoi *triangulation = new PatCC_Delaunay_Voronoi(expanded_CoR_grid->get_grid_size(), expanded_bound_min_lon, expanded_bound_max_lon, expanded_bound_min_lat, expanded_bound_max_lat, 
																		   (double*)expanded_CoR_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf, (double*)expanded_CoR_grid->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf, 
																		   (bool*) expanded_CoR_grid->get_grid_mask_field()->get_grid_data_field()->data_buf, words_are_the_same(expanded_CoR_grid->get_sphere_grid_coord_unit(), COORD_UNIT_DEGREES), expanded_CoR_grid_global_cell_indexes,
				   															original_cor_grid->get_boundary_min_lon(), original_cor_grid->get_boundary_max_lon(), original_cor_grid->get_boundary_min_lat(), original_cor_grid->get_boundary_max_lat(), lon_bounds.first, lon_bounds.second, 
                                                                            lat_bounds.first, lat_bounds.second, subdomain_CoR_grid == NULL? 0 : subdomain_CoR_grid->get_grid_size(), subdomain_type);
		current_triangulation = triangulation;
		triangulation->generate_all_result_triangles();

		result_triangles = new Triangle_inline [triangulation->num_result_triangles];
		subdomain_triangulation_info_index = subdomain_ID / comp_comm_group_mgt_mgr->search_global_node(comp_id)->get_num_procs();
		subdomain_triangulation_info[subdomain_triangulation_info_index*3+0] = subdomain_ID;
		subdomain_triangulation_info[subdomain_triangulation_info_index*3+1] = triangulation->num_result_triangles;
		subdomain_triangulation_info[subdomain_triangulation_info_index*3+2] = triangulation_comm_data_buf_content_size;

		for(int triangle_index = 0; triangle_index < triangulation->num_result_triangles; triangle_index ++) {
			result_triangles[triangle_index] = (triangulation->result_triangles)[triangle_index];
			for (int point_index = 0; point_index < 3; point_index ++) {
				int tmp_int = result_triangles[triangle_index].v[point_index].id;
				write_data_into_array_buffer(&tmp_int, sizeof(int), &triangulation_comm_data_buf, triangulation_comm_data_buf_max_size, triangulation_comm_data_buf_content_size);
				write_data_into_array_buffer(&(result_triangles[triangle_index].v[point_index].x), sizeof(double), &triangulation_comm_data_buf, triangulation_comm_data_buf_max_size, triangulation_comm_data_buf_content_size);
				write_data_into_array_buffer(&(result_triangles[triangle_index].v[point_index].y), sizeof(double), &triangulation_comm_data_buf, triangulation_comm_data_buf_max_size, triangulation_comm_data_buf_content_size);
			}
		}

		subdomains_final_triangles[subdomain_triangulation_info_index] = result_triangles;
		num_subdomains_final_triangles[subdomain_triangulation_info_index] = triangulation->num_result_triangles;

		if (subdomain_CoR_grid != NULL) 
			triangulation->generate_Voronoi_diagram(new_max_num_voronoi_diagram_vertex, &new_vertex_lon_values, &new_vertex_lat_values);
		delete triangulation;
		if (subdomain_CoR_grid != NULL) {
			subdomain_CoR_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, DATA_TYPE_DOUBLE, (char*)new_vertex_lon_values, new_max_num_voronoi_diagram_vertex);
			subdomain_CoR_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, DATA_TYPE_DOUBLE, (char*)new_vertex_lat_values, new_max_num_voronoi_diagram_vertex);
			subdomain_CoR_grid->set_vertex_values_generated_in_default();
			delete [] new_vertex_lat_values;
			delete [] new_vertex_lon_values;
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish generating vertexes of the grid \"%s\" at the subdomain %d", subdomain_CoR_grid->get_grid_name(), subdomain_ID);
		}

		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, local_subdomains_index[subdomain_triangulation_info_index] == subdomain_ID, "Software error after triangulation, %d %d.", local_subdomains_index[subdomain_triangulation_info_index], subdomain_ID);
	}


	wtime(&time7);

	EXECUTION_REPORT(REPORT_LOG, -1, true, "generate_subdomain_expanded_CoR_grid_with_triangulation times: %lf  %lf  %lf  %lf  %lf  %lf  %lf", time2-time1, time3-time2, time4-time3, time5-time4, time6-time5, time7-time6);
	src_subdomains_expanded_CoR_grids[subdomain_ID/num_total_procs] = expanded_CoR_grid;
	if (expanded_CoR_grid_global_cell_indexes != NULL)
		delete [] expanded_CoR_grid_global_cell_indexes;
}


Remap_grid_class *Remapping_grid_domain_decomp_engine::get_src_expanded_subdomain(int subdomain_index)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_index >= 0 && subdomain_index <= src_subdomains_expanded_CoR_grids.size(), "Software error in Remapping_grid_domain_decomp_engine::get_src_expanded_subdomain");
	return src_subdomains_expanded_CoR_grids[subdomain_index];
}


Remap_grid_class *Remapping_grid_domain_decomp_engine::get_dst_subdomain(int subdomain_index)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_index >= 0 && subdomain_index <= dst_subdomains_CoR_grids.size(), "Software error in Remapping_grid_domain_decomp_engine::get_src_expanded_subdomain");
	return dst_subdomains_CoR_grids[subdomain_index].second;
}


Distributed_H2D_grid_engine *Distributed_H2D_grid_mgt::search_distributed_H2D_grid(const char *comp_full_name, Original_grid_info *original_grid)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid->get_H2D_sub_CoR_grid() != NULL, "Software error in Distributed_H2D_grid_mgt::search_distributed_H2D_grid");
	return search_distributed_H2D_grid(comp_full_name, original_grid->get_H2D_sub_CoR_grid());
}


Distributed_H2D_grid_engine *Distributed_H2D_grid_mgt::search_distributed_H2D_grid(const char *comp_full_name, Remap_grid_class *CoR_grid_initialization)
{
	for (int i = 0; i < distributed_H2D_grids.size(); i ++)
		if ((comp_full_name == NULL || words_are_the_same(comp_full_name, distributed_H2D_grids[i]->get_comp_full_name())) && CoR_grid_initialization == distributed_H2D_grids[i]->get_CoR_grid_initialization())
			return distributed_H2D_grids[i];

	return NULL;
}



Distributed_H2D_grid_mgt::~Distributed_H2D_grid_mgt()
{
	for (int i = 0; i < distributed_H2D_grids.size(); i ++)
		delete distributed_H2D_grids[i];
}


Distributed_H2D_grid_engine *Distributed_H2D_grid_mgt::generate_distributed_H2D_grid_engine(const char *comp_full_name, const char *grid_name, Remap_grid_class *CoR_grid_initialization, long global_grid_size, int num_lons, int num_lats, Original_grid_info *original_grid, bool generate_basic_decomp_grid)
{
	if (search_distributed_H2D_grid(comp_full_name, CoR_grid_initialization) != NULL)
		return search_distributed_H2D_grid(comp_full_name, CoR_grid_initialization);

	EXECUTION_REPORT_LOG(REPORT_LOG, comp_comm_group_mgt_mgr->search_global_node(comp_full_name)->get_comp_id(), true, "Generate basic distributed grid of original grid \"%s\"  \"%s\"", comp_full_name, CoR_grid_initialization->get_grid_name());
	Distributed_H2D_grid_engine *distributed_H2D_grid = new Distributed_H2D_grid_engine(comp_full_name, grid_name, CoR_grid_initialization, global_grid_size, num_lons, num_lats, original_grid, generate_basic_decomp_grid);
	distributed_H2D_grids.push_back(distributed_H2D_grid);
	return distributed_H2D_grid;
}


void Distributed_H2D_grid_mgt::calculate_common_grid_domain_for_remapping(int comp_id, Remap_grid_class *src_CoR_grid, Remap_grid_class *dst_CoR_grid, double &common_min_lon, double &common_max_lon, double &common_min_lat, double &common_max_lat, bool extrapolation_enabled)
{
	double src_min_lon, src_max_lon, src_min_lat, src_max_lat, src_lon_diff, src_lat_diff, src_point_dens;
	double dst_min_lon, dst_max_lon, dst_min_lat, dst_max_lat, dst_lon_diff, num_expanded_degrees, expansion_factor=4.0;


	src_min_lon = src_CoR_grid->get_boundary_min_lon();
	src_max_lon = src_CoR_grid->get_boundary_max_lon();
	src_min_lat = src_CoR_grid->get_boundary_min_lat();
	src_max_lat = src_CoR_grid->get_boundary_max_lat();
	dst_min_lon = dst_CoR_grid->get_boundary_min_lon();
	dst_max_lon = dst_CoR_grid->get_boundary_max_lon();
	dst_min_lat = dst_CoR_grid->get_boundary_min_lat();
	dst_max_lat = dst_CoR_grid->get_boundary_max_lat();

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_max_lon != NULL_COORD_VALUE && dst_min_lon != NULL_COORD_VALUE && src_min_lat < src_max_lat && dst_min_lat < dst_max_lat, "Software error in Distributed_H2D_grid_mgt::common_grid_domain_for_remapping");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(src_CoR_grid->get_coord_unit(), dst_CoR_grid->get_coord_unit()), "Software error in Distributed_H2D_grid_mgt::common_grid_domain_for_remapping");

	src_lat_diff = src_max_lat - src_min_lat;
	if (src_min_lon < src_max_lon)
		src_lon_diff = src_max_lon - src_min_lon;
	else {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(dst_CoR_grid->get_coord_unit(), COORD_UNIT_DEGREES), "Software error in Distributed_H2D_grid_mgt::calculate_common_grid_domain_for_remapping");
		src_lon_diff = src_max_lon - src_min_lon + 360.0;
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_lon_diff > 0 && src_lat_diff > 0, "Software error in Distributed_H2D_grid_mgt::calculate_common_grid_domain_for_remapping");
	src_point_dens = sqrt(((double)src_CoR_grid->get_grid_size())/(src_lon_diff*src_lat_diff));
	num_expanded_degrees = expansion_factor / src_point_dens;
	num_expanded_degrees = 0.0;
	if (words_are_the_same(src_CoR_grid->get_coord_unit(), COORD_UNIT_DEGREES)) {
		dst_min_lat = MAX(-90.0, dst_min_lat-num_expanded_degrees);
		dst_max_lat = MIN(90.0, dst_max_lat+num_expanded_degrees);
	}
	if (dst_min_lon < dst_max_lon)
		dst_lon_diff = dst_max_lon - dst_min_lon;
	else {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(dst_CoR_grid->get_coord_unit(), COORD_UNIT_DEGREES), "Software error in Distributed_H2D_grid_mgt::calculate_common_grid_domain_for_remapping");
		dst_lon_diff = dst_max_lon - dst_min_lon + 360.0;	
	}
	if (words_are_the_same(src_CoR_grid->get_coord_unit(), COORD_UNIT_DEGREES))
		if (dst_lon_diff + 2*num_expanded_degrees >= 360.0) {
			dst_min_lon = 0.0;
			dst_max_lon = 360.0;
		}
		else {
			dst_min_lon = dst_min_lon - num_expanded_degrees;
			dst_max_lon = dst_max_lon + num_expanded_degrees;
			if (dst_min_lon < 0.0) 
				dst_min_lon += 360.0;
			if (dst_max_lon > 360.0)
				dst_max_lon -= 360.0;
		}

	common_min_lon = NULL_COORD_VALUE; common_max_lon = NULL_COORD_VALUE; 
	common_min_lat = NULL_COORD_VALUE; common_max_lat = NULL_COORD_VALUE;

	common_min_lat = MIN(src_min_lat, dst_min_lat);
	common_max_lat = MAX(src_max_lat, dst_max_lat);

	if (src_min_lon < src_max_lon) {
		if (dst_min_lon < dst_max_lon) {
			common_min_lon = MIN(src_min_lon, dst_min_lon);
			common_max_lon = MAX(src_max_lon, dst_max_lon);				
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, common_min_lon < common_max_lon, "C-Coupler encounters a special case of no common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
		}
		else {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(dst_min_lon >= src_max_lon && src_min_lon >= dst_max_lon), "C-Coupler encounters a special case of no common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
			if (dst_min_lon >= src_max_lon) {
				common_min_lon = dst_min_lon;
				common_max_lon = MAX(dst_max_lon, src_max_lon);
			}
			else if (src_min_lon >= dst_max_lon) {
				common_min_lon = MIN(src_min_lon, dst_min_lon);
				common_max_lon = dst_max_lon;
			}
			else if (dst_max_lon >= src_max_lon || dst_min_lon <= src_min_lon) {
				common_min_lon = dst_min_lon;
				common_max_lon = dst_max_lon;
			}
			else {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_min_lon == 0.0 && src_max_lon == 360.0 && words_are_the_same(src_CoR_grid->get_coord_unit(), COORD_UNIT_DEGREES), "C-Coupler encounters a special case of two separate segments in the common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
				common_min_lon = 0.0;
				common_max_lon = 360.0;
			}
		}
	}
	else {
		if (dst_min_lon < dst_max_lon) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(dst_min_lon >= src_max_lon && src_min_lon >= dst_max_lon), "C-Coupler encounters a special case of no common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
			if (src_min_lon >= dst_max_lon) {
				common_min_lon = src_min_lon;
				common_max_lon = MAX(src_max_lon, dst_max_lon);
			}
			else if (src_max_lon <= dst_min_lon) {
				common_min_lon = MIN(dst_min_lon, src_min_lon);
				common_max_lon = src_max_lon;
			}
			else if (src_max_lon >= dst_max_lon || src_min_lon <= dst_min_lon) {
				common_min_lon = src_min_lon;
				common_max_lon = src_max_lon;
			}
			else {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_min_lon == 0.0 && dst_max_lon == 360.0 && words_are_the_same(src_CoR_grid->get_coord_unit(), COORD_UNIT_DEGREES), "C-Coupler encounters a special case of two separate segments in the common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
				common_min_lon = 0.0;
				common_max_lon = 360.0;
			}
		}
		else {
			common_min_lon = MIN(src_min_lon, dst_min_lon);
			common_max_lon = MAX(src_max_lon, dst_max_lon);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, common_max_lon < common_min_lon, "C-Coupler encounters a special case of two separate segments in the common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
			common_min_lon = 0.0;
		}
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, common_min_lat < common_max_lat, "C-Coupler encounters a special case of no common domain between source and target grids in interpolation. Please verify or ask help from Dr. Li Liu");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, common_min_lon != common_max_lon, "Software error in Distributed_H2D_grid_mgt::common_grid_domain_for_remapping");
}


void Distributed_H2D_grid_mgt::generate_vertexes_via_distributed_triangulation(Distributed_H2D_grid_engine *distributed_H2D_grid)
{
	double grid_min_lon = distributed_H2D_grid->get_CoR_grid_initialization()->get_boundary_min_lon();
    double grid_max_lon = distributed_H2D_grid->get_CoR_grid_initialization()->get_boundary_max_lon();
	double grid_min_lat = distributed_H2D_grid->get_CoR_grid_initialization()->get_boundary_min_lat();
    double grid_max_lat = distributed_H2D_grid->get_CoR_grid_initialization()->get_boundary_max_lat();
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(distributed_H2D_grid->get_original_grid()->get_comp_id());
	int  max_num_local_vertex = 0, max_num_global_vertex, current_num_vertex;
	int i, j, k, iter = 0;
	Remap_grid_class *current_subdomain_grid, *subdomains_decomp_grid;


	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_min_lon != NULL_COORD_VALUE && grid_max_lon != NULL_COORD_VALUE && grid_min_lat != NULL_COORD_VALUE && grid_max_lat != NULL_COORD_VALUE, "Software error in Distributed_H2D_grid_mgt::generate_vertexes_via_distributed_triangulation");

	if (distributed_H2D_grid->get_CoR_grid_initialization()->get_num_vertexes() > 0 || distributed_H2D_grid->get_CoR_grid_initialization()->are_all_vertex_fields_specified_by_user())
		return;

	EXECUTION_REPORT_LOG(REPORT_LOG, distributed_H2D_grid->get_original_grid()->get_comp_id(), true, "Start to generate vertexes for H2D grid \"%s\" via distributed triangulation", distributed_H2D_grid->get_original_grid()->get_grid_name());
	Remapping_grid_domain_decomp_engine *remapping_grid_domain_decomp_engine = new Remapping_grid_domain_decomp_engine(distributed_H2D_grid->get_original_grid()->get_grid_id(), distributed_H2D_grid->get_original_grid()->get_comp_id(), distributed_H2D_grid->get_CoR_grid_initialization(), NULL, grid_min_lon, grid_max_lon, grid_min_lat, grid_max_lat);
	subdomains_decomp_grid = remapping_grid_domain_decomp_engine->get_src_subdomain_decomp_grid()->get_decomp_grid();
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_decomp_grid->get_num_vertexes() == 0, "Software error in Distributed_H2D_grid_mgt::generate_vertexes_via_distributed_triangulation");
	for (int i = 0; i < remapping_grid_domain_decomp_engine->src_subdomains_CoR_grids.size(); i ++)
		if (remapping_grid_domain_decomp_engine->src_subdomains_CoR_grids[i].second != NULL) {
			current_num_vertex = remapping_grid_domain_decomp_engine->src_subdomains_CoR_grids[i].second->get_num_vertexes();
			if (current_num_vertex > max_num_local_vertex)
				max_num_local_vertex = current_num_vertex;
		}
	MPI_Allreduce(&max_num_local_vertex, &max_num_global_vertex, 1, MPI_INT, MPI_MAX, comp_node->get_comm_group());
	double *new_vertex_lon_buf = new double [(subdomains_decomp_grid->get_grid_size()>0? subdomains_decomp_grid->get_grid_size(): 1) * max_num_global_vertex];
	double *new_vertex_lat_buf = new double [(subdomains_decomp_grid->get_grid_size()>0? subdomains_decomp_grid->get_grid_size(): 1) * max_num_global_vertex];
	const int *subdomains_cell_global_indexes = subdomains_decomp_grid->get_local_cell_global_indexes();
	for (int i = 0; i < remapping_grid_domain_decomp_engine->src_subdomains_CoR_grids.size(); i ++) {
		current_subdomain_grid = remapping_grid_domain_decomp_engine->src_subdomains_CoR_grids[i].second;
		if (current_subdomain_grid == NULL)
			continue;
		if (current_subdomain_grid->get_grid_size() > 0) {
			double *original_vertex_lon_buf = (double*) (current_subdomain_grid->get_grid_vertex_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf);
			double *original_vertex_lat_buf = (double*) (current_subdomain_grid->get_grid_vertex_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf);
			const int *current_cell_global_indexes = current_subdomain_grid->get_local_cell_global_indexes();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_cell_global_indexes != NULL, "Software error in Distributed_H2D_grid_mgt::generate_vertexes_via_distributed_triangulation");
			for (j = 0; j < current_subdomain_grid->get_grid_size(); j ++) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomains_cell_global_indexes[iter] == current_cell_global_indexes[j], "Software error in Distributed_H2D_grid_mgt::generate_vertexes_via_distributed_triangulation");
				for (k = 0; k < current_subdomain_grid->get_num_vertexes(); k ++) {
					new_vertex_lon_buf[iter*max_num_global_vertex+k] = original_vertex_lon_buf[j*current_subdomain_grid->get_num_vertexes()+k];
					new_vertex_lat_buf[iter*max_num_global_vertex+k] = original_vertex_lat_buf[j*current_subdomain_grid->get_num_vertexes()+k];
				}
				for (; k < max_num_global_vertex; k ++) {
					new_vertex_lon_buf[iter*max_num_global_vertex+k] = NULL_COORD_VALUE;
					new_vertex_lat_buf[iter*max_num_global_vertex+k] = NULL_COORD_VALUE;
				}
				iter ++;
			}
		}	
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, iter == subdomains_decomp_grid->get_grid_size(), "Software error in Distributed_H2D_grid_mgt::generate_vertexes_via_distributed_triangulation");
	subdomains_decomp_grid->read_grid_data_from_array("vertex", COORD_LABEL_LON, DATA_TYPE_DOUBLE, (char*)new_vertex_lon_buf, max_num_global_vertex);
	subdomains_decomp_grid->read_grid_data_from_array("vertex", COORD_LABEL_LAT, DATA_TYPE_DOUBLE, (char*)new_vertex_lat_buf, max_num_global_vertex);	
	subdomains_decomp_grid->set_vertex_values_generated_in_default();
	delete [] new_vertex_lat_buf;
	delete [] new_vertex_lon_buf;
	if (report_error_enabled) {
		is_coord_unit_degree[0] = words_are_the_same(distributed_H2D_grid->get_CoR_grid_initialization()->get_coord_unit(), COORD_UNIT_DEGREES);
		is_coord_unit_degree[1] = is_coord_unit_degree[0];
		subdomains_decomp_grid->check_center_vertex_values_consistency_2D(true);
	}

	distributed_H2D_grid->get_CoR_grid_initialization()->read_grid_data_from_array("vertex", COORD_LABEL_LON, DATA_TYPE_DOUBLE, NULL, max_num_global_vertex);
	distributed_H2D_grid->get_CoR_grid_initialization()->read_grid_data_from_array("vertex", COORD_LABEL_LAT, DATA_TYPE_DOUBLE, NULL, max_num_global_vertex);
	distributed_H2D_grid->get_CoR_grid_initialization()->set_vertex_values_generated_in_default();
	distributed_H2D_grid->update_basic_decomp_grid(remapping_grid_domain_decomp_engine->get_src_subdomain_decomp_grid());
	delete remapping_grid_domain_decomp_engine;

	EXECUTION_REPORT_LOG(REPORT_LOG, distributed_H2D_grid->get_original_grid()->get_comp_id(), true, "Finish generate vertexes for H2D grid \"%s\" via distributed triangulation", distributed_H2D_grid->get_original_grid()->get_grid_name());
}


void Distributed_H2D_grid_mgt::generate_full_grid_data_from_distributed_grid(int comp_id, Remap_grid_class *global_grid_without_data, Remap_grid_class *distributed_basic_grid, bool do_bcast)
{
	double *full_center_lon_buf = NULL, *full_center_lat_buf = NULL, *full_vertex_lon_buf = NULL, *full_vertex_lat_buf = NULL, *full_area_buf = NULL;
	bool *full_mask_buf = NULL;
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	Remap_data_field *grid_data_field;
	long global_size;
	int already_global = 0;

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "generate full grid data for %s (%lx) with %lx", global_grid_without_data->get_grid_name(), global_grid_without_data, distributed_basic_grid);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !do_bcast, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, distributed_basic_grid->get_grid_size() > 0, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
	if (report_error_enabled) {
		distributed_basic_grid->end_grid_definition_stage(NULL);
		gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size(), sizeof(double), NULL, (void **)(&full_center_lon_buf), global_size, comp_node->get_comm_group());
	}

	if (comp_node->get_current_proc_local_id() == 0)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, global_size == global_grid_without_data->get_grid_size(), "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");

	if (do_bcast || comp_node->get_current_proc_local_id() == 0) {
		grid_data_field = global_grid_without_data->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field();
		if (grid_data_field->data_buf != NULL) {
			if (grid_data_field->required_data_size == global_grid_without_data->get_grid_size())
				for (int i = 0; i < grid_data_field->required_data_size; i ++)
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, ((double*)grid_data_field->data_buf)[i] == full_center_lon_buf[i], "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid: %s (%lx) %d: %lf vs %lf", global_grid_without_data->get_grid_name(), global_grid_without_data, i, ((double*)grid_data_field->data_buf)[i], full_center_lon_buf[i]);
			already_global = 1;
		}
	}
	MPI_Bcast(&already_global, 1, MPI_INT, 0, comp_node->get_comm_group());
	if (full_center_lon_buf != NULL) {
		delete [] full_center_lon_buf;
		full_center_lon_buf = NULL;
	}

	if (already_global == 1)
		return;

	gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size(), sizeof(double), NULL, (void **)(&full_center_lon_buf), global_size, comp_node->get_comm_group());
	gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size(), sizeof(double), NULL, (void **)(&full_center_lat_buf), global_size, comp_node->get_comm_group());
	if (distributed_basic_grid->get_grid_vertex_field(COORD_LABEL_LON) != NULL) {
		gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_vertex_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size()*global_grid_without_data->get_num_vertexes(), sizeof(double), NULL, (void **)(&full_vertex_lon_buf), global_size, comp_node->get_comm_group());
		gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_vertex_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size()*global_grid_without_data->get_num_vertexes(), sizeof(double), NULL, (void **)(&full_vertex_lat_buf), global_size, comp_node->get_comm_group());
		if (comp_node->get_current_proc_local_id() == 0)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, global_size == global_grid_without_data->get_grid_size() * global_grid_without_data->get_num_vertexes(), "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
	}
	if (distributed_basic_grid->get_grid_imported_area() != NULL)
		gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_imported_area()->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size(), sizeof(double), NULL, (void **)(&full_area_buf), global_size, comp_node->get_comm_group());

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, distributed_basic_grid->get_grid_mask_field() != NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
	gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), distributed_basic_grid->get_grid_mask_field()->get_grid_data_field()->data_buf, distributed_basic_grid->get_grid_size(), sizeof(bool), NULL, (void **)(&full_mask_buf), global_size, comp_node->get_comm_group());

	if (do_bcast || comp_node->get_current_proc_local_id() == 0) {
		grid_data_field = global_grid_without_data->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, full_center_lon_buf != NULL && grid_data_field->data_buf == NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid: %s  %lx   %d", global_grid_without_data->get_grid_name(), global_grid_without_data, grid_data_field->required_data_size);
		grid_data_field->data_buf = full_center_lon_buf;
		grid_data_field->read_data_size = grid_data_field->required_data_size = global_grid_without_data->get_grid_size();
		grid_data_field = global_grid_without_data->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, full_center_lat_buf != NULL && grid_data_field->data_buf == NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
		grid_data_field->data_buf = full_center_lat_buf;
		grid_data_field->read_data_size = grid_data_field->required_data_size = global_grid_without_data->get_grid_size();
		if (full_vertex_lat_buf != NULL) {
			grid_data_field = global_grid_without_data->get_grid_vertex_field(COORD_LABEL_LON)->get_grid_data_field();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field->data_buf == NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
			grid_data_field->data_buf = full_vertex_lon_buf;
			grid_data_field->read_data_size = grid_data_field->required_data_size = global_grid_without_data->get_grid_size() * global_grid_without_data->get_num_vertexes();
			grid_data_field = global_grid_without_data->get_grid_vertex_field(COORD_LABEL_LAT)->get_grid_data_field();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field->data_buf == NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
			grid_data_field->data_buf = full_vertex_lat_buf;
			grid_data_field->read_data_size = grid_data_field->required_data_size = global_grid_without_data->get_grid_size() * global_grid_without_data->get_num_vertexes();
		}
		if (full_area_buf != NULL) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, global_grid_without_data->get_grid_imported_area() != NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
			grid_data_field = global_grid_without_data->get_grid_imported_area()->get_grid_data_field();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_data_field->data_buf == NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
			grid_data_field->data_buf = full_area_buf;
			grid_data_field->required_data_size = global_grid_without_data->get_grid_size();
		}
			
		grid_data_field = global_grid_without_data->get_grid_mask_field()->get_grid_data_field();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, full_mask_buf != NULL && grid_data_field->data_buf == NULL, "Software error in Distributed_H2D_grid_mgt::generate_full_data_grid_from_distributed_grid");
		grid_data_field->data_buf = full_mask_buf;
		grid_data_field->read_data_size = grid_data_field->required_data_size = global_grid_without_data->get_grid_size();
		if (report_error_enabled)
			global_grid_without_data->end_grid_definition_stage(NULL);
	}
}

