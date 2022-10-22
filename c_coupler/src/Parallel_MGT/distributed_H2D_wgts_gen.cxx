/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "distributed_H2D_wgts_gen.h"
#include "remap_utils_nearest_points.h"


void sort_normal_distributed_remapping_weight_elements_locally(common_sort_struct<Normal_distributed_wgt_element> *normal_distributed_wgt_map, int num_weights)
{
	for (int i = 0; i < num_weights; i ++)
		normal_distributed_wgt_map[i].key = normal_distributed_wgt_map[i].content.dst_cell_index;
	do_quick_sort(normal_distributed_wgt_map, (int*) NULL, 0, num_weights-1);
	int weight_group_begin_pos = 0;
	for (int i = 0; i < num_weights; i ++) {
		normal_distributed_wgt_map[i].key = normal_distributed_wgt_map[i].content.src_cell_index;
		if (i+1 == num_weights || normal_distributed_wgt_map[i].content.dst_cell_index != normal_distributed_wgt_map[i+1].content.dst_cell_index) {
			for (int j = weight_group_begin_pos; j < i; j ++)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, normal_distributed_wgt_map[j].content.dst_cell_index == normal_distributed_wgt_map[i].content.dst_cell_index);
			do_quick_sort(normal_distributed_wgt_map, (int*) NULL, weight_group_begin_pos, i);
			weight_group_begin_pos = i+1;
		}
	}
}


void sort_normal_remapping_weights_in_sparse_matrix_locally(Remap_weight_sparse_matrix *weight_matrix)
{
	if (weight_matrix->get_num_weights() == 0)
		return;

	common_sort_struct<Normal_distributed_wgt_element> *normal_distributed_wgt_map = new common_sort_struct<Normal_distributed_wgt_element> [weight_matrix->get_num_weights()];
	for (int i = 0; i < weight_matrix->get_num_weights(); i ++) {
		normal_distributed_wgt_map[i].content.src_cell_index = weight_matrix->get_indexes_src_grid()[i];
		normal_distributed_wgt_map[i].content.dst_cell_index = weight_matrix->get_indexes_dst_grid()[i];
		normal_distributed_wgt_map[i].content.wgt_value = weight_matrix->get_weight_values()[i];
	}

	sort_normal_distributed_remapping_weight_elements_locally(normal_distributed_wgt_map,weight_matrix->get_num_weights());

	for (int i = 0; i < weight_matrix->get_num_weights(); i ++) {
		weight_matrix->get_indexes_src_grid()[i] = normal_distributed_wgt_map[i].content.src_cell_index;
		weight_matrix->get_indexes_dst_grid()[i] = normal_distributed_wgt_map[i].content.dst_cell_index;
		weight_matrix->get_weight_values()[i] = normal_distributed_wgt_map[i].content.wgt_value;
	}

	delete [] normal_distributed_wgt_map;
}


Remap_weight_sparse_matrix *rearrange_for_local_H2D_parallel_weights(Decomp_info *dst_decomp_info, Remap_operator_basis *entire_remap_operator, common_sort_struct<Normal_distributed_wgt_element> *normal_distributed_wgt_map, int num_wgts_after_redistribution, const char *full_default_wgt_file_name, Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid)
{
	Comp_comm_group_mgt_node *dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_decomp_info->get_comp_id());	
	int num_cells_per_proc = (dst_decomp_info->get_num_global_cells()+dst_comp_node->get_num_procs()-1)/dst_comp_node->get_num_procs();
	int num_cells_after_redistribution;
	common_sort_struct<Grid_cell_rearrange_map_element> *dst_grid_cell_rearrange_map = NULL;
	std::vector<Normal_distributed_wgt_element> normal_distributed_wgt_elements_vect;
	long i, j, k;
 	Remap_weight_sparse_matrix *normal_remap_weights = NULL;

	comp_comm_group_mgt_mgr->get_root_component_model()->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "rearrange_for_local_H2D_parallel_weights");

	for (int i = 0; i < num_wgts_after_redistribution; i ++) {
		normal_distributed_wgt_map[i].target_proc_id = normal_distributed_wgt_map[i].content.dst_cell_index / num_cells_per_proc;
		normal_distributed_wgt_map[i].content.owner_process_id = normal_distributed_wgt_map[i].target_proc_id;
	}
	Distribute_merge_sort<Normal_distributed_wgt_element> *wgt_rearrange_distribute_sorting = new Distribute_merge_sort<Normal_distributed_wgt_element>(dst_comp_node->get_current_proc_local_id(), dst_comp_node->get_current_proc_local_id(), dst_comp_node, dst_comp_node);

	wgt_rearrange_distribute_sorting->do_data_sorting_with_target_process_id(&normal_distributed_wgt_map, dst_decomp_info, dst_comp_node, &num_wgts_after_redistribution);
	delete wgt_rearrange_distribute_sorting;

	if (full_default_wgt_file_name != NULL) {
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Generate the default remapping_weight file is \"%s\"", full_default_wgt_file_name);
		sort_normal_distributed_remapping_weight_elements_locally(normal_distributed_wgt_map, num_wgts_after_redistribution);
		normal_remap_weights = new Remap_weight_sparse_matrix(entire_remap_operator);
		for (int i = 0; i < num_wgts_after_redistribution; i ++)
			normal_remap_weights->add_weights(&(normal_distributed_wgt_map[i].content.src_cell_index), normal_distributed_wgt_map[i].content.dst_cell_index, &(normal_distributed_wgt_map[i].content.wgt_value), 1, false);
		Remap_weight_sparse_matrix *overall_remap_weight_sparse_matrix = normal_remap_weights->gather(dst_comp_node->get_comp_id());
		delete normal_remap_weights;
		if (dst_comp_node->get_current_proc_local_id() == 0) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid != NULL && dst_original_grid != NULL, "Software error in rearrange_for_local_H2D_parallel_weights");
			IO_netcdf *io_netcdf = new IO_netcdf(full_default_wgt_file_name, full_default_wgt_file_name, "w", true);
			int last_execution_phase_number = execution_phase_number;
			io_netcdf->write_remap_weights(entire_remap_operator, overall_remap_weight_sparse_matrix);
			long checksum = src_original_grid->get_checksum_H2D_mask();
			io_netcdf->put_global_attr("mask_checksum_a", &checksum, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
			checksum = src_original_grid->get_checksum_H2D_center_lon();
			io_netcdf->put_global_attr("xc_checksum_a", &checksum, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
			checksum = src_original_grid->get_checksum_H2D_center_lat();
			io_netcdf->put_global_attr("yc_checksum_a", &checksum, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
			checksum = dst_original_grid->get_checksum_H2D_mask();
			io_netcdf->put_global_attr("mask_checksum_b", &checksum, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
			checksum = dst_original_grid->get_checksum_H2D_center_lon();
			io_netcdf->put_global_attr("xc_checksum_b", &checksum, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
			checksum = dst_original_grid->get_checksum_H2D_center_lat();
			io_netcdf->put_global_attr("yc_checksum_b", &checksum, DATA_TYPE_LONG, DATA_TYPE_LONG, 1);
			int num_procs = comp_comm_group_mgt_mgr->search_global_node(dst_original_grid->get_comp_id())->get_num_procs();
			io_netcdf->put_global_attr("number of processes for generation", &num_procs, DATA_TYPE_INT, DATA_TYPE_INT, 1);
			delete io_netcdf;
		}
		if (overall_remap_weight_sparse_matrix != NULL)
			delete overall_remap_weight_sparse_matrix;
	}

	for (int i = 0; i < num_wgts_after_redistribution; i ++) {
		normal_distributed_wgt_map[i].key = normal_distributed_wgt_map[i].content.dst_cell_index;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, normal_distributed_wgt_map[i].content.owner_process_id == dst_comp_node->get_current_proc_local_id(), "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
	}
	do_quick_sort(normal_distributed_wgt_map, (int*) NULL, 0, num_wgts_after_redistribution-1);

	num_cells_after_redistribution = 0;
	if (dst_decomp_info->get_num_local_cells() > 0) {
		dst_grid_cell_rearrange_map = new common_sort_struct<Grid_cell_rearrange_map_element> [dst_decomp_info->get_num_local_cells()];
		for (int i = 0; i < dst_decomp_info->get_num_local_cells(); i ++) {
			if (dst_decomp_info->get_local_cell_global_indx()[i] == CCPL_NULL_INT)
				continue;
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_decomp_info->get_local_cell_global_indx()[i] >= 0 && dst_decomp_info->get_local_cell_global_indx()[i] < dst_decomp_info->get_num_global_cells(), "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
			new (&dst_grid_cell_rearrange_map[num_cells_after_redistribution].content) Grid_cell_rearrange_map_element(dst_decomp_info->get_local_cell_global_indx()[i]/num_cells_per_proc, dst_comp_node->get_current_proc_local_id(), dst_decomp_info->get_local_cell_global_indx()[i], i);
			dst_grid_cell_rearrange_map[num_cells_after_redistribution].target_proc_id = dst_grid_cell_rearrange_map[num_cells_after_redistribution].content.domain_decomp_process_id;
			num_cells_after_redistribution ++;
		}
	}
	Distribute_merge_sort<Grid_cell_rearrange_map_element> *cell_rearrange_distribute_sorting = new Distribute_merge_sort<Grid_cell_rearrange_map_element>(dst_comp_node->get_current_proc_local_id(), dst_comp_node->get_current_proc_local_id(), dst_comp_node, dst_comp_node);
	cell_rearrange_distribute_sorting->do_data_sorting_with_target_process_id(&dst_grid_cell_rearrange_map, dst_decomp_info, dst_comp_node, &num_cells_after_redistribution);
	for (int i = 0; i < num_cells_after_redistribution; i ++)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_grid_cell_rearrange_map[i].content.domain_decomp_process_id == dst_comp_node->get_current_proc_local_id(), "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
	for (int i = 0; i < num_cells_after_redistribution; i ++)
		dst_grid_cell_rearrange_map[i].key = dst_grid_cell_rearrange_map[i].content.local_cell_global_index;
	do_quick_sort(dst_grid_cell_rearrange_map, (int*) NULL, 0, num_cells_after_redistribution-1);

	if (num_wgts_after_redistribution > 0) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_cells_after_redistribution > 0, "Distributed_H2D_weights_generator::read_normal_remap_weights");
		j = 0;
		for (i = 0, j = 0; i < num_cells_after_redistribution; i ++) {
			if (dst_grid_cell_rearrange_map[i].content.local_cell_global_index < normal_distributed_wgt_map[j].content.dst_cell_index)
				continue;
			for (; j < num_wgts_after_redistribution && normal_distributed_wgt_map[j].content.dst_cell_index < dst_grid_cell_rearrange_map[i].content.local_cell_global_index; j ++);
			if (j >= num_wgts_after_redistribution)
				break;
			for (k = j; k < num_wgts_after_redistribution && dst_grid_cell_rearrange_map[i].content.local_cell_global_index == normal_distributed_wgt_map[k].content.dst_cell_index; k ++) {
				normal_distributed_wgt_map[k].content.owner_process_id = dst_grid_cell_rearrange_map[i].content.owner_process_id;
				normal_distributed_wgt_map[k].content.dst_local_index = dst_grid_cell_rearrange_map[i].content.owner_local_index;
				normal_distributed_wgt_elements_vect.push_back(normal_distributed_wgt_map[k].content);
			}
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, normal_distributed_wgt_map != NULL && normal_distributed_wgt_elements_vect.size() > 0, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
		delete [] normal_distributed_wgt_map;
		normal_distributed_wgt_map = new common_sort_struct<Normal_distributed_wgt_element> [normal_distributed_wgt_elements_vect.size()];
		for (int i = 0; i < normal_distributed_wgt_elements_vect.size(); i ++) {
			normal_distributed_wgt_map[i].content = normal_distributed_wgt_elements_vect[i];
			normal_distributed_wgt_map[i].target_proc_id = normal_distributed_wgt_elements_vect[i].owner_process_id;
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, k == num_wgts_after_redistribution, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights: %d vs %d", k, num_wgts_after_redistribution);
	}
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, normal_distributed_wgt_map == NULL, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
	num_wgts_after_redistribution = normal_distributed_wgt_elements_vect.size();
	normal_distributed_wgt_elements_vect.clear();

	Distribute_merge_sort<Normal_distributed_wgt_element> *wgts_rearrange_distribute_sorting = new Distribute_merge_sort<Normal_distributed_wgt_element>(dst_comp_node->get_current_proc_local_id(), dst_comp_node->get_current_proc_local_id(), dst_comp_node, dst_comp_node);
	wgts_rearrange_distribute_sorting->do_data_sorting_with_target_process_id(&normal_distributed_wgt_map, dst_decomp_info, dst_comp_node, &num_wgts_after_redistribution);
	delete wgts_rearrange_distribute_sorting;

	if (num_wgts_after_redistribution > 0)
		current_remap_local_cell_global_indexes = new int [num_wgts_after_redistribution];
	else current_remap_local_cell_global_indexes = new int [1];
	for (i = 0; i < num_wgts_after_redistribution; i ++)
		normal_distributed_wgt_map[i].key = normal_distributed_wgt_map[i].content.src_cell_index;
	do_quick_sort(normal_distributed_wgt_map, (int*) NULL, 0, num_wgts_after_redistribution-1);
	for (i = 0, j = 0; i < num_wgts_after_redistribution; i ++) {
		if (i > 0 && normal_distributed_wgt_map[i-1].key != normal_distributed_wgt_map[i].key)
			j ++;
		current_remap_local_cell_global_indexes[j] = normal_distributed_wgt_map[i].content.src_cell_index + 1;
		normal_distributed_wgt_map[i].content.src_cell_index = j;
	}
	num_current_remap_local_cell_global_indexes = 0;
	if (num_wgts_after_redistribution > 0)
		num_current_remap_local_cell_global_indexes = j + 1;
	sort_normal_distributed_remapping_weight_elements_locally(normal_distributed_wgt_map, num_wgts_after_redistribution);
	normal_remap_weights = new Remap_weight_sparse_matrix(entire_remap_operator);
	for (int i = 0; i < num_wgts_after_redistribution; i ++)
		normal_remap_weights->add_weights(&(normal_distributed_wgt_map[i].content.src_cell_index), normal_distributed_wgt_map[i].content.dst_local_index, &(normal_distributed_wgt_map[i].content.wgt_value), 1, false);

	if (dst_grid_cell_rearrange_map != NULL)
		delete [] dst_grid_cell_rearrange_map;
	delete cell_rearrange_distribute_sorting;
	if (normal_distributed_wgt_map != NULL)
		delete [] normal_distributed_wgt_map;

	comp_comm_group_mgt_mgr->get_root_component_model()->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "rearrange_for_local_H2D_parallel_weights");

	return normal_remap_weights;
}


Grid_cell_rearrange_map_element::Grid_cell_rearrange_map_element(int target_proc_id, int owner_proc_id, int cell_global_index, int local_index)
{
	this->domain_decomp_process_id = target_proc_id;
	this->owner_process_id = owner_proc_id;
	this->local_cell_global_index = cell_global_index;
	this->owner_local_index = local_index;
}


Normal_distributed_wgt_element::Normal_distributed_wgt_element(long src_cell_index, long dst_cell_index, double wgt_value)
{
	this->src_cell_index = src_cell_index;
	this->dst_cell_index = dst_cell_index;
	this->wgt_value = wgt_value;
}


void Distributed_H2D_weights_generator::append_wgt_sparse_matrix(std::vector<Normal_distributed_wgt_element> &normal_distributed_wgt_elements, Remap_weight_sparse_matrix *wgt_matrix)
{
	long *indexes_src = wgt_matrix->get_indexes_src_grid();
	long *indexes_dst = wgt_matrix->get_indexes_dst_grid();
	double *wgt_values = wgt_matrix->get_weight_values();

	for (int i = 0; i < wgt_matrix->get_num_weights(); i ++)
		normal_distributed_wgt_elements.push_back(Normal_distributed_wgt_element(indexes_src[i], indexes_dst[i], wgt_values[i]));
}


Distributed_H2D_weights_generator::Distributed_H2D_weights_generator(int comp_id, int src_original_grid_id, int dst_original_grid_id, int dst_decomp_id, Remap_operator_basis *entire_remap_operator)
{
	double common_min_lon, common_max_lon, common_min_lat, common_max_lat;
	double src_subdomain_min_lon, src_subdomain_max_lon, src_subdomain_min_lat, src_subdomain_max_lat;
	Remap_operator_basis *subdomain_remap_operator, *backup_remap_operator;
	Remap_operator_grid *backup_operator_grid_src, *backup_operator_grid_dst;
	std::vector<Normal_distributed_wgt_element> normal_distributed_wgt_elements_vect;
	Remap_grid_class *dst_decomp_grid;
	common_sort_struct<Normal_distributed_wgt_element> *normal_distributed_wgt_map = NULL;
	int num_active_dst_decomp_grid_local_cells = 0, send_recv_mark, proc_id_send_to, proc_id_recv_from, j, prev_j, next_j;
	int num_original_normal_distributed_wgt_elements, num_dst_decomp_grid_local_wgts, weight_group_begin_pos;
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	double time1, time2, time3, time4;
	char src_H2D_sub_grid_name[NAME_STR_SIZE], dst_H2D_sub_grid_name[NAME_STR_SIZE], full_default_wgt_file_name[NAME_STR_SIZE];
	Original_grid_info *src_original_grid, *dst_original_grid;
	Distributed_H2D_grid_engine *distributed_grid_src, *distributed_grid_dst;


	EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "begin to generate remapping distributed weights for operator \"%s\" from \%s\" to \"%s\": %d", entire_remap_operator->get_operator_name(), entire_remap_operator->get_src_grid()->get_grid_name(), entire_remap_operator->get_dst_grid()->get_grid_name(), entire_remap_operator->get_extrapolate_enabled()?1:0);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, decomps_info_mgr->is_decomp_id_legal(dst_decomp_id), "Software error in Distributed_H2D_weights_generator::Distributed_H2D_weights_generator: %lx", dst_decomp_id);

	this->comp_id = comp_id;
	this->dst_decomp_id = dst_decomp_id;
	this->entire_remap_operator = entire_remap_operator;
	backup_operator_grid_src = current_runtime_remap_operator_grid_src;
	backup_operator_grid_dst = current_runtime_remap_operator_grid_dst;
	backup_remap_operator = current_runtime_remap_operator;
	current_distributed_H2D_weights_generator = this;

	wtime(&time1);
	distributed_H2D_grid_mgr->calculate_common_grid_domain_for_remapping(comp_id, entire_remap_operator->get_src_grid(), entire_remap_operator->get_dst_grid(), common_min_lon, common_max_lon, common_min_lat, common_max_lat, false);
	remapping_grid_domain_decomp_engine = new Remapping_grid_domain_decomp_engine(src_original_grid_id, comp_id, entire_remap_operator->get_src_grid(), entire_remap_operator->get_dst_grid(),common_min_lon, common_max_lon, common_min_lat, common_max_lat);
	wtime(&time2);
	MPI_Barrier(comp_node->get_comm_group());
	EXECUTION_REPORT(REPORT_LOG, -1, true, "Time for generating Remapping_grid_domain_decomp_engine is %lf", time2-time1);
	for (int i = 0; i < remapping_grid_domain_decomp_engine->get_num_subdomains(); i ++) {
		Remap_grid_class *current_dst_subdomain_grid = remapping_grid_domain_decomp_engine->get_dst_subdomain(i);
		if (current_dst_subdomain_grid == NULL || current_dst_subdomain_grid->get_num_active_cells() == 0)
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_dst_subdomain_grid->get_num_vertexes(), "Software error in Distributed_H2D_weights_generator::Distributed_H2D_weights_generator");
		Remap_grid_class *current_src_subdomain_grid = remapping_grid_domain_decomp_engine->get_src_expanded_subdomain(i);
		remapping_grid_domain_decomp_engine->get_src_current_expanded_subdomain_boundaries(i, src_subdomain_min_lon, src_subdomain_max_lon, src_subdomain_min_lat, src_subdomain_max_lat);
		while (!current_dst_subdomain_grid->are_all_active_cells_in_subdomain(src_subdomain_min_lon, src_subdomain_max_lon, src_subdomain_min_lat, src_subdomain_max_lat)) {
			current_src_subdomain_grid = remapping_grid_domain_decomp_engine->expand_src_subdomain_halo_grid(i);
			remapping_grid_domain_decomp_engine->get_src_current_expanded_subdomain_boundaries(i, src_subdomain_min_lon, src_subdomain_max_lon, src_subdomain_min_lat, src_subdomain_max_lat);
		}
		if (!entire_remap_operator->get_extrapolate_enabled() && current_src_subdomain_grid == NULL)
			continue;
		while (current_src_subdomain_grid == NULL) {
			current_src_subdomain_grid = remapping_grid_domain_decomp_engine->expand_src_subdomain_halo_grid(i);			
			remapping_grid_domain_decomp_engine->get_src_current_expanded_subdomain_boundaries(i, src_subdomain_min_lon, src_subdomain_max_lon, src_subdomain_min_lat, src_subdomain_max_lat);
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_src_subdomain_grid->get_grid_size() > 0 && current_src_subdomain_grid->get_num_vertexes() > 0, "Software error in Distributed_H2D_weights_generator::Distributed_H2D_weights_generator");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_src_subdomain_grid->get_boundary_min_lon() == src_subdomain_min_lon && current_src_subdomain_grid->get_boundary_max_lon() == src_subdomain_max_lon && current_src_subdomain_grid->get_boundary_min_lat() == src_subdomain_min_lat && current_src_subdomain_grid->get_boundary_max_lat() == src_subdomain_max_lat, "Software error in Distributed_H2D_weights_generator::Distributed_H2D_weights_generator: %lx, %lx  vs  %lx, %lx", current_src_subdomain_grid->get_boundary_min_lon(), src_subdomain_min_lon && current_src_subdomain_grid->get_boundary_max_lon(), src_subdomain_max_lon);
		subdomain_remap_operator = entire_remap_operator->duplicate_remap_operator(false);
		subdomain_remap_operator->initialize();
		subdomain_remap_operator->set_src_grid(current_src_subdomain_grid);
		subdomain_remap_operator->set_dst_grid(current_dst_subdomain_grid);
		current_subdomain_index = i;
	    current_runtime_remap_operator_grid_src = new Remap_operator_grid(current_src_subdomain_grid, subdomain_remap_operator, true, false);
    	current_runtime_remap_operator_grid_dst = new Remap_operator_grid(current_dst_subdomain_grid, subdomain_remap_operator, false, false);
		current_runtime_remap_operator = subdomain_remap_operator;
		current_runtime_remap_operator_grid_src->update_operator_grid_data();
		current_runtime_remap_operator_grid_dst->update_operator_grid_data();
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "use %s to calculate remapping weight from \"%s\"(%lf %lf %lf %lf) to \"%s\"(%lf %lf %lf %lf) (%d grid cells) : %d %d :  %d", subdomain_remap_operator->get_operator_name(), entire_remap_operator->get_src_grid()->get_grid_name(), src_subdomain_min_lon, src_subdomain_max_lon, src_subdomain_min_lat, src_subdomain_max_lat, entire_remap_operator->get_dst_grid()->get_grid_name(), current_dst_subdomain_grid->get_boundary_min_lon(), current_dst_subdomain_grid->get_boundary_max_lon(), current_dst_subdomain_grid->get_boundary_min_lat(), current_dst_subdomain_grid->get_boundary_max_lat(), current_dst_subdomain_grid->get_grid_size(), current_src_subdomain_grid->get_grid_size(), current_src_subdomain_grid->get_num_active_cells(), entire_remap_operator->get_extrapolate_enabled()?1:0);
		subdomain_remap_operator->calculate_remap_weights();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_remap_operator->get_num_remap_weights_groups() == 1, "Software error in Distributed_H2D_weights_generator::Distributed_H2D_weights_generator");
		if (words_are_the_same(subdomain_remap_operator->get_operator_name(), REMAP_OPERATOR_NAME_BILINEAR) || words_are_the_same(subdomain_remap_operator->get_operator_name(), REMAP_OPERATOR_NAME_DISTWGT) || words_are_the_same(subdomain_remap_operator->get_operator_name(), REMAP_OPERATOR_NAME_CONSERV_2D))
			append_wgt_sparse_matrix(normal_distributed_wgt_elements_vect, subdomain_remap_operator->get_remap_weights_group(0));
		delete current_runtime_remap_operator_grid_src;
		delete current_runtime_remap_operator_grid_dst;
		delete subdomain_remap_operator;
	}

	current_runtime_remap_operator_grid_src = backup_operator_grid_src;
	current_runtime_remap_operator_grid_dst = backup_operator_grid_dst;
	current_runtime_remap_operator = backup_remap_operator;

	if (normal_distributed_wgt_elements_vect.size() > 0) {
		normal_distributed_wgt_map = new common_sort_struct<Normal_distributed_wgt_element> [normal_distributed_wgt_elements_vect.size()];
		for (int i = 0; i < normal_distributed_wgt_elements_vect.size(); i ++) {
			normal_distributed_wgt_map[i].content = normal_distributed_wgt_elements_vect[i];
		}
	}
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Distributed_H2D_weights_generator has %d weights locally", normal_distributed_wgt_elements_vect.size());
	if (num_active_dst_decomp_grid_local_cells > 0 && entire_remap_operator->get_extrapolate_enabled())
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, normal_distributed_wgt_elements_vect.size() > 0, "Software error Distributed_H2D_weights_generator::Distributed_H2D_weights_generator");

	original_grid_mgr->get_original_grid(src_original_grid_id)->get_H2D_sub_grid_full_name(src_H2D_sub_grid_name);
	original_grid_mgr->get_original_grid(dst_original_grid_id)->get_H2D_sub_grid_full_name(dst_H2D_sub_grid_name);
	sprintf(full_default_wgt_file_name, "%s/DEFAULT_WGT__%s__FROM__%s__TO__%s.nc", comp_comm_group_mgt_mgr->get_internal_remapping_weights_dir(), entire_remap_operator->get_operator_name(), src_H2D_sub_grid_name, dst_H2D_sub_grid_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid_id != -1 && dst_original_grid_id != -1, "Software error in rearrange_for_local_H2D_parallel_weights");
	if (src_original_grid_id != -1) {
		src_original_grid = original_grid_mgr->get_original_grid(src_original_grid_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid != NULL && src_original_grid->get_H2D_sub_CoR_grid(), "Software error in rearrange_for_local_H2D_parallel_weights");
		distributed_grid_src = distributed_H2D_grid_mgr->search_distributed_H2D_grid(comp_node->get_full_name(), src_original_grid->get_H2D_sub_CoR_grid());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, distributed_grid_src != NULL && distributed_grid_src->get_basic_decomp_grid() != NULL, "Software error in rearrange_for_local_H2D_parallel_weights");
		distributed_H2D_grid_mgr->generate_full_grid_data_from_distributed_grid(comp_id, src_original_grid->get_H2D_sub_CoR_grid(), distributed_grid_src->get_basic_decomp_grid()->get_decomp_grid(), false);
	}
	if (dst_original_grid_id != -1) {
		dst_original_grid = original_grid_mgr->get_original_grid(dst_original_grid_id);		
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid != NULL, "Software error in rearrange_for_local_H2D_parallel_weights");
		distributed_grid_dst = distributed_H2D_grid_mgr->search_distributed_H2D_grid(comp_node->get_full_name(), dst_original_grid->get_H2D_sub_CoR_grid());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, distributed_grid_dst != NULL && distributed_grid_dst->get_basic_decomp_grid() != NULL, "Software error in rearrange_for_local_H2D_parallel_weights");
		distributed_H2D_grid_mgr->generate_full_grid_data_from_distributed_grid(comp_id, dst_original_grid->get_H2D_sub_CoR_grid(), distributed_grid_dst->get_basic_decomp_grid()->get_decomp_grid(), false);
	}
	
	normal_remap_weights = rearrange_for_local_H2D_parallel_weights(decomps_info_mgr->get_decomp_info(dst_decomp_id), entire_remap_operator, normal_distributed_wgt_map, normal_distributed_wgt_elements_vect.size(), full_default_wgt_file_name, src_original_grid, dst_original_grid);

	delete remapping_grid_domain_decomp_engine;
	remapping_grid_domain_decomp_engine = NULL;

	wtime(&time4);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "Time for distributed weight generation is %lf vs %lf", time3-time2, time4-time3);

	EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "Finish generating remapping distributed weights for operator \"%s\" from \%s\" to \"%s\" with %d local weights", entire_remap_operator->get_operator_name(), entire_remap_operator->get_src_grid()->get_grid_name(), entire_remap_operator->get_dst_grid()->get_grid_name(), num_dst_decomp_grid_local_wgts);
}


void Distributed_H2D_weights_generator::check_consistency_of_normal_remap_weights(Remap_weight_sparse_matrix *another_normal_remap_weights)
{
	if (!report_error_enabled)
		return;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, normal_remap_weights->get_num_weights() == another_normal_remap_weights->get_num_weights(), "Software error in Distributed_H2D_weights_generator::check_consistency_of_normal_remap_weights: %d vs %d", normal_remap_weights->get_num_weights(), another_normal_remap_weights->get_num_weights());

	long *this_indexes_src = normal_remap_weights->get_indexes_src_grid();
	long *this_indexes_dst = normal_remap_weights->get_indexes_dst_grid();
	double *this_wgt_values = normal_remap_weights->get_weight_values();
	long *another_indexes_src = another_normal_remap_weights->get_indexes_src_grid();
	long *another_indexes_dst = another_normal_remap_weights->get_indexes_dst_grid();
	double *another_wgt_values = another_normal_remap_weights->get_weight_values();
	
	for (int i = 0; i < normal_remap_weights->get_num_weights(); i ++) {
		if (!(this_indexes_dst[i] == another_indexes_dst[i] && this_wgt_values[i] == another_wgt_values[i])) {
			int start_j = i - 4;
			if (start_j < 0)
				start_j = 0;
			int end_j = i + 4;
			if (end_j > normal_remap_weights->get_num_weights())
				end_j = normal_remap_weights->get_num_weights();
			for (int j = start_j; j < end_j; j ++) 
				EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Software error in Distributed_H2D_weights_generator::check_consistency_of_normal_remap_weights: %d (%d %d %lf) vs (%d %d %lf)", j, this_indexes_src[j], this_indexes_dst[j], this_wgt_values[j], another_indexes_src[j], another_indexes_dst[j], another_wgt_values[j]);			
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, this_indexes_dst[i] == another_indexes_dst[i] && this_wgt_values[i] == another_wgt_values[i], "Software error in Distributed_H2D_weights_generator::check_consistency_of_normal_remap_weights: %d (%d %d %lf) vs (%d %d %lf)", i, this_indexes_src[i], this_indexes_dst[i], this_wgt_values[i], another_indexes_src[i], another_indexes_dst[i], another_wgt_values[i]);
		}
	}
}


bool Distributed_H2D_weights_generator::should_enlarge_src_subdomain_grid_for_remapping(long dst_cell_index, double radius)
{
	double subdomain_grid_min_lon, subdomain_grid_max_lon, subdomain_grid_min_lat, subdomain_grid_max_lat;
	double overall_grid_min_lon, overall_grid_max_lon, overall_grid_min_lat, overall_grid_max_lat;
	double common_grid_min_lon, common_grid_max_lon, common_grid_min_lat, common_grid_max_lat;
	bool is_sphere_grid = words_are_the_same(current_runtime_remap_operator->get_src_grid()->get_a_leaf_grid(COORD_LABEL_LON)->get_coord_unit(), COORD_UNIT_DEGREES);
	double distances_to_boundaries[4], dst_center_values[2];
	int iter = 1;
	bool check_necessity[4];

	dst_center_values[0] = ((double*)current_runtime_remap_operator->get_dst_grid()->get_grid_center_field(COORD_LABEL_LON)->get_grid_data_field()->data_buf)[dst_cell_index];
	dst_center_values[1] = ((double*)current_runtime_remap_operator->get_dst_grid()->get_grid_center_field(COORD_LABEL_LAT)->get_grid_data_field()->data_buf)[dst_cell_index];
	get_cell_center_coord_values_of_grid(current_runtime_remap_operator_grid_dst, dst_cell_index, dst_center_values);
	subdomain_grid_min_lon = current_runtime_remap_operator->get_src_grid()->get_boundary_min_lon();
	subdomain_grid_max_lon = current_runtime_remap_operator->get_src_grid()->get_boundary_max_lon();
	subdomain_grid_min_lat = current_runtime_remap_operator->get_src_grid()->get_boundary_min_lat();
	subdomain_grid_max_lat = current_runtime_remap_operator->get_src_grid()->get_boundary_max_lat();	
	overall_grid_min_lon = entire_remap_operator->get_src_grid()->get_boundary_min_lon();
	overall_grid_max_lon = entire_remap_operator->get_src_grid()->get_boundary_max_lon();
	overall_grid_min_lat = entire_remap_operator->get_src_grid()->get_boundary_min_lat();
	overall_grid_max_lat = entire_remap_operator->get_src_grid()->get_boundary_max_lat();	

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, radius > 0, "Software error in Remapping_grid_domain_decomp_engine::should_enlarge_subdomain_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, subdomain_grid_min_lon != NULL_COORD_VALUE && subdomain_grid_min_lat != NULL_COORD_VALUE && overall_grid_min_lon != NULL_COORD_VALUE && overall_grid_min_lat != NULL_COORD_VALUE, "Software error in Remapping_grid_domain_decomp_engine::should_enlarge_subdomain_grid");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_center_values[1] <= subdomain_grid_max_lat && dst_center_values[1] >= subdomain_grid_min_lat, "Software error in Remapping_grid_domain_decomp_engine::should_enlarge_subdomain_grid");
	if (subdomain_grid_min_lon < subdomain_grid_max_lon) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_center_values[0] <= subdomain_grid_max_lon && dst_center_values[0] >= subdomain_grid_min_lon, "Software error in Remapping_grid_domain_decomp_engine::should_enlarge_subdomain_grid");
	}
	else {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_sphere_grid, "Software error in Remapping_grid_domain_decomp_engine::should_enlarge_subdomain_grid");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_center_values[0] <= subdomain_grid_max_lon || dst_center_values[0] >= subdomain_grid_min_lon);
	}

	((Distributed_H2D_grid_mgt*) NULL)->calculate_common_grid_domain_for_remapping(comp_id, current_runtime_remap_operator->get_src_grid(), entire_remap_operator->get_src_grid(), common_grid_min_lon, common_grid_max_lon, common_grid_min_lat, common_grid_max_lat, true);

	check_necessity[0] = (subdomain_grid_min_lon != common_grid_min_lon);
	check_necessity[1] = (subdomain_grid_max_lon != common_grid_max_lon);
	if (entire_remap_operator->get_src_grid()->get_grid_cyclic() && (subdomain_grid_min_lon != common_grid_min_lon || subdomain_grid_max_lon != common_grid_max_lon)) {
		check_necessity[0] = true;
		check_necessity[1] = true;
	}
	check_necessity[2] = (subdomain_grid_min_lat != common_grid_min_lat);
	check_necessity[3] = (subdomain_grid_max_lat != common_grid_max_lat);
	distances_to_boundaries[0] = calculate_distance_of_two_points_2D(subdomain_grid_min_lon, dst_center_values[1], dst_center_values[0], dst_center_values[1], is_sphere_grid);
	distances_to_boundaries[1] = calculate_distance_of_two_points_2D(subdomain_grid_max_lon, dst_center_values[1], dst_center_values[0], dst_center_values[1], is_sphere_grid);
	distances_to_boundaries[2] = calculate_distance_of_two_points_2D(dst_center_values[0], subdomain_grid_min_lat, dst_center_values[0], dst_center_values[1], is_sphere_grid);
	distances_to_boundaries[3] = calculate_distance_of_two_points_2D(dst_center_values[0], subdomain_grid_max_lat, dst_center_values[0], dst_center_values[1], is_sphere_grid);
	for (int i = 0; i < 4; i ++)
		if (check_necessity[i] && (distances_to_boundaries[i] < radius || relative_eq(distances_to_boundaries[i], radius))) {
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Should enlarge src subdomain grid %d with distance %lf (vs %lf) from center %lf %lf", current_subdomain_index, radius, distances_to_boundaries[i], dst_center_values[0], dst_center_values[1]);
			return true;
		}

	return false;
}


bool Distributed_H2D_weights_generator::confirm_or_enlarge_current_src_subdomain_grid_for_remapping(long dst_cell_index, double radius)
{	
	Remap_grid_class *current_src_subdomain_grid = NULL;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remapping_grid_domain_decomp_engine != NULL, "Software error in Distributed_H2D_weights_generator::enlarge_current_src_subdomain_grid_for_remapping");

	while (should_enlarge_src_subdomain_grid_for_remapping(dst_cell_index, radius)) {
		if (current_src_subdomain_grid == NULL)
			delete current_runtime_remap_operator_grid_src;
		current_src_subdomain_grid = remapping_grid_domain_decomp_engine->expand_src_subdomain_halo_grid(current_subdomain_index);
		current_runtime_remap_operator->set_src_grid(current_src_subdomain_grid);
	}

	if (current_src_subdomain_grid != NULL) {
    	current_runtime_remap_operator_grid_src = new Remap_operator_grid(current_src_subdomain_grid, current_runtime_remap_operator, true, false);
		current_runtime_remap_operator_grid_src->update_operator_grid_data();
		return true;
	}

	return false;
}


Remap_weight_sparse_matrix *Distributed_H2D_weights_generator::extract_normal_remap_weights()
{
	Remap_weight_sparse_matrix *current_normal_remap_weights = normal_remap_weights;
	normal_remap_weights = NULL;
	return current_normal_remap_weights;
}


Remap_weight_sparse_matrix *Distributed_H2D_weights_generator::read_normal_remap_weights(Decomp_info *dst_decomp_info, const char *wgt_file_name, Remap_operator_basis *entire_remap_operator)
{
	int tile_size = 100000, max_num_processes_for_reading = 400, IO_process_status = 0, IO_proc_id;
	size_t start = 0, count = 0;
	int ncfile_id, var_id;
	common_sort_struct<Grid_cell_rearrange_map_element> *dst_grid_cell_rearrange_map = NULL;
	common_sort_struct<Normal_distributed_wgt_element> *normal_distributed_wgt_map = NULL;
	int num_cells_after_redistribution, num_wgts_after_redistribution;
	std::vector<Normal_distributed_wgt_element> normal_distributed_wgt_elements_vect;
	long i, j, k;

	Comp_comm_group_mgt_node *dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_decomp_info->get_comp_id());	
    IO_netcdf *netcdf_file_object = new IO_netcdf("remapping weights file for H2D interpolation", wgt_file_name, "r", false);
	long num_total_weights = netcdf_file_object->get_dimension_size("n_s", dst_comp_node->get_comm_group(), dst_comp_node->get_current_proc_local_id() == 0);
	int num_processes_for_reading = (num_total_weights+tile_size-1) / tile_size;
	int num_cells_per_proc = (dst_decomp_info->get_num_global_cells()+dst_comp_node->get_num_procs()-1)/dst_comp_node->get_num_procs();

	delete netcdf_file_object;
	if (num_processes_for_reading > dst_comp_node->get_num_procs())
		num_processes_for_reading = dst_comp_node->get_num_procs();
	if (num_processes_for_reading > max_num_processes_for_reading)
		num_processes_for_reading = max_num_processes_for_reading;
	int proc_id_interval = dst_comp_node->get_num_procs() / num_processes_for_reading;
	if (dst_comp_node->get_current_proc_local_id() % proc_id_interval == 0 && dst_comp_node->get_current_proc_local_id() < proc_id_interval * num_processes_for_reading)
		IO_process_status = 1;

	if (report_error_enabled) {
		int num_IO_processes;
		MPI_Allreduce(&IO_process_status, &num_IO_processes, 1, MPI_INT, MPI_SUM, dst_comp_node->get_comm_group());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_IO_processes == num_processes_for_reading, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
	}

	if (IO_process_status == 1) {
		IO_proc_id = dst_comp_node->get_current_proc_local_id() / proc_id_interval;
		count = num_total_weights / num_processes_for_reading;
		start = IO_proc_id * count;
		if (IO_proc_id < num_total_weights % num_processes_for_reading) {
			count ++;
			start += IO_proc_id;
		}
		else start += num_total_weights % num_processes_for_reading;
	}

	if (report_error_enabled) {
		if (IO_process_status == 1 && IO_proc_id == num_processes_for_reading -1)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, start + count == num_total_weights, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
		int total_counts;
		MPI_Allreduce(&count, &total_counts, 1, MPI_INT, MPI_SUM, dst_comp_node->get_comm_group());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, total_counts == num_total_weights, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
	}

	if (IO_process_status == 1) {
		comp_comm_group_mgt_mgr->get_root_component_model()->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "read_normal_remap_weights from NC file");
		int rcode = nc_open(wgt_file_name, NC_NOWRITE, &ncfile_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, count > 0, "Software error in Distributed_H2D_weights_generator::read_normal_remap_weights");
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		int *wgts_row = new int [count];
		int *wgts_col = new int [count];
		double *wgts_values = new double [count];
		rcode = nc_inq_varid(ncfile_id, "col", &var_id);
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		rcode = nc_get_vara_int(ncfile_id, var_id, &start, &count, wgts_col);		
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		rcode = nc_inq_varid(ncfile_id, "row", &var_id);
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		rcode = nc_get_vara_int(ncfile_id, var_id, &start, &count, wgts_row);		
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		rcode = nc_inq_varid(ncfile_id, "S", &var_id);
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		rcode = nc_get_vara_double(ncfile_id, var_id, &start, &count, wgts_values);		
		EXECUTION_REPORT(REPORT_ERROR, -1, rcode == NC_NOERR, "Netcdf error: %s for file %s\n", nc_strerror(rcode), wgt_file_name);
		comp_comm_group_mgt_mgr->get_root_component_model()->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "read_normal_remap_weights from NC file");
		normal_distributed_wgt_map = new common_sort_struct<Normal_distributed_wgt_element> [count];
		for (int i = 0; i < count; i ++)
			new (&normal_distributed_wgt_map[i].content) Normal_distributed_wgt_element(wgts_col[i]-1, wgts_row[i]-1, wgts_values[i]);
		delete [] wgts_col;
		delete [] wgts_row;
		delete [] wgts_values;
	}
	num_wgts_after_redistribution = count;
	return rearrange_for_local_H2D_parallel_weights(dst_decomp_info, entire_remap_operator, normal_distributed_wgt_map, num_wgts_after_redistribution, NULL, NULL, NULL);
}

