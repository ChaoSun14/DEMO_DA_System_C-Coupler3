/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "runtime_remapping_weights_mgt.h"
#include "remap_operator_bilinear.h"
#include "remap_operator_linear.h"
#include "remap_operator_distwgt.h"
#include "remap_operator_conserv_2D.h"
#include "remap_operator_spline_1D.h"
#include "global_data.h"



Runtime_remapping_weights::Runtime_remapping_weights()
{
    src_comp_full_name = NULL;
    dst_comp_full_name = NULL;
    src_original_grid = NULL;
    dst_original_grid = NULL;
    src_decomp_info = NULL;
    dst_decomp_info = NULL;
    remapping_setting = NULL;
    remapping_strategy = NULL;
    sequential_remapping_weights = NULL;
    parallel_remapping_weights = NULL;
    intermediate_V3D_grid_bottom_field = NULL;
    dynamic_V1D_remap_weight_of_operator = NULL;
    runtime_V1D_remap_grid_src = NULL;
    runtime_V1D_remap_grid_dst = NULL;
    src_H2D_grid_area = NULL;
    dst_H2D_grid_area = NULL;
}


Runtime_remapping_weights::Runtime_remapping_weights(const char *src_comp_full_name, const char *dst_comp_full_name, Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, Remapping_setting *remapping_setting, Decomp_info *dst_decomp_info, H2D_remapping_wgt_file_info *H2D_remapping_weight_file)
{
    Remap_operator_basis *remap_operator_H2D = NULL;
    Remap_operator_basis *remap_operator_V1D = NULL;
    Remap_operator_basis *remap_operator_T1D = NULL;
    Remap_operator_basis *remap_operators[3];
    Remap_grid_class *remap_grids[2];
    Remapping_setting *cloned_remapping_setting = remapping_setting->clone();
    char parameter_name[NAME_STR_SIZE], parameter_value[NAME_STR_SIZE], remap_weight_name[NAME_STR_SIZE];
    int num_remap_operators = 0;


    this->src_comp_full_name = strdup(src_comp_full_name);
    this->dst_comp_full_name = strdup(dst_comp_full_name);
    this->src_original_grid = src_original_grid;
    this->dst_original_grid = dst_original_grid;
    this->remapping_setting = cloned_remapping_setting;
    this->dst_decomp_info = dst_decomp_info;
    this->src_decomp_info = NULL;
    this->sequential_remapping_weights = NULL;
    this->parallel_remapping_weights = NULL;
    this->intermediate_V3D_grid_bottom_field = NULL;
    this->dynamic_V1D_remap_weight_of_operator = NULL;
    this->runtime_V1D_remap_grid_src = NULL;
    this->runtime_V1D_remap_grid_dst = NULL;
    this->src_H2D_grid_area = NULL;
    this->dst_H2D_grid_area = NULL;


    if (src_original_grid->get_H2D_sub_CoR_grid() != NULL) {
        remap_grids[0] = src_original_grid->get_H2D_sub_CoR_grid();
        remap_grids[1] = dst_original_grid->get_H2D_sub_CoR_grid();
        if (words_are_the_same(cloned_remapping_setting->get_H2D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_BILINEAR))
            remap_operator_H2D = new Remap_operator_bilinear("H2D_algorithm", 2, remap_grids);
        else if (words_are_the_same(cloned_remapping_setting->get_H2D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_CONSERV_2D)) 
            remap_operator_H2D = new Remap_operator_conserv_2D("H2D_algorithm", 2,  remap_grids);
        else if (words_are_the_same(cloned_remapping_setting->get_H2D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_DISTWGT))
            remap_operator_H2D = new Remap_operator_distwgt("H2D_algorithm", 2,  remap_grids);
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Runtime_remapping_weights::Runtime_remapping_weights: wrong H2D algorithm");
        for (int i = 0; i < cloned_remapping_setting->get_H2D_remapping_algorithm()->get_num_parameters(); i ++) {
            cloned_remapping_setting->get_H2D_remapping_algorithm()->get_parameter(i, parameter_name, parameter_value);
            remap_operator_H2D->set_parameter(parameter_name, parameter_value);
        }
        remap_operators[num_remap_operators++] = remap_operator_H2D;
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "new remap operator \"%s\" from grid \"%s\" to \"%s\": %d", remap_operator_H2D->get_operator_name(), src_original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), dst_original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), remap_operator_H2D->get_extrapolate_enabled()?1:0);
    }
    if (src_original_grid->get_V1D_sub_CoR_grid() != NULL) {
        remap_grids[0] = src_original_grid->get_V1D_sub_CoR_grid();
        remap_grids[1] = dst_original_grid->get_V1D_sub_CoR_grid();
        if (words_are_the_same(cloned_remapping_setting->get_V1D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_LINEAR))
            remap_operator_V1D = new Remap_operator_linear("V1D_algorithm", 2, remap_grids);
        else if (words_are_the_same(cloned_remapping_setting->get_V1D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_SPLINE_1D))
            remap_operator_V1D = new Remap_operator_spline_1D("V1D_algorithm", 2, remap_grids);
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Runtime_remapping_weights::Runtime_remapping_weights: wrong V1D algorithm");
        for (int i = 0; i < cloned_remapping_setting->get_V1D_remapping_algorithm()->get_num_parameters(); i ++) {
            cloned_remapping_setting->get_V1D_remapping_algorithm()->get_parameter(i, parameter_name, parameter_value);
            remap_operator_V1D->set_parameter(parameter_name, parameter_value);
        }
        remap_operators[num_remap_operators++] = remap_operator_V1D;
		if (!src_original_grid->get_original_CoR_grid()->is_sigma_grid() && !src_original_grid->get_original_CoR_grid()->does_use_V3D_level_coord())
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, dst_original_grid->get_comp_id(), src_original_grid->get_V1D_sub_CoR_grid()->get_super_grid_of_setting_coord_values() != NULL, "Error happens when generating remapping weights from the grid \"%s\" to \"%s\": the coordinate values of the vertical sub grid \"%s\" of the source grid has not been specified. Please verify.", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name(), src_original_grid->get_V1D_sub_CoR_grid()->get_grid_name());
		if (!dst_original_grid->get_original_CoR_grid()->is_sigma_grid() && !dst_original_grid->get_original_CoR_grid()->does_use_V3D_level_coord())
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, dst_original_grid->get_comp_id(), dst_original_grid->get_V1D_sub_CoR_grid()->get_super_grid_of_setting_coord_values() != NULL, "Error happens when generating remapping weights from the grid \"%s\" to \"%s\": the coordinate values of the vertical sub grid \"%s\" of the target grid has not been specified. Please verify.", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name(), dst_original_grid->get_V1D_sub_CoR_grid()->get_grid_name());
    }
    if (src_original_grid->get_Time1D_sub_CoR_grid() != NULL) {
        remap_grids[0] = src_original_grid->get_Time1D_sub_CoR_grid();
        remap_grids[1] = dst_original_grid->get_Time1D_sub_CoR_grid();
        if (words_are_the_same(cloned_remapping_setting->get_T1D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_LINEAR))
            remap_operator_T1D = new Remap_operator_linear("T1D_algorithm", 2, remap_grids);
        else if (words_are_the_same(cloned_remapping_setting->get_T1D_remapping_algorithm()->get_algorithm_name(), REMAP_OPERATOR_NAME_SPLINE_1D))
            remap_operator_T1D = new Remap_operator_spline_1D("T1D_algorithm", 2, remap_grids);
        else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Runtime_remapping_weights::Runtime_remapping_weights: wrong T1D algorithm");        
        for (int i = 0; i < cloned_remapping_setting->get_T1D_remapping_algorithm()->get_num_parameters(); i ++) {
            cloned_remapping_setting->get_T1D_remapping_algorithm()->get_parameter(i, parameter_name, parameter_value);
            remap_operator_T1D->set_parameter(parameter_name, parameter_value);
        }
        remap_operators[num_remap_operators++] = remap_operator_T1D;
    }
    execution_phase_number = 1;
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_remap_operators > 0, "Software error in Runtime_remapping_weights::Runtime_remapping_weights: no remapping operator");
    remapping_strategy = new Remap_strategy_class("runtime_remapping_strategy", num_remap_operators, remap_operators);
    EXECUTION_REPORT_LOG(REPORT_LOG, dst_original_grid->get_comp_id(), true, "before generating sequential_remapping_weights from original grid %s to %s", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name());    
    sprintf(remap_weight_name, "weights_%lx_%s(%s)_to_%s(%s)", remapping_setting->calculate_checksum(), src_original_grid->get_grid_name(), src_comp_full_name, dst_original_grid->get_grid_name(), dst_comp_full_name);
    if (H2D_remapping_weight_file != NULL) {
        EXECUTION_REPORT(REPORT_PROGRESS, dst_original_grid->get_comp_id(), true, "The remapping weight file \"%s\" will be used for data remapping from the horizontal grid \"%s\" (of the component model \"%s\") to the horizontal grid \"%s\" (of the component model \"%s\").", H2D_remapping_weight_file->get_wgt_file_name(), src_original_grid->get_grid_name(), src_comp_full_name, dst_original_grid->get_grid_name(), dst_comp_full_name);
        sequential_remapping_weights = new Remap_weight_of_strategy_class(remap_weight_name, remapping_strategy, src_original_grid->get_original_CoR_grid()->get_ordered_similar_grid_under_V3D(), dst_original_grid->get_original_CoR_grid()->get_ordered_similar_grid_under_V3D(), H2D_remapping_weight_file->get_wgt_file_name(), true, comp_comm_group_mgt_mgr->search_global_node(dst_comp_full_name)->get_comp_id(), dst_decomp_info->get_decomp_id(), src_original_grid->get_grid_id(), dst_original_grid->get_grid_id());
        if (src_original_grid->is_H2D_grid()) 
            set_H2D_grids_area(H2D_remapping_weight_file->get_src_area(), H2D_remapping_weight_file->get_dst_area(), src_original_grid->get_original_CoR_grid()->get_grid_size(), dst_original_grid->get_original_CoR_grid()->get_grid_size());
    }    
    else if (dst_original_grid->get_H2D_sub_CoR_grid() != NULL) {    
        EXECUTION_REPORT(REPORT_PROGRESS, dst_original_grid->get_comp_id(), true, "No remapping weight file has been specified for data remapping from the horizontal sub grid of \"%s\" (of the component model \"%s\") to the horizontal sub grid of \"%s\" (of the component model \"%s\"). So the remapping weights will be generated by C-Coupler", src_original_grid->get_grid_name(), src_comp_full_name, dst_original_grid->get_grid_name(), dst_comp_full_name);
        sequential_remapping_weights = new Remap_weight_of_strategy_class(remap_weight_name, remapping_strategy, src_original_grid->get_original_CoR_grid()->get_ordered_similar_grid_under_V3D(), dst_original_grid->get_original_CoR_grid()->get_ordered_similar_grid_under_V3D(), NULL, true, comp_comm_group_mgt_mgr->search_global_node(dst_comp_full_name)->get_comp_id(), dst_decomp_info->get_decomp_id(), src_original_grid->get_grid_id(), dst_original_grid->get_grid_id());
        if (src_original_grid->is_H2D_grid() && src_original_grid->get_original_CoR_grid()->get_area_or_volumn() != NULL)
            set_H2D_grids_area(src_original_grid->get_original_CoR_grid()->get_area_or_volumn(), dst_original_grid->get_original_CoR_grid()->get_area_or_volumn(), src_original_grid->get_original_CoR_grid()->get_grid_size(), dst_original_grid->get_original_CoR_grid()->get_grid_size());
    }    
    EXECUTION_REPORT_LOG(REPORT_LOG, dst_original_grid->get_comp_id(), true, "after generating sequential_remapping_weights from original grid %s to %s", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name());    
    execution_phase_number = 2;

    if (dst_original_grid->get_H2D_sub_CoR_grid() == NULL || dst_decomp_info == NULL) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_original_grid->get_H2D_sub_CoR_grid() == NULL && dst_decomp_info == NULL, "Software error in Coupling_connection::generate_interpolation: conflict between grid and decomp");
        parallel_remapping_weights = sequential_remapping_weights;
    }    
    else {
        generate_parallel_remapping_weights();
        delete sequential_remapping_weights;
        sequential_remapping_weights = NULL;
    }
}


Runtime_remapping_weights::~Runtime_remapping_weights()
{
    if (src_comp_full_name != NULL)
        delete [] src_comp_full_name;
    if (dst_comp_full_name != NULL)
        delete [] dst_comp_full_name;
    if (remapping_setting != NULL)
        delete remapping_setting;
    if (remapping_strategy != NULL)
        delete remapping_strategy;
    if (parallel_remapping_weights != NULL)
        delete parallel_remapping_weights;
    if (runtime_V1D_remap_grid_src != NULL)
        delete runtime_V1D_remap_grid_src;
    if (runtime_V1D_remap_grid_dst != NULL)
        delete runtime_V1D_remap_grid_dst;
    if (src_H2D_grid_area != NULL)
        delete [] src_H2D_grid_area;
    if (dst_H2D_grid_area != NULL)
        delete [] dst_H2D_grid_area;
}


Field_mem_info *Runtime_remapping_weights::allocate_intermediate_V3D_grid_bottom_field()
{
    if (intermediate_V3D_grid_bottom_field == NULL) {	
		if (src_original_grid->get_original_CoR_grid()->is_sigma_grid())
	        intermediate_V3D_grid_bottom_field = memory_manager->alloc_mem("V3D_grid_bottom_field", dst_decomp_info->get_decomp_id(), decomps_info_mgr->get_decomp_info(dst_decomp_info->get_decomp_id())->get_grid_id(), -dst_original_grid->get_grid_id(), DATA_TYPE_DOUBLE, "unitless", "Runtime_remapping_weights::allocate_intermediate_V3D_grid_bottom_field", false, false);
		else if (src_original_grid->get_original_CoR_grid()->does_use_V3D_level_coord()) {
			Remap_grid_class *full_original_grid = decomp_grids_mgr->search_decomp_grid_original_grid(dst_decomp_info->get_decomp_id(), get_parallel_remapping_weights()->get_dynamic_V1D_remap_weight_of_operator()->get_field_data_grid_src());
			if (full_original_grid == NULL)
				full_original_grid = get_parallel_remapping_weights()->get_dynamic_V1D_remap_weight_of_operator()->get_field_data_grid_src();
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, full_original_grid->get_num_dimensions() == 3 && full_original_grid->does_use_V3D_level_coord(), "Software error in Runtime_remapping_weights::allocate_intermediate_V3D_grid_bottom_field");
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->search_grid_info(dst_decomp_info->get_comp_id(), full_original_grid) != NULL, "Software error in Runtime_remapping_weights::allocate_intermediate_V3D_grid_bottom_field");
			int grid_id = original_grid_mgr->search_grid_info(dst_decomp_info->get_comp_id(), full_original_grid)->get_grid_id();
			intermediate_V3D_grid_bottom_field = memory_manager->alloc_mem(V3D_GRID_3D_LEVEL_FIELD_NAME, dst_decomp_info->get_decomp_id(), grid_id, -grid_id, DATA_TYPE_DOUBLE, "unitless", "Runtime_remapping_weights::allocate_intermediate_V3D_grid_bottom_field", false, false);
		}
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, true, "Software error in Runtime_remapping_weights::allocate_intermediate_V3D_grid_bottom_field");
    }
    return intermediate_V3D_grid_bottom_field;
}


bool Runtime_remapping_weights::match_requirements(const char *src_comp_full_name, const char *dst_comp_full_name, Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, Remapping_setting *remapping_setting, Decomp_info *dst_decomp_info)
{
    return words_are_the_same(this->src_comp_full_name,src_comp_full_name)&& words_are_the_same(this->dst_comp_full_name, dst_comp_full_name) && 
           this->src_original_grid == src_original_grid && this->dst_original_grid == dst_original_grid && 
           this->remapping_setting->is_the_same_as_another(remapping_setting) && this->dst_decomp_info == dst_decomp_info;
}


void Runtime_remapping_weights::generate_parallel_remapping_weights()
{
    Remap_grid_class **remap_related_decomp_grids;
	std::vector<std::pair<Remap_grid_class *, bool> > remap_related_grids;
    Remap_grid_class *decomp_original_grids[256];
    int i, j;


    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR,-1, sequential_remapping_weights != NULL, "C-Coupler software error remap weights is not found\n");
    cpl_check_remap_weights_format(sequential_remapping_weights);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR,-1, src_original_grid->get_H2D_sub_CoR_grid()->is_subset_of_grid(sequential_remapping_weights->get_data_grid_src()) && dst_original_grid->get_H2D_sub_CoR_grid()->is_subset_of_grid(sequential_remapping_weights->get_data_grid_dst()),
                     "Software error in Runtime_remapping_weights::generate_parallel_remapping_weights: grid inconsistency");

    EXECUTION_REPORT_LOG(REPORT_LOG, dst_decomp_info->get_comp_id(), true, "before generating remap_weights_src_decomp");
    src_decomp_info = decomps_info_mgr->generate_remap_weights_src_decomp(dst_decomp_info, src_original_grid, dst_original_grid, sequential_remapping_weights);
    EXECUTION_REPORT_LOG(REPORT_LOG, dst_decomp_info->get_comp_id(), true, "after generating remap_weights_src_decomp");
    EXECUTION_REPORT_LOG(REPORT_LOG, dst_decomp_info->get_comp_id(), true, "before generating parallel remap weights for runtime_remap_algorithm");

    decomp_original_grids[0] = src_original_grid->get_H2D_sub_CoR_grid();
    decomp_original_grids[1] = dst_original_grid->get_H2D_sub_CoR_grid();

	sequential_remapping_weights->get_remap_related_grids(remap_related_grids);
	remap_related_decomp_grids = new Remap_grid_class *[remap_related_grids.size()];
    for (i = 0; i < remap_related_grids.size(); i ++) {
        remap_related_decomp_grids[i] = remap_related_grids[i].first;
		if (!(remap_related_decomp_grids[i]->has_grid_coord_label(COORD_LABEL_LON) || remap_related_decomp_grids[i]->has_grid_coord_label(COORD_LABEL_LAT)))
			continue;
		if (!remap_related_grids[i].second) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, decomp_original_grids[0]->is_subset_of_grid(remap_related_grids[i].first),  "Software error in Runtime_remapping_weights::generate_parallel_remapping_weights");
            remap_related_decomp_grids[i] = decomp_grids_mgr->search_decomp_grid_info(src_decomp_info->get_decomp_id(), remap_related_grids[i].first, false)->get_decomp_grid();
		}
		else {
			if (remap_related_decomp_grids[i]->has_grid_coord_label(COORD_LABEL_LON) || remap_related_decomp_grids[i]->has_grid_coord_label(COORD_LABEL_LAT))
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, decomp_original_grids[1]->is_subset_of_grid(remap_related_grids[i].first), "Software error in Runtime_remapping_weights::generate_parallel_remapping_weights");
            remap_related_decomp_grids[i] = decomp_grids_mgr->search_decomp_grid_info(dst_decomp_info->get_decomp_id(), remap_related_grids[i].first, false)->get_decomp_grid();
		}
		if (remap_related_decomp_grids[i] == NULL)
			remap_related_decomp_grids[i] = remap_related_grids[i].first;
    }
 
    parallel_remapping_weights = sequential_remapping_weights->generate_parallel_remap_weights(remap_related_decomp_grids, decomp_original_grids, dst_decomp_info->get_decomp_id());
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, parallel_remapping_weights != NULL);
    dynamic_V1D_remap_weight_of_operator = parallel_remapping_weights->get_dynamic_V1D_remap_weight_of_operator();
    if (dynamic_V1D_remap_weight_of_operator != NULL) {
        runtime_V1D_remap_grid_src = dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->generate_remap_operator_runtime_grid(dynamic_V1D_remap_weight_of_operator->get_original_remap_operator()->get_src_grid(), dynamic_V1D_remap_weight_of_operator->get_original_remap_operator(), NULL);
        runtime_V1D_remap_grid_dst = dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->generate_remap_operator_runtime_grid(dynamic_V1D_remap_weight_of_operator->get_original_remap_operator()->get_dst_grid(), dynamic_V1D_remap_weight_of_operator->get_original_remap_operator(), NULL);
    }

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "find dynamic_V1D_remap_weight_of_operator corresponding to the grid \"%s\", with surface type %d", get_dst_original_grid()->get_grid_name(), get_dst_original_grid()->get_bottom_field_variation_type());

    if (dynamic_V1D_remap_weight_of_operator != NULL && get_dst_original_grid()->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_EXTERNAL && get_dst_original_grid()->get_bottom_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET) {
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "try to set level_V3D_coord_dynamic_trigger_field for the grid \"%s\" at %lx", dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_grid_name(), dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst());
        if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field() != NULL) {
            EXECUTION_REPORT(REPORT_ERROR, -1, dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field() == memory_manager->get_field_instance(get_dst_original_grid()->get_bottom_field_id())->get_field_data(), "Software error in Coupling_connection::add_bottom_field_coupling_info: the surface field of the same grid has been set to different data fields");
        }	
        else {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, memory_manager->check_is_legal_field_instance_id(get_dst_original_grid()->get_bottom_field_id()), "software error in Runtime_remapping_weights::generate_parallel_remapping_weights");
            printf("bottome_field:%s\n", decomps_info_mgr->get_decomp_info(memory_manager->get_field_instance(get_dst_original_grid()->get_bottom_field_id())->get_decomp_id())->get_decomp_name());
			EXECUTION_REPORT(REPORT_ERROR, -1, dst_decomp_info == decomps_info_mgr->get_decomp_info(memory_manager->get_field_instance(get_dst_original_grid()->get_bottom_field_id())->get_decomp_id()), "Encounter a special problem of parallel decomposition of data interpolation related to surface field of 3D grid: please contact Dr. Li Liu for support");
			dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->set_level_V3D_coord_dynamic_trigger_field(memory_manager->get_field_instance(get_dst_original_grid()->get_bottom_field_id())->get_field_data());
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Set level_V3D_coord_dynamic_trigger_field for the grid \"%s\" to surface field at decomp %s", dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_grid_name(), decomps_info_mgr->get_decomp_info(memory_manager->get_field_instance(get_dst_original_grid()->get_bottom_field_id())->get_decomp_id())->get_decomp_name());
        }
    }
	if (dynamic_V1D_remap_weight_of_operator != NULL && get_dst_original_grid()->get_V3D_lev_field_variation_type() != BOTTOM_FIELD_VARIATION_UNSET) {
        if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field() != NULL)
            EXECUTION_REPORT(REPORT_ERROR, -1, dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field() == memory_manager->get_field_instance(get_dst_original_grid()->get_V3D_lev_field_id())->get_field_data(), "Software error in Coupling_connection::add_bottom_field_coupling_info: the V3D field of the same grid has been set to different data fields");
        else dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->set_level_V3D_coord_dynamic_trigger_field(memory_manager->get_field_instance(get_dst_original_grid()->get_V3D_lev_field_id())->get_field_data());
	}
    
    EXECUTION_REPORT_LOG(REPORT_LOG, dst_decomp_info->get_comp_id(), true, "after generating parallel remap weights for runtime_remap_algorithm");

    delete [] remap_related_decomp_grids;
    remap_related_grids.clear();
}


void Runtime_remapping_weights::set_H2D_grids_area(const double *src_area, const double *dst_area, long src_grid_size, long dst_grid_size)
{
	if (src_area != NULL) {
	    src_H2D_grid_area = new double [src_grid_size];
    	memcpy(src_H2D_grid_area, src_area, src_grid_size*sizeof(double));
	    size_src_H2D_grid_area = src_grid_size;
	}
	if (dst_area != NULL) {
    	dst_H2D_grid_area = new double [dst_grid_size];
	    memcpy(dst_H2D_grid_area, dst_area, dst_grid_size*sizeof(double));
	    size_dst_H2D_grid_area = dst_grid_size;
	}
}


void Runtime_remapping_weights::renew_dynamic_V1D_remapping_weights()
{
    bool src_bottom_value_updated = false, dst_bottom_value_updated = false;
    bool src_bottom_value_specified = false, dst_bottom_value_specified = false;

    
    if (dynamic_V1D_remap_weight_of_operator == NULL)
        return;
	
    if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->is_sigma_grid() || dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->does_use_V3D_level_coord()) {
		if (intermediate_V3D_grid_bottom_field != NULL)
			intermediate_V3D_grid_bottom_field->transformation_between_chunks_array(true);
        src_bottom_value_specified = dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->is_level_V3D_coord_trigger_field_specified();
        src_bottom_value_updated = dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->is_level_V3D_coord_trigger_field_updated();
        if (src_original_grid->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_STATIC)
            EXECUTION_REPORT(REPORT_WARNING, src_original_grid->get_comp_id(), !src_bottom_value_updated || !src_bottom_value_specified, "the surface field of the 3-D grid \"%s\" (registered in the component \"%s\") is updated while the surface field has been specified as a static one. Please check.", src_original_grid->get_grid_name(), src_original_grid->get_comp_full_name());
        if (src_original_grid->get_V3D_lev_field_variation_type() == BOTTOM_FIELD_VARIATION_STATIC)
            EXECUTION_REPORT(REPORT_WARNING, src_original_grid->get_comp_id(), !src_bottom_value_updated || !src_bottom_value_specified, "the 3-D level field of the 3-D grid \"%s\" (registered in the component \"%s\") is updated while the 3-D level field has been specified as a static one. Please check.", src_original_grid->get_grid_name(), src_original_grid->get_comp_full_name());
    }
    if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->is_sigma_grid() || dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->does_use_V3D_level_coord()) {
		if (dst_original_grid->get_bottom_field_id() != -1)
			memory_manager->get_field_instance(dst_original_grid->get_bottom_field_id())->transformation_between_chunks_array(true);
		if (dst_original_grid->get_V3D_lev_field_id() != -1)
			memory_manager->get_field_instance(dst_original_grid->get_V3D_lev_field_id())->transformation_between_chunks_array(true);
        dst_bottom_value_specified = dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->is_level_V3D_coord_trigger_field_specified();
        dst_bottom_value_updated = dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->is_level_V3D_coord_trigger_field_updated();
        if (dst_original_grid->get_bottom_field_variation_type() == BOTTOM_FIELD_VARIATION_STATIC)
            EXECUTION_REPORT(REPORT_ERROR, dst_original_grid->get_comp_id(), !dst_bottom_value_updated || !dst_bottom_value_specified, "the surface field of the 3-D grid \"%s\" is updated while the surface field has been specified as a static one. Please verify", dst_original_grid->get_grid_name());
        if (dst_original_grid->get_V3D_lev_field_variation_type() == BOTTOM_FIELD_VARIATION_STATIC)
            EXECUTION_REPORT(REPORT_ERROR, dst_original_grid->get_comp_id(), !dst_bottom_value_updated || !dst_bottom_value_specified, "the 3-D level field of the 3-D grid \"%s\" (registered in the component \"%s\") is updated while the 3-D level field has been specified as a static one. Please verify", dst_original_grid->get_grid_name(), dst_original_grid->get_comp_full_name());
    }
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "src_bottom_value_updated is %lx, dst_bottom_value_updated is %lx", dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->get_level_V3D_coord_dynamic_trigger_field(), dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->get_level_V3D_coord_dynamic_trigger_field());
	comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "vertical coord update");
    if (src_bottom_value_updated)
		if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->is_sigma_grid())
	        dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->calculate_lev_sigma_values();
		else dynamic_V1D_remap_weight_of_operator->get_field_data_grid_src()->update_grid_center_3D_level_field_from_external();
    if (dst_bottom_value_updated)
		if (dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->is_sigma_grid())
	        dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->calculate_lev_sigma_values();
		else dynamic_V1D_remap_weight_of_operator->get_field_data_grid_dst()->update_grid_center_3D_level_field_from_external();
	comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "vertical coord update");

	comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "dyn v1d wgt");
    if (src_bottom_value_updated || dst_bottom_value_updated)
        dynamic_V1D_remap_weight_of_operator->renew_vertical_remap_weights(runtime_V1D_remap_grid_src, runtime_V1D_remap_grid_dst);
	comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "dyn v1d wgt");
}


Runtime_remapping_weights_container::Runtime_remapping_weights_container(const char *src_comp_full_name, const char *dst_comp_full_name, Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, Remapping_setting *remapping_setting, Decomp_info *dst_decomp_info, H2D_remapping_wgt_file_info *H2D_remapping_weight_file)
{
    this->src_comp_full_name = strdup(src_comp_full_name);
    this->dst_comp_full_name = strdup(dst_comp_full_name);
    this->src_original_grid = src_original_grid;
    this->dst_original_grid = dst_original_grid;
    this->dst_decomp_info = dst_decomp_info;
    this->src_decomp_info = NULL;
	runtime_remapping_weights_under_V3D = NULL;
	runtime_remapping_weights_Time1D = NULL;


	EXECUTION_REPORT_LOG(REPORT_LOG, dst_original_grid->get_comp_id(), true, "Generate remapping weights from grid \"%s\" to grid \"%s\" ", src_original_grid->get_grid_name(), dst_original_grid->get_grid_name());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(src_original_grid->get_grid_id()), "Error happens when generating the remapping weights from the grid \"%s\" of the component model \"%s\" to the grid \"%s\" of the component model \"%s\": these two grids do not have the same dimensions. Please verify.", src_original_grid->get_grid_name(), src_comp_full_name, dst_original_grid->get_grid_name(), dst_comp_full_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, dst_original_grid->get_comp_id(), src_original_grid->get_original_CoR_grid()->have_the_same_dimensions_with(dst_original_grid->get_original_CoR_grid()), "Error happens when generating the remapping weights from the grid \"%s\" of the component model \"%s\" to the grid \"%s\" of the component model \"%s\": these two grids do not have the same dimensions. Please verify.", src_original_grid->get_grid_name(), src_comp_full_name, dst_original_grid->get_grid_name(), dst_comp_full_name);
	if (dst_original_grid->get_Tracer1D_sub_grid() != NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid->get_Tracer1D_sub_grid() != NULL, "Software error in Runtime_remapping_weights_container::Runtime_remapping_weights_container");		
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, dst_original_grid->get_comp_id(), src_original_grid->get_Tracer1D_sub_grid()->get_original_CoR_grid()->get_grid_size() == dst_original_grid->get_Tracer1D_sub_grid()->get_original_CoR_grid()->get_grid_size(), "Error happens when generating the remapping weights from the grid \"%s\" of the component model \"%s\" to the grid \"%s\" of the component model \"%s\": these two grids have tracer sub grids but with different number of tracers. Please verify.", src_original_grid->get_grid_name(), src_comp_full_name, dst_original_grid->get_grid_name(), dst_comp_full_name);
	}	
	
	if (dst_original_grid->get_max_sub_grid_under_V3D() != NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid->get_max_sub_grid_under_V3D() != NULL, "Software error in Runtime_remapping_weights_container::Runtime_remapping_weights_container");
		EXECUTION_REPORT_LOG(REPORT_LOG, dst_original_grid->get_comp_id(), true, "Generate remapping weights under 3-D from grid \"%s\" to grid \"%s\" ", src_original_grid->get_max_sub_grid_under_V3D()->get_grid_name(), dst_original_grid->get_max_sub_grid_under_V3D()->get_grid_name());
		runtime_remapping_weights_under_V3D = runtime_remapping_weights_mgr->search_or_generate_runtime_remapping_weights(src_comp_full_name, dst_comp_full_name, src_original_grid->get_max_sub_grid_under_V3D(), dst_original_grid->get_max_sub_grid_under_V3D(), remapping_setting, dst_decomp_info, H2D_remapping_weight_file);
		src_decomp_info = runtime_remapping_weights_under_V3D->get_src_decomp_info();
	}

	if (dst_original_grid->get_Time1D_sub_grid() != NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid->get_Time1D_sub_grid() != NULL, "Software error in Runtime_remapping_weights_container::Runtime_remapping_weights_container");
		EXECUTION_REPORT_LOG(REPORT_LOG, dst_original_grid->get_comp_id(), true, "Generate remapping weights of Time 1D from grid \"%s\" to grid \"%s\" ", src_original_grid->get_Time1D_sub_grid()->get_grid_name(), dst_original_grid->get_Time1D_sub_grid()->get_grid_name());
		runtime_remapping_weights_Time1D = runtime_remapping_weights_mgr->search_or_generate_runtime_remapping_weights(src_comp_full_name, dst_comp_full_name, src_original_grid->get_Time1D_sub_grid(), dst_original_grid->get_Time1D_sub_grid(), remapping_setting, dst_decomp_info, NULL);
	}	
	comp_id = dst_original_grid->get_comp_id();
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remapping_weights_Time1D != NULL || runtime_remapping_weights_under_V3D != NULL, "Software error in Runtime_remapping_weights_container::Runtime_remapping_weights_container");
	temp_remap_data_after_remapping_under_V3D = NULL;
	interchange_grid_before_remapping_under_V3D = NULL;
	interchange_grid_before_remapping_Time1D = NULL;
	interchange_grid_after_remapping_under_V3D = NULL;
	interchange_grid_after_remapping_Time1D = NULL;
	partial_data_field_for_src_remapping_under_V3D = NULL;
	partial_data_field_for_dst_remapping_under_V3D = NULL;
	partial_data_field_for_src_remapping_Time1D = NULL;
	partial_data_field_for_dst_remapping_Time1D = NULL;
}


Runtime_remapping_weights_container::~Runtime_remapping_weights_container()
{
	if (temp_remap_data_after_remapping_under_V3D != NULL)
		delete temp_remap_data_after_remapping_under_V3D;
	if (partial_data_field_for_src_remapping_under_V3D != NULL) {
		partial_data_field_for_src_remapping_under_V3D->get_grid_data_field()->data_buf = NULL;
		partial_data_field_for_dst_remapping_under_V3D->get_grid_data_field()->data_buf = NULL;
		delete partial_data_field_for_src_remapping_under_V3D;
		delete partial_data_field_for_dst_remapping_under_V3D;
	}
	if (partial_data_field_for_src_remapping_Time1D != NULL) {
		partial_data_field_for_src_remapping_Time1D->get_grid_data_field()->data_buf = NULL;
		partial_data_field_for_dst_remapping_Time1D->get_grid_data_field()->data_buf = NULL;
		delete partial_data_field_for_src_remapping_Time1D;
		delete partial_data_field_for_dst_remapping_Time1D;
	}
}


bool Runtime_remapping_weights_container::is_empty()
{
	if (runtime_remapping_weights_under_V3D != NULL && runtime_remapping_weights_under_V3D->get_parallel_remapping_weights() != NULL)
		return false;

	return runtime_remapping_weights_Time1D == NULL;
}


void Runtime_remapping_weights_container::do_remap(Remap_grid_data_class *src_field_data, Remap_grid_data_class *dst_field_data)
{
	Remap_grid_data_class *current_src_field_data, *current_dst_field_data;
	Remap_grid_data_class *current_partial_data_field_for_src_remapping_under_V3D, *current_partial_data_field_for_dst_remapping_under_V3D;
	Remap_grid_data_class *current_partial_data_field_for_src_remapping_Time1D, *current_partial_data_field_for_dst_remapping_Time1D;


	if (runtime_remapping_weights_under_V3D != NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remapping_weights_under_V3D->get_parallel_remapping_weights() != NULL, "Software error in Runtime_remapping_weights_container::do_remap");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_src()->is_subset_of_grid(src_field_data->get_coord_value_grid()), "Software error in Runtime_remapping_weights_container::do_remap");
		if (!words_are_the_same(src_field_data->get_grid_data_field()->field_name_in_application, V3D_GRID_3D_LEVEL_FIELD_NAME))
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_dst()->is_subset_of_grid(dst_field_data->get_coord_value_grid()), "Software error in Runtime_remapping_weights_container::do_remap");
		if (src_field_data->get_coord_value_grid()->get_num_dimensions() > runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_src()->get_num_dimensions() && interchange_grid_before_remapping_under_V3D == NULL) {
			src_field_data->get_coord_value_grid()->generate_interchange_grids(runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_src(), &interchange_grid_before_remapping_under_V3D, NULL, 0);
			interchange_grid_after_remapping_under_V3D = new Remap_grid_class(interchange_grid_before_remapping_under_V3D, runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_src(), runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_dst(), true);
			remap_grid_manager->add_temp_grid(interchange_grid_after_remapping_under_V3D);
			if (runtime_remapping_weights_Time1D != NULL) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, temp_remap_data_after_remapping_under_V3D == NULL, "Software error in Runtime_remapping_weights_container::do_remap");
				temp_remap_data_after_remapping_under_V3D = src_field_data->duplicate_grid_data_field(interchange_grid_after_remapping_under_V3D, 1, false, false);
			}	
			partial_data_field_for_src_remapping_under_V3D = src_field_data->duplicate_grid_data_field(runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_src(), 1, false, false);
			partial_data_field_for_dst_remapping_under_V3D = dst_field_data->duplicate_grid_data_field(runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_dst(), 1, false, false);
			delete [] partial_data_field_for_src_remapping_under_V3D->get_grid_data_field()->data_buf;
			delete [] partial_data_field_for_dst_remapping_under_V3D->get_grid_data_field()->data_buf;
		}
	}

	if (runtime_remapping_weights_Time1D != NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src()->is_subset_of_grid(src_field_data->get_coord_value_grid()), "Software error in Runtime_remapping_weights_container::do_remap");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_dst()->is_subset_of_grid(dst_field_data->get_coord_value_grid()), "Software error in Runtime_remapping_weights_container::do_remap");
		if (src_field_data->get_coord_value_grid()->get_num_dimensions() > runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src()->get_num_dimensions() && interchange_grid_before_remapping_Time1D == NULL) {
			if (interchange_grid_after_remapping_under_V3D == NULL)
				src_field_data->get_coord_value_grid()->generate_interchange_grids(runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src(), &interchange_grid_before_remapping_Time1D, NULL, 0);
			else interchange_grid_after_remapping_under_V3D->generate_interchange_grids(runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src(), &interchange_grid_before_remapping_Time1D, NULL, 0);
			interchange_grid_after_remapping_Time1D = new Remap_grid_class(interchange_grid_before_remapping_Time1D, runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src(), runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_dst(), true);
			remap_grid_manager->add_temp_grid(interchange_grid_before_remapping_Time1D);
			remap_grid_manager->add_temp_grid(interchange_grid_after_remapping_Time1D);
			partial_data_field_for_src_remapping_Time1D = src_field_data->duplicate_grid_data_field(runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src(), 1, false, false);
			partial_data_field_for_dst_remapping_Time1D = dst_field_data->duplicate_grid_data_field(runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_dst(), 1, false, false);
			delete [] partial_data_field_for_src_remapping_Time1D->get_grid_data_field()->data_buf;
			delete [] partial_data_field_for_dst_remapping_Time1D->get_grid_data_field()->data_buf;
		}
	}

	if (runtime_remapping_weights_under_V3D != NULL) {
		current_src_field_data = src_field_data;
		current_dst_field_data = temp_remap_data_after_remapping_under_V3D != NULL? temp_remap_data_after_remapping_under_V3D : dst_field_data;
		current_partial_data_field_for_src_remapping_under_V3D = partial_data_field_for_src_remapping_under_V3D != NULL? partial_data_field_for_src_remapping_under_V3D : current_src_field_data;
		current_partial_data_field_for_dst_remapping_under_V3D = partial_data_field_for_dst_remapping_under_V3D != NULL? partial_data_field_for_dst_remapping_under_V3D : current_dst_field_data;
		if (interchange_grid_before_remapping_under_V3D != NULL) {
			current_src_field_data->interchange_grid_data(interchange_grid_before_remapping_under_V3D);
			current_dst_field_data->interchange_grid_data(interchange_grid_after_remapping_under_V3D);
		}
		long src_size_sub_grid_under_V3D = runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_src()->get_grid_size();
		long dst_size_sub_grid_under_V3D = runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->get_data_grid_dst()->get_grid_size();
		int num_iterations = current_src_field_data->get_grid_data_field()->required_data_size / src_size_sub_grid_under_V3D;
		if (!words_are_the_same(src_field_data->get_grid_data_field()->field_name_in_application, V3D_GRID_3D_LEVEL_FIELD_NAME))
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (current_src_field_data->get_grid_data_field()->required_data_size % src_size_sub_grid_under_V3D) == 0 && (current_dst_field_data->get_grid_data_field()->required_data_size % dst_size_sub_grid_under_V3D) == 0 && num_iterations == current_dst_field_data->get_grid_data_field()->required_data_size / dst_size_sub_grid_under_V3D, "Software error in Runtime_remapping_weights_container::do_remap");
		if (!words_are_the_same(src_field_data->get_grid_data_field()->field_name_in_application,V3D_GRID_3D_LEVEL_FIELD_NAME))
			runtime_remapping_weights_under_V3D->renew_dynamic_V1D_remapping_weights();
		comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "remapping cal");
		for (int i = 0; i < num_iterations; i ++) {
			current_partial_data_field_for_src_remapping_under_V3D->get_grid_data_field()->data_buf = (double*)(current_src_field_data->get_grid_data_field()->data_buf) + i * src_size_sub_grid_under_V3D;
			current_partial_data_field_for_dst_remapping_under_V3D->get_grid_data_field()->data_buf = (double*)(current_dst_field_data->get_grid_data_field()->data_buf) + i * dst_size_sub_grid_under_V3D;
			runtime_remapping_weights_under_V3D->get_parallel_remapping_weights()->do_remap(dst_original_grid->get_comp_id(), current_partial_data_field_for_src_remapping_under_V3D, current_partial_data_field_for_dst_remapping_under_V3D);
		}
		comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "remapping cal");	
	}

	if (runtime_remapping_weights_Time1D != NULL) {
		current_src_field_data = temp_remap_data_after_remapping_under_V3D != NULL? temp_remap_data_after_remapping_under_V3D : src_field_data;
		current_dst_field_data = dst_field_data;
		current_partial_data_field_for_src_remapping_Time1D = partial_data_field_for_src_remapping_Time1D != NULL? partial_data_field_for_src_remapping_Time1D : current_src_field_data;
		current_partial_data_field_for_dst_remapping_Time1D = partial_data_field_for_dst_remapping_Time1D != NULL? partial_data_field_for_dst_remapping_Time1D : current_dst_field_data;
		if (interchange_grid_before_remapping_Time1D != NULL) {
			current_src_field_data->interchange_grid_data(interchange_grid_before_remapping_Time1D);
			current_dst_field_data->interchange_grid_data(interchange_grid_after_remapping_Time1D);
		}	
		long src_size_sub_grid_Time1D = runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_src()->get_grid_size();
		long dst_size_sub_grid_Time1D = runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->get_data_grid_dst()->get_grid_size();
		int num_iterations = current_src_field_data->get_grid_data_field()->required_data_size / src_size_sub_grid_Time1D;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (current_src_field_data->get_grid_data_field()->required_data_size % src_size_sub_grid_Time1D) == 0 && (current_dst_field_data->get_grid_data_field()->required_data_size % dst_size_sub_grid_Time1D) == 0 && num_iterations == current_dst_field_data->get_grid_data_field()->required_data_size / dst_size_sub_grid_Time1D, "Software error in Runtime_remapping_weights_container::do_remap");
		comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "remapping cal");
		for (int i = 0; i < num_iterations; i ++) {
			current_partial_data_field_for_src_remapping_Time1D->get_grid_data_field()->data_buf = (double*)(current_src_field_data->get_grid_data_field()->data_buf) + i * src_size_sub_grid_Time1D;
			current_partial_data_field_for_dst_remapping_Time1D->get_grid_data_field()->data_buf = (double*)(current_dst_field_data->get_grid_data_field()->data_buf) + i * dst_size_sub_grid_Time1D;
			runtime_remapping_weights_Time1D->get_parallel_remapping_weights()->do_remap(dst_original_grid->get_comp_id(), current_partial_data_field_for_src_remapping_Time1D, current_partial_data_field_for_dst_remapping_Time1D);
		}
		comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_original_grid->get_comp_id(),false,"")->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "remapping cal");	
	}

	src_field_data->interchange_grid_data(src_field_data->get_coord_value_grid());
	dst_field_data->interchange_grid_data(dst_field_data->get_coord_value_grid());
}


Runtime_remapping_weights_mgt::~Runtime_remapping_weights_mgt()
{
    for (int i = 0; i < runtime_remapping_weights.size(); i ++)
        delete runtime_remapping_weights[i];
}


Runtime_remapping_weights *Runtime_remapping_weights_mgt::search_or_generate_runtime_remapping_weights(const char *src_comp_full_name, const char *dst_comp_full_name, Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, Remapping_setting *remapping_setting, Decomp_info *dst_decomp_info, H2D_remapping_wgt_file_info *H2D_remapping_weight_file)
{
    remapping_setting->shrink(src_original_grid, dst_original_grid);
    
    for (int i = 0; i < runtime_remapping_weights.size(); i ++)
        if (runtime_remapping_weights[i]->match_requirements(src_comp_full_name, dst_comp_full_name, src_original_grid, dst_original_grid, remapping_setting, dst_decomp_info))
            return runtime_remapping_weights[i];

    runtime_remapping_weights.push_back(new Runtime_remapping_weights(src_comp_full_name, dst_comp_full_name, src_original_grid, dst_original_grid, remapping_setting, dst_decomp_info, H2D_remapping_weight_file));
    return runtime_remapping_weights[runtime_remapping_weights.size()-1];
}


void Runtime_remapping_weights_mgt::transfer_runtime_remapping_weights(Runtime_remapping_weights *remapping_weights_from, Runtime_remapping_weights **remapping_weights_to, Comp_comm_group_mgt_node *comp_node_from, Comp_comm_group_mgt_node *comp_node_to)
{
    double *temp_src_H2D_grid_area = NULL, *temp_dst_H2D_grid_area = NULL;
    long temp_src_H2D_grid_size = 0, temp_dst_H2D_grid_size = 0;


    if (comp_node_from->get_current_proc_local_id() != -1) {
        if (remapping_weights_from->get_src_H2D_grid_area() != NULL) {
            temp_src_H2D_grid_size = remapping_weights_from->get_src_original_grid()->get_original_CoR_grid()->get_grid_size()*sizeof(double);
            temp_src_H2D_grid_area = new double [remapping_weights_from->get_src_original_grid()->get_original_CoR_grid()->get_grid_size()];
            memcpy(temp_src_H2D_grid_area, remapping_weights_from->get_src_H2D_grid_area(), temp_src_H2D_grid_size);
        }
        if (remapping_weights_from->get_dst_H2D_grid_area() != NULL) {
            temp_dst_H2D_grid_size = remapping_weights_from->get_dst_original_grid()->get_original_CoR_grid()->get_grid_size()*sizeof(double);
            temp_dst_H2D_grid_area = new double [remapping_weights_from->get_dst_original_grid()->get_original_CoR_grid()->get_grid_size()];
            memcpy(temp_dst_H2D_grid_area, remapping_weights_from->get_dst_H2D_grid_area(), temp_dst_H2D_grid_size);
        }
    }
    transfer_array_from_one_comp_to_another(comp_node_from->get_current_proc_local_id(), comp_node_from->get_local_proc_global_id(0), comp_node_to->get_current_proc_local_id(), comp_node_to->get_local_proc_global_id(0), comp_node_to->get_comm_group(), (char**)(&temp_src_H2D_grid_area), temp_src_H2D_grid_size);
    transfer_array_from_one_comp_to_another(comp_node_from->get_current_proc_local_id(), comp_node_from->get_local_proc_global_id(0), comp_node_to->get_current_proc_local_id(), comp_node_to->get_local_proc_global_id(0), comp_node_to->get_comm_group(), (char**)(&temp_dst_H2D_grid_area), temp_dst_H2D_grid_size);
    
    if (comp_node_to->get_current_proc_local_id() != -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, *remapping_weights_to == NULL, "Software error in Runtime_remapping_weights_mgt::transfer_runtime_remapping_weights");
        *remapping_weights_to = new Runtime_remapping_weights();
        if (temp_src_H2D_grid_size != 0)
            (*remapping_weights_to)->set_H2D_grids_area(temp_src_H2D_grid_area, temp_dst_H2D_grid_area, temp_src_H2D_grid_size/sizeof(double), temp_dst_H2D_grid_size/sizeof(double));
        runtime_remapping_weights.push_back(*remapping_weights_to);
    }

    if (temp_src_H2D_grid_area != NULL) {
        delete [] temp_src_H2D_grid_area;
        delete [] temp_dst_H2D_grid_area;
    }
}


