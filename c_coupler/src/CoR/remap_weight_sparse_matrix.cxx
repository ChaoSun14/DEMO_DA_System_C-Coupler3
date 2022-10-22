/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file is initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "cor_global_data.h"
#include "global_data.h"
#include "remap_weight_sparse_matrix.h"
#include "remap_operator_basis.h"
#include <stdio.h>
#include <string.h>
#include <math.h>


Remap_weight_sparse_matrix::Remap_weight_sparse_matrix(Remap_operator_basis *remap_operator, 
                                                       long num_weights, long *cells_indexes_src, long *cells_indexes_dst, double *weight_values, 
                                                       long num_remaped_dst_cells_indexes, long *remaped_dst_cells_indexes)
{
    this->remap_operator = remap_operator;
    this->num_weights = num_weights;
    this->weight_arrays_size = num_weights;
    this->cells_indexes_src = cells_indexes_src;
    this->cells_indexes_dst = cells_indexes_dst;
    this->weight_values = weight_values;
    this->remaped_dst_cells_indexes = remaped_dst_cells_indexes;
    this->num_remaped_dst_cells_indexes = num_remaped_dst_cells_indexes;
    
    if (remaped_dst_cells_indexes == NULL) {
        int *mask = new int [remap_operator->get_dst_grid()->get_grid_size()];
        for (int i = 0; i < remap_operator->get_dst_grid()->get_grid_size(); i ++)
            mask[i] = 0;
        for (int i = 0; i < num_weights; i ++)
            mask[cells_indexes_dst[i]] = 1;
        this->num_remaped_dst_cells_indexes = 0;
        for (int i = 0; i < remap_operator->get_dst_grid()->get_grid_size(); i ++)
            if (mask[i] == 1)
                this->num_remaped_dst_cells_indexes ++;
        this->remaped_dst_cells_indexes = new long [this->num_remaped_dst_cells_indexes];
        this->num_remaped_dst_cells_indexes = 0;
        for (int i = 0; i < remap_operator->get_dst_grid()->get_grid_size(); i ++)
            if (mask[i] == 1)
                this->remaped_dst_cells_indexes[this->num_remaped_dst_cells_indexes ++] = i;
        delete [] mask;
    }

    this->remaped_dst_cells_indexes_array_size = this->num_remaped_dst_cells_indexes;
}


Remap_weight_sparse_matrix::Remap_weight_sparse_matrix(Remap_operator_basis *remap_operator)
{
    this->remap_operator = remap_operator;
    num_weights = 0;
    weight_arrays_size = 1;
    num_remaped_dst_cells_indexes = 0;
    remaped_dst_cells_indexes_array_size = weight_arrays_size;
    cells_indexes_src = new long [weight_arrays_size];
    cells_indexes_dst = new long [weight_arrays_size];
    weight_values = new double [weight_arrays_size];
    remaped_dst_cells_indexes = new long [remaped_dst_cells_indexes_array_size];
}


Remap_weight_sparse_matrix::~Remap_weight_sparse_matrix()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cells_indexes_src != NULL && cells_indexes_dst != NULL && weight_values != NULL, "Software error in ~Remap_weight_sparse_matrix");

    delete [] cells_indexes_src;
    delete [] cells_indexes_dst;
    if (remaped_dst_cells_indexes != NULL)
        delete [] remaped_dst_cells_indexes;
    delete [] weight_values;
}


void Remap_weight_sparse_matrix::clear_weights_info()
{
    num_weights = 0; 
    num_remaped_dst_cells_indexes = 0;
}


void Remap_weight_sparse_matrix::add_weights(long *indexes_src, long index_dst, double *added_weight_values, int num_added_weights, bool is_real_weight)
{
    long *new_indexes_src_grid, *new_indexes_dst_grid, *new_remaped_dst_cells_indexes;
    double *new_weight_values;
    long i, new_array_size;
	const int *local_cells_global_indexes_src, *local_cells_global_indexes_dst;


    if (is_real_weight) {
        for (i = 0; i < num_added_weights; i ++)
            EXECUTION_REPORT(REPORT_ERROR, -1, indexes_src[i] >= 0 && indexes_src[i] < remap_operator->get_src_grid()->get_grid_size(), "C-Coupler error1 in add_weights of Remap_weight_sparse_matrix");
        EXECUTION_REPORT(REPORT_ERROR, -1, remap_operator->get_dst_grid() != NULL && index_dst >= 0 && index_dst < remap_operator->get_dst_grid()->get_grid_size(), "C-Coupler error2 in add_weights of Remap_weight_sparse_matrix");
    }
    
    if (num_weights + num_added_weights > weight_arrays_size) {
        new_array_size = 2 * (num_weights+num_added_weights);
        new_indexes_src_grid = new long [new_array_size];
        new_indexes_dst_grid = new long [new_array_size];
        new_weight_values = new double [new_array_size];
        for (i = 0; i < weight_arrays_size; i ++) {
            new_indexes_src_grid[i] = cells_indexes_src[i];
            new_indexes_dst_grid[i] = cells_indexes_dst[i];
            new_weight_values[i] = weight_values[i];
        }
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cells_indexes_src != NULL && cells_indexes_dst != NULL && weight_values != NULL, "Software error in Remap_weight_sparse_matrix");
        delete [] cells_indexes_src;
        delete [] cells_indexes_dst;
        delete [] weight_values;
        cells_indexes_src = new_indexes_src_grid;
        cells_indexes_dst = new_indexes_dst_grid;
        weight_values = new_weight_values;
        weight_arrays_size = new_array_size;
    }

    if (num_remaped_dst_cells_indexes >= remaped_dst_cells_indexes_array_size) {
        new_remaped_dst_cells_indexes = new long [remaped_dst_cells_indexes_array_size*2];
        for (i = 0; i < num_remaped_dst_cells_indexes; i ++)
            new_remaped_dst_cells_indexes[i] = remaped_dst_cells_indexes[i];
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remaped_dst_cells_indexes != NULL && remaped_dst_cells_indexes_array_size > 0, "Software error in Remap_weight_sparse_matrix");
        delete [] remaped_dst_cells_indexes;
        remaped_dst_cells_indexes = new_remaped_dst_cells_indexes;
        remaped_dst_cells_indexes_array_size *= 2;
    }

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remap_operator != NULL && remap_operator->get_src_grid() != NULL && remap_operator->get_dst_grid() != NULL , "Software error in Remap_weight_sparse_matrix::add_weights");
	local_cells_global_indexes_src = remap_operator->get_src_grid()->get_local_cell_global_indexes();
	local_cells_global_indexes_dst = remap_operator->get_dst_grid()->get_local_cell_global_indexes();
	index_dst = local_cells_global_indexes_dst != NULL? local_cells_global_indexes_dst[index_dst] : index_dst;
    for (i = 0; i < num_added_weights; i ++) 
		indexes_src[i] = local_cells_global_indexes_src != NULL? local_cells_global_indexes_src[indexes_src[i]] : indexes_src[i];

	do_quick_sort(indexes_src, added_weight_values, 0, num_added_weights-1);
    for (i = 0; i < num_added_weights; i ++) {
        cells_indexes_src[num_weights] = indexes_src[i];
        cells_indexes_dst[num_weights] = index_dst;
        weight_values[num_weights] = added_weight_values[i];
        num_weights ++;
    }

	if (num_remaped_dst_cells_indexes == 0 || remaped_dst_cells_indexes[num_remaped_dst_cells_indexes-1] != index_dst)
	    remaped_dst_cells_indexes[num_remaped_dst_cells_indexes++] = index_dst;
}


void Remap_weight_sparse_matrix::get_weight(long *index_src, long *index_dst, double *weight_value, int index_weight)
{
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, index_weight >= 0 && index_weight < num_weights, "software error when get remapping weight of sparse matrix\n");
    *index_src = cells_indexes_src[index_weight];
    *index_dst = cells_indexes_dst[index_weight];
    *weight_value = weight_values[index_weight];
}


void Remap_weight_sparse_matrix::remap_values(double *data_values_src, double *data_values_dst, int dst_array_size)
{
    for (long i = 0; i < num_remaped_dst_cells_indexes; i ++)
        data_values_dst[remaped_dst_cells_indexes[i]] = 0.0;

    for (long i = 0; i < num_weights; i ++)
        data_values_dst[cells_indexes_dst[i]] += data_values_src[cells_indexes_src[i]] * weight_values[i];
}


void Remap_weight_sparse_matrix::calc_src_decomp(long *decomp_map_src, const long *decomp_map_dst)
{
    for (long i = 0; i < num_weights; i ++)
        decomp_map_src[cells_indexes_src[i]] = (decomp_map_src[cells_indexes_src[i]] | decomp_map_dst[cells_indexes_dst[i]]);
}


Remap_weight_sparse_matrix *Remap_weight_sparse_matrix::duplicate_remap_weight_of_sparse_matrix()
{
    Remap_weight_sparse_matrix *duplicated_remap_weight_of_sparse_matrix;


    duplicated_remap_weight_of_sparse_matrix = new Remap_weight_sparse_matrix(remap_operator);
    duplicated_remap_weight_of_sparse_matrix->weight_arrays_size = this->weight_arrays_size;
    duplicated_remap_weight_of_sparse_matrix->num_weights = this->num_weights;
    duplicated_remap_weight_of_sparse_matrix->remaped_dst_cells_indexes_array_size = this->remaped_dst_cells_indexes_array_size;
    duplicated_remap_weight_of_sparse_matrix->num_remaped_dst_cells_indexes = this->num_remaped_dst_cells_indexes;
    delete [] duplicated_remap_weight_of_sparse_matrix->cells_indexes_src;
    delete [] duplicated_remap_weight_of_sparse_matrix->cells_indexes_dst;
    delete [] duplicated_remap_weight_of_sparse_matrix->weight_values;
    delete [] duplicated_remap_weight_of_sparse_matrix->remaped_dst_cells_indexes;
    duplicated_remap_weight_of_sparse_matrix->cells_indexes_src = new long [this->weight_arrays_size];
    duplicated_remap_weight_of_sparse_matrix->cells_indexes_dst = new long [this->weight_arrays_size];
    duplicated_remap_weight_of_sparse_matrix->weight_values = new double [this->weight_arrays_size];
    duplicated_remap_weight_of_sparse_matrix->remaped_dst_cells_indexes = new long [this->remaped_dst_cells_indexes_array_size];    
    memcpy(duplicated_remap_weight_of_sparse_matrix->cells_indexes_src, this->cells_indexes_src, this->weight_arrays_size*sizeof(long));
    memcpy(duplicated_remap_weight_of_sparse_matrix->cells_indexes_dst, this->cells_indexes_dst, this->weight_arrays_size*sizeof(long));
    memcpy(duplicated_remap_weight_of_sparse_matrix->weight_values, this->weight_values, this->weight_arrays_size*sizeof(double));
    memcpy(duplicated_remap_weight_of_sparse_matrix->remaped_dst_cells_indexes, this->remaped_dst_cells_indexes, this->remaped_dst_cells_indexes_array_size*sizeof(long));

    return duplicated_remap_weight_of_sparse_matrix;
}


void Remap_weight_sparse_matrix::compare_to_another_sparse_matrix(Remap_weight_sparse_matrix *another_sparse_matrix)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_weights == another_sparse_matrix->num_weights, "C-Coupler error1 in compare_to_another_sparse_matrix");
    EXECUTION_REPORT(REPORT_ERROR, -1, this->num_remaped_dst_cells_indexes == another_sparse_matrix->num_remaped_dst_cells_indexes, "C-Coupler error2 in compare_to_another_sparse_matrix");

    for (long i = 0; i < num_weights; i ++) {
        EXECUTION_REPORT(REPORT_ERROR, -1, this->cells_indexes_src[i] == another_sparse_matrix->cells_indexes_src[i], "C-Coupler error3 in compare_to_another_sparse_matrix");
        EXECUTION_REPORT(REPORT_ERROR, -1, this->cells_indexes_dst[i] == another_sparse_matrix->cells_indexes_dst[i], "C-Coupler error4 in compare_to_another_sparse_matrix");
        EXECUTION_REPORT(REPORT_ERROR, -1, this->weight_values[i] == another_sparse_matrix->weight_values[i], "C-Coupler error5 in compare_to_another_sparse_matrix");        
    }
    
    for (long i = 0; i < num_remaped_dst_cells_indexes; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, this->remaped_dst_cells_indexes[i] == another_sparse_matrix->remaped_dst_cells_indexes[i], "C-Coupler error6 in compare_to_another_sparse_matrix");
}


void Remap_weight_sparse_matrix::print()
{
    for (int i = 0; i < num_weights; i ++)
        printf("remapping weight (%d): src_index=%d, dst_index=%d, weight_value=%lf\n", i, cells_indexes_src[i], cells_indexes_dst[i], weight_values[i]);
}


Remap_weight_sparse_matrix *Remap_weight_sparse_matrix::gather(int comp_id)
{
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);
	long *overall_cells_indexes_src = NULL;
	long *overall_cells_indexes_dst = NULL;
	double *overall_wgt_values = NULL;
	long num_overall_wgts, true_num_overall_wgts, i, j, offset;
	int *all_array_size = new int [comp_node->get_num_procs()];
	Remap_weight_sparse_matrix *overall_sparse_matrix = NULL;
	
	gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), (void*)cells_indexes_src, num_weights, sizeof(long), all_array_size, (void**)(&overall_cells_indexes_src), num_overall_wgts, comp_node->get_comm_group());
	gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), (void*)cells_indexes_dst, num_weights, sizeof(long), all_array_size, (void**)(&overall_cells_indexes_dst), num_overall_wgts, comp_node->get_comm_group());
	gather_array_in_one_comp(comp_node->get_num_procs(), comp_node->get_current_proc_local_id(), (void*)weight_values, num_weights, sizeof(long), all_array_size, (void**)(&overall_wgt_values), num_overall_wgts, comp_node->get_comm_group());

	if (comp_node->get_current_proc_local_id() == 0) {
		overall_sparse_matrix = new Remap_weight_sparse_matrix(remap_operator, num_overall_wgts, overall_cells_indexes_src, overall_cells_indexes_dst, overall_wgt_values, 0, NULL);
		if (report_error_enabled) {
			bool *checking_mask = new bool [remap_operator->get_dst_grid()->get_grid_size()];
			for (i = 0; i < remap_operator->get_dst_grid()->get_grid_size(); i ++)
				checking_mask[i] = false;		
			for (i = 0, offset = 0; i < comp_node->get_num_procs(); i ++) {
				for (j = 0; j < all_array_size[i]; j ++) {
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, overall_cells_indexes_dst[offset+j] >= 0 && overall_cells_indexes_dst[offset+j] < remap_operator->get_dst_grid()->get_grid_size(), "Software error in Remap_weight_sparse_matrix::gather: %ld vs %ld", overall_cells_indexes_dst[offset+j], remap_operator->get_dst_grid()->get_grid_size());
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !checking_mask[overall_cells_indexes_dst[offset+j]], "Software error in Remap_weight_sparse_matrix::gather");
				}
				for (j = 0; j < all_array_size[i]; j ++)
					if (overall_cells_indexes_dst[offset+j] != -1)
						checking_mask[overall_cells_indexes_dst[offset+j]] = true;
				offset += all_array_size[i];
			}
			delete [] checking_mask;
		}				
		EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "The overall remapping sparse matrix have %ld weights", num_overall_wgts);
	}	

	delete [] all_array_size;
	return overall_sparse_matrix;
}


