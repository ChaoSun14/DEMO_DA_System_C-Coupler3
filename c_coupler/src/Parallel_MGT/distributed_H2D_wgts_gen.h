/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DISTRIBUTED_H2D_WGT
#define DISTRIBUTED_H2D_WGT


#include "remap_operator_basis.h"
#include "distributed_H2D_grid.h"


class Distributed_H2D_weights_generator;

class Grid_cell_rearrange_map_element
{
	public:
		int domain_decomp_process_id;
		int owner_process_id;
		int owner_local_index;
		int local_cell_global_index;
		
		Grid_cell_rearrange_map_element(int, int, int, int);
		Grid_cell_rearrange_map_element() {}
		~Grid_cell_rearrange_map_element() {}
};


class Normal_distributed_wgt_element
{
	public:
		int owner_process_id;
		int dst_local_index;
		long src_cell_index;
		long dst_cell_index;
		double wgt_value;

		Normal_distributed_wgt_element(long, long, double);
		Normal_distributed_wgt_element(){}
		~Normal_distributed_wgt_element() {}
};


class Distributed_H2D_weights_generator
{
	private:
		int comp_id;
		int dst_decomp_id;
		int current_subdomain_index;
		Remap_operator_basis *entire_remap_operator;
		Remap_weight_sparse_matrix *normal_remap_weights;
		Remapping_grid_domain_decomp_engine *remapping_grid_domain_decomp_engine;

		void append_wgt_sparse_matrix(std::vector<Normal_distributed_wgt_element> &, Remap_weight_sparse_matrix *);

	public:
		Distributed_H2D_weights_generator(int, int, int, int, Remap_operator_basis *);
		~Distributed_H2D_weights_generator() { if (normal_remap_weights != NULL) delete normal_remap_weights; }
		bool should_enlarge_src_subdomain_grid_for_remapping(long, double);
		void check_consistency_of_normal_remap_weights(Remap_weight_sparse_matrix *);
		bool confirm_or_enlarge_current_src_subdomain_grid_for_remapping(long dst_cell_index, double radius);
		Remap_weight_sparse_matrix *extract_normal_remap_weights();		
		Remap_weight_sparse_matrix *read_normal_remap_weights(Decomp_info *, const char *, Remap_operator_basis *);
};


extern void sort_normal_remapping_weights_in_sparse_matrix_locally(Remap_weight_sparse_matrix *);


#endif


