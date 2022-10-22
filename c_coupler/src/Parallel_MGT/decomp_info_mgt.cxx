/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "decomp_info_mgt.h"
#include "global_data.h"
#include "memory_mgt.h"
#include "cor_global_data.h"
#include "cor_cpl_interface.h"
#include "datamodel_mgt.h"
#include <string.h>
#include <math.h>
      

Decomp_info::Decomp_info(const char *decomp_name, int decomp_id, int host_comp_id, int grid_id, int num_local_cells, const int *cell_indexes_in_decomp, int num_chunks, const int *chunks_size, const char *annotation, bool registered)
{
    Remap_grid_class *CoR_grid;
    int i, iter;
    

    EXECUTION_REPORT(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(grid_id), "Software error in Decomp_info::Decomp_info: wrong grid_id");
    
    this->decomp_id = decomp_id;
    this->grid_id = grid_id;
    this->halo_host_decomp = NULL;
    this->comp_id = original_grid_mgr->get_comp_id_of_grid(grid_id);
	this->chunks_size = NULL;
	this->chunks_start = NULL;
	this->local_cell_chunk_id = NULL;
	this->num_chunks = num_chunks;
    if (host_comp_id == -1)
        host_comp_id = this->comp_id;
    this->host_comp_id = host_comp_id;
    this->num_global_cells = original_grid_mgr->get_original_CoR_grid(grid_id)->get_grid_size();
    this->num_local_cells = num_local_cells;
    this->local_cell_global_indx = NULL;
    strcpy(this->decomp_name, decomp_name);
    strcpy(this->grid_name, original_grid_mgr->search_grid_info(grid_id)->get_grid_name());
    annotation_mgr->add_annotation(decomp_id, "register decomposition", annotation);
    is_registered = registered;
    synchronize_comp_processes_for_API(host_comp_id, API_ID_DECOMP_MGT_REG_NORMAL_DECOMP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id,"C-Coupler code in Decomp_info::Decomp_info for getting comm group"), "for register a parallel decomposition of a grid", annotation);
    check_API_parameter_string(host_comp_id, API_ID_DECOMP_MGT_REG_NORMAL_DECOMP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id,"C-Coupler code in Decomp_info::Decomp_info for getting comm group"), "for register a parallel decomposition of a grid", decomp_name, "decomp_name", annotation);
	if (is_registered)
	    check_API_parameter_string(host_comp_id, API_ID_DECOMP_MGT_REG_NORMAL_DECOMP, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id,"C-Coupler code in Decomp_info::Decomp_info for getting comm group"), "for register a parallel decomposition of a grid", original_grid_mgr->get_name_of_grid(grid_id), "grid_id (the corresponding grid name)", annotation);    

    CoR_grid = original_grid_mgr->search_grid_info(grid_id)->get_original_CoR_grid();
    EXECUTION_REPORT(REPORT_ERROR, -1, CoR_grid->get_is_H2D_grid(), "Software error in Decomp_info::Decomp_info: not on a horizontal grid");
    EXECUTION_REPORT(REPORT_ERROR, -1, num_local_cells >= 0, "Software error in Decomp_info::Decomp_info: wrong num_local_cells");

	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, 
					 "parallel decomposition \"%s\" on the grid \"%s\" has %d local grid cells", 
					 decomp_name, original_grid_mgr->search_grid_info(grid_id)->get_grid_name(), num_local_cells);
    if (num_local_cells == 0)
        local_cell_global_indx = NULL;
    else {
        local_cell_global_indx = new int [num_local_cells];
        for (i = 0; i < num_local_cells; i ++) {
            if (cell_indexes_in_decomp[i] == CCPL_NULL_INT)
                local_cell_global_indx[i] = CCPL_NULL_INT;
            else {
                EXECUTION_REPORT(REPORT_ERROR, -1, cell_indexes_in_decomp[i] > 0 && cell_indexes_in_decomp[i] <= CoR_grid->get_grid_size(), "Software error in Decomp_info::Decomp_info: wrong cell index: %d vs [1,%d]", cell_indexes_in_decomp[i], CoR_grid->get_grid_size());
                local_cell_global_indx[i] = cell_indexes_in_decomp[i] - 1;  // -1 because fortran array index starts from 1 but c/c++ starts from 0
            }
        }
		if (num_chunks > 0) {
			this->chunks_size = new int [num_chunks];
			this->chunks_start = new int [num_chunks];
			this->local_cell_chunk_id = new int [num_local_cells];
			for (i = 0, iter = 0; i < num_chunks; i ++) {
				this->chunks_start[i] = iter;
				this->chunks_size[i] = chunks_size[i];
				iter += chunks_size[i];
				for (int j = 0; j < this->chunks_size[i]; j ++)
					this->local_cell_chunk_id[this->chunks_start[i]+j] = i;
			}
			EXECUTION_REPORT(REPORT_ERROR, -1, num_local_cells == iter, "Software error in Decomp_info::Decomp_info");
		}
    }
}


Decomp_info::~Decomp_info()
{
    if (local_cell_global_indx != NULL)
        delete [] local_cell_global_indx;
	if (chunks_size != NULL) {
		delete [] chunks_size;
		delete [] chunks_start;
		delete [] local_cell_chunk_id;
	}
}


void Decomp_info::check_local_cell_global_indx()
{
    for (int i = 0; i < num_local_cells; i ++)
        if (local_cell_global_indx[i] != CCPL_NULL_INT)
            EXECUTION_REPORT(REPORT_ERROR,-1, local_cell_global_indx[i] >= 0 && local_cell_global_indx[i] < num_global_cells, "C-Coupler error in check_local_cell_global_indx\n");
}


void Decomp_info::set_halo_host_decomp(Decomp_info *halo_host_decomp)
{
    this->halo_host_decomp = halo_host_decomp;
	if (halo_host_decomp->num_chunks > 0) {
		this->num_chunks = halo_host_decomp->num_chunks;
		this->chunks_size = new int [this->num_chunks];
		this->chunks_start = new int [this->num_chunks];
		for (int i = 0; i < this->num_chunks; i ++) {
			this->chunks_size[i] = halo_host_decomp->chunks_size[i];
			this->chunks_start[i] = halo_host_decomp->chunks_start[i];
		}
		this->local_cell_chunk_id = new int [num_local_cells];
		for (int i = 0; i < num_local_cells; i ++)
			this->local_cell_chunk_id[i] = halo_host_decomp->local_cell_chunk_id[i];
	}
}


Decomp_info_mgt::~Decomp_info_mgt()
{
    for (int i = 0; i < decomps_info.size(); i ++)
        delete decomps_info[i];
}

int Decomp_info_mgt::generate_empty_decomp(int original_decomp_id)
{
    char empty_decomp_name[NAME_STR_SIZE];
    Decomp_info *empty_decomp;

    if (empty_decomps_map.find(original_decomp_id) != empty_decomps_map.end())
        return empty_decomps_map[original_decomp_id];

    sprintf(empty_decomp_name, "empty_decomp_for_%s", get_decomp_info(original_decomp_id)->get_decomp_name());
    empty_decomp = search_decomp_info(empty_decomp_name, get_decomp_info(original_decomp_id)->get_comp_id());
    if (empty_decomp != NULL)
        return empty_decomp->get_decomp_id();

    empty_decomp = new Decomp_info(empty_decomp_name, (TYPE_DECOMP_ID_PREFIX|decomps_info.size()), get_decomp_info(original_decomp_id)->get_host_comp_id(), get_decomp_info(original_decomp_id)->get_grid_id(), 0, NULL, 0, NULL, "empty decomp", false);
    decomps_info.push_back(empty_decomp);
    empty_decomps_map[original_decomp_id] = empty_decomp->get_decomp_id();
    return empty_decomp->get_decomp_id();
}

int Decomp_info_mgt::generate_fully_decomp(int original_decomp_id)
{
    char fully_decomp_name[NAME_STR_SIZE];
    Decomp_info *fully_decomp;
    int *local_cells_global_indexes, num_global_cells;


    if (fully_decomps_map.find(original_decomp_id) != fully_decomps_map.end())
        return fully_decomps_map[original_decomp_id];

    sprintf(fully_decomp_name, "fully_decomp_for_%s", get_decomp_info(original_decomp_id)->get_decomp_name());
    fully_decomp = search_decomp_info(fully_decomp_name, get_decomp_info(original_decomp_id)->get_comp_id());
    if (fully_decomp != NULL)
        return fully_decomp->get_decomp_id();

    if (comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(get_decomp_info(original_decomp_id)->get_host_comp_id(), "in Decomp_info_mgt::generate_fully_decomp") == 0) {
        num_global_cells = get_decomp_info(original_decomp_id)->get_num_global_cells();
        local_cells_global_indexes = new int [num_global_cells];
        for (int i = 0; i < num_global_cells; i ++)
            local_cells_global_indexes[i] = i + 1;
        fully_decomp = new Decomp_info(fully_decomp_name, (TYPE_DECOMP_ID_PREFIX|decomps_info.size()), get_decomp_info(original_decomp_id)->get_host_comp_id(), get_decomp_info(original_decomp_id)->get_grid_id(), num_global_cells, local_cells_global_indexes, 0, NULL, "fully decomp", false);
        delete [] local_cells_global_indexes;    
    }
    else fully_decomp = new Decomp_info(fully_decomp_name, (TYPE_DECOMP_ID_PREFIX|decomps_info.size()), get_decomp_info(original_decomp_id)->get_host_comp_id(), get_decomp_info(original_decomp_id)->get_grid_id(), 0, NULL, 0, NULL, "fully decomp", false);
    decomps_info.push_back(fully_decomp);

    fully_decomps_map[original_decomp_id] = fully_decomp->get_decomp_id();

    return fully_decomp->get_decomp_id();
}


Decomp_info *Decomp_info_mgt::generate_remap_weights_src_decomp(Decomp_info *dst_decomp_info, Original_grid_info *src_original_grid, Original_grid_info *dst_original_grid, Remap_weight_of_strategy_class *remap_weights)
{
    Decomp_info *decomp_for_remap;    
    char decomp_name_remap[NAME_STR_SIZE];
    int i, j, src_H2D_original_grid_id;


    sprintf(decomp_name_remap, "src_decomp_for_%s_%s_%s", remap_weights->get_object_name(), dst_decomp_info->get_decomp_name(), comp_comm_group_mgt_mgr->get_global_node_of_local_comp(dst_decomp_info->get_comp_id(),true,"in Decomp_info_mgt::generate_remap_weights_src_decomp")->get_comp_full_name());
	if (report_error_enabled)
	    dst_decomp_info->check_local_cell_global_indx();
    Original_grid_info *existing_src_H2D_original_grid = original_grid_mgr->search_grid_info(src_original_grid->get_comp_id(), src_original_grid->get_H2D_sub_CoR_grid());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, existing_src_H2D_original_grid != NULL && existing_src_H2D_original_grid->get_comp_id() == src_original_grid->get_comp_id(), "Software error in Decomp_info_mgt::generate_remap_weights_src_decomp: %s", src_original_grid->get_H2D_sub_CoR_grid()->get_grid_name());
    if (existing_src_H2D_original_grid != NULL)
        src_H2D_original_grid_id = existing_src_H2D_original_grid->get_grid_id();
    else src_H2D_original_grid_id = original_grid_mgr->add_original_grid(src_original_grid->get_comp_id(), src_original_grid->get_H2D_sub_CoR_grid()->get_grid_name(), src_original_grid->get_H2D_sub_CoR_grid());
    decomp_for_remap = search_decomp_info(decomp_name_remap, src_original_grid->get_comp_id());
    if (decomp_for_remap == NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current_remap_local_cell_global_indexes != NULL , "Software error in Decomp_info_mgt::generate_remap_weights_src_decomp");
        decomp_for_remap = new Decomp_info(decomp_name_remap, (TYPE_DECOMP_ID_PREFIX|decomps_info.size()), dst_original_grid->get_comp_id(), src_H2D_original_grid_id, num_current_remap_local_cell_global_indexes, current_remap_local_cell_global_indexes, 0, NULL, "in Decomp_info_mgt::generate_remap_weights_src_decomp", false);
        decomps_info.push_back(decomp_for_remap);
    }
	if (current_remap_local_cell_global_indexes != NULL) {
		delete [] current_remap_local_cell_global_indexes;
		current_remap_local_cell_global_indexes = NULL;
	}

    return decomp_for_remap;
}


Decomp_info *Decomp_info_mgt::search_decomp_info(const char *decomp_name, int comp_id)
{
    for (int i = 0; i < decomps_info.size(); i ++)
        if (words_are_the_same(decomps_info[i]->get_decomp_name(), decomp_name) && decomps_info[i]->get_comp_id() == comp_id) {
            return decomps_info[i];
        }

    return NULL;
}


int Decomp_info_mgt::register_H2D_parallel_decomposition(const char *decomp_name, int grid_id, int host_comp_id, int num_local_cells, const int *cell_indexes_in_decomp, int num_chunks, const int *chunks_size, const char *annotation)
{
    Decomp_info *new_decomp = new Decomp_info(decomp_name, (TYPE_DECOMP_ID_PREFIX|decomps_info.size()), host_comp_id, grid_id, num_local_cells, cell_indexes_in_decomp, num_chunks, chunks_size, annotation, true);

    if (search_decomp_info(decomp_name, original_grid_mgr->get_comp_id_of_grid(grid_id)) != NULL)
        EXECUTION_REPORT(REPORT_ERROR, new_decomp->get_comp_id(), false, 
                         "Error happens when calling the API \"CCPL_register_normal_parallel_decomp\" to register a parallel decomposition \"%s\" at the model code with the annotation \"%s\": a parallel decomposition with the same name has already been registered before at the model code with the annotations \"%s\". Please verify.",
                         decomp_name, annotation, annotation_mgr->get_annotation(search_decomp_info(decomp_name, original_grid_mgr->get_comp_id_of_grid(grid_id))->get_decomp_id(), "register decomposition"), annotation);

    decomps_info.push_back(new_decomp);

    return new_decomp->get_decomp_id();
} 
    
    
bool Decomp_info_mgt::is_decomp_id_legal(int decomp_id) const
{
    int true_decomp_id = decomp_id & TYPE_ID_SUFFIX_MASK;

    
    if ((decomp_id & TYPE_ID_PREFIX_MASK) != TYPE_DECOMP_ID_PREFIX)
        return false;

    if (true_decomp_id < 0 || true_decomp_id >= decomps_info.size())
        return false;

    return true;
}


int Decomp_info_mgt::get_comp_id_of_decomp(int decomp_id) const
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_decomp_id_legal(decomp_id), "Software error when calling Decomp_info_mgt::get_comp_id_of_decomp");
    return decomps_info[decomp_id&TYPE_ID_SUFFIX_MASK]->get_comp_id();
}


Remap_grid_class *Decomp_info_mgt::get_CoR_grid_of_decomp(int decomp_id) const
{    
    EXECUTION_REPORT(REPORT_ERROR, -1, is_decomp_id_legal(decomp_id), "Software error when calling Decomp_info_mgt::get_comp_id_of_decomp");
    return original_grid_mgr->get_original_CoR_grid(decomps_info[decomp_id&TYPE_ID_SUFFIX_MASK]->get_grid_id());
}


Decomp_info *Decomp_info_mgt::get_decomp_info(int decomp_id)
{
    EXECUTION_REPORT(REPORT_ERROR, -1, is_decomp_id_legal(decomp_id), "Software error when calling Decomp_info_mgt::get_decomp_info");
    return decomps_info[decomp_id&TYPE_ID_SUFFIX_MASK];
}


int Decomp_info_mgt::register_halo_parallel_decomposition(const char *halo_name, int decomp_id, int num_local_cells, int *local_cells_local_indexes, int *local_cells_global_indexes, const char *annotation)
{
    char halo_full_name[NAME_STR_SIZE];
    Decomp_info *host_decomp = decomps_info_mgr->get_decomp_info(decomp_id);
    int comp_id = host_decomp->get_comp_id();
    int *halo_decomp_local_cells_global_indexes = NULL;


    EXECUTION_REPORT(REPORT_ERROR, host_decomp->get_host_comp_id(), num_local_cells <= host_decomp->get_num_local_cells(), "Error happens when calling API \"CCPL_register_halo_region\" to register a halo region named \"%s\" for a parallel decomposition named \"%s\": the value of \"num_local_cells\" (%d) cannot not be larger than the number of local cells in the parallel decomposition (the number is %d). Please check the model code with the annotation \"%s\"", halo_name, host_decomp->get_decomp_name(), num_local_cells, host_decomp->get_num_local_cells(), annotation);
    synchronize_comp_processes_for_API(host_decomp->get_host_comp_id(), API_ID_DECOMP_MGT_REG_HALO, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_decomp->get_host_comp_id(),"register_halo_parallel_decomposition"), "for register a halo region corresponding to a parallel decomposition", annotation);
    check_API_parameter_string(host_decomp->get_host_comp_id(), API_ID_DECOMP_MGT_REG_HALO, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_decomp->get_host_comp_id(),"register_halo_parallel_decomposition"), "for register a halo region corresponding to a parallel decomposition", halo_name, "halo_name", annotation);
    check_API_parameter_string(host_decomp->get_host_comp_id(), API_ID_DECOMP_MGT_REG_HALO, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_decomp->get_host_comp_id(),"register_halo_parallel_decomposition"), "for register a halo region corresponding to a parallel decomposition", host_decomp->get_decomp_name(), "decomp_id (the corresponding parallel decomposition name)", annotation);
    sprintf(halo_full_name, "%s@%s", host_decomp->get_decomp_name(), halo_name);
    for (int i = 0; i < decomps_info.size(); i ++)
        if (decomps_info[i]->get_comp_id() == comp_id && words_are_the_same(halo_full_name, decomps_info[i]->get_decomp_name()))
            EXECUTION_REPORT(REPORT_ERROR, host_decomp->get_host_comp_id(), false, "Error happens when calling API \"CCPL_register_halo_region\" at the model code with the annotation \"%s\" to register a halo region named \"%s\" for a parallel decomposition named \"%s\": a halo region with the same name and the same parallel decomposition has already been registered at the model code with the annotation \"%s\". Please verify.", annotation, halo_name, host_decomp->get_decomp_name(), annotation_mgr->get_annotation(decomps_info[i]->get_decomp_id(), "register decomposition"));
    if (host_decomp->get_num_local_cells() > 0) {
        halo_decomp_local_cells_global_indexes = new int [host_decomp->get_num_local_cells()];
        for (int i = 0; i < host_decomp->get_num_local_cells(); i ++)
            halo_decomp_local_cells_global_indexes[i] = CCPL_NULL_INT;
        for (int i = 0; i < num_local_cells; i ++) {
            EXECUTION_REPORT(REPORT_ERROR, host_decomp->get_host_comp_id(), local_cells_local_indexes[i] >= 1 && local_cells_local_indexes[i] <= host_decomp->get_num_local_cells(), "Error happens when calling API \"CCPL_register_halo_region\" to register a halo region named \"%s\" for a parallel decomposition named \"%s\": the number %d element of \"local_cells_local_indexes\" (the value is %d) is out of bounds determined by the number of local cells in the parallel decomposition (the value is %d). Please check the model code with the annotation \"%s\"", halo_name, host_decomp->get_decomp_name(), i+1, local_cells_local_indexes[i], host_decomp->get_num_local_cells(), annotation);
            EXECUTION_REPORT(REPORT_ERROR, host_decomp->get_host_comp_id(), local_cells_global_indexes[i] >= 1 && local_cells_global_indexes[i] <= host_decomp->get_num_global_cells(), "Error happens when calling API \"CCPL_register_halo_region\" to register a halo region named \"%s\" for a parallel decomposition named \"%s\": the number %d element of \"local_cells_global_indexes\" (the value is %d) is out of bounds determined by the grid size of the corresponding horizontal grid (the grid size is %d). Please check the model code with the annotation \"%s\"", halo_name, host_decomp->get_decomp_name(), i+1, local_cells_global_indexes[i], host_decomp->get_num_global_cells(), annotation);
            EXECUTION_REPORT(REPORT_ERROR, host_decomp->get_host_comp_id(), host_decomp->get_local_cell_global_indx()[local_cells_local_indexes[i]-1] == CCPL_NULL_INT, "Error happens when calling API \"CCPL_register_halo_region\" to register a halo region named \"%s\" for a parallel decomposition named \"%s\": regarding the No. %d element of \"local_cells_local_indexes\" (value is %d), the corresponding global grid index in the parallel decomposition is not CCPL_NULL_INT (the global index is %d). Please note that the halo region cannot contain valid cells in the parallel decomposition. Please check the model code with the annotation \"%s\"", halo_name, host_decomp->get_decomp_name(), i+1, local_cells_local_indexes[i], host_decomp->get_local_cell_global_indx()[local_cells_local_indexes[i]-1]+1, annotation);
            EXECUTION_REPORT(REPORT_ERROR, host_decomp->get_host_comp_id(), halo_decomp_local_cells_global_indexes[local_cells_local_indexes[i]-1] == CCPL_NULL_INT || halo_decomp_local_cells_global_indexes[local_cells_local_indexes[i]-1] == local_cells_global_indexes[i], "Error happens when calling API \"CCPL_register_halo_region\" to register a halo region named \"%s\" for a parallel decomposition named \"%s\": there are more than one element in \"local_cells_local_indexes\" corresponding to the same local cell in the parallel decomposition (the local index of the local cell is %d). Please check the model code with the annotation \"%s\"", halo_name, host_decomp->get_decomp_name(), local_cells_local_indexes[i], annotation);
            halo_decomp_local_cells_global_indexes[local_cells_local_indexes[i]-1] = local_cells_global_indexes[i];
        }
    }
    Decomp_info *halo_decomp = new Decomp_info(halo_full_name, (TYPE_DECOMP_ID_PREFIX|decomps_info.size()), -1, host_decomp->get_grid_id(), host_decomp->get_num_local_cells(), halo_decomp_local_cells_global_indexes, 0, NULL, annotation, false);
    halo_decomp->set_halo_host_decomp(host_decomp);
    decomps_info.push_back(halo_decomp);
    if (halo_decomp_local_cells_global_indexes != NULL)
        delete [] halo_decomp_local_cells_global_indexes;

    return halo_decomp->get_decomp_id();
}


Decomp_info *Decomp_info_mgt::generate_parallel_decomp_for_parallel_IO(Original_grid_info *original_grid, int io_proc_stride, int num_io_procs, int io_proc_mark) {
    Original_grid_info *H2D_sub_grid = original_grid->get_H2D_sub_grid();
    char parallel_IO_decomp_name[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grid->get_comp_id());

    if (H2D_sub_grid == NULL)
        return NULL;
    sprintf(parallel_IO_decomp_name, "parallel_IO_decomp_name_of_%s", H2D_sub_grid->get_grid_name());

    Decomp_info *existing_decomp_info = search_decomp_info(parallel_IO_decomp_name, original_grid->get_comp_id());
    if (existing_decomp_info != NULL)
        return existing_decomp_info;

    existing_decomp_info = generate_consecutive_parallel_decomp_for_parallel_IO(parallel_IO_decomp_name, original_grid, io_proc_stride, num_io_procs, io_proc_mark);

    return existing_decomp_info;
}


Decomp_info *Decomp_info_mgt::set_rectangle_decomp_parameter(int &num_columns, int &num_local_lons, int &num_rows, int &num_local_lats, int &num_io_procs, int &io_proc_mark, int &global_index_start_lon, int &global_index_start_lat, int num_total_procs, Original_grid_info *original_grid) {

    int io_proc_stride, num_sized_sub_grids;
    Remap_grid_class *sized_sub_grids[256];
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grid->get_comp_id());
    original_grid->get_original_CoR_grid()->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);

    int CoR_lon_grid_size = sized_sub_grids[0]->get_grid_size();
    int CoR_lat_grid_size = sized_sub_grids[1]->get_grid_size();
    int num_cell_each_dim = sqrt(CoR_lon_grid_size * CoR_lat_grid_size / num_io_procs);
    num_columns = CoR_lon_grid_size / num_cell_each_dim;//unchagable
    num_rows = CoR_lat_grid_size / num_cell_each_dim;//unchagable
    num_local_lons = CoR_lon_grid_size / num_columns;
    num_local_lats = CoR_lat_grid_size / num_rows;
    num_io_procs = num_columns * num_rows;

    io_proc_stride = num_total_procs / num_io_procs;
    if (io_proc_stride == 1)
        io_proc_mark = comp_node->get_current_proc_local_id() < num_io_procs? 1: 0;
    else io_proc_mark = (comp_node->get_current_proc_local_id()) % io_proc_stride == 0? 1: 0;
    int io_proc_index = comp_node->get_current_proc_local_id() / io_proc_stride;
    global_index_start_lon = num_local_lons * io_proc_index;
    global_index_start_lat = num_local_lats * io_proc_index;

    if (io_proc_mark == 1) {
        if (CoR_lon_grid_size % num_columns > io_proc_index) {
            num_local_lons++;
            global_index_start_lon += io_proc_index;
        }
        else global_index_start_lon += CoR_lon_grid_size % num_columns;
        if (CoR_lat_grid_size % num_rows > io_proc_index) {
            num_local_lats ++;
            global_index_start_lat += io_proc_index;
        }
        else global_index_start_lat += CoR_lat_grid_size % num_rows;
    }
    else {
        num_local_lons = 0;
        num_local_lats = 0;
    }
}

Decomp_info *Decomp_info_mgt::generate_rectangle_parallel_decomp_for_parallel_IO(char *parallel_IO_decomp_name, Original_grid_info *original_grid, int num_total_procs, int &num_io_procs, MPI_Comm &io_comm, int &io_proc_mark) {
    int num_local_cells, num_columns, num_rows, num_local_lons, num_local_lats;
    int color = 0, global_index_start_lon, global_index_start_lat, num_sized_sub_grids;
    int *local_cells_global_index;
    Remap_grid_class *sized_sub_grids[256];

    original_grid->get_original_CoR_grid()->get_sized_sub_grids(&num_sized_sub_grids, sized_sub_grids);
    Original_grid_info *H2D_sub_grid = original_grid->get_H2D_sub_grid();
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grid->get_comp_id());

    set_rectangle_decomp_parameter(num_columns, num_local_lons, num_rows, num_local_lats, num_io_procs, io_proc_mark, global_index_start_lon, global_index_start_lat, num_total_procs, original_grid);

    num_local_cells = num_local_lons * num_local_lats;

    local_cells_global_index = new int [num_local_cells];
    for (int i = 0; i < num_local_lons; i++) {
        for (int j = 0; j < num_local_lats; j++) {
            int index_lon = global_index_start_lon + i;
            int index_lat = global_index_start_lat + j;
            local_cells_global_index[i * num_local_lats + j] = index_lon * num_rows + index_lat + 1;
        }
    }
    Decomp_info *existing_decomp_info = get_decomp_info(register_H2D_parallel_decomposition(parallel_IO_decomp_name, H2D_sub_grid->get_grid_id(), original_grid->get_comp_id(), num_local_cells, local_cells_global_index, 0, NULL, "in Decomp_info_mgt::generate_parallel_decomp_for_parallel_IO"));
    EXECUTION_REPORT(REPORT_LOG, -1, true, "Decomp_info_mgt::generate_parallel_decomp_for_parallel_IO lon:%d, lat:%d", sized_sub_grids[0]->get_grid_size(), sized_sub_grids[1]->get_grid_size());
    delete [] local_cells_global_index;
    if (io_proc_mark == 1)
        color = 1;
    MPI_Comm_split(comp_node->get_comm_group(), color, 0, &io_comm);
    return existing_decomp_info;
}


Decomp_info *Decomp_info_mgt::generate_consecutive_parallel_decomp_for_parallel_IO(char *parallel_IO_decomp_name, Original_grid_info *original_grid, int io_proc_stride, int num_io_procs, int io_proc_mark) {
    Original_grid_info *H2D_sub_grid = original_grid->get_H2D_sub_grid();
    Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grid->get_comp_id());
    int num_local_cells, *local_cells_global_index = NULL, global_index_start, color = 0;

    EXECUTION_REPORT(REPORT_ERROR, -1, num_io_procs != 0, "Software error in Decomp_info_mgt::generate_consecutive_parallel_decomp_for_parallel_IO");

    if (io_proc_mark == 1) {
		int io_proc_index = comp_node->get_current_proc_local_id() / io_proc_stride;
		num_local_cells = H2D_sub_grid->get_original_CoR_grid()->get_grid_size() / num_io_procs;
		global_index_start = num_local_cells * io_proc_index;
        if ((H2D_sub_grid->get_original_CoR_grid()->get_grid_size() % num_io_procs) > io_proc_index) {
            num_local_cells ++;
            global_index_start += io_proc_index;
        }
        else global_index_start += H2D_sub_grid->get_original_CoR_grid()->get_grid_size() % num_io_procs;
		local_cells_global_index = new int [num_local_cells + 1];
		for (int i = 0; i < num_local_cells; i++)
			local_cells_global_index[i] = global_index_start + i + 1;
    }
    else num_local_cells = 0;

    Decomp_info *existing_decomp_info = get_decomp_info(register_H2D_parallel_decomposition(parallel_IO_decomp_name, H2D_sub_grid->get_grid_id(), H2D_sub_grid->get_comp_id(), num_local_cells, local_cells_global_index, 0, NULL, "in Decomp_info_mgt::generate_parallel_decomp_for_parallel_IO"));

	if (io_proc_mark == 1)
	    delete [] local_cells_global_index;

    return existing_decomp_info;
}


Decomp_info *Decomp_info_mgt::generate_default_parallel_decomp_serial(Original_grid_info *original_grid) 
{
    Original_grid_info *H2D_sub_grid = original_grid->get_H2D_sub_grid();
    char default_parallel_decomp_name[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *comp_node;
    int num_local_cells, *local_cells_global_index, global_index_start;

    if (H2D_sub_grid == NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid->get_H2D_sub_CoR_grid() == NULL, "Software error in Decomp_info_mgt::generate_default_parallel_decomp_serial: %s", original_grid->get_grid_name());
        return NULL;
	}

    sprintf(default_parallel_decomp_name, "default_parallel_decomp_of_%s", H2D_sub_grid->get_grid_name());
    Decomp_info *existing_decomp_info = search_decomp_info(default_parallel_decomp_name, original_grid->get_comp_id());
    if (existing_decomp_info != NULL)
        return existing_decomp_info;

    comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grid->get_comp_id());
    int local_proc_id = comp_node->get_current_proc_local_id();
    num_local_cells = local_proc_id == 0? H2D_sub_grid->get_original_CoR_grid()->get_grid_size(): 0;
    global_index_start = 0;
    local_cells_global_index = new int [num_local_cells + 1];
    for (int i = 0; i < num_local_cells; i ++)
        local_cells_global_index[i] = global_index_start+i+1;
    existing_decomp_info = get_decomp_info(register_H2D_parallel_decomposition(default_parallel_decomp_name, H2D_sub_grid->get_grid_id(), H2D_sub_grid->get_comp_id(), num_local_cells, local_cells_global_index, 0, NULL, "in Decomp_info_mgt::generate_default_parallel_decomp_serial"));
    delete [] local_cells_global_index;

    return existing_decomp_info;
}

Decomp_info *Decomp_info_mgt::generate_default_parallel_decomp(Original_grid_info *original_grid, int host_comp_id)
{
	Original_grid_info *H2D_sub_grid = original_grid->get_H2D_sub_grid();
	char default_parallel_decomp_name[NAME_STR_SIZE];
	Comp_comm_group_mgt_node *comp_node;
	int num_local_cells, *local_cells_global_index, global_index_start;

	if (H2D_sub_grid == NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid->get_H2D_sub_CoR_grid() == NULL, "Software error in Decomp_info_mgt::generate_default_parallel_decomp");
		return NULL;
	}

	comp_node = comp_comm_group_mgt_mgr->search_global_node(original_grid->get_comp_id());
	if (host_comp_id != -1)
		comp_node = comp_comm_group_mgt_mgr->search_global_node(host_comp_id);	
	sprintf(default_parallel_decomp_name, "default_parallel_decomp_of_%s_at_COMP_%s", H2D_sub_grid->get_grid_name(), comp_node->get_full_name());
	Decomp_info *existing_decomp_info = search_decomp_info(default_parallel_decomp_name, original_grid->get_comp_id());
	if (existing_decomp_info != NULL)
		return existing_decomp_info;
	num_local_cells = H2D_sub_grid->get_original_CoR_grid()->get_grid_size() / comp_node->get_num_procs();
	global_index_start = num_local_cells * comp_node->get_current_proc_local_id();
	if (H2D_sub_grid->get_original_CoR_grid()->get_grid_size() % comp_node->get_num_procs() > comp_node->get_current_proc_local_id()) {
		num_local_cells ++;
		global_index_start += comp_node->get_current_proc_local_id();
	}
	else global_index_start += (H2D_sub_grid->get_original_CoR_grid()->get_grid_size() % comp_node->get_num_procs());
	local_cells_global_index = new int [num_local_cells + 1];
	for (int i = 0; i < num_local_cells; i ++)
		local_cells_global_index[i] = global_index_start+i+1;
	existing_decomp_info = get_decomp_info(register_H2D_parallel_decomposition(default_parallel_decomp_name, H2D_sub_grid->get_grid_id(), host_comp_id, num_local_cells, local_cells_global_index, 0, NULL, "in Decomp_info_mgt::generate_default_parallel_decomp"));
	delete [] local_cells_global_index;

	return existing_decomp_info;
}


