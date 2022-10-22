/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DECOMP_INFO
#define DECOMP_INFO

#include "common_utils.h"
#include "remap_grid_class.h"
#include "remap_weight_of_strategy_class.h"
#include "original_grid_mgt.h"
#include <vector>
#include <map>


class Original_grid_info;

class Decomp_info
{
    private:
        char decomp_name[NAME_STR_SIZE];
        char grid_name[NAME_STR_SIZE];
        int num_global_cells;
        int num_local_cells;
        int *local_cell_global_indx;
		int *local_cell_chunk_id;
        bool is_registered;
        int comp_id;
        int host_comp_id;
        int grid_id;
        int decomp_id;
		int num_chunks;
		int *chunks_size;
		int *chunks_start;
        Decomp_info *halo_host_decomp;

    public:
        Decomp_info(const char *, int, int, int, int, const int *, int, const int *, const char *, bool);
        ~Decomp_info();
        int get_num_local_cells() { return num_local_cells; }
        const int *get_local_cell_global_indx() { return local_cell_global_indx; }
        const char *get_grid_name() { return grid_name; }
        const char *get_decomp_name() { return decomp_name; }
        int get_grid_id() { return grid_id; }
        int get_comp_id() { return comp_id; }
        int get_host_comp_id() { return host_comp_id; }
        int get_decomp_id() { return decomp_id; }
        int get_num_global_cells() { return num_global_cells; }
        bool get_is_fully_grid_decomp() { return num_local_cells == num_global_cells; }
        bool is_registered_decomp() { return is_registered; }
        void check_local_cell_global_indx();
        void set_halo_host_decomp(Decomp_info*);
        Decomp_info *get_halo_host_decomp() { return halo_host_decomp; }
		int get_num_chunks() { return num_chunks; }
		int get_chunk_size(int i) { return chunks_size[i]; }
		int get_chunk_start(int i) { return chunks_start[i]; }
		const int *get_chunks_start() { return chunks_start; }
		const int *get_local_cell_chunk_id() { return local_cell_chunk_id; }
};


class Decomp_info_mgt
{
    private:
        std::vector<Decomp_info *> decomps_info;
        std::map<int, int> fully_decomps_map;
        std::map<int, int> empty_decomps_map;
        
    public:
        Decomp_info_mgt() {}
        ~Decomp_info_mgt();
        int generate_fully_decomp(int);
        int generate_empty_decomp(int);
        Decomp_info *generate_remap_weights_src_decomp(Decomp_info*, Original_grid_info*, Original_grid_info*, Remap_weight_of_strategy_class*);
        Decomp_info *search_decomp_info(const char*, int);
        int register_H2D_parallel_decomposition(const char *, int, int, int, const int *, int, const int *, const char *);
        bool is_decomp_id_legal(int) const;
        int get_comp_id_of_decomp(int) const;
        Remap_grid_class *get_CoR_grid_of_decomp(int) const;
        Decomp_info *get_decomp_info(int);    
        int register_halo_parallel_decomposition(const char *, int, int, int *, int *, const char *);
		Decomp_info *generate_default_parallel_decomp(Original_grid_info *, int);
        Decomp_info *generate_default_parallel_decomp_serial(Original_grid_info *);
        Decomp_info *generate_consecutive_parallel_decomp_for_parallel_IO(char*, Original_grid_info*, int, int, int);
        Decomp_info *generate_parallel_decomp_for_parallel_IO(Original_grid_info*, int, int, int);
        Decomp_info *generate_rectangle_parallel_decomp_for_parallel_IO(char*, Original_grid_info*, int, int&, MPI_Comm&, int &);
        Decomp_info *set_rectangle_decomp_parameter(int&, int&, int&, int&, int&, int&, int&, int&, int, Original_grid_info*);
};


#endif
