/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DISTRIBUTED_H2D_GRID
#define DISTRIBUTED_H2D_GRID


#include "remap_grid_class.h"
#include "remap_grid_data_class.h"
#include "decomp_grid_mgt.h"
#include "memory_mgt.h"
#include "triangle.h"
#include "triangulation.h"
#include "distributed_merge_sort.h"
#include <vector>
#include <mpi.h>

#define EARTH_RADIUS 6371.393

class Original_grid_info;
class Grid_cell_rearrange_map_element;


class Distributed_H2D_grid_engine
{
	private:
		const char *comp_full_name;
		const char *grid_name;
		int num_lons;
		int num_lats;
		long global_grid_size;
		Remap_grid_class *CoR_grid_initialization;
		Original_grid_info *original_grid;
		Decomp_grid_info *basic_decomp_grid;
		std::vector<Field_mem_info*> grid_fields_insts;
		
	public:
		Distributed_H2D_grid_engine(const char *, const char *, Remap_grid_class *, long, int, int, Original_grid_info *, bool);
		~Distributed_H2D_grid_engine();
		void initialize_data(const char *, const char *, long, int, int, Remap_grid_class *, Original_grid_info *);
		Decomp_grid_info *get_basic_decomp_grid();
		const char *get_comp_full_name() { return comp_full_name; }
		Original_grid_info *get_original_grid() { return original_grid; }
		Remap_grid_class *get_CoR_grid_initialization() { return CoR_grid_initialization; }
		Decomp_grid_info *generate_basic_decomp_grid_from_decomp_grid(const char *, Decomp_grid_info *, const char *);
		Decomp_grid_info *generate_decomp_grid_via_data_transfer(const char *, const char *, Decomp_grid_info *, const char *, const char *, int, Remap_grid_class *);
		void allocate_grid_field_inst(Remap_grid_data_class *, int, int, bool, Field_mem_info **, int &, int);
		void update_basic_decomp_grid(Decomp_grid_info *);
		bool is_the_same_as_another_distributed_H2D_grid(Distributed_H2D_grid_engine *);
};


struct Grid_point_subdomain_mapping_info
{
	int global_grid_cell_index;
	int subdomain_id;
	int owner_proc_id;
};


class Distributed_H2D_grid_mgt;


class Remapping_grid_domain_decomp_engine
{
	private:
		friend class Distributed_H2D_grid_mgt;
		int comp_id;
		int num_total_procs;
		int src_original_grid_id;
		Remap_grid_class *src_original_cor_grid;
		Remap_grid_class *dst_original_cor_grid;
		int src_grid_original_num_vertex;
		int dst_grid_original_num_vertex;
		double whole_domain_min_lon;
		double whole_domain_max_lon;
		double whole_domain_min_lat;
		double whole_domain_max_lat;
		double nonPole_subdomain_x_side_length;
		double nonPole_subdomain_y_side_length;
		double nonPole_region_min_lat;
		double nonPole_region_max_lat;
		double current_subdomain_min_lon;
		double current_subdomain_max_lon;
		double current_subdomain_min_lat;
		double current_subdomain_max_lat;
		int num_northPole_subdomain;
		int num_southPole_subdomain;
		int num_nonPole_region_x_subdomains;
		int num_nonPole_region_y_subdomains;
		int num_total_subdomains;
		Decomp_grid_info *src_subdomain_decomp_grid;
		Decomp_grid_info *dst_subdomain_decomp_grid;
		std::vector<std::pair<int, Remap_grid_class*> > src_subdomains_CoR_grids;
		std::vector<std::pair<int, Remap_grid_class*> > dst_subdomains_CoR_grids;
		std::vector<std::vector<Remap_grid_class *> > subdomains_halo_subdomains_CoR_grids;
		std::vector<std::pair<double, double> > subdomains_halo_lon_bounds;
		std::vector<std::pair<double, double> > subdomains_halo_lat_bounds;
		std::vector<Remap_grid_class *> src_subdomains_expanded_CoR_grids;
		std::vector<Remap_grid_class *> remote_subdomain_CoR_grids;
		std::vector<std::vector<int> > subdomains_halo_subdomains_IDs;
		std::vector<int> subdomains_halo_num_levels;
		std::vector<int> local_subdomains_index;
		std::vector<Triangle_inline*> subdomains_final_triangles;
		std::vector<int> num_subdomains_final_triangles;

		char *src_grid_one_sided_comm_data_buf, *dst_grid_one_sided_comm_data_buf;
		char *triangulation_comm_data_buf;
		long triangulation_comm_data_buf_max_size, triangulation_comm_data_buf_content_size;

		MPI_Win src_grid_MPI_Win, dst_grid_MPI_Win;
		MPI_Win triangulation_comm_win;
		int *subdomain_triangulation_info;

	public:
		Remapping_grid_domain_decomp_engine(int, int, Remap_grid_class*, Remap_grid_class*, double, double, double, double);
		~Remapping_grid_domain_decomp_engine();
		bool is_north_pole_subdomain(int subdomain_ID) { return subdomain_ID < num_northPole_subdomain; } 
		bool is_south_pole_subdomain(int subdomain_ID) { return subdomain_ID == num_northPole_subdomain && num_southPole_subdomain > 0; }
		int get_subdomain_owner_proc_id(int subdomain_ID) { return subdomain_ID % num_total_procs; }
		int calculate_grid_point_subdomain_ID(double, double);
		Decomp_grid_info *generate_subdomains_decomp_grid(Remap_grid_class *, std::vector<std::pair<int, Remap_grid_class*> > &);
		void generate_grid_one_sided_comm_info(std::vector<std::pair<int, Remap_grid_class*> >&, char**, MPI_Win*);
		Remap_grid_class *get_subdomain_grid(int, std::vector<std::pair<int, Remap_grid_class*> >&, bool);
		void get_subdomain_halo_subdomains_IDs(Remap_grid_class*, int, int, std::vector<int> &);
		void expand_subdomain_halo_grid(int, Remap_grid_class*, std::vector<std::pair<int, Remap_grid_class*> >&, bool);
		int get_subdomain_ID_from_XY_index(int, int);
		int calculate_two_subdomains_distance(Remap_grid_class*, int, int);
		void get_subdomain_XY_index(int, int &, int &);
		void calculate_subdomain_halo_bounds(Remap_grid_class*, int, int, std::pair<double, double>&, std::pair<double, double>&);
		void generate_subdomain_expanded_CoR_grid_with_triangulation(int, std::vector<std::pair<int, Remap_grid_class*> >&, Remap_grid_class*, bool);
		bool is_point_in_a_domain(double, double, double, double, double, double);
		void *generate_expanded_grid_field_data(int, std::vector<Remap_grid_class*> &, std::vector<const void*> &, bool *, int, std::vector<int> &);
        unsigned long cal_triangles_checksum_on_boundary(Triangle_inline*, int, PatCC_Point, PatCC_Point, int, double = 0);
		void get_halo_subdomain_triangulation_voronoi(Remap_grid_class*, int, bool, Comp_comm_group_mgt_node*);
		void check_boundary_triangles_checksum(Remap_grid_class*, int, int, Triangle_inline*, int);
		void check_halo_subdomain_triangulation_checksum(Remap_grid_class*, int, bool, Comp_comm_group_mgt_node*);
		void confirm_or_generate_grid_vertexes(Remap_grid_class*, std::vector<std::pair<int, Remap_grid_class*> >&, Comp_comm_group_mgt_node*, bool);
		void destroy_grid_one_sided_comm_data_buf(bool);
		void destroy_remote_subdomain_CoR_grids();		
		void initialize_subdomains_halo_info(Remap_grid_class *, std::vector<std::pair<int, Remap_grid_class*> >&, Comp_comm_group_mgt_node *);
		void clear_subdomains_halo_info();
		int get_num_subdomains() { return src_subdomains_CoR_grids.size(); }
		Remap_grid_class *get_src_expanded_subdomain(int);
		Remap_grid_class *get_dst_subdomain(int);	
		void get_src_current_expanded_subdomain_boundaries(int, double &, double &, double &, double &);
		Decomp_grid_info *get_src_subdomain_decomp_grid() { return src_subdomain_decomp_grid; }		
		Remap_grid_class *expand_src_subdomain_halo_grid(int);
};


class Distributed_H2D_grid_mgt
{
	private:
		std::vector<Distributed_H2D_grid_engine	*> distributed_H2D_grids;

	public:
		Distributed_H2D_grid_mgt() {}
		~Distributed_H2D_grid_mgt();
		Distributed_H2D_grid_engine	*search_distributed_H2D_grid(const char *, Original_grid_info *);
		Distributed_H2D_grid_engine	*search_distributed_H2D_grid(const char *, Remap_grid_class *);
		Distributed_H2D_grid_engine	*generate_distributed_H2D_grid_engine(const char *, const char *, Remap_grid_class *, long, int, int, Original_grid_info *, bool);
		void calculate_common_grid_domain_for_remapping(int, Remap_grid_class *, Remap_grid_class *, double &, double &, double &, double &, bool);
		void generate_vertexes_via_distributed_triangulation(Distributed_H2D_grid_engine *);
		void generate_full_grid_data_from_distributed_grid(int, Remap_grid_class *, Remap_grid_class *, bool);
};

#endif

