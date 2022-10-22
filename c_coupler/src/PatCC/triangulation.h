/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef __PATCC_DELAUNAY_VORONOI_H__
#define __PATCC_DELAUNAY_VORONOI_H__


#include <vector>
#include <map>
#include <list>
#include <cmath>
#include <iostream>
#include "patcc_common_utils.h"
#include "memory_pool.h"
#include "triangle.h"

#ifdef UNITTEST
#include "gtest/gtest_prod.h"
#endif

#define PDLN_CHECKSUM_FALSE (0xFFFFFFFF)

#define PDLN_UP     3
#define PDLN_LEFT   0
#define PDLN_DOWN   1
#define PDLN_RIGHT  2

#define PDLN_EXPECTED_EXPANDING_TIMES (3)

#define PAT_NUM_LOCAL_VPOINTS (4)
#define PAT_CYCLIC_EDGE_THRESHOLD (180.)

using std::vector;
using std::pair;

void sort_points_in_triangle(Triangle_inline*, int);
void sort_triangles(Triangle_inline*, int);

bool is_triangle_intersecting_with_segment(Triangle_inline*, PatCC_Point, PatCC_Point, double);
unsigned long hash_triangle_by_id(Triangle_inline);
bool on_a_line_projected(const PatCC_Point*, const PatCC_Point*, const PatCC_Point*);

struct Bound;

struct PatCC_Cell {
    PatCC_Point center;
    vector<double> vertex_lons;
    vector<double> vertex_lats;
};

class PatCC_Delaunay_Voronoi
{
    public:
		PatCC_Delaunay_Voronoi(int, double, double, double, double, double*, double*, bool *, bool, int *, double, double, double, double, double, double, double, double, long, int);
        ~PatCC_Delaunay_Voronoi();

        void add_points(double**, double**, const bool*, int, int&);
        void triangulate();
        bool try_fast_triangulate(double, double, double, double);
        void set_regional(bool);
        bool get_is_spheric_grid() { return is_spheric_grid; }
        void get_triangles_in_region(double, double, double, double, Triangle_inline *, int *, int);
        void add_point_in_rectangle(std::pair<double, double>*, double**, double**, int&, double, double, double, double, int, int, int&, bool);
        void extend_coordinate_buffer(double**, double**, int&);
        void make_final_triangle();
        void delete_irregular_triangles_recursively(PatCC_Triangle*);
        void delete_irregular_triangles();

		void generate_all_result_triangles();
        void generate_Voronoi_diagram(int &, double **, double **);
        void load_polars_info();
        bool on_a_line_lon_lat(const PatCC_Point*, const PatCC_Point*, const PatCC_Point*);
        void convert_cyclic_triangle_point(double*, int);
		void check_and_add_vertex_to_map(int, double, double, char **, long &, long &, int *);
        void delete_point_on_common_line(double*, double*, long*, int&);
        void delete_overlap_point_in_cell_vertexes(double*, double*, long*, int&);
		bool is_real_point(int);
		void rotate_sphere_coordinate(double*, double*, int, double, double, int);
		void build_check_twin_edge_relationship_between_triangles(PatCC_Triangle*, PatCC_Triangle*, bool);

        /* distributed H2D grid */
        int num_result_triangles; 
        int origin_num_points;
        Triangle_inline* result_triangles;
		double *origin_lon_values, *origin_lat_values;
        int num_polar_point_in_origin_coordinate;
		std::vector<std::pair<int, int> > overlapping_points;
		std::vector<int> candicate_points_indexes;
		bool do_projection;

#ifdef OPENCV
        void plot_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_projection_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_original_points_into_file(const char*, double min_x=0.0, double max_x=0.0, double min_y=0.0, double max_y=0.0);
        void plot_current_step_into_file(const char*);
#endif
    private:

        unsigned triangulating_process(PatCC_Triangle*, unsigned);
        void fast_triangulate(int, int, bool);
        void map_buffer_index_to_point_index();
        void push(unsigned *, PatCC_Triangle*);
		bool are_two_points_the_same(int, int);
		void add_redistribute_head_tail(PatCC_Triangle *);

        /* preparing function */
        void initialize(int);
        void extend_points_buffer(int);
        void distribute_initial_points(const double*, const double*, int, int**);
        int enlarge_super_rectangle(double**, double**, int, int&);
        Bound* make_bounding_box();
        bool point_in_triangle(double, double, PatCC_Triangle*);
        bool point_in_bound(double, double, Bound*);

        void distribute_points_into_triangles(int, int, unsigned, unsigned);
        void link_remained_list(unsigned, unsigned, int*, int*);
        void swap_points(int, int);

        void mark_special_triangles();
        bool check_uniqueness(int, const PatCC_Edge*);
        int  get_lowest_point_of_four(const PatCC_Edge*);
        PatCC_Edge* make_twins_edge(PatCC_Edge*);

        bool is_edge_legal(int, const PatCC_Edge *PatCC_Edge);
        bool is_triangle_legal(const PatCC_Triangle *);
		int circum_circle_contains_reliably_3D(const PatCC_Edge *, PatCC_Point *);
        int  circum_circle_contains_reliably(const PatCC_Edge*, PatCC_Point*, double);
		int  circum_circle_contains_reliably_new(const PatCC_Edge*, PatCC_Point*, double);
        int  get_index_in_array(const PatCC_Point*);

        void legalize_triangles(int, PatCC_Edge *PatCC_Edge, unsigned, unsigned*);

        PatCC_Edge* allocate_edge(int, int);
        PatCC_Triangle* allocate_triangle(PatCC_Edge*, PatCC_Edge*, PatCC_Edge*, bool = false);
        PatCC_Triangle* allocate_triangle(int, int, int, bool = false);
        void initialize_triangle_with_edges(PatCC_Triangle*, PatCC_Edge*, PatCC_Edge*, PatCC_Edge*, bool = false);
        void initialize_edge(PatCC_Edge* e, int head, int tail);

        inline void ref_inc(PatCC_Edge* e) {e->ref_count++;};
        inline void ref_dec(PatCC_Edge* e) {
            e->ref_count--;
            if(e->ref_count <= 0) {
                if (e->twin_edge)
                    e->twin_edge->twin_edge = NULL;
                edge_allocator->deleteElement(e);
				num_allocate_edges --;
            }
        };
        void clean_triangle(PatCC_Triangle*);

        inline PatCC_Point* vertex(const PatCC_Triangle* t, int i);
        inline PatCC_Point* head(const PatCC_Edge* e);
        inline PatCC_Point* tail(const PatCC_Edge* e);
		inline void diag_triangle(PatCC_Triangle *);

        /* Storage */
        PatCC_Point*            all_points;
        vector<PatCC_Triangle*> all_leaf_triangles;
        int               max_points;

        /* Memory management */
        Triangle_pool *triangle_allocator;
        Edge_pool *edge_allocator;
        PatCC_Triangle** triangle_stack;
        unsigned   stack_size;

        /* Property */
        bool   is_regional;
        bool   polar_mode;
        bool   fast_mode;
        double tolerance;

        /* Grid info */
        int num_points;
        int vpolar_local_index;
        double *x_ref;
        double *y_ref;
        int *global_index;
        double original_lon_center;
		double original_lat_center;
		double original_min_lon;
		double original_max_lon;
		double original_min_lat;
		double original_max_lat;
		double grid_min_lon;
		double grid_max_lon;
		double grid_min_lat;
		double grid_max_lat;
		double kernel_min_lon;
        double kernel_max_lon;
        double kernel_min_lat;
        double kernel_max_lat;
		long kernel_grid_size;
		
		bool is_spheric_grid;

        /* Triangulating stuff */
        PatCC_Point *virtual_point[4];
        int *point_idx_to_buf_idx;
        unsigned dirty_triangles_count;

        /* Consistency checking boundary */
        PatCC_Point  avoiding_line_head[2];
        PatCC_Point  avoiding_line_tail[2];
        PatCC_Point  avoiding_circle_center[2];
        double avoiding_circle_radius[2];
		int *redistribute_head_tails;
		long redistribute_head_tails_max_size;
		long redistribute_head_tails_content_size;
		int num_allocate_triangles;
		int num_allocate_edges;
};

#ifdef OPENCV
void plot_triangles_into_file(const char *filename, Triangle_inline *t, int num, bool plot_cyclic_triangles=true);
void plot_triangles_into_file(const char *filename, std::vector<PatCC_Triangle*>, PatCC_Point*);
#endif
void save_triangles_info_file(const char *filename, Triangle_inline *t, int num);

#endif
