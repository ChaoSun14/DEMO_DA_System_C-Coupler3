/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef PDLN_TRIANGLE_H
#define PDLN_TRIANGLE_H

#include "patcc_common_utils.h"
#ifdef OPENCV
#include "opencv2/opencv.hpp"
#endif



class PatCC_Triangle;
class Triangle_inline;

class PatCC_Point
{
    public:
        double x;
        double y;
        double xx;
        double yy;
		double zz;
        int    id;
        int    mask:1;
        int    next;
        int    prev;
		bool is_spheric_grid;

        PatCC_Point();
        PatCC_Point(double, double, bool);
        PatCC_Point(double, double, int, bool, int = -1, int = -1);
        PatCC_Point(double, double, int, bool, bool, int = -1, int = -1);
        PatCC_Point(PAT_REAL, PAT_REAL, int, bool, bool, int = -1, int = -1);
        ~PatCC_Point();
        double calculate_distance(const PatCC_Point*) const;
        double calculate_distance_MD(double, double) const;
        double calculate_distance(double, double, double, double) const;
		int position_to_edge(const PatCC_Point*, const PatCC_Point*) const;
        int position_to_triangle(const PatCC_Point*, const PatCC_Point*, const PatCC_Point*, bool check_on_a_line = false) const;
        int position_to_triangle(const Triangle_inline*) const;
        int is_in_region(double min_x, double max_x, double min_y, double max_y) const;
		bool is_the_same_with_another(const PatCC_Point *) const;
		void transform_coord_from_sphere_to_3D();
};


class PatCC_Edge
{
    private:
        int   head;
        int   tail;                      /* the tail of this PatCC_Edge, constant */
        PatCC_Edge* twin_edge;                 /* the twin_edge PatCC_Edge, whose tail is the head of this PatCC_Edge and head is the tail of this PatCC_Edge */
        PatCC_Edge* next_edge_in_triangle;     /* the next_edge_in_triangle PatCC_Edge, whose tail is the head of this PatCC_Edge but head isn't the tail of this PatCC_Edge */
        PatCC_Edge* prev_edge_in_triangle;     /* the prev_edge_in_triangle PatCC_Edge, whose head is the tail of this PatCC_Edge but tail isn't the head of this PatCC_Edge */
        int   ref_count;
        PatCC_Triangle* triangle;               /* the PatCC_Triangle which is composed by this PatCC_Edge and its next_edge_in_triangle and prev_edge_in_triangle */

    public:
        PatCC_Edge();
        ~PatCC_Edge();
        PatCC_Edge *generate_twins_edge();

        inline void ref_inc() {ref_count++;};
#ifdef OPENCV
        friend void draw_line(cv::Mat, PatCC_Point*, PatCC_Edge*, double, double, double, double, cv::Scalar);
#endif
        friend class PatCC_Triangle;
        friend class PatCC_Delaunay_Voronoi;
};


class PatCC_Triangle
{
    private:
        PatCC_Edge*    edge[3];
        double circum_center[3];
        double circum_radius2;
        int      v[3];    /* index of vertexes */
        int      remained_points_head;
        int      remained_points_tail;
        int      stack_ref_count;
        bool is_leaf;
        bool is_cyclic;
        bool is_virtual;

        int  circum_circle_position_to(PatCC_Point*, double=PDLN_ABS_TOLERANCE);

    public:
        PatCC_Triangle();
        ~PatCC_Triangle();
        int find_best_candidate_point(PatCC_Point*, std::vector<int> &);
        bool contain_vertex(int);
        void calulate_circum_circle(const PatCC_Point*, const PatCC_Point*, const PatCC_Point*);
		void calulate_circum_circle_2D(const PatCC_Point*, const PatCC_Point*, const PatCC_Point*);
        void calulate_circum_circle_3D(const PatCC_Point*, const PatCC_Point*, const PatCC_Point*);
        int find_dividing_point(PatCC_Point*);
        void set_remained_points(int, int);
        int pop_tail(PatCC_Point*, int);
		PAT_REAL get_circum_radius2() { return circum_radius2; }

        friend class PatCC_Delaunay_Voronoi;
        friend class PatCC_Point;
        friend void  plot_triangles_into_file(const char *filename, std::vector<PatCC_Triangle*>, PatCC_Point*);
        friend class Triangle_pool;
};


class Triangle_inline
{
    public:
        PatCC_Point v[3];
        bool is_cyclic;
        Triangle_inline() {};
        Triangle_inline(PatCC_Point, PatCC_Point, PatCC_Point, bool = false);
        void check_cyclic();
        friend bool operator == (Triangle_inline, Triangle_inline);
};

#endif
