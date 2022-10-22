/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu, 
  *  Haoyu Yang and Hao Yu. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "mpi.h"
#include "triangulation.h"
#include "patcc_common_utils.h"
#include "merge_sort.h"
#include "coordinate_hash.h"
#include "global_data.h"
#include "projection.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sys/time.h>
#include <tr1/unordered_map>
#include <list>
#include <utility>

#include "remap_common_utils.h"
#include "remap_utils_nearest_points.h"

#include <vector>
#ifdef OPENCV
#include "opencv_utils.h"
#endif


double calculate_angle(double x1, double y1, double x2, double y2, double x3, double y3) 
{
	double theta = atan2(x3 - x1, y3 - y1) - atan2(x2 - x1, y2 - y1);
	double result = theta * 180.0 / M_PI;
	if (result < 0)
		result += 360;
	return result;
}
			

/*
 *           o Center
 *           ^
 *          /
 *         /
 *     V1 o -----> o V2
 */
double compute_three_2D_points_cross_product(double center_x, double center_y,
											 double v1_x, double v1_y,
											 double v2_x, double v2_y)
{
	double delta1_x = v2_x - v1_x;
	double delta1_y = v2_y - v1_y;
	double delta2_x = center_x - v1_x;
	double delta2_y = center_y - v1_y;

	return delta1_x*delta2_y - delta2_x*delta1_y;
}


double det(const PatCC_Point *pt1, const PatCC_Point *pt2, const PatCC_Point *pt3)
{
	if (pt1->is_spheric_grid)
		return compute_three_3D_points_cross_product(pt3->xx, pt3->yy, pt3->zz, pt1->xx, pt1->yy, pt1->zz, pt2->xx, pt2->yy, pt2->zz);
	else return compute_three_2D_points_cross_product(pt3->x, pt3->y, pt1->x, pt1->y, pt2->x, pt2->y);
}


static inline void swap(PatCC_Point *p1, PatCC_Point *p2)
{
	PatCC_Point tmp = *p1;
	*p1 = *p2;
	*p2 = tmp;
}


static int compare_v2(const void* a, const void* b)
{
	Triangle_inline t1 = *(const Triangle_inline*)a;
	Triangle_inline t2 = *(const Triangle_inline*)b;

	if(t1.v[2].id < t2.v[2].id) return -1;
	if(t1.v[2].id > t2.v[2].id) return  1;
	return 0;
}


static int compare_v1(const void* a, const void* b)
{
	Triangle_inline t1 = *(const Triangle_inline*)a;
	Triangle_inline t2 = *(const Triangle_inline*)b;

	if(t1.v[1].id < t2.v[1].id) return -1;
	if(t1.v[1].id > t2.v[1].id) return  1;
	return 0;
}


static int compare_v0(const void* a, const void* b)
{
	Triangle_inline t1 = *(const Triangle_inline*)a;
	Triangle_inline t2 = *(const Triangle_inline*)b;

	if(t1.v[0].id < t2.v[0].id) return -1;
	if(t1.v[0].id > t2.v[0].id) return  1;
	if(t1.v[1].id < t2.v[1].id) return -1;
	if(t1.v[1].id > t2.v[1].id) return  1;	
	if(t1.v[2].id < t2.v[2].id) return -1;
	if(t1.v[2].id > t2.v[2].id) return  1;

	return 0;
}


static inline void radix_sort(Triangle_inline *triangles, int num_triangles)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, sizeof(Triangle_inline) > sizeof(void *)/2, "Software error in radix_sort");
	merge_sort(triangles, num_triangles, sizeof(Triangle_inline), compare_v2);
	merge_sort(triangles, num_triangles, sizeof(Triangle_inline), compare_v1);
	merge_sort(triangles, num_triangles, sizeof(Triangle_inline), compare_v0);
}


void sort_points_in_triangle(Triangle_inline *triangles, int num_triangles)
{
	for(int i = 0; i < num_triangles; i++) {
		if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
		if(triangles[i].v[1].id > triangles[i].v[2].id) swap(&triangles[i].v[1], &triangles[i].v[2]);
		if(triangles[i].v[0].id > triangles[i].v[1].id) swap(&triangles[i].v[0], &triangles[i].v[1]);
	}
}


inline void sort_points_in_triangle(Triangle_inline& triangle)
{
	if(triangle.v[0].id > triangle.v[1].id) std::swap(triangle.v[0], triangle.v[1]);
	if(triangle.v[1].id > triangle.v[2].id) std::swap(triangle.v[1], triangle.v[2]);
	if(triangle.v[0].id > triangle.v[1].id) std::swap(triangle.v[0], triangle.v[1]);
}


void sort_triangles(Triangle_inline *triangles, int num_triangles)
{
	radix_sort(triangles, num_triangles);
}


PatCC_Point::PatCC_Point()
{
	is_spheric_grid = false;
}


PatCC_Point::PatCC_Point(double x, double y, bool is_spheric)
	: x(x)
	, y(y)
	, mask(true)
	, is_spheric_grid(is_spheric)
{
	transform_coord_from_sphere_to_3D();
}


PatCC_Point::PatCC_Point(double x, double y, int id, bool is_spheric, int fd, int bk)
	: x(x)
	, y(y)
	, id(id)
	, mask(true)
	, next(fd)
	, prev(bk)
	, is_spheric_grid(is_spheric)
{
	transform_coord_from_sphere_to_3D();
}


PatCC_Point::PatCC_Point(double x, double y, int id, bool msk, bool is_spheric, int fd, int bk)
	: x(x)
	, y(y)
	, id(id)
	, mask(msk)
	, next(fd)
	, prev(bk)
	, is_spheric_grid(is_spheric)
{
	transform_coord_from_sphere_to_3D();
}


double PatCC_Point::calculate_distance(const PatCC_Point *pt) const
{
	double dx = pt->x - x;
	double dy = pt->y - y;
	return sqrt(dx * dx + dy * dy);
}


double PatCC_Point::calculate_distance_MD(double pt_x, double pt_y) const
{
	return calculate_distance_of_two_points_2D(x, y, pt_x, pt_y, is_spheric_grid);
}


static inline double calculate_distance(double x0, double y0, double x1, double y1)
{
	double dx = x0 - x1;
	double dy = y0 - y1;
	return sqrt(dx * dx + dy * dy);
}


PatCC_Point::~PatCC_Point()
{
}


PatCC_Point operator - (PatCC_Point p1, PatCC_Point p2)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "software error in PatCC_Point operator -");
	return PatCC_Point(p1.x - p2.x, p1.y - p2.y, -1, p1.is_spheric_grid);
}


bool operator == (PatCC_Point p1, PatCC_Point p2)
{
	if(p1.id != -1 && p2.id != -1)
		return p1.id == p2.id;

	if(float_eq(p1.x, p2.x) && float_eq(p1.y, p2.y))
		return true;
	else if(current_triangulation->get_is_spheric_grid() && float_eq(fabs(p1.x - p2.x), 360.0) && float_eq(p1.y, p2.y))
		return true;
	else
		return false;
}


bool operator != (PatCC_Point p1, PatCC_Point p2)
{
	return !(p1 == p2);
}


/**
 * Check point's position relative to an PatCC_Edge<pt1, pt2>
 * Points should be distinct
 * @param  pt1    the head of the PatCC_Edge
 * @param  pt2    the head of the PatCC_Edge 
 * @return    1    left
 *            0    on
 *            -1    right
 */
int PatCC_Point::position_to_edge(const PatCC_Point *pt1, const PatCC_Point *pt2) const
{
	double res;
	int result1, result2;

    res = det(pt1, pt2, this);

	if (float_eq_hi(res, 0))
		return 0;
	else if (res > 0)
		return 1;
	else
		return -1;
}


/**
 * Check point's position relative to a triangle
 * This point and points of the triangle should be distinct
 * @return     0    inside
 *            -1    outside
 *             1    lies on the PatCC_Edge <pt1, pt2>
 *             2    lies on the PatCC_Edge <pt2, pt3>
 *             3    lies on the PatCC_Edge <pt3, pt1>
 */
int PatCC_Point::position_to_triangle(const PatCC_Point* v0, const PatCC_Point* v1, const PatCC_Point* v2, bool check_on_a_line) const
{
	if (is_the_same_with_another(v0) || is_the_same_with_another(v1) || is_the_same_with_another(v2)) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in PatCC_Point::position_to_triangle");
		return 4;
	}

	if (report_error_enabled) {
		bool on1 = position_to_edge(v0, v1) == 0;
		bool on2 = position_to_edge(v1, v2) == 0;
		bool on3 = position_to_edge(v2, v0) == 0;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(on1 && on2), "Software error in PatCC_Point::position_to_triangle");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(on2 && on3), "Software error in PatCC_Point::position_to_triangle");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(on3 && on1), "Software error in PatCC_Point::position_to_triangle");
	}

	if (check_on_a_line && false) {
		if (on_a_line_projected(this, v0, v1) || on_a_line_projected(this, v1, v2) || on_a_line_projected(this, v0, v2))
			return 0;
	}

	int ret = 0;
	int pos = position_to_edge(v0, v1);
	if (pos == -1)
		return -1;
	else if (pos == 0)
		ret = 1;
	pos = position_to_edge(v1, v2);
	if (pos == -1)
		return -1;
	else if (pos == 0)
		ret = 2;
	pos = position_to_edge(v2, v0);
	if (pos == -1)
		return -1;
	else if (pos == 0)
		ret = 3;
	return ret;
}


int PatCC_Point::position_to_triangle(const Triangle_inline *triangle) const
{
	if (is_the_same_with_another(&(triangle->v[0])) || is_the_same_with_another(&(triangle->v[1])) || is_the_same_with_another(&(triangle->v[2])))
		return 4;

	if (report_error_enabled) {
		bool on1 = position_to_edge(&triangle->v[0], &triangle->v[1]) == 0;
		bool on2 = position_to_edge(&triangle->v[1], &triangle->v[2]) == 0;
		bool on3 = position_to_edge(&triangle->v[2], &triangle->v[0]) == 0;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(on1 && on2), "Software error in PatCC_Point::position_to_triangle: =(%lf, %lf)  (%lf, %lf)  (%lf, %lf)  (%lf, %lf, %lf)  ", triangle->v[0].x, triangle->v[0].y, triangle->v[1].x, triangle->v[1].y, triangle->v[2].x, triangle->v[2].y, triangle->v[0].xx, triangle->v[0].yy, triangle->v[0].zz);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(on2 && on3), "Software error in PatCC_Point::position_to_triangle");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(on3 && on1), "Software error in PatCC_Point::position_to_triangle");
	}

	int ret = 0;
	int pos = position_to_edge(&triangle->v[0], &triangle->v[1]);
	if (pos == -1)
		return -1;
	else if (pos == 0)
		ret = 1;
	pos = position_to_edge(&triangle->v[1], &triangle->v[2]);
	if (pos == -1)
		return -1;
	else if (pos == 0)
		ret = 2;
	pos = position_to_edge(&triangle->v[2], &triangle->v[0]);
	if (pos == -1)
		return -1;
	else if (pos == 0)
		ret = 3;
	return ret;
}


inline int PatCC_Point::is_in_region(double min_x, double max_x, double min_y, double max_y) const
{
	return x > min_x && x < max_x && y > min_y && y < max_y;
}

PatCC_Edge::PatCC_Edge()
{
}


PatCC_Edge::~PatCC_Edge()
{
}


PatCC_Edge* PatCC_Delaunay_Voronoi::make_twins_edge(PatCC_Edge *e)
{
	PatCC_Edge *twins_edge = allocate_edge(e->tail, e->head);
	twins_edge->twin_edge = e;
	e->twin_edge = twins_edge;

	return twins_edge;
}


int PatCC_Delaunay_Voronoi::get_lowest_point_of_four(const PatCC_Edge* shared_edge)
{
	PatCC_Point* points[4];

	points[0] = &all_points[shared_edge->prev_edge_in_triangle->head];
	points[1] = &all_points[shared_edge->head];
	points[2] = &all_points[shared_edge->tail];
	points[3] = &all_points[shared_edge->twin_edge->prev_edge_in_triangle->head];

	for (int i = 0; i < 4; i++)
		if (points[i]->id == -1 || global_index[points[i]->id] < 0)
			return get_index_in_array(points[0]);

	double x_fixed[4], y_fixed[4];

	for (int i = 0; i < 4; i++) {
		x_fixed[i] = x_ref[points[i]->id];
		y_fixed[i] = y_ref[points[i]->id];
	}

	if (polar_mode && (shared_edge->triangle->is_cyclic || (shared_edge->twin_edge != NULL && shared_edge->twin_edge->triangle->is_cyclic))) {
		for (int i = 0; i < 4; i++)
			if (x_fixed[i] > 180) x_fixed[i] -= 360;
	}

	int lowest = 0;

	for (int i = 1; i < 4; i++)
		if (x_fixed[i] <  x_fixed[lowest] ||
		   (x_fixed[i] == x_fixed[lowest] && y_fixed[i] < y_fixed[lowest]) ) {
			lowest = i;
		}
	return get_index_in_array(points[lowest]);
}


int PatCC_Delaunay_Voronoi::get_index_in_array(const PatCC_Point* p)
{
	return p - all_points;
}


bool PatCC_Delaunay_Voronoi::is_edge_legal(int p_idx, const PatCC_Edge *PatCC_Edge)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, p_idx == PatCC_Edge->prev_edge_in_triangle->head, "Software error in check_uniqueness: p_idx=%d, PatCC_Edge->prev_edge_in_triangle->head=%d", p_idx, PatCC_Edge->prev_edge_in_triangle->head);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, PatCC_Edge->triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::is_edge_legal");
	if (!PatCC_Edge->twin_edge) {
		return true;
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, PatCC_Edge->head == PatCC_Edge->twin_edge->tail && PatCC_Edge->tail == PatCC_Edge->twin_edge->head, "Software error in PatCC_Delaunay_Voronoi::is_edge_legal");

	if(!PatCC_Edge->twin_edge->triangle) {
		return true;
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, PatCC_Edge->twin_edge->triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::is_edge_legal");

	int ret = circum_circle_contains_reliably(PatCC_Edge, head(PatCC_Edge->twin_edge->prev_edge_in_triangle), tolerance);

	if (ret == -1) {
		return true;
	}
	
	if (ret == 0) {
		return check_uniqueness(p_idx, PatCC_Edge);
	}

	return false;
}

bool PatCC_Delaunay_Voronoi::check_uniqueness(int p_idx, const PatCC_Edge *PatCC_Edge)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, p_idx == PatCC_Edge->prev_edge_in_triangle->head, "Software error in check_uniqueness: p_idx=%d, PatCC_Edge->prev_edge_in_triangle->head=%d", p_idx, PatCC_Edge->prev_edge_in_triangle->head);

	int lowest = get_lowest_point_of_four(PatCC_Edge);
	bool is_lowest = (p_idx == lowest || PatCC_Edge->twin_edge->prev_edge_in_triangle->head == lowest);

	return is_lowest;
}


bool PatCC_Delaunay_Voronoi::is_triangle_legal(const PatCC_Triangle *t)
{
	for(int i = 0; i < 3; i++) {
		if(!is_edge_legal(t->edge[i]->prev_edge_in_triangle->head, t->edge[i])) {
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "[%d] +illegal triangle: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", 1, vertex(t, 0)->x, vertex(t, 0)->y, vertex(t, 1)->x, vertex(t, 1)->y, vertex(t, 2)->x, vertex(t, 2)->y);
			PatCC_Triangle *tt = t->edge[i]->twin_edge->triangle;
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "[%d] -illegal triangle: (%lf, %lf), (%lf, %lf), (%lf, %lf)\n", 1, vertex(tt, 0)->x, vertex(tt, 0)->y, vertex(tt, 1)->x, vertex(tt, 1)->y, vertex(tt, 2)->x, vertex(tt, 2)->y);
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "[%d] +: %d, -: %d\n", 1, t->is_leaf, tt->is_leaf);
			return false;
		}
	}
	return true;
}

/*
 *
 *     Vi________ Vk
 *      /\      /
 *     /  \ 2  /
 *    / 1  \  /
 *   /______\/
 *  Vr      Vj
 *
 *  1: triangle
 *  2: twin triangle
 */
void PatCC_Delaunay_Voronoi::legalize_triangles(int vr_idx, PatCC_Edge *edge, unsigned stack_base, unsigned *stack_top)
{
	if (is_edge_legal(vr_idx, edge))
		return;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, edge->triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::legalize_triangles");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, edge->twin_edge->triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::legalize_triangles");

	int vk_idx = edge->twin_edge->prev_edge_in_triangle->head;
	PatCC_Edge *eij = edge;
	PatCC_Edge *ejr = eij->next_edge_in_triangle;
	PatCC_Edge *eri = ejr->next_edge_in_triangle;
	PatCC_Edge *eji = eij->twin_edge;
	PatCC_Edge *eik = eji->next_edge_in_triangle;
	PatCC_Edge *ekj = eik->next_edge_in_triangle;
	PatCC_Edge *erk = allocate_edge(vr_idx, vk_idx);
	PatCC_Edge *ekr = make_twins_edge(erk);
	
//	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "the another point is (%lf, %lf)(%lf, %lf)", all_points[edge->head].x, all_points[edge->head].y, all_points[edge->tail].x, all_points[edge->tail].y);
	PatCC_Triangle* tikr = allocate_triangle(eik, ekr, eri);
	PatCC_Triangle* tjrk = allocate_triangle(ejr, erk, ekj);
	push(stack_top, tikr);
	int stack_pos_mark = *stack_top;
	push(stack_top, tjrk);
	add_redistribute_head_tail(edge->triangle);
	add_redistribute_head_tail(edge->twin_edge->triangle);

	if (tikr->stack_ref_count == stack_pos_mark) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tikr->is_leaf, "Software error in PatCC_Delaunay_Voronoi::legalize_triangles");
		legalize_triangles(vr_idx, eik, stack_base, stack_top);
	}
	if (tikr->stack_ref_count == stack_pos_mark) { 
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tikr->is_leaf, "Software error in PatCC_Delaunay_Voronoi::legalize_triangles");
		legalize_triangles(vk_idx, eri, stack_base, stack_top);
	}
	if (tjrk->stack_ref_count == stack_pos_mark+1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tjrk->is_leaf, "Software error in PatCC_Delaunay_Voronoi::legalize_triangles");
		legalize_triangles(vr_idx, ekj, stack_base, stack_top);
	}
	if (tjrk->stack_ref_count == stack_pos_mark+1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tjrk->is_leaf, "Software error in PatCC_Delaunay_Voronoi::legalize_triangles");
		legalize_triangles(vk_idx, ejr, stack_base, stack_top); 
	}
}


PatCC_Triangle::PatCC_Triangle()
	: stack_ref_count(-1)
{
	edge[0] = NULL;
	edge[1] = NULL;
	edge[2] = NULL;
}


PatCC_Triangle::~PatCC_Triangle()
{
}


void PatCC_Triangle::calulate_circum_circle_3D(const PatCC_Point *a, const PatCC_Point *b, const PatCC_Point *c) 
{
    double bax, bay, baz, cax, cay, caz;
    double balength, calength;
    double crossbcx, crossbcy, crossbcz;
    double denominator;

    bax = b->xx - a->xx;
    bay = b->yy - a->yy;
    baz = b->zz - a->zz;
    cax = c->xx - a->xx;
    cay = c->yy - a->yy;
    caz = c->zz - a->zz;

    balength = bax * bax + bay * bay + baz * baz;
    calength = cax * cax + cay * cay + caz * caz;

    crossbcx = bay * caz - cay * baz;
    crossbcy = baz * cax - caz * bax;
    crossbcz = bax * cay - cax * bay;

    denominator = 2 * (crossbcx * crossbcx + crossbcy * crossbcy + crossbcz * crossbcz);
    circum_center[0] = a->xx + ((balength * cay - calength * bay) * crossbcz - (balength * caz - calength * baz) * crossbcy) / denominator;
    circum_center[1] = a->yy + ((balength * caz - calength * baz) * crossbcx - (balength * cax - calength * bax) * crossbcz) / denominator;
    circum_center[2] = a->zz + ((balength * cax - calength * bax) * crossbcy - (balength * cay - calength * bay) * crossbcx) / denominator;
}



void PatCC_Triangle::calulate_circum_circle(const PatCC_Point* v0, const PatCC_Point* v1, const PatCC_Point* v2)
{
	PAT_REAL ab, cd, ef;

	if (v0->is_spheric_grid)
		calulate_circum_circle_3D(v0, v1, v2);
	else calulate_circum_circle_2D(v0, v1, v2);
}


void PatCC_Triangle::calulate_circum_circle_2D(const PatCC_Point* v0, const PatCC_Point* v1, const PatCC_Point* v2)
{
	PAT_REAL ab, cd, ef;

	ab = ((PAT_REAL)v0->x * (PAT_REAL)v0->x) + ((PAT_REAL)v0->y * (PAT_REAL)v0->y);
	cd = ((PAT_REAL)v1->x * (PAT_REAL)v1->x) + ((PAT_REAL)v1->y * (PAT_REAL)v1->y);
	ef = ((PAT_REAL)v2->x * (PAT_REAL)v2->x) + ((PAT_REAL)v2->y * (PAT_REAL)v2->y);

	circum_center[0] = (ab * ((PAT_REAL)v2->y - (PAT_REAL)v1->y) + cd * ((PAT_REAL)v0->y - (PAT_REAL)v2->y) + ef * ((PAT_REAL)v1->y - (PAT_REAL)v0->y)) /
					   ((PAT_REAL)v0->x * ((PAT_REAL)v2->y - (PAT_REAL)v1->y) + (PAT_REAL)v1->x * ((PAT_REAL)v0->y - (PAT_REAL)v2->y) + (PAT_REAL)v2->x * ((PAT_REAL)v1->y - (PAT_REAL)v0->y)) * 0.5;
	circum_center[1] = (ab * ((PAT_REAL)v2->x - (PAT_REAL)v1->x) + cd * ((PAT_REAL)v0->x - (PAT_REAL)v2->x) + ef * ((PAT_REAL)v1->x - (PAT_REAL)v0->x)) /
					   ((PAT_REAL)v0->y * ((PAT_REAL)v2->x - (PAT_REAL)v1->x) + (PAT_REAL)v1->y * ((PAT_REAL)v0->x - (PAT_REAL)v2->x) + (PAT_REAL)v2->y * ((PAT_REAL)v1->x - (PAT_REAL)v0->x)) * 0.5;
	circum_radius2 = ((v0->x - circum_center[0]) * (v0->x - circum_center[0])) + ((v0->y - circum_center[1]) * (v0->y - circum_center[1]));
}


void PatCC_Delaunay_Voronoi::initialize_triangle_with_edges(PatCC_Triangle* t, PatCC_Edge *edge1, PatCC_Edge *edge2, PatCC_Edge *edge3, bool force)
{
	PatCC_Point *pt1, *pt2, *pt3;

	t->is_leaf = true;
	t->is_cyclic = false;
	t->is_virtual = false;
	if(force) {
		t->v[0] = edge1->head;
		t->v[1] = edge2->head;
		t->v[2] = edge3->head;
		t->edge[0] = edge1;
		t->edge[1] = edge2;
		t->edge[2] = edge3;
		pt1 = vertex(t, 0);
		pt2 = vertex(t, 1);
		pt3 = vertex(t, 2);
	} else {
		pt1 = head(edge1);
		pt2 = head(edge2);
		pt3 = head(edge3);

		if (report_error_enabled) {
			if(float_eq_hi(det(pt1, pt2, pt3), 0) || float_eq_hi(det(pt2, pt3, pt1), 0) || float_eq_hi(det(pt3, pt1, pt2), 0)) {
				EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "(%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)", pt1->x, pt1->y, pt2->x, pt2->y, pt3->x, pt3->y);
				EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "std::fabs(det(pt1, pt2, pt3)): %.20lf std::fabs(det(pt2, pt3, pt1)): %.20lf std::fabs(det(pt3, pt1, pt2)): %.20lf", std::fabs(det(pt1, pt2, pt3)), std::fabs(det(pt2, pt3, pt1)), std::fabs(det(pt3, pt1, pt2)));
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in initialize_triangle_with_edges");
			}
		}

		/* if there are unmarked redundant points, the PDASSERTion may fail */
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !float_eq_hi(det(pt1, pt2, pt3), 0) && !float_eq_hi(det(pt2, pt3, pt1), 0) && !float_eq_hi(det(pt3, pt1, pt2), 0), "Software error in initialize_triangle_with_edges");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, edge1->tail==edge2->head && edge2->tail==edge3->head && edge3->tail==edge1->head, "Software error in initialize_triangle_with_edges");
		   
		t->v[0] = edge1->head;
		if (fast_mode || pt1->position_to_edge(pt2, pt3) == 1) {
			t->v[1] = edge2->head;
			t->v[2] = edge3->head;
			t->edge[0] = edge1;
			t->edge[1] = edge2;
			t->edge[2] = edge3;
		}
		else {
			EXECUTION_REPORT(REPORT_LOG, -1, true, "not counterclockwise (%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)\n", x_ref[pt1->id], y_ref[pt1->id], x_ref[pt2->id], y_ref[pt2->id], x_ref[pt3->id], y_ref[pt3->id]);
			EXECUTION_REPORT(REPORT_LOG, -1, true, "not counterclockwise (%.20lf, %.20lf), (%.20lf, %.20lf), (%.20lf, %.20lf)\n", pt1->x, pt1->y, pt2->x, pt2->y, pt3->x, pt3->y);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in initialize_triangle_with_edges");
			t->v[1] = edge3->head;
			t->v[2] = edge2->head;
			if (edge1->twin_edge == NULL)
				edge1->twin_edge = make_twins_edge(edge1);
			if (edge2->twin_edge == NULL)
				edge2->twin_edge = make_twins_edge(edge2);
			if (edge3->twin_edge == NULL)
				edge3->twin_edge = make_twins_edge(edge3);
			t->edge[0] = edge3->twin_edge;
			t->edge[1] = edge2->twin_edge;
			t->edge[2] = edge1->twin_edge;
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, edge3->twin_edge != NULL && edge2->twin_edge != NULL && edge1->twin_edge != NULL, "Software error in initialize_triangle_with_edges");
		}

	}

	t->remained_points_head = -1;
	t->remained_points_tail = -1;
	
	t->edge[0]->next_edge_in_triangle = t->edge[1];
	t->edge[1]->next_edge_in_triangle = t->edge[2];
	t->edge[2]->next_edge_in_triangle = t->edge[0];
	t->edge[0]->prev_edge_in_triangle = t->edge[2];
	t->edge[1]->prev_edge_in_triangle = t->edge[0];
	t->edge[2]->prev_edge_in_triangle = t->edge[1];

	t->edge[0]->triangle = t;
	t->edge[1]->triangle = t;
	t->edge[2]->triangle = t;

	t->edge[0]->ref_inc();
	t->edge[1]->ref_inc();
	t->edge[2]->ref_inc();
	t->calulate_circum_circle(pt1, pt2, pt3);

	if (!is_regional) {
		int id[3];
		for (int j = 0; j < 3; j++)
			id[j] = vertex(t, j)->id;

		if (id[0] != -1 && id[1] != -1 && id[2] != -1 && global_index[id[0]] > -1 && global_index[id[1]] > -1 && global_index[id[2]] > -1)
			for (int j = 0; j < 3; j++) {
				if (is_spheric_grid && calculate_distance(x_ref[id[j]], y_ref[id[j]], x_ref[id[(j+1)%3]], y_ref[id[(j+1)%3]]) > PAT_CYCLIC_EDGE_THRESHOLD) {
					t->is_cyclic = true;
					break;
				}
			}
	}
}


void PatCC_Delaunay_Voronoi::initialize_edge(PatCC_Edge* e, int head, int tail)
{
	e->head = head;
	e->tail = tail;
	e->twin_edge = NULL;
	e->next_edge_in_triangle = NULL;
	e->prev_edge_in_triangle = NULL;
	e->ref_count = 0;
	e->triangle = NULL;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, head >= 0 && head <= num_points && tail >= 0 && tail <= num_points, "Software error in PatCC_Delaunay_Voronoi::initialize_edge: %d %d vs %d", head, tail, num_points);
}


/*
 * Input : PatCC_Point to be checked
 * Return:  1    point is in circum circle
 *          0    point is on circum circle
 *         -1    point is out of circum circle
 */
int PatCC_Delaunay_Voronoi::circum_circle_contains_reliably_3D(const PatCC_Edge *PatCC_Edge, PatCC_Point *p)
{
	PatCC_Point pnts[4];
	PatCC_Point tmp_pnts[3];
	double e = 1.0e-12;


	pnts[0] = all_points[PatCC_Edge->prev_edge_in_triangle->head];
	pnts[1] = all_points[PatCC_Edge->head];
    pnts[2] = all_points[PatCC_Edge->tail];
	pnts[3] = *p;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, pnts[2].position_to_edge(&(pnts[0]), &(pnts[1])) == 1, "Software error in PatCC_Delaunay_Voronoi::circum_circle_contains_reliably_3D");
	for (int i = 0; i < 3; i ++) {
		tmp_pnts[i].xx = pnts[i].xx - pnts[3].xx;
		tmp_pnts[i].yy = pnts[i].yy - pnts[3].yy;
		tmp_pnts[i].zz = pnts[i].zz - pnts[3].zz;
	}

	double res = compute_three_3D_points_cross_product(tmp_pnts[2].xx, tmp_pnts[2].yy, tmp_pnts[2].zz, tmp_pnts[0].xx, tmp_pnts[0].yy, tmp_pnts[0].zz, tmp_pnts[1].xx, tmp_pnts[1].yy, tmp_pnts[1].zz);
	if (fabs(res) <= e)
		return 0;
	else if (res > 0)
		return -1;
	else return 1;
}



/*
 * Input : PatCC_Point to be checked
 * Return:  1    point is in circum circle
 *          0    point is on circum circle
 *         -1    point is out of circum circle
 */
int PatCC_Delaunay_Voronoi::circum_circle_contains_reliably_new(const PatCC_Edge *PatCC_Edge, PatCC_Point *p, double tolerance)
{
	PatCC_Point pnts[4];
	bool redo_projection = do_projection;
	PatCC_Triangle tmp_triangles[4];
	int indie_position_values[4], cross_position_values[4];
	PAT_REAL current_max_radius, counterpart_max_radius;
	int current_result, counterpart_result, overall_result;


	pnts[0] = all_points[PatCC_Edge->prev_edge_in_triangle->head];
	pnts[1] = all_points[PatCC_Edge->head];
	pnts[2] = *p;
    pnts[3] = all_points[PatCC_Edge->tail];
	if (do_projection) {
		for (int i = 0; i < 4; i ++)
			if (pnts[i].id < 0 || pnts[i].id >= num_points)
				redo_projection = false;
		if (redo_projection)
			for (int i = 0; i < 4; i ++)
				for (int j = i+1; j < 4; j ++)
					if (fabs(pnts[i].x - pnts[j].x) >= 5.0 || fabs(pnts[i].y - pnts[j].y) >= 5.0)
						redo_projection = false;
	}
	if (redo_projection) {
		int center_point_index = 0;
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special trace-1 (%lf, %lf) : %d  (%lf, %lf) %d  (%lf, %lf) %d  (%lf, %lf) %d", pnts[0].x, pnts[0].y, pnts[0].id, pnts[1].x, pnts[1].y, pnts[1].id, pnts[2].x, pnts[2].y, pnts[2].id, pnts[3].x, pnts[3].y, pnts[3].id);
		for (int i = 0; i < 4; i ++) {
			tmp_triangles[i].calulate_circum_circle(&pnts[(i+3)%4], &pnts[i], &pnts[(i+1)%4]);
			indie_position_values[i] = tmp_triangles[i].circum_circle_position_to(&pnts[(i+2)%4], tolerance);
			EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special trace0: %d: (%lf, %lf): %lf : %d\n", i, tmp_triangles[i].circum_center[0], tmp_triangles[i].circum_center[1], tmp_triangles[i].circum_radius2, indie_position_values[i]);
		}
		for (int i = 0; i < 4; i ++) {
			if (global_index[pnts[center_point_index].id] < global_index[pnts[i].id])
				center_point_index = i;
			pnts[i].x = x_ref[pnts[i].id];
			pnts[i].y = y_ref[pnts[i].id];
		}
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special trace1 (%lf, %lf)  (%lf, %lf)  (%lf, %lf)  (%lf, %lf)", pnts[0].x, pnts[0].y, pnts[1].x, pnts[1].y, pnts[2].x, pnts[2].y, pnts[3].x, pnts[3].y);
		PAT_REAL uv1_x, uv1_y, uv1_z, uv2_x, uv2_y, uv2_z, center_x, center_y, center_z;
		lonlat2xyz(pnts[center_point_index].x, pnts[center_point_index].y, &center_x, &center_y, &center_z);
		calculate_unit_vectors(pnts[center_point_index].x, pnts[center_point_index].y, &uv1_x, &uv1_y, &uv1_z, &uv2_x, &uv2_y, &uv2_z);		
		for (int i = 0; i < 4; i ++)
			fast_stereographic_projection(pnts[i].x, pnts[i].y, center_x, center_y, center_z, uv1_x, uv1_y, uv1_z,
										  uv2_x, uv2_y, uv2_z, pnts[i].x, pnts[i].y);
		pnts[center_point_index].x = 0.0;
		pnts[center_point_index].y = 0.0;
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special trace2 (%lf, %lf)  (%lf, %lf)  (%lf, %lf)  (%lf, %lf)", pnts[0].x, pnts[0].y, pnts[1].x, pnts[1].y, pnts[2].x, pnts[2].y, pnts[3].x, pnts[3].y);
	}
	for (int i = 0; i < 4; i ++) {
		tmp_triangles[i].calulate_circum_circle(&pnts[(i+3)%4], &pnts[i], &pnts[(i+1)%4]);
		indie_position_values[i] = tmp_triangles[i].circum_circle_position_to(&pnts[(i+2)%4], tolerance);
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special trace3: %d: (%lf, %lf): %lf : %d\n", i, tmp_triangles[i].circum_center[0], tmp_triangles[i].circum_center[1], tmp_triangles[i].circum_radius2, indie_position_values[i]);
	}
	cross_position_values[0] = tmp_triangles[0].circum_radius2 > tmp_triangles[2].circum_radius2 ? indie_position_values[0] : indie_position_values[2];
	cross_position_values[2] = tmp_triangles[0].circum_radius2 <= tmp_triangles[2].circum_radius2 ? indie_position_values[0] : indie_position_values[2];
	cross_position_values[1] = tmp_triangles[1].circum_radius2 > tmp_triangles[1].circum_radius2 ? indie_position_values[1] : indie_position_values[1];
	cross_position_values[3] = tmp_triangles[3].circum_radius2 <= tmp_triangles[3].circum_radius2 ? indie_position_values[3] : indie_position_values[3];
	current_max_radius = std::max(tmp_triangles[0].circum_radius2, tmp_triangles[2].circum_radius2);
	counterpart_max_radius = std::max(tmp_triangles[1].circum_radius2, tmp_triangles[3].circum_radius2);

	if (indie_position_values[0] > 0 && indie_position_values[2] > 0)
		current_result = 0;
	else if (cross_position_values[0] == 0 || cross_position_values[0] == 2 || cross_position_values[2] == 0)
	    current_result = 1;
	else current_result = -1;

	if (indie_position_values[1] > 0 && indie_position_values[3] > 0)
		counterpart_result = 0;
	else if (cross_position_values[1] == 0 || cross_position_values[1] == 2 || cross_position_values[3] == 0)
	    counterpart_result = 1;
	else counterpart_result = -1;	

	if (current_result * counterpart_result < 0)
		overall_result = current_result;
	else if (current_result == counterpart_result) {
	    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(current_result == 1 && counterpart_result == 1), "Software error4 in circum_circle_contains_reliably");	
		overall_result = 0;
	}
	else {
		if (relative_eq_int(counterpart_max_radius, current_max_radius, tolerance))
			overall_result = 0;
		else if (current_max_radius < counterpart_max_radius && current_result == 0 ||
				 counterpart_max_radius < current_max_radius && counterpart_result == 0)
			overall_result = 0;
		else {
			if (current_max_radius < counterpart_max_radius)
				overall_result = current_result;
			else overall_result = counterpart_result * -1;
		}
	}

	return overall_result;
}


/*
 * Input : PatCC_Point to be checked
 * Return:  1    point is in circum circle
 *          0    point is on circum circle
 *         -1    point is out of circum circle
 */
int PatCC_Delaunay_Voronoi::circum_circle_contains_reliably(const PatCC_Edge *PatCC_Edge, PatCC_Point *p, double tolerance)
{

	int current_result, counterpart_result, overall_result;
	PatCC_Triangle tmp_tri1, tmp_tri2;
	int ret1 = PatCC_Edge->triangle->circum_circle_position_to(p, tolerance);
	int ret2 = PatCC_Edge->twin_edge->triangle->circum_circle_position_to(&all_points[PatCC_Edge->prev_edge_in_triangle->head], tolerance);
	int ret3 = PatCC_Edge->triangle->circum_radius2 > PatCC_Edge->twin_edge->triangle->circum_radius2 ? ret1 : ret2;
	int ret4 = PatCC_Edge->triangle->circum_radius2 <= PatCC_Edge->twin_edge->triangle->circum_radius2 ? ret1 : ret2;
	PAT_REAL current_max_radius = std::max(PatCC_Edge->triangle->circum_radius2, PatCC_Edge->twin_edge->triangle->circum_radius2);


	if (p->is_spheric_grid) {
		counterpart_result = circum_circle_contains_reliably_3D(PatCC_Edge, p);
		if (counterpart_result != overall_result) {
			PatCC_Point pnts[4];
			pnts[0] = all_points[PatCC_Edge->prev_edge_in_triangle->head];
			pnts[1] = all_points[PatCC_Edge->head];
			pnts[3] = *p;
			pnts[2] = all_points[PatCC_Edge->tail];
		}
		return counterpart_result;
	}

	if (ret1 > 0 && ret2 > 0)
		current_result = 0;
	else if (ret3 == 0 || ret3 == 2 || ret4 == 0)
	    current_result = 1;
	else current_result = -1;

	tmp_tri1.calulate_circum_circle(p, &all_points[PatCC_Edge->tail], &all_points[PatCC_Edge->next_edge_in_triangle->tail]);
	tmp_tri2.calulate_circum_circle(p, &all_points[PatCC_Edge->head], &all_points[PatCC_Edge->next_edge_in_triangle->tail]);
	PAT_REAL counterpart_max_radius = std::max(tmp_tri1.circum_radius2, tmp_tri2.circum_radius2);
	ret1 = tmp_tri1.circum_circle_position_to(&all_points[PatCC_Edge->head], tolerance);
	ret2 = tmp_tri2.circum_circle_position_to(&all_points[PatCC_Edge->tail], tolerance);
	ret3 = tmp_tri1.circum_radius2 > tmp_tri2.circum_radius2 ? ret1 : ret2;
	ret4 = tmp_tri1.circum_radius2 <= tmp_tri2.circum_radius2 ? ret1 : ret2;

	if (ret1 > 0 && ret2 > 0)
		counterpart_result = 0;
	else if (ret3 == 0 || ret3 == 2 || ret4 == 0)
	    counterpart_result = 1;
	else counterpart_result = -1;

	if (current_result * counterpart_result < 0)
		overall_result = current_result;
	else if (current_result == counterpart_result) {
		if (current_result == 1 && counterpart_result == 1) {
			counterpart_result = circum_circle_contains_reliably_3D(PatCC_Edge, p);
			PatCC_Point pnts[4];
			pnts[0] = all_points[PatCC_Edge->prev_edge_in_triangle->head];
			pnts[1] = all_points[PatCC_Edge->head];
			pnts[3] = *p;
			pnts[2] = all_points[PatCC_Edge->tail];			
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(current_result == 1 && counterpart_result == 1), "Software error4 in circum_circle_contains_reliably: (%lf, %lf) (%lf, %lf) (%lf, %lf) (%lf, %lf): %d", pnts[0].x, pnts[0].y, pnts[1].x, pnts[1].y, pnts[2].x, pnts[2].y, pnts[3].x, pnts[3].y, counterpart_result);	
		}
		overall_result = 0;
	}
	else {
		if (relative_eq_int(counterpart_max_radius, current_max_radius, tolerance))
			overall_result = 0;
		else if (current_max_radius < counterpart_max_radius && current_result == 0 ||
				 counterpart_max_radius < current_max_radius && counterpart_result == 0)
			overall_result = 0;
		else {
			if (current_max_radius < counterpart_max_radius)
				overall_result = current_result;
			else overall_result = counterpart_result * -1;
		}
	}

	if (report_error_enabled && false) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, overall_result == circum_circle_contains_reliably_new(PatCC_Edge, p, tolerance), "Software erorr in PatCC_Delaunay_Voronoi::circum_circle_contains_reliably");
	}

	return overall_result;
}


/*
 * Input : PatCC_Point to be checked
 * Return:  2  on/in
 *          1  on/out
 *          0  in
 *         -1  out
 */
int PatCC_Triangle::circum_circle_position_to(PatCC_Point *p, double tolerance)
{
	double dist2 = ((p->x - circum_center[0]) * (p->x - circum_center[0])) + ((p->y - circum_center[1]) * (p->y - circum_center[1]));

	int ret = dist2 > circum_radius2 ? -1 : 0;

	if(relative_eq_int(dist2, circum_radius2, tolerance))
		ret += 2;

	return ret;
}




int PatCC_Triangle::find_best_candidate_point(PatCC_Point* buf, std::vector<int> &candicate_points_indexes)
{
	if (remained_points_tail == -1)
		return -1;

	double min_dist=1e10, dist;
	double center_x = (buf[v[0]].x+buf[v[1]].x+buf[v[2]].x) * 0.3333333333333333333333333333;
	double center_y = (buf[v[0]].y+buf[v[1]].y+buf[v[2]].y) * 0.3333333333333333333333333333;
	int best_candidate_id = -1;
	bool no_more_mask = true;


	for (int i = remained_points_head; i > -1; i = buf[i].next) {
		if (buf[i].mask)
			no_more_mask = false;
		dist = buf[i].calculate_distance_MD(center_x, center_y);
		if (best_candidate_id == -1 || min_dist > dist) {
			min_dist = dist;
			best_candidate_id = i;
		} 
	}

	candicate_points_indexes.clear();
    for (int i = remained_points_head; i > -1; i = buf[i].next) {
		if (best_candidate_id == i || (&(buf[best_candidate_id]))->is_the_same_with_another(&(buf[i]))) 
			candicate_points_indexes.push_back(i);
	}
	
	return best_candidate_id;
}


bool PatCC_Triangle::contain_vertex(int pt)
{
	for(int i = 0; i < 3; i++)
		if(v[i] == pt)
			return true;
	return false;
}


void PatCC_Delaunay_Voronoi::swap_points(int idx1, int idx2)
{
	double tmp_x    = all_points[idx1].x;
	double tmp_y    = all_points[idx1].y;
	double tmp_xx   = all_points[idx1].xx;
	double tmp_yy   = all_points[idx1].yy;
	double tmp_zz   = all_points[idx1].zz;
	int    tmp_id   = all_points[idx1].id;
	bool   tmp_mask = all_points[idx1].mask;

	all_points[idx1].x    = all_points[idx2].x;
	all_points[idx1].y    = all_points[idx2].y;
	all_points[idx1].xx	  = all_points[idx2].xx;
	all_points[idx1].yy   = all_points[idx2].yy;
	all_points[idx1].zz   = all_points[idx2].zz;
	all_points[idx1].id   = all_points[idx2].id;
	all_points[idx1].mask = all_points[idx2].mask;

	all_points[idx2].x    = tmp_x;
	all_points[idx2].y    = tmp_y;
	all_points[idx2].xx   = tmp_xx;
	all_points[idx2].yy   = tmp_yy;
	all_points[idx2].zz   = tmp_zz;
	all_points[idx2].id   = tmp_id;
	all_points[idx2].mask = tmp_mask;
}


/*
 *   head:  head of linked list of points.
 *   tail:  tail of linked list of points.
 *   base:  base of the triangle stack.
 *   top:   top of the triangle stack.
 */
void PatCC_Delaunay_Voronoi::distribute_points_into_triangles(int head, int tail, unsigned base, unsigned top)
{
	int start = head;
	int end, j, num_total_points = 0, num_distributed_points = 0;

	if (tail == -1)
		return;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, head != -1 && all_points[tail].next == -1, "Software error in PatCC_Delaunay_Voronoi::distribute_points_into_triangles");

	for (j = start; j > -1; j = all_points[j].next)
		num_total_points ++;

	for (unsigned i = base+1; i <= top; i++) {
		if (!triangle_stack[i] || !triangle_stack[i]->is_leaf)
			continue;

		end = tail;

		/* put points in current triangle together */
		for (j = start; j != end && j > -1;) {
			int* v_idx = triangle_stack[i]->v;
			if (all_points[j].position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) >= 0) { // in triangles or on the PatCC_Edge
				j = all_points[j].next;
				num_distributed_points ++;
			} else {
				swap_points(j, end);
				end = all_points[end].prev;
			}
		}

		int* v_idx = triangle_stack[i]->v;
		if (j > -1 && all_points[end].position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) < 0) // Case "j == end"
			end = all_points[end].prev;
		else if (j > -1)
			num_distributed_points ++;

		/* set triangle remained points info */
		triangle_stack[i]->set_remained_points(start, end);

		/* go on */
		triangle_stack[i]->set_remained_points(start, end);
		if (end > -1) {
			start = all_points[end].next;
			all_points[end].next = -1;
			if (start > -1)
				all_points[start].prev = -1;
			else
				break;
		}
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, start == -1 && num_total_points == num_distributed_points, "Software error in PatCC_Delaunay_Voronoi::distribute_points_into_triangles");
}


void PatCC_Triangle::set_remained_points(int head, int tail)
{
	if (tail > -1 && head > -1) {
		remained_points_head = head;
		remained_points_tail = tail;
	} else {
		remained_points_head = -1;
		remained_points_tail = -1;
	}
}


int PatCC_Triangle::pop_tail(PatCC_Point* buf, int num_points)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remained_points_tail != -1, "Software error in PatCC_Triangle::pop_tail");
	int old_tail = remained_points_tail;

	remained_points_tail = buf[remained_points_tail].prev;
	if (remained_points_tail > -1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, remained_points_tail < num_points, "Software error in PatCC_Triangle::pop_tail");
		buf[remained_points_tail].next = -1;
	}
	else
		remained_points_head = -1;

	return old_tail;
}


void PatCC_Delaunay_Voronoi::push(unsigned *stack_top, PatCC_Triangle* t)
{
	diag_triangle(t);

	*stack_top = *stack_top + 1;
	if (*stack_top >= stack_size) {
		stack_size *= 2;
		PatCC_Triangle** tmp = new PatCC_Triangle*[stack_size];
		memcpy(tmp, triangle_stack, sizeof(PatCC_Triangle*) * *stack_top);
		delete [] triangle_stack;
		triangle_stack = tmp;
	}

	triangle_stack[*stack_top] = t;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t->stack_ref_count == -1, "Software error in PatCC_Delaunay_Voronoi::push");
	t->stack_ref_count = *stack_top;
}


static inline double point_distence_to_line(double px, double py, double x1, double y1, double x2, double y2)
{
	//TODO: may be slow
	if(x1 == x2)
		return fabs(px - x1);
	else if(y1 == y2)
		return fabs(py - y1);
	else {
		double A=(y1-y2)/(x1-x2);
		double B=y1-A*x1;
		return fabs((A*px+B-py)/sqrt(A*A+1));
	}
}

static inline bool point_distence_in_threshold(PatCC_Point v, PatCC_Point p1, PatCC_Point p2, double threshold) {
	return point_distence_to_line(v.x, v.y, p1.x, p1.y, p2.x, p2.y) <= threshold;
}

static inline bool all_distence_in_threshold(Triangle_inline* triangle, PatCC_Point p1, PatCC_Point p2, double threshold) {
	return point_distence_in_threshold(triangle->v[0], p1, p2, threshold) &&
		   point_distence_in_threshold(triangle->v[1], p1, p2, threshold) &&
		   point_distence_in_threshold(triangle->v[2], p1, p2, threshold);
}

static inline bool is_segment_in_triangle(Triangle_inline* triangle, PatCC_Point p1, PatCC_Point p2) {
	return p1.position_to_triangle(triangle) >= 0 && p2.position_to_triangle(triangle) >= 0;
}

static inline bool is_segment_intersected_with_edge(PatCC_Point p_e1, PatCC_Point p_e2, PatCC_Point p1, PatCC_Point p2) {
	return (p1.position_to_edge(&p_e1, &p_e2) * p2.position_to_edge(&p_e1, &p_e2) <= 0) &&
		   (p_e1.position_to_edge(&p1, &p2) * p_e2.position_to_edge(&p1, &p2) <= 0);
}

static inline bool is_segment_intersected_with_any_edges(Triangle_inline* triangle, PatCC_Point p1, PatCC_Point p2) {
	return is_segment_intersected_with_edge(triangle->v[0], triangle->v[1], p1, p2) ||
		   is_segment_intersected_with_edge(triangle->v[1], triangle->v[2], p1, p2) ||
		   is_segment_intersected_with_edge(triangle->v[2], triangle->v[0], p1, p2);
}

inline void get_bounding_box(Triangle_inline* t, double& x_min, double& x_max, double& y_min, double& y_max) {
	PatCC_Point* v = t->v;
	x_min = v[0].x;
	x_max = v[0].x;
	y_min = v[0].y;
	y_max = v[0].y;
	for (int i = 1; i < 3; i++) {
		if (v[i].x < x_min) x_min = v[i].x;
		if (v[i].x > x_max) x_max = v[i].x;
		if (v[i].y < y_min) y_min = v[i].y;
		if (v[i].y > y_max) y_max = v[i].y;
	}
}

bool is_triangle_intersecting_with_segment(Triangle_inline* triangle, PatCC_Point p1, PatCC_Point p2, double threshold) 
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, p1.x == p2.x || p1.y == p2.y, "Software error in is_triangle_intersecting_with_segment");
	double x_min, x_max, y_min, y_max;
	get_bounding_box(triangle, x_min, x_max, y_min, y_max);

	double seg_min, seg_max;
	if (p1.x == p2.x) {
		if (p1.y > p2.y) {
			seg_min = p2.y;
			seg_max = p1.y;
		} else {
			seg_min = p1.y;
			seg_max = p2.y;
		}

		if(!(x_min <= p1.x && x_max >= p1.x && !(y_max < seg_min || y_min > seg_max)))
			return false;
	} else {
		if (p1.x > p2.x) {
			seg_min = p2.x;
			seg_max = p1.x;
		} else {
			seg_min = p1.x;
			seg_max = p2.x;
		}

		if(!(y_min <= p1.y && y_max >= p1.y && !(x_max < seg_min || x_min > seg_max))) {
			return false;
		}
	}
	return (is_segment_intersected_with_any_edges(triangle, p1, p2) || is_segment_in_triangle(triangle, p1, p2));
}


unsigned PatCC_Delaunay_Voronoi::triangulating_process(PatCC_Triangle *triangle, unsigned stack_base)
{
	unsigned stack_top = stack_base, stack_pos_mark;
	PatCC_Triangle *tilr = NULL, *tjrl = NULL;


	redistribute_head_tails_content_size = 0;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
	diag_triangle(triangle);

	int candidate_id = triangle->find_best_candidate_point(&all_points[0], candicate_points_indexes);

	if (candidate_id == -1) 
		return stack_top;

	int* v_idx = triangle->v;
	int dividing_idx;
	candidate_id = -1;
	for (int i = candicate_points_indexes.size()-1; i >=0; i --) {
		if (candidate_id == -1)
			candidate_id = candicate_points_indexes[i];
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, all_points[candicate_points_indexes[i]].id < 0 && all_points[candidate_id].id < 0 || all_points[candicate_points_indexes[i]].id >= 0 && all_points[candidate_id].id >= 0, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
		if (all_points[candidate_id].id >= 0 && global_index[all_points[candidate_id].id] > global_index[all_points[candicate_points_indexes[i]].id])
			candidate_id = candicate_points_indexes[i];
	}
	int candidate_point_id = all_points[candidate_id].id;
	for (int i = candicate_points_indexes.size()-1; i >=0; i --) {
		if (candicate_points_indexes[i] != candidate_id) 
	        overlapping_points.push_back(std::make_pair(candidate_point_id, all_points[candicate_points_indexes[i]].id));
		swap_points(candicate_points_indexes[i], triangle->remained_points_tail);
		if (candicate_points_indexes[i] != candidate_id)
			triangle->pop_tail(&all_points[0], num_points);
		else dividing_idx = triangle->pop_tail(&all_points[0], num_points);
	}
	PatCC_Point* dividing_point = &all_points[dividing_idx];
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, candidate_point_id == dividing_point->id, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");

	for (int i = 0; i < 3; i ++)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !are_two_points_the_same(v_idx[i],dividing_idx), "Software error in PatCC_Delaunay_Voronoi::triangulating_process");

	if (dividing_point->position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]) == 0) { // inside
		PatCC_Edge *e_v1_can = allocate_edge(triangle->v[0], dividing_idx);
		PatCC_Edge *e_can_v1 = make_twins_edge(e_v1_can);
		PatCC_Edge *e_v2_can = allocate_edge(triangle->v[1], dividing_idx);
		PatCC_Edge *e_can_v2 = make_twins_edge(e_v2_can);
		PatCC_Edge *e_v3_can = allocate_edge(triangle->v[2], dividing_idx);
		PatCC_Edge *e_can_v3 = make_twins_edge(e_v3_can);
		PatCC_Triangle *t_can_v1_v2 = allocate_triangle(e_can_v1, triangle->edge[0], e_v2_can);
		PatCC_Triangle *t_can_v2_v3 = allocate_triangle(e_can_v2, triangle->edge[1], e_v3_can);
		PatCC_Triangle *t_can_v3_v1 = allocate_triangle(e_can_v3, triangle->edge[2], e_v1_can);
		push(&stack_top, t_can_v1_v2);
		stack_pos_mark = stack_top;
		push(&stack_top, t_can_v2_v3);
		push(&stack_top, t_can_v3_v1);
		
		/* Actually, vertex(t_can_v1_v2, 0) is dividing_point */
		if(t_can_v1_v2->stack_ref_count == stack_pos_mark) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t_can_v1_v2->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(t_can_v1_v2->v[0], t_can_v1_v2->edge[1], stack_base, &stack_top);
		}
		if (t_can_v2_v3->stack_ref_count == stack_pos_mark+1) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t_can_v2_v3->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(t_can_v2_v3->v[0], t_can_v2_v3->edge[1], stack_base, &stack_top);
		}
		if (t_can_v3_v1->stack_ref_count == stack_pos_mark+2) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t_can_v3_v1->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(t_can_v3_v1->v[0], t_can_v3_v1->edge[1], stack_base, &stack_top);
		}
	} else { // on the side
		int idx_i, idx_j, idx_k;
		PatCC_Edge *eij;
		switch (dividing_point->position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]])) {
			case 1:
				idx_i = triangle->v[0];
				idx_j = triangle->v[1];
				idx_k = triangle->v[2];
				eij = triangle->edge[0];
				break;
			case 2:
				idx_i = triangle->v[1];
				idx_j = triangle->v[2];
				idx_k = triangle->v[0];
				eij = triangle->edge[1];
				break;
			case 3:
				idx_i = triangle->v[2];
				idx_j = triangle->v[0];
				idx_k = triangle->v[1];
				eij = triangle->edge[2];
				break;
			default:
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "point, which should be found in triangle, is outside of triangle\n");
				break;
		}
		PatCC_Point* vi = &all_points[idx_i];
		PatCC_Point* vj = &all_points[idx_j];
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dividing_point->position_to_edge(vi, vj) == 0, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, eij->twin_edge == NULL || eij->twin_edge->triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");

		PatCC_Edge* ejk = eij->next_edge_in_triangle;
		PatCC_Edge* eki = ejk->next_edge_in_triangle;
		PatCC_Edge *eil, *elj, *eji;

		PatCC_Edge *eir = allocate_edge(idx_i, dividing_idx);
		PatCC_Edge *erk = allocate_edge(dividing_idx, idx_k);
		PatCC_Edge *ekr = make_twins_edge(erk);
		PatCC_Edge *erj = allocate_edge(dividing_idx, idx_j);
		PatCC_Triangle* tirk = allocate_triangle(eir, erk, eki);
		PatCC_Triangle* tjkr = allocate_triangle(ejk, ekr, erj);

		push(&stack_top, tjkr);
		stack_pos_mark = stack_top;
		push(&stack_top, tirk);
		if (eij->twin_edge != NULL) {
			eji = eij->twin_edge;
			eil = eji->next_edge_in_triangle;
			elj = eil->next_edge_in_triangle;
			int idx_l = elj->head;
			PatCC_Edge *eri = make_twins_edge(eir);
			PatCC_Edge *ejr = make_twins_edge(erj);
			PatCC_Edge *erl = allocate_edge(dividing_idx, idx_l);
			PatCC_Edge *elr = make_twins_edge(erl);
			tilr = allocate_triangle(eil, elr, eri);
			tjrl = allocate_triangle(ejr, erl, elj);
			push(&stack_top, tilr);
			push(&stack_top, tjrl);
			add_redistribute_head_tail(eij->twin_edge->triangle);
		}
		if (tjkr->stack_ref_count == stack_pos_mark) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tjkr->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(dividing_idx, ejk, stack_base, &stack_top);
		}
		if (tirk->stack_ref_count == stack_pos_mark+1) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tirk->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(dividing_idx, eki, stack_base, &stack_top); 
		}
		if (tilr != NULL && tilr->stack_ref_count == stack_pos_mark+2) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tilr->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(dividing_idx, eil, stack_base, &stack_top);
		}
		if (tjrl != NULL && tjrl->stack_ref_count == stack_pos_mark+3) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, tjrl->is_leaf, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
			legalize_triangles(dividing_idx, elj, stack_base, &stack_top);
		}
	}

	add_redistribute_head_tail(triangle);

	if (report_error_enabled) {
		for (unsigned i = stack_base+1; i <= stack_top; i ++)
			if (triangle_stack[i] && triangle_stack[i]->is_leaf)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, triangle_stack[i]->remained_points_head == -1 && triangle_stack[i]->remained_points_tail == -1, "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
	}

	int list_head, list_tail;
	link_remained_list(stack_base, stack_top, &list_head, &list_tail);

	distribute_points_into_triangles(list_head, list_tail, stack_base, stack_top);

	for (unsigned i = stack_base+1; i <= stack_top; i ++)
		if(triangle_stack[i] != NULL && triangle_stack[i]->is_leaf) {
			for (int j = i + 1; j <= stack_top; j ++)
				diag_triangle(triangle_stack[i]);
			triangulating_process(triangle_stack[i], stack_top);
		}

	return stack_top;
}

void PatCC_Delaunay_Voronoi::clean_triangle(PatCC_Triangle* t)
{
	for (int i = 0; i < 3; i++)
		ref_dec(t->edge[i]);
	triangle_allocator->deleteElement(t);
	num_allocate_triangles --;
}


static int compare_node_index(const void* a, const void* b)
{
	if (*(const int*)a < *(const int*)b) 
		return -1;
	if (*(const int*)a > *(const int*)b) 
		return 1;
	return 0;
}


void PatCC_Delaunay_Voronoi::add_redistribute_head_tail(PatCC_Triangle *triangle)
{
	if (triangle->remained_points_tail > -1) {
		write_data_into_array_buffer(&(triangle->remained_points_head), sizeof(int), (char**)(&redistribute_head_tails), redistribute_head_tails_max_size, redistribute_head_tails_content_size);
		write_data_into_array_buffer(&(triangle->remained_points_tail), sizeof(int), (char**)(&redistribute_head_tails), redistribute_head_tails_max_size, redistribute_head_tails_content_size);
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, triangle->is_leaf, "Software error in PatCC_Delaunay_Voronoi::add_redistribute_head_tail");
	if (triangle_stack[triangle->stack_ref_count] == triangle)
		triangle_stack[triangle->stack_ref_count] = NULL;
	triangle->remained_points_head = -1;
	triangle->remained_points_tail = -1;
	triangle->stack_ref_count = -1;
	triangle->is_leaf = false;
	clean_triangle(triangle);
}


void PatCC_Delaunay_Voronoi::link_remained_list(unsigned base, unsigned top, int* head, int* tail)
{
	unsigned i;


	if (report_error_enabled) {
		bool* map = new bool[num_points]();
		for (i = base+1; i <= top; i ++)
			if (triangle_stack[i] && !triangle_stack[i]->is_leaf && triangle_stack[i]->remained_points_tail > -1) {
				for (int j = triangle_stack[i]->remained_points_head; j > -1; j = all_points[j].next) {
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !map[j], "Software error in PatCC_Delaunay_Voronoi::triangulating_process");
					map[j] = true;
				}
			}
		delete[] map;
	}

	redistribute_head_tails_content_size = redistribute_head_tails_content_size/2/sizeof(int);

	unsigned point_count = 0;
	if (redistribute_head_tails_content_size > 0) {
		merge_sort(redistribute_head_tails, redistribute_head_tails_content_size, sizeof(int)*2, compare_node_index);
		if (report_error_enabled) {
			for (unsigned i = 0; i < redistribute_head_tails_content_size; i++)
				for (int j = redistribute_head_tails[i*2]; j > -1; j = all_points[j].next)
					point_count++;
		}

		for (unsigned i = 0; i < redistribute_head_tails_content_size - 1; i++) {
			int cur_tail_id = redistribute_head_tails[i*2+1];
			int nxt_head_id = redistribute_head_tails[i*2+2];
			all_points[cur_tail_id].next = nxt_head_id;
			all_points[nxt_head_id].prev = cur_tail_id;
		}
		*head = redistribute_head_tails[0];
		*tail = redistribute_head_tails[redistribute_head_tails_content_size*2-1];
	} else {
		*head = -1;
		*tail = -1;
	}

	if (report_error_enabled) {
		unsigned point_count2 = 0;
		for (int i = *head; i > -1; i = all_points[i].next)
			point_count2 ++;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, point_count == point_count2, "Software error PatCC_Delaunay_Voronoi::link_remained_list");
	}
}


void PatCC_Delaunay_Voronoi::map_buffer_index_to_point_index()
{
	if (point_idx_to_buf_idx != NULL)
		delete[] point_idx_to_buf_idx;
	point_idx_to_buf_idx = new int[num_points - PAT_NUM_LOCAL_VPOINTS]();
	for (int i = PAT_NUM_LOCAL_VPOINTS; i < num_points; i++) {
		if (all_points[i].id == -1)
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, point_idx_to_buf_idx[all_points[i].id] == 0, "Software error in PatCC_Delaunay_Voronoi::map_buffer_index_to_point_index");
		point_idx_to_buf_idx[all_points[i].id] = i;
	}
}


void PatCC_Delaunay_Voronoi::mark_special_triangles()
{
	for (unsigned i = 0; i < all_leaf_triangles.size(); i++) {
		if(!all_leaf_triangles[i]->is_leaf)
			continue;

		for (int j = 0; j < 3; j ++) {
			if (!is_real_point(all_leaf_triangles[i]->v[j])) {
				all_leaf_triangles[i]->is_virtual = true;
				break;
			}
		}
	}
}

PatCC_Delaunay_Voronoi::PatCC_Delaunay_Voronoi(int num_grid_points, double min_lon_x, double max_lon_x, double min_lat_y, double max_lat_y, double *lon_x_values, double *lat_y_values, bool *mask, bool is_spheric_grid, int *cells_global_index, double grid_min_lon, double grid_max_lon, double grid_min_lat, double grid_max_lat, double kernel_min_lon, double kernel_max_lon, double kernel_min_lat, double kernel_max_lat, long kernel_grid_size, int subdomain_type)
	: triangle_stack(NULL)
	, stack_size(0)
	, is_regional(false)
	, polar_mode(false)
	, fast_mode(false)
	, tolerance(1e-9)
	, vpolar_local_index(-1)
	, original_lon_center(-1e10)
	, point_idx_to_buf_idx(NULL)
	, do_projection(false)
{
	overlapping_points.clear();
	num_polar_point_in_origin_coordinate = 0;
	int max_extend_coordinate_buffer_size;
	bool *mask_ref;


	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "conduct triangulation for kernal subdomain (%lf, %lf, %lf, %lf) in expanded subdomain (%lf, %lf, %lf, %lf) ", kernel_min_lon, kernel_max_lon, kernel_min_lat, kernel_max_lat, min_lon_x, max_lon_x, min_lat_y, max_lat_y);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_grid_points > 0, "Software error in PatCC_Delaunay_Voronoi::PatCC_Delaunay_Voronoi");
	if (min_lon_x > max_lon_x)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_spheric_grid, "Software error in PatCC_Delaunay_Voronoi::PatCC_Delaunay_Voronoi");
	num_points = num_grid_points;
	polar_mode = subdomain_type == 0 ? false : true;
	triangle_allocator = new Triangle_pool();
	edge_allocator = new Edge_pool();
	result_triangles = NULL;

	max_extend_coordinate_buffer_size = num_points + 5;
	x_ref = new double [max_extend_coordinate_buffer_size];
	y_ref = new double [max_extend_coordinate_buffer_size];
	global_index = new int [max_extend_coordinate_buffer_size];
	mask_ref = new bool [max_extend_coordinate_buffer_size];
	memcpy(mask_ref, mask, sizeof(bool)*num_points);
	for(int i = num_points; i < max_extend_coordinate_buffer_size; i ++)
		mask_ref[i] = true;

	redistribute_head_tails = NULL;
	original_min_lat = min_lat_y;
	original_max_lat = max_lat_y;
	original_min_lon = min_lon_x;
	original_max_lon = max_lon_x;

	origin_lon_values = lon_x_values;
	origin_lat_values = lat_y_values;
	origin_num_points = num_points;

	this->grid_min_lon = grid_min_lon;
	this->grid_max_lon = grid_max_lon;
	this->grid_min_lat = grid_min_lat;
	this->grid_max_lat = grid_max_lat;

	this->kernel_min_lon = kernel_min_lon;
	this->kernel_max_lon = kernel_max_lon;
	this->kernel_min_lat = kernel_min_lat;
	this->kernel_max_lat = kernel_max_lat;
	this->kernel_grid_size = kernel_grid_size;
	this->is_spheric_grid = is_spheric_grid;
	this->num_allocate_triangles = 0;
	this->num_allocate_edges = 0;

	if (min_lon_x > max_lon_x) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_spheric_grid, "Software error in PatCC_Delaunay_Voronoi::PatCC_Delaunay_Voronoi, is_spheric_grid wrong.");
		original_min_lon -= 360.0;
	}

	if (is_spheric_grid) {
		if (min_lat_y == -90.0 && polar_mode) {
			original_lon_center = 0.0;
			original_lat_center = -90.0;
		}
		else if (max_lat_y == 90.0 && polar_mode) {
			original_lon_center = 0.0;
			original_lat_center = 90.0;
		}
		else {
			original_lat_center = (min_lat_y+max_lat_y) / 2.0;
			original_lon_center = (original_min_lon+original_max_lon) / 2.0;
		}
	}

	for (int i = 0; i < num_points; i ++) {
		if (min_lon_x > max_lon_x && lon_x_values[i] >= min_lon_x)
			x_ref[i] = lon_x_values[i] - 360.0;
		else x_ref[i] = lon_x_values[i];
		y_ref[i] = lat_y_values[i];
		global_index[i] = cells_global_index[i];
	}

	avoiding_line_head[0].x = avoiding_line_head[0].y = 0;
	avoiding_line_head[1].x = avoiding_line_head[1].y = 0;
	avoiding_line_tail[0].x = avoiding_line_tail[0].y = 0;
	avoiding_line_tail[1].x = avoiding_line_tail[1].y = 0;
	avoiding_line_head[0].id = avoiding_line_head[1].id = 0;

	avoiding_circle_center[0].x = avoiding_circle_center[0].y = 0;
	avoiding_circle_center[1].x = avoiding_circle_center[1].y = 0;
	avoiding_circle_center[0].id = avoiding_circle_center[1].id = 0;
	avoiding_circle_radius[0] = avoiding_circle_radius[1] = 0.0;

	// check whether have polar point
	if (subdomain_type != 0) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_spheric_grid, "Software error in PatCC_Delaunay_Voronoi::PatCC_Delaunay_Voronoi");
		for(int i = 0; i < origin_num_points; i ++)
			if(relative_eq_int(fabs(origin_lat_values[i]), 90, 0.0000002))
				num_polar_point_in_origin_coordinate ++;
		if (subdomain_type == 1 && num_polar_point_in_origin_coordinate == 0) { //north pole
			vpolar_local_index = num_points;
			x_ref[num_points] = 0.0;
			y_ref[num_points] = 90.0;
			global_index[num_points++] = -2;
		}
		if (subdomain_type == -1 && num_polar_point_in_origin_coordinate == 0) {
			vpolar_local_index = num_points;
			x_ref[num_points] = 0.0;
			y_ref[num_points] = -90.0;
			global_index[num_points++] = -2;
		}
	}

	if (try_fast_triangulate(original_min_lon, original_max_lon, original_min_lat, original_max_lat)) {
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "fast triangulation for lon-lat grid is enough");
		mark_special_triangles();
	} else if (is_spheric_grid) {
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "fast triangulation for lon-lat grid is not enough, num_points=%d: projection will be conducted with center (%lf %lf)", num_points, original_lon_center, original_lat_center);
		do_projection = true;

		load_polars_info();

		double *projected_lon_values, *projected_lat_values;
		PAT_REAL uv1_x, uv1_y, uv1_z, uv2_x, uv2_y, uv2_z, center_x, center_y, center_z;
		if (is_spheric_grid && !(max_lon_x == 360.0 && min_lon_x == 0.0))
			set_regional(true);
/*
		lonlat2xyz(original_lon_center, original_lat_center, &center_x, &center_y, &center_z);
		calculate_unit_vectors(original_lon_center, original_lat_center, &uv1_x, &uv1_y, &uv1_z, &uv2_x, &uv2_y, &uv2_z);
		projected_lon_values = new double [num_points*2];
		projected_lat_values = new double [num_points*2];
		for (int i = 0; i < num_points; i ++)
			fast_stereographic_projection(x_ref[i], y_ref[i], center_x, center_y, center_z, uv1_x, uv1_y, uv1_z,
										  uv2_x, uv2_y, uv2_z, projected_lon_values[i], projected_lat_values[i]);
*/
		add_points(&x_ref, &y_ref, mask_ref, num_points, max_extend_coordinate_buffer_size);
		triangulate();
//		delete [] projected_lon_values;
//		delete [] projected_lat_values;
//		fast_triangulate = false;
	}
	else {
		set_regional(true);
		add_points(&x_ref, &y_ref, mask_ref, num_points, max_extend_coordinate_buffer_size);
		triangulate();
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "fast triangulation for lon-lat grid is not enough");
	}
	delete [] mask_ref;
}


PatCC_Delaunay_Voronoi::~PatCC_Delaunay_Voronoi()
{
	delete[] global_index;
	if (point_idx_to_buf_idx != NULL)
		delete[] point_idx_to_buf_idx;
	if (triangle_stack != NULL)
		delete[] triangle_stack;
	delete[] x_ref;
	delete[] y_ref;
	delete all_points;
	if (result_triangles != NULL)
		delete [] result_triangles;
	delete triangle_allocator;
	delete edge_allocator;
}


inline bool PatCC_Delaunay_Voronoi::point_in_triangle(double x, double y, PatCC_Triangle* t)
{
	int* v_idx = t->v;

	PatCC_Point p(x, y, all_points[v_idx[0]].is_spheric_grid);
	int ret = p.position_to_triangle(&all_points[v_idx[0]], &all_points[v_idx[1]], &all_points[v_idx[2]]);
	return ret != -1;
}


struct Bound {
	double min_x;
	double max_x;
	double min_y;
	double max_y;
};


inline bool PatCC_Delaunay_Voronoi::point_in_bound(double x, double y, Bound* b)
{
	return x >= b->min_x && x <= b->max_x && y >= b->min_y && y <= b->max_y;
}


Bound* PatCC_Delaunay_Voronoi::make_bounding_box()
{
	unsigned num_triangles = all_leaf_triangles.size();
	Bound* bound = new Bound[num_triangles];

	for (unsigned i = 0; i < num_triangles; i++){
		PatCC_Triangle* lf = all_leaf_triangles[i];
		double x[3], y[3];
		for (int j = 0; j < 3; j++) {
			x[j] = vertex(lf, j)->x;
			y[j] = vertex(lf, j)->y;
		}
		if (x[0] > x[1]) std::swap(x[0], x[1]);
		if (x[1] > x[2]) std::swap(x[1], x[2]);
		if (x[0] > x[1]) std::swap(x[0], x[1]);
		if (y[0] > y[1]) std::swap(y[0], y[1]);
		if (y[1] > y[2]) std::swap(y[1], y[2]);
		if (y[0] > y[1]) std::swap(y[0], y[1]);

		bound[i].min_x = x[0];
		bound[i].max_x = x[2];
		bound[i].min_y = y[0];
		bound[i].max_y = y[2];
	}

	return bound;
}


static inline unsigned hash(double x, double y, double block_size, double min_x, double len_x, double min_y, double len_y)
{
	unsigned x_idx = (x-min_x) / len_x * block_size;
	unsigned y_idx = (y-min_y) / len_y * block_size;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, x_idx >= 0 && x_idx < block_size && y_idx >= 0 && y_idx < block_size);
	return x_idx+y_idx*block_size;
}


static inline void hash(const Bound* bound, double block_size, double min_x, double len_x, double min_y, double len_y,
				   unsigned* x_idx_bgn, unsigned* x_idx_end, unsigned* y_idx_bgn, unsigned* y_idx_end)
{
	*x_idx_bgn = (bound->min_x-min_x) / len_x * block_size; // FIXME: optimise
	*y_idx_bgn = (bound->min_y-min_y) / len_y * block_size;
	*x_idx_end = (bound->max_x-min_x) / len_x * block_size;
	*y_idx_end = (bound->max_y-min_y) / len_y * block_size;
	if (*x_idx_end == block_size) (*x_idx_end)--;
	if (*y_idx_end == block_size) (*y_idx_end)--;
}


#define PAT_TRIANGLES_PER_BLOCK (50)
#define PAT_MAX_BLOCK (10000)
void PatCC_Delaunay_Voronoi::distribute_initial_points(const double* x, const double* y, int num, int** output_nexts)
{
	int j;
	double min_x = all_points[0].x;
	double max_y = all_points[0].y;
	double max_x = all_points[2].x;
	double min_y = all_points[2].y;

	double len_x = max_x - min_x;
	double len_y = max_y - min_y;


	unsigned num_triangles = all_leaf_triangles.size();
	int* nexts = new int[num];

	memset(nexts, -1, num*sizeof(int));

	unsigned block_size = std::sqrt(std::min(num_triangles / PAT_TRIANGLES_PER_BLOCK, (unsigned)PAT_MAX_BLOCK));
	
	if (block_size*block_size > 2 && false) {
		Bound* bound = make_bounding_box();

		for (unsigned i = 0; i < num_triangles; i++) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, bound[i].min_x >= min_x);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, bound[i].max_x <= max_x);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, bound[i].min_y >= min_y);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, bound[i].max_y <= max_y);
		}

		unsigned* num_including_points = new unsigned[block_size*block_size]();
		unsigned* num_including_triangles = new unsigned[block_size*block_size]();

		/* first scan: counting */

		for (int i = 0; i < num; i++) {
			unsigned idx = hash(x[i], y[i], block_size, min_x, len_x, min_y, len_y);
			num_including_points[idx]++;
		}

		for (unsigned i = 0; i < num_triangles; i++) {
			unsigned x_idx_bgn, y_idx_bgn, x_idx_end, y_idx_end;
			hash(&bound[i], block_size, min_x, len_x, min_y, len_y, &x_idx_bgn, &x_idx_end, &y_idx_bgn, &y_idx_end);
			for (unsigned j = y_idx_bgn; j <= y_idx_end; j++)
				for (unsigned k = x_idx_bgn; k <= x_idx_end; k++) {
					unsigned idx = k + j * block_size;
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, j < block_size);
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, k < block_size);
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, idx < block_size*block_size);
					num_including_triangles[idx]++;
				}
		}

		/* allocing memory */
		unsigned** including_points = new unsigned*[block_size*block_size];
		unsigned** including_triangles = new unsigned*[block_size*block_size];
		for (unsigned i = 0; i < block_size*block_size; i++) {
			if (num_including_points[i] > 0) {
				including_points[i] = new unsigned[num_including_points[i]];
				including_triangles[i] = new unsigned[num_including_triangles[i]];
			} else {
				including_points[i] = NULL;
				including_triangles[i] = NULL;
			}
		}

		/* second scan: distributing points and triangles into mesh */
		memset(num_including_points, 0, sizeof(unsigned)*block_size*block_size);
		memset(num_including_triangles, 0, sizeof(unsigned)*block_size*block_size);

		for (int i = 0; i < num; i++) {
			unsigned idx = hash(x[i], y[i], block_size, min_x, len_x, min_y, len_y);
			including_points[idx][num_including_points[idx]] = i;
			num_including_points[idx]++;
		}

		for (unsigned i = 0; i < num_triangles; i++) {
			unsigned x_idx_bgn, y_idx_bgn, x_idx_end, y_idx_end;
			hash(&bound[i], block_size, min_x, len_x, min_y, len_y, &x_idx_bgn, &x_idx_end, &y_idx_bgn, &y_idx_end);
			for (unsigned j = y_idx_bgn; j <= y_idx_end; j++)
				for (unsigned k = x_idx_bgn; k <= x_idx_end; k++) {
					unsigned idx = k + j * block_size;
					if (including_triangles[idx]) {
						including_triangles[idx][num_including_triangles[idx]] = i;
						num_including_triangles[idx]++;
					}
				}
		}

		/* distributing points into triangles in the same mesh */
		int distributed = 0;
		for (unsigned m = 0; m < block_size*block_size; m++) {
			if (including_points[m] == NULL)
				continue;

			unsigned* p_idxs = including_points[m];
			unsigned* t_idxs = including_triangles[m];
			unsigned  p_num = num_including_points[m];
			unsigned  t_num = num_including_triangles[m];

			for (unsigned i = 0; i < p_num; i++) {
				unsigned p_idx = p_idxs[i];
				for (unsigned j = 0; j < t_num; j++) {
					unsigned t_idx = t_idxs[j];
					if (!all_leaf_triangles[t_idx]->is_leaf)
						continue;

					if (!point_in_bound(x[p_idx], y[p_idx], &bound[t_idx]))
						continue;

					if (point_in_triangle(x[p_idx], y[p_idx], all_leaf_triangles[t_idx])) {
						if (all_leaf_triangles[t_idx]->remained_points_head == -1)
							all_leaf_triangles[t_idx]->remained_points_head = all_leaf_triangles[t_idx]->remained_points_tail = p_idx;
						else
							all_leaf_triangles[t_idx]->remained_points_tail = nexts[all_leaf_triangles[t_idx]->remained_points_tail] = p_idx;
						distributed++;
						break;
					}
				}
			}
		}

		if (report_error_enabled) {
			unsigned total = 0;
			for (unsigned m = 0; m < block_size*block_size; m++) {
				total += num_including_points[m];
			}
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num == total, "Software error in PatCC_Delaunay_Voronoi::distribute_initial_points");
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num == distributed, "Software error in PatCC_Delaunay_Voronoi::distribute_initial_points");

		/* freeing memory */
		for (unsigned i = 0; i < block_size*block_size; i++) {
			if (including_points[i] != NULL)
				delete[] including_points[i];
			if (including_triangles[i] != NULL)
				delete[] including_triangles[i];
		}
		delete[] including_points;
		delete[] including_triangles;

		delete[] num_including_points;
		delete[] num_including_triangles;


		delete[] bound;
	} else {	
		int distributed = 0;
		for (int i = 0; i < num; i++) {
			while (true) {
				for (j = 0; j < num_triangles; j++) {
					if (!all_leaf_triangles[j]->is_leaf)
						continue;

					if (point_in_triangle(x[i], y[i], all_leaf_triangles[j])) {
						if (all_leaf_triangles[j]->remained_points_head == -1)
							all_leaf_triangles[j]->remained_points_head = all_leaf_triangles[j]->remained_points_tail = i;
						else
							all_leaf_triangles[j]->remained_points_tail = nexts[all_leaf_triangles[j]->remained_points_tail] = i;
						distributed++;
						break;
					} 
				}	
				if (j < num_triangles)
					break;
				EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "will enlarge the triangulation box (%lf %lf, %lf %lf, %lf %lf, %lf %lf) for the point (%lf %lf)", all_points[0].x, all_points[0].y, all_points[1].x, all_points[1].y, all_points[2].x, all_points[2].y, all_points[3].x, all_points[3].y, x[i], y[i]);
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_spheric_grid, "Software error in in \"PatCC_Delaunay_Voronoi::distribute_initial_points\"");
				if (original_max_lat == 90.0) {
					all_points[0].y -= 1;
					all_points[1].y -= 1;
					all_points[2].y -= 1;
					all_points[3].y -= 1;	
					EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "further enlarge grid domain at a north pole case");
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, all_points[0].y > 10, "Software error in in \"PatCC_Delaunay_Voronoi::distribute_initial_points\"");
				}
				else if (original_min_lat == -90.0) {
					all_points[0].y += 1;
					all_points[1].y += 1;
					all_points[2].y += 1;
					all_points[3].y += 1;
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, all_points[0].y < -10, "Software error in in \"PatCC_Delaunay_Voronoi::distribute_initial_points\"");
					EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "further enlarge grid domain at a south pole case");
				}
				else if (all_points[0].y < 0) {  // max lat
					all_points[0].y += 1;
					all_points[3].y += 1;
					EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "further enlarge grid domain at a north boundary case");
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, all_points[0].y <= 89.7, "Software error in in \"PatCC_Delaunay_Voronoi::distribute_initial_points\"");
				}
				else if (all_points[1].y > 0) { // min lat
					all_points[1].y -= 1;
					all_points[2].y -= 1;
					EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "further enlarge grid domain at a south boundary case");
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, all_points[2].y >= -89.7, "Software error in in \"PatCC_Delaunay_Voronoi::distribute_initial_points\"");
				}
				else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in in \"PatCC_Delaunay_Voronoi::distribute_initial_points\"");
				for (int k = 0; k < 4; k ++)
					all_points[k].transform_coord_from_sphere_to_3D();
			}				
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num == distributed, "in \"PatCC_Delaunay_Voronoi::distribute_initial_points\" distributed %d, expected %d\n", distributed, num);
	}

	*output_nexts = nexts;
}


void PatCC_Delaunay_Voronoi::extend_coordinate_buffer(double** x, double** y, int& max_coordinate_buffer_size)
{
	int extend_max_coordinate_buffer_size = max_coordinate_buffer_size * 2;
	double* x_tmp = new double [extend_max_coordinate_buffer_size];
	double* y_tmp = new double [extend_max_coordinate_buffer_size];
	memcpy(x_tmp, *x, max_coordinate_buffer_size*sizeof(double));
	memcpy(y_tmp, *y, max_coordinate_buffer_size*sizeof(double));
	if (max_coordinate_buffer_size > 0) {
		delete [] *x;
		delete [] *y;
	}
	*x = x_tmp;
	*y = y_tmp;
	max_coordinate_buffer_size = extend_max_coordinate_buffer_size;
}


void PatCC_Delaunay_Voronoi::add_point_in_rectangle(std::pair<double, double> *buckets, double** x, double** y, int& max_coordinate_buffer_size, double bounding_min_x, double bounding_max_x, double bounding_min_y, double bounding_max_y, int num_x_virtual_bound_points, int num_y_virtual_bound_points, int& num, bool lon_first) {
	
	double x_density = (bounding_max_x - bounding_min_x) / (double) num_x_virtual_bound_points;
	double y_density = (bounding_max_y - bounding_min_y) / (double) num_y_virtual_bound_points;

	int buf_idx_cur = 0;
	for (int i = 0; i < num_x_virtual_bound_points; i ++) {
		(buckets[i]).first = bounding_max_y;
		(buckets[i]).second = bounding_min_y;
	}

	for (int i = 0; i < num; i++) {
		int index = (int)(((*x)[i] - bounding_min_x) / x_density);
		index = std::min(std::max(0, index), num_x_virtual_bound_points-1);
		buckets[index].first = std::min(buckets[index].first, (*y)[i]);
		buckets[index].second = std::max(buckets[index].second, (*y)[i]);
	}

	for (int i = 0; i < num_x_virtual_bound_points; i ++) {
		if (buckets[i].first == bounding_max_y && buckets[i].second == bounding_min_y)
			continue;
//		double num_min_added_points = (buckets[i].first - bounding_min_y)/y_density;
		double num_min_added_points = 1;
		double new_y_density = (buckets[i].first - bounding_min_y) / (double)num_min_added_points;
		for (int j = 1; j <= (int)num_min_added_points; j ++) {
			if (is_spheric_grid && original_max_lat == 90.0 && (lon_first && buckets[i].first-new_y_density*j > original_min_lat || !lon_first && bounding_min_x+i*x_density+x_density/2 > original_min_lat))
				continue;
			if (is_spheric_grid && original_min_lat == -90.0 && (lon_first && buckets[i].first-new_y_density*j < original_max_lat || !lon_first && bounding_min_x+i*x_density+x_density/2 < original_max_lat))
				continue;
			if (is_spheric_grid && (lon_first && fabs(buckets[i].first - new_y_density*j) > 89.5 || !lon_first && fabs(bounding_min_x + i*x_density + x_density/2) > 89.5))
				continue;
			if (num+buf_idx_cur >= max_coordinate_buffer_size)
				extend_coordinate_buffer(x, y, max_coordinate_buffer_size);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num+buf_idx_cur < max_coordinate_buffer_size, "Software error in PatCC_Delaunay_Voronoi::add_point_in_rectangle");
			(*x)[num+buf_idx_cur] = bounding_min_x + i*x_density + x_density/2;
			(*y)[num+buf_idx_cur] = buckets[i].first - new_y_density*j;
			buf_idx_cur ++;
		}

//		double num_max_added_points = (bounding_max_y - buckets[i].second)/y_density;
		double num_max_added_points = 1;
		new_y_density = (bounding_max_y - buckets[i].second) / (double)num_max_added_points;
		for (int j = 1; j <= (int)num_max_added_points; j ++) {
			if (is_spheric_grid && original_max_lat == 90.0 && (lon_first && buckets[i].second + new_y_density*j > original_min_lat || !lon_first && bounding_min_x + i*x_density + x_density/2 > original_min_lat))
				continue;
			if (is_spheric_grid && original_min_lat == -90.0 && (lon_first && buckets[i].second + new_y_density*j < original_max_lat || !lon_first && bounding_min_x + i*x_density + x_density/2 < original_max_lat))
				continue;
			if (is_spheric_grid && (lon_first && fabs(buckets[i].second + new_y_density*j) > 89.5 || !lon_first && fabs(bounding_min_x + i*x_density + x_density/2) > 89.5))
				continue;
			if (num+buf_idx_cur >= max_coordinate_buffer_size)
				extend_coordinate_buffer(x, y, max_coordinate_buffer_size);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num+buf_idx_cur < max_coordinate_buffer_size, "Software error in PatCC_Delaunay_Voronoi::add_point_in_rectangle");
			(*x)[num+buf_idx_cur] = bounding_min_x + i*x_density + x_density/2;
			(*y)[num+buf_idx_cur] = buckets[i].second + new_y_density*j;
			buf_idx_cur ++;
		}
	}
	num += buf_idx_cur;
}


int PatCC_Delaunay_Voronoi::enlarge_super_rectangle(double** x, double** y, int num, int& max_coordinate_buffer_size)
{
	double min_x = 1e10;
	double max_x = -1e10;
	double min_y = 1e10;
	double max_y = -1e10;
	for (int i = 0; i < num; i++) {
		if ((*x)[i] < min_x) min_x = (*x)[i];
		if ((*x)[i] > max_x) max_x = (*x)[i];
		if ((*y)[i] < min_y) min_y = (*y)[i];
		if ((*y)[i] > max_y) max_y = (*y)[i];
	}
	if (is_spheric_grid && (original_min_lat == -90.0 || original_max_lat == 90.0)) {
		min_x = original_min_lon;
		max_x = original_max_lon;
		min_y = original_min_lat;
		max_y = original_max_lat;
	}

	double x_expected_points = std::sqrt((max_x - min_x)/(max_y - min_y)*(double)num);
	double y_expected_points = std::sqrt((max_y - min_y)/(max_x - min_x)*(double)num);
	double density = (max_x - min_x) / x_expected_points;
	int num_expanded_level_factor = 3;
	int inner_virtual_bound_factor = 2;
	double delta = density * num_expanded_level_factor * inner_virtual_bound_factor;
	int min_num_virtual_bound_points = 3;
	x_expected_points = (max_x - min_x + delta*2) / density;
	y_expected_points = (max_y - min_y + delta*2) / density;
	int num_x_virtual_bound_points = std::max((int)(x_expected_points/inner_virtual_bound_factor), min_num_virtual_bound_points);
	int num_y_virtual_bound_points = std::max((int)(y_expected_points/inner_virtual_bound_factor), min_num_virtual_bound_points);
	std::pair<double, double> *x_buckets = new std::pair<double, double> [num_x_virtual_bound_points];
	std::pair<double, double> *y_buckets = new std::pair<double, double> [num_y_virtual_bound_points];

	if (is_spheric_grid && (min_y == -90.0 || max_y == 90.0)) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, min_x == 0 && max_x == 360.0, "Software error in PatCC_Delaunay_Voronoi::enlarge_super_rectangle");
		if (min_y == -90.0) {
			all_points[0].y = max_y + delta;
			all_points[1].y = max_y + delta;
			all_points[2].y = max_y + delta;
			all_points[3].y = max_y + delta;
			all_points[0].x = 180;
			all_points[1].x = 90;
			all_points[2].x = 0;
			all_points[3].x = 270;
		}
		else {
			all_points[0].y = min_y-delta;
			all_points[1].y = min_y-delta;
			all_points[2].y = min_y-delta;
			all_points[3].y = min_y-delta;
			all_points[0].x = 0;
			all_points[1].x = 90;
			all_points[2].x = 180;
			all_points[3].x = 270;
		}
	}
	else {
		double enlarged_min_y = min_y - delta;
		double enlarged_max_y = max_y + delta;
		if (is_spheric_grid) {
			enlarged_min_y = std::min(std::max(enlarged_min_y, -89.5), min_y);
			enlarged_max_y = std::max(std::min(enlarged_max_y, 89.5), max_y);
		}
		all_points[0].x = min_x-delta;
		all_points[0].y = enlarged_max_y;
		all_points[1].x = min_x-delta;
		all_points[1].y = enlarged_min_y;
		all_points[2].x = max_x+delta;
		all_points[2].y = enlarged_min_y;
		all_points[3].x = max_x+delta;
		all_points[3].y = enlarged_max_y;
	}
	for (int i = 0; i < 4; i ++)
		all_points[i].transform_coord_from_sphere_to_3D();

	add_point_in_rectangle(x_buckets, x, y, max_coordinate_buffer_size, min_x-0.95*delta, max_x+0.95*delta, min_y-0.95*delta, max_y+0.95*delta, num_x_virtual_bound_points, num_y_virtual_bound_points, num, true);
	add_point_in_rectangle(y_buckets, y, x, max_coordinate_buffer_size, min_y-0.95*delta, max_y+0.95*delta, min_x-0.95*delta, max_x+0.95*delta, num_y_virtual_bound_points, num_x_virtual_bound_points, num, false);

	delete [] x_buckets;
	delete [] y_buckets;
	return num;
}


/*
 *      0 ---- 3
 *      |      |
 *      |      |
 *      1 ---- 2
 */
void PatCC_Delaunay_Voronoi::initialize(int num){
	int triangles_count_estimate = 2*(num+PAT_NUM_LOCAL_VPOINTS);
	all_leaf_triangles.reserve(triangles_count_estimate);

	stack_size = triangles_count_estimate * 2;
	triangle_stack = new PatCC_Triangle*[stack_size];

	max_points = num*3/2 + PAT_NUM_LOCAL_VPOINTS;
	all_points = (PatCC_Point *)::operator new(max_points * sizeof(PatCC_Point));

	new(&all_points[0]) PatCC_Point(-0.1,  0.1, -1, is_spheric_grid, -1, -1);
	new(&all_points[1]) PatCC_Point(-0.1, -0.1, -1, is_spheric_grid, -1, -1);
	new(&all_points[2]) PatCC_Point( 0.1, -0.1, -1, is_spheric_grid, -1, -1);
	new(&all_points[3]) PatCC_Point( 0.1,  0.1, -1, is_spheric_grid, -1, -1);

	num_points = PAT_NUM_LOCAL_VPOINTS;

	all_leaf_triangles.push_back(allocate_triangle(allocate_edge(0, 1), allocate_edge(1, 2), allocate_edge(2, 0)));
	all_leaf_triangles.push_back(allocate_triangle(allocate_edge(0, 2), allocate_edge(2, 3), allocate_edge(3, 0)));

	PatCC_Edge* e1 = all_leaf_triangles[0]->edge[2];
	PatCC_Edge* e2 = all_leaf_triangles[1]->edge[0];
	e1->twin_edge = e2;
	e2->twin_edge = e1;
}


void PatCC_Delaunay_Voronoi::extend_points_buffer(int introduced_points)
{
	int full_size = num_points + introduced_points;
	int realloc_size = full_size + introduced_points * PDLN_EXPECTED_EXPANDING_TIMES;

	void* tmp_buf = ::operator new(realloc_size * sizeof(PatCC_Point));
	memcpy(tmp_buf, all_points, num_points * sizeof(PatCC_Point));

	::operator delete(all_points);
	all_points = (PatCC_Point *)tmp_buf;
	max_points = realloc_size;
}


void PatCC_Delaunay_Voronoi::add_points(double** x, double** y, const bool* mask, int num, int& max_coordinate_buffer_size)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, *x != NULL, "Software error in PatCC_Delaunay_Voronoi::add_points");
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, *y != NULL, "Software error in PatCC_Delaunay_Voronoi::add_points");

	if(stack_size == 0)
		initialize(num);
	else if (num_points + num > max_points)
		extend_points_buffer(num);
	int num_enlarged_points = enlarge_super_rectangle(x, y, num, max_coordinate_buffer_size);
	if (num_points + num_enlarged_points > max_points)
		extend_points_buffer(num_enlarged_points);

	int *nexts;

	distribute_initial_points(*x, *y, num_enlarged_points, &nexts);
	dirty_triangles_count = 0;
	int buf_idx_cur = num_points;
	int local_idx_start = num_points - PAT_NUM_LOCAL_VPOINTS;

	for (unsigned i = 0; i < all_leaf_triangles.size(); i++) {
		int head = buf_idx_cur;
		for (int p = all_leaf_triangles[i]->remained_points_head; p != -1; p = nexts[p]) {
			if (p < num)
				new(&all_points[buf_idx_cur]) PatCC_Point((*x)[p], (*y)[p], local_idx_start+p, mask ? mask[p] : true, is_spheric_grid,  buf_idx_cur+1, buf_idx_cur-1);
			else
				new(&all_points[buf_idx_cur]) PatCC_Point((*x)[p], (*y)[p], -1, false, is_spheric_grid, buf_idx_cur+1, buf_idx_cur-1);
			buf_idx_cur++;
		}
		if (head != buf_idx_cur) {
			all_points[head].prev = -1;
			all_points[buf_idx_cur-1].next = -1;
			all_leaf_triangles[i]->set_remained_points(head, buf_idx_cur-1);
			push(&dirty_triangles_count, all_leaf_triangles[i]);
		}
	}
	for (unsigned i = 0; i < all_leaf_triangles.size(); i++)
		all_leaf_triangles[i]->calulate_circum_circle(&all_points[all_leaf_triangles[i]->v[0]], &all_points[all_leaf_triangles[i]->v[1]], &all_points[all_leaf_triangles[i]->v[2]]);

	num_points += num_enlarged_points;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, buf_idx_cur == num_points, "Software error in PatCC_Delaunay_Voronoi::add_points");
	
	if (report_error_enabled) {
		bool* check = new bool[num_points]();
		for (int i = 4; i < num_points; i++) {
			if (!is_real_point(i))
				continue;
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, check[all_points[i].id] == 0, "Software error in add_points, %d all_points[%d].id=%d", check[all_points[i].id], i, all_points[i].id);
			check[all_points[i].id] = 1;
		}
		delete [] check;
	}

	delete[] nexts;
}


// only check coordinate value after projection, whether coordinate value before projection cannot be check here
bool on_a_line_projected(const PatCC_Point* pt1, const PatCC_Point* pt2, const PatCC_Point* pt3)
{
	double distance[3];
	
	distance[0] = pt1->calculate_distance_MD(pt2->x, pt2->y);
	distance[1] = pt2->calculate_distance_MD(pt3->x, pt3->y);
	distance[2] = pt3->calculate_distance_MD(pt1->x, pt1->y);

	return (relative_eq((distance[0]+distance[1]), distance[2]) || relative_eq((distance[0]+distance[2]), distance[1]) || relative_eq((distance[1]+distance[2]), distance[0]));
}

bool PatCC_Delaunay_Voronoi::on_a_line_lon_lat(const PatCC_Point* pt1, const PatCC_Point* pt2, const PatCC_Point* pt3)
{
	return ((float_eq(x_ref[pt1->id], x_ref[pt2->id]) && float_eq(x_ref[pt2->id], x_ref[pt3->id])) || (float_eq(y_ref[pt1->id], y_ref[pt2->id]) && float_eq(y_ref[pt2->id], y_ref[pt3->id])));

}


void PatCC_Delaunay_Voronoi::generate_all_result_triangles()
{
	result_triangles = new Triangle_inline[all_leaf_triangles.size()];

	get_triangles_in_region(original_min_lon, original_max_lon, original_min_lat, original_max_lat, result_triangles, &num_result_triangles, all_leaf_triangles.size()); 
	sort_points_in_triangle(result_triangles, num_result_triangles);
//	sort_triangles(result_triangles, num_result_triangles);

}


void PatCC_Delaunay_Voronoi::check_and_add_vertex_to_map(int idx, double value_x, double value_y, char **temp_array_buffer,long &buffer_max_size,long &buffer_content_size, int *num_vertexes_of_kernel_points)
{
	int point_id = all_points[idx].id;


	if (!is_real_point(idx) || !(point_id < kernel_grid_size && point_id >= 0))
		return;

	num_vertexes_of_kernel_points[point_id] ++;
	write_data_into_array_buffer(&point_id, sizeof(int), temp_array_buffer, buffer_max_size, buffer_content_size);
	write_data_into_array_buffer(&value_x, sizeof(double), temp_array_buffer, buffer_max_size, buffer_content_size);
	write_data_into_array_buffer(&value_y, sizeof(double), temp_array_buffer, buffer_max_size, buffer_content_size);
}


bool PatCC_Delaunay_Voronoi::is_real_point(int idx) 
{
	return all_points[idx].id >= 0 && global_index[all_points[idx].id] >= 0;
}

void PatCC_Delaunay_Voronoi::delete_overlap_point_in_cell_vertexes(double* vertex_lons, double* vertex_lats, long *vertex_index, int& num_vertexes_points) {
	int i, i_prev;
	int num_active_vertexes = 0;

	for (i = 0; i < num_vertexes_points; i ++) {
		i_prev = (i - 1 + num_vertexes_points) % num_vertexes_points;
		if (!(std::fabs(vertex_lons[i] - vertex_lons[i_prev]) < 1e-5 && std::fabs(vertex_lats[i] - vertex_lats[i_prev]) < 1e-5)) {
			vertex_lons[num_active_vertexes] = vertex_lons[i];
			vertex_lats[num_active_vertexes] = vertex_lats[i];
			vertex_index[num_active_vertexes ++] = vertex_index[i];
		}
	}
	num_vertexes_points = num_active_vertexes;
}


void PatCC_Delaunay_Voronoi::delete_point_on_common_line(double* vertex_lons, double* vertex_lats, long *vertex_index, int& num_vertexes_points) {
	int i = 0, j, i_prev, i_next;
	double distance12, distance13, distance23;
	int num_active_vertexes = num_vertexes_points;

	while (i < num_active_vertexes) {
		i_prev = (i - 1 + num_active_vertexes) % num_active_vertexes;
		i_next = (i + 1) % num_active_vertexes;
		distance12 = calculate_distance(vertex_lons[i], vertex_lats[i], vertex_lons[i_prev], vertex_lats[i_prev]);
		distance13 = calculate_distance(vertex_lons[i_prev], vertex_lats[i_prev], vertex_lons[i_next], vertex_lats[i_next]);
		distance23 = calculate_distance(vertex_lons[i], vertex_lats[i], vertex_lons[i_next], vertex_lats[i_next]);
		if (relative_eq((distance12 + distance23), distance13)) {
			for (j = i+1; j < num_active_vertexes; j ++) {
				vertex_lons[j-1] = vertex_lons[j];
				vertex_lats[j-1] = vertex_lats[j];
				vertex_index[j-1] = vertex_index[j];
			}
			vertex_index[--num_active_vertexes] = -1;
		}
		else i ++;
	}

	for (int i = num_active_vertexes; i < num_vertexes_points; i ++)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, vertex_index[i] == -1, "Software error in PatCC_Delaunay_Voronoi::delete_point_on_common_line");
	
	num_vertexes_points = num_active_vertexes;
}


//https://gis.stackexchange.com/questions/10808/manually-transforming-rotated-lat-lon-to-regular-lat-lon
void PatCC_Delaunay_Voronoi::rotate_sphere_coordinate(double* lon_original, double* lat_original, int num_vertexes_points, double theta_double, double phi_double, int type)
{
	double theta = DEGREE_TO_RADIAN(theta_double);
	double cosl_theta, sinl_theta, cosl_lat;
	double x_new, z_new, x, y, z, lon_new, lat_new;

	if (!is_spheric_grid || num_vertexes_points == 0 || fast_mode)
		return;

	cosl_theta = cos(theta);
    sinl_theta = sin(theta);
	if (type == 1) {
		for (int i = 0; i < num_vertexes_points; i ++)
			lon_original[i] -= phi_double;
	}

	for (int i = 0; i < num_vertexes_points; i ++) {
		lon_original[i] = DEGREE_TO_RADIAN(lon_original[i]);
		lat_original[i] = DEGREE_TO_RADIAN(lat_original[i]);
		cosl_lat = cos(lat_original[i]);
		x = cos(lon_original[i])*cosl_lat;
		y = sin(lon_original[i])*cosl_lat;
		z = sin(lat_original[i]);
		x_new = cosl_theta*x + sinl_theta*z;
		z_new = -sinl_theta*x + cosl_theta*z;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(float_eq(x_new, 0.0) && float_eq(y, 0.0)), "Software error in rotate sphere coordinate");
		lon_original[i] = RADIAN_TO_DEGREE(atan2(y, x_new));
		lat_original[i] = RADIAN_TO_DEGREE(asin(z_new));
	}

	if (type == 2) {
		for (int i = 0; i < num_vertexes_points; i ++) {
			lon_original[i] += phi_double;
			if (lon_original[i] < 0)
				lon_original[i] += 360.0;
			if (lon_original[i] >= 360.0)
				lon_original[i] -= 360.0;
		}
	}
}


void PatCC_Delaunay_Voronoi::generate_Voronoi_diagram(int &max_num_voronoi_diagram_vertex, double **vertex_lon_values, double **vertex_lat_values)
{
	double vertex_value_x, vertex_value_y, max_lon, difference1, difference2;
	int num_real_vertexes_index, real_vertexes_index[3], i, j, k, max_global_index, max_global_index_local_index, point_id, num_kernel_vertexes;
	double temp_lon[3], temp_lat[3];
	PatCC_Point pt[3], temp_pt;
	PatCC_Triangle* temp_triangle = new PatCC_Triangle();
	char *vertex_info_buf = NULL;
	long vertex_info_buf_iter, vertex_info_buf_size;
	double *temp_vertex_lon_values, *temp_vertex_lat_values;


	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, kernel_grid_size > 0, "Software error in PatCC_Delaunay_Voronoi::generate_Voronoi_diagram");

	int *num_vertexes_of_kernel_points = new int [kernel_grid_size];
	for (int i = 0; i < kernel_grid_size; i ++)
		num_vertexes_of_kernel_points[i] = 0;

	for (i = 0; i < all_leaf_triangles.size(); i ++) {
		if (!all_leaf_triangles[i]->is_leaf)
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, fast_mode || is_triangle_legal(all_leaf_triangles[i]), "Software error in generate_Voronoi_diagram");

		for (j = 0, num_kernel_vertexes = 0; j < 3; j ++)
			if (all_points[all_leaf_triangles[i]->v[j]].id >= 0 && all_points[all_leaf_triangles[i]->v[j]].id < kernel_grid_size)
				num_kernel_vertexes ++;
		if (num_kernel_vertexes == 0)
			continue;

		for (j = 0, num_real_vertexes_index = 0; j < 3; j ++)
			if (is_real_point(all_leaf_triangles[i]->v[j]))
				real_vertexes_index[num_real_vertexes_index ++] = j;

		max_global_index = -999999;
		if (num_real_vertexes_index == 3) {
			for (j = 0; j < 3; j ++) {
				point_id = all_points[all_leaf_triangles[i]->v[j]].id;
				temp_lon[j] = origin_lon_values[point_id];
				temp_lat[j] = origin_lat_values[point_id];
				if (global_index[point_id] > max_global_index) { 
					max_global_index = global_index[point_id];
					max_global_index_local_index = j;
				}
			}

			if (num_polar_point_in_origin_coordinate == 1)
				for (j = 0; j < 3; j ++)
					if (relative_eq_int(fabs(temp_lat[j]), 90, 0.0000002)) {
						match_degree_values(temp_lon[(j+1)%3], temp_lon[(j+2)%3]);
						temp_lon[j] = (temp_lon[(j+1)%3] + temp_lon[(j+2)%3]) / 2;
					}
			convert_cyclic_triangle_point(temp_lon, 3);
			difference1 = temp_lat[max_global_index_local_index];
			difference2 = temp_lon[max_global_index_local_index];
			rotate_sphere_coordinate(temp_lon, temp_lat, 3, difference1, difference2, 1);

			for (j = 0; j < 3; j ++) {
				pt[j].x = temp_lon[j];
				pt[j].y = temp_lat[j];
			}
			temp_triangle->calulate_circum_circle(&pt[0], &pt[1], &pt[2]);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !pt[0].is_spheric_grid, "Software error in calculate circum circle, %d has %d real point");
			vertex_value_x = temp_triangle->circum_center[0];
			vertex_value_y = temp_triangle->circum_center[1];
			rotate_sphere_coordinate(&vertex_value_x, &vertex_value_y, 1, -difference1, difference2, 2);

			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !isnan(vertex_value_x) && !isnan(vertex_value_y), "Software error in calculate circum circle, %d has %d real point", i, num_real_vertexes_index);
		}
		else {
			int virtual_point_idx = 3 - real_vertexes_index[1] - real_vertexes_index[0];
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_spheric_grid && relative_eq_int(fabs(y_ref[all_points[all_leaf_triangles[i]->v[virtual_point_idx]].id]), 90.0, 0.0000002) && num_real_vertexes_index == 2, "Software error in generate_Voronoi_diagram, all_points[all_leaf_triangles[i]->v[virtual_point_idx]].y = %lf, num_real_vertexes_index = %d", y_ref[all_points[all_leaf_triangles[i]->v[virtual_point_idx]].id], num_real_vertexes_index);

			for(j = 0; j < 2; j ++) {
 				temp_lon[j] = origin_lon_values[all_points[all_leaf_triangles[i]->v[real_vertexes_index[j]]].id];
				temp_lat[j] = origin_lat_values[all_points[all_leaf_triangles[i]->v[real_vertexes_index[j]]].id];
			}
			convert_cyclic_triangle_point(temp_lon, 2);
			vertex_value_x = (temp_lon[0] + temp_lon[1]) / 2;
			vertex_value_y = y_ref[all_points[all_leaf_triangles[i]->v[virtual_point_idx]].id];
		}

		for (j = 0; j < 3; j ++)
			check_and_add_vertex_to_map(all_leaf_triangles[i]->v[j], vertex_value_x, vertex_value_y, &vertex_info_buf, vertex_info_buf_size, vertex_info_buf_iter, num_vertexes_of_kernel_points);

		for (j = 0; j < 3; j ++) {
			PatCC_Edge *edge = all_leaf_triangles[i]->edge[j];
			if (edge->twin_edge == NULL) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_real_point(edge->head) && is_real_point(edge->tail), "Software error in generate_Voronoi_diagram, the edge whose twin edge is null has a virtual point.");
				temp_lon[0] = origin_lon_values[all_points[edge->head].id];
				temp_lon[1] = origin_lon_values[all_points[edge->tail].id];
				convert_cyclic_triangle_point(temp_lon, 2);

				vertex_value_x = (temp_lon[0] + temp_lon[1]) / 2;
				vertex_value_y = (origin_lat_values[all_points[edge->head].id] + origin_lat_values[all_points[edge->tail].id]) / 2;
				check_and_add_vertex_to_map(edge->head, vertex_value_x, vertex_value_y, &vertex_info_buf, vertex_info_buf_size, vertex_info_buf_iter, num_vertexes_of_kernel_points);
				check_and_add_vertex_to_map(edge->tail, vertex_value_x, vertex_value_y, &vertex_info_buf, vertex_info_buf_size, vertex_info_buf_iter, num_vertexes_of_kernel_points);
			}
		}
	}
	delete temp_triangle;

	max_num_voronoi_diagram_vertex = 0;
	for (int i = 0; i < kernel_grid_size; i ++) {
		if (num_vertexes_of_kernel_points[i] > max_num_voronoi_diagram_vertex)
			max_num_voronoi_diagram_vertex = num_vertexes_of_kernel_points[i];
		num_vertexes_of_kernel_points[i] = 0;
	}
	max_num_voronoi_diagram_vertex ++;
	temp_vertex_lon_values = new double [kernel_grid_size*max_num_voronoi_diagram_vertex];
	temp_vertex_lat_values = new double [kernel_grid_size*max_num_voronoi_diagram_vertex];
	for (int i = 0; i < vertex_info_buf_iter; i += sizeof(int)+2*sizeof(double)) {
		point_id = ((int*)(vertex_info_buf+i))[0];
		vertex_value_x = ((double*)(vertex_info_buf+i+sizeof(int)))[0];
		vertex_value_y = ((double*)(vertex_info_buf+i+sizeof(int)))[1];
		temp_vertex_lon_values[point_id*max_num_voronoi_diagram_vertex+num_vertexes_of_kernel_points[point_id]] = vertex_value_x;
		temp_vertex_lat_values[point_id*max_num_voronoi_diagram_vertex+num_vertexes_of_kernel_points[point_id]] = vertex_value_y;
		num_vertexes_of_kernel_points[point_id] ++;
	}
	double *vertex_lon_check = new double [max_num_voronoi_diagram_vertex];
	double *vertex_lat_check = new double [max_num_voronoi_diagram_vertex];
	double *vertex_lon = new double [max_num_voronoi_diagram_vertex];
	double *vertex_lat = new double [max_num_voronoi_diagram_vertex];
	long *vertex_index_check = new long [max_num_voronoi_diagram_vertex];
	for (i = 0; i < kernel_grid_size; i ++) {
		int num_vertexes_points = num_vertexes_of_kernel_points[i];
		int num_polar_points;
		double center_x, center_y;
		if (num_vertexes_points == 0)
			continue;
		for (int j = 0; j < num_vertexes_points; j ++) {
		    vertex_lon[j] = temp_vertex_lon_values[i*max_num_voronoi_diagram_vertex+j];
		    vertex_lat[j] = temp_vertex_lat_values[i*max_num_voronoi_diagram_vertex+j];
		}
		vertex_lon[num_vertexes_points] = origin_lon_values[i];
		vertex_lat[num_vertexes_points] = origin_lat_values[i];
		convert_cyclic_triangle_point(vertex_lon, num_vertexes_points+1);
		center_x = vertex_lon[num_vertexes_points];
		center_y = vertex_lat[num_vertexes_points];

		num_polar_points = 0;
		if (is_spheric_grid) 
			for (j = 0; j < num_vertexes_points; j ++)
				if (relative_eq(fabs(vertex_lat[j]), 90.0))
					num_polar_points ++;
		
		for (j = 0; j < num_vertexes_points; j ++) {
			vertex_lon_check[j] = vertex_lon[j];
			vertex_lat_check[j] = vertex_lat[j];
			vertex_index_check[j] = j;
		}

		if ((num_polar_points < 2 || num_polar_point_in_origin_coordinate <= 1) && is_spheric_grid && !fast_mode) {
			rotate_sphere_coordinate(vertex_lon_check, vertex_lat_check, num_vertexes_points, center_y, center_x, 1);
			center_x = 0.0;
			center_y = 0.0;
		}
		sort_polygon_vertexes(center_x, center_y, vertex_lon_check, vertex_lat_check, vertex_index_check, num_vertexes_points);

		delete_overlap_point_in_cell_vertexes(vertex_lon_check, vertex_lat_check, vertex_index_check, num_vertexes_points);
		
		bool center_point_in_cell;
		if (num_polar_points < 2 || num_polar_point_in_origin_coordinate <= 1) {
			center_point_in_cell = is_point_in_2D_cell(center_x, center_y, vertex_lon_check, vertex_lat_check, num_vertexes_points, is_spheric_grid, is_spheric_grid, is_spheric_grid);
		}
		else {
			center_point_in_cell = is_point_in_2D_cell(center_x, center_y, vertex_lon_check, vertex_lat_check, num_vertexes_points, is_spheric_grid, is_spheric_grid, false);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, float_eq(fabs(center_y), 90) || center_point_in_cell, "Software error in generate_Voronoi_diagram, the center point is wrong %lf, %lf", center_x, center_y);
		}
		if (!center_point_in_cell) {
			vertex_lon[num_vertexes_of_kernel_points[i]] = origin_lon_values[i];
			vertex_lon_check[num_vertexes_points] = center_x;
			vertex_lat[num_vertexes_of_kernel_points[i]] = origin_lat_values[i];
			vertex_lat_check[num_vertexes_points] = center_y;
			vertex_index_check[num_vertexes_points] = num_vertexes_of_kernel_points[i];
			num_vertexes_points ++;
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_point_in_2D_cell(center_x, center_y, vertex_lon_check, vertex_lat_check, num_vertexes_points, is_spheric_grid, is_spheric_grid, false), "Software error in generate_Voronoi_diagram, the center point is wrong %lf, %lf", center_x, center_y);

		int last_num_vertexes_points = num_vertexes_points;
		delete_point_on_common_line(vertex_lon_check, vertex_lat_check, vertex_index_check, num_vertexes_points);

		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_vertexes_points >= 3, "Software error in PatCC_Delaunay_Voronoi::delete_point_on_common_line");

		for (j = 0; j < num_vertexes_points; j ++) {
			double lon_value = vertex_lon[vertex_index_check[j]];
			if (is_spheric_grid) {
				if (lon_value >= 360.0)
					lon_value -= 360.0;
				if (lon_value < 0)
					lon_value += 360.0;
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lon_value >= 0 && lon_value < 360.0, "Software error in PatCC_Delaunay_Voronoi::generate_Voronoi_diagram: %lf", lon_value);
			}
		    temp_vertex_lon_values[i*max_num_voronoi_diagram_vertex+j] = lon_value;
		    temp_vertex_lat_values[i*max_num_voronoi_diagram_vertex+j] = vertex_lat[vertex_index_check[j]];
		}
		num_vertexes_of_kernel_points[i] = num_vertexes_points;
	}	

	int num_overlapping_point = overlapping_points.size();
	for (j = 0; j < num_overlapping_point; j ++) {
		int point_id2 = overlapping_points[j].second;
		int point_id1 = overlapping_points[j].first;
		if (point_id1 >= 0 && point_id1 < kernel_grid_size) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, point_id2 >= 0 && point_id2 < kernel_grid_size, "Software error in PatCC_Delaunay_Voronoi::generate_Voronoi_diagram");
		}
		else {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, point_id2 < 0 || point_id2 >= kernel_grid_size, "Software error in PatCC_Delaunay_Voronoi::generate_Voronoi_diagram");
			continue;
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, float_eq(origin_lon_values[point_id1], origin_lon_values[point_id2]) && float_eq(origin_lat_values[point_id1], origin_lat_values[point_id2]), "Software ERROR3 in generate_Voronoi_diagram in generating overlapping cell vertexes, %lf, %lf, %lf, %lf", origin_lon_values[point_id], origin_lon_values[overlapping_points[j].first], origin_lat_values[point_id], origin_lat_values[overlapping_points[j].first]);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_vertexes_of_kernel_points[point_id1] > 0 && num_vertexes_of_kernel_points[point_id2] == 0, "Software error in generate_Voronoi_diagram in generating overlapping cell vertexes: %d vs %d", num_vertexes_of_kernel_points[point_id1], num_vertexes_of_kernel_points[point_id2]);
		for (int i = 0; i < num_vertexes_of_kernel_points[point_id1]; i ++) {
			temp_vertex_lon_values[point_id2*max_num_voronoi_diagram_vertex+i] = temp_vertex_lon_values[point_id1*max_num_voronoi_diagram_vertex+i];
			temp_vertex_lat_values[point_id2*max_num_voronoi_diagram_vertex+i] = temp_vertex_lat_values[point_id1*max_num_voronoi_diagram_vertex+i];
		}
		num_vertexes_of_kernel_points[point_id2] = num_vertexes_of_kernel_points[point_id1];
	}

	int new_max_num_voronoi_diagram_vertex = 0;
	for (int i = 0; i < kernel_grid_size; i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_vertexes_of_kernel_points[i] >= 3, "Software error in PatCC_Delaunay_Voronoi::generate_Voronoi_diagram");
		if (num_vertexes_of_kernel_points[i] > new_max_num_voronoi_diagram_vertex)
			new_max_num_voronoi_diagram_vertex = num_vertexes_of_kernel_points[i];
	}

	*vertex_lon_values = new double [kernel_grid_size*new_max_num_voronoi_diagram_vertex];
	*vertex_lat_values = new double [kernel_grid_size*new_max_num_voronoi_diagram_vertex];
	for (int i = 0; i < kernel_grid_size; i ++) {
		for (j = 0; j < num_vertexes_of_kernel_points[i]; j ++) {
			(*vertex_lon_values)[new_max_num_voronoi_diagram_vertex*i+j] = temp_vertex_lon_values[max_num_voronoi_diagram_vertex*i+j];
			(*vertex_lat_values)[new_max_num_voronoi_diagram_vertex*i+j] = temp_vertex_lat_values[max_num_voronoi_diagram_vertex*i+j];
		}	
		for (; j < new_max_num_voronoi_diagram_vertex; j ++) {
			(*vertex_lon_values)[new_max_num_voronoi_diagram_vertex*i+j] = NULL_COORD_VALUE;
			(*vertex_lat_values)[new_max_num_voronoi_diagram_vertex*i+j] = NULL_COORD_VALUE;
		}
	}	
	max_num_voronoi_diagram_vertex = new_max_num_voronoi_diagram_vertex;

	delete [] vertex_lon_check;
	delete [] vertex_lat_check;
	delete [] vertex_lon;
	delete [] vertex_lat;
	delete [] vertex_index_check;
	delete [] vertex_info_buf;
	delete [] temp_vertex_lon_values;
	delete [] temp_vertex_lat_values;
	delete [] num_vertexes_of_kernel_points;
}


void PatCC_Delaunay_Voronoi::build_check_twin_edge_relationship_between_triangles(PatCC_Triangle *t1, PatCC_Triangle *t2, bool do_check)
{
	bool have_twin_edge_relationship = false;


	for (int i = 0; i < 3; i ++)
		for (int j = 0; j < 3; j ++) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !(t1->edge[i]->head == t2->edge[j]->head && t1->edge[i]->tail == t2->edge[j]->tail), "Software error in PatCC_Delaunay_Voronoi::build_check_twin_edge_relationship_between_triangles");
			if (t1->edge[i]->head == t2->edge[j]->tail && t1->edge[i]->tail == t2->edge[j]->head) {
				if (do_check) {
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t1->edge[i]->twin_edge == t2->edge[j] && t2->edge[j]->twin_edge == t1->edge[i], "Software error in PatCC_Delaunay_Voronoi::build_check_twin_edge_relationship_between_triangles: fail to pass check, %d %d %d %d, %d %d %d %d", t1->edge[i]->head, t1->edge[i]->tail, t2->edge[j]->head, t2->edge[j]->tail, t1->edge[i]->twin_edge, t2->edge[j], t2->edge[j]->twin_edge, t1->edge[i]);
				}
				else {
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !have_twin_edge_relationship && t1->edge[i]->twin_edge == NULL && t2->edge[j]->twin_edge == NULL, "Software error in PatCC_Delaunay_Voronoi::build_check_twin_edge_relationship_between_triangles");
					have_twin_edge_relationship = true;
					t1->edge[i]->twin_edge = t2->edge[j];
					t2->edge[j]->twin_edge = t1->edge[i];
				}
			}
		}
}


void PatCC_Delaunay_Voronoi::load_polars_info() {
	vector<int> polars_local_index;
	double nearest_point_lat, shifted_polar_lat, polar_lat;

	if (!is_spheric_grid || original_max_lat != 90.0 && original_min_lat != -90.0)
		return;

	polar_lat = original_max_lat == 90.0? original_max_lat : original_min_lat;
	nearest_point_lat = 0.0;
	for (int i = 0; i < origin_num_points; i ++) {
		if(relative_eq_int(fabs(y_ref[i]), 90.0, 0.0000002))
			polars_local_index.push_back(i);
		else if (fabs(polar_lat - nearest_point_lat) > fabs(polar_lat - y_ref[i]))
			nearest_point_lat = y_ref[i];
	}

	if (polars_local_index.size() > 1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, nearest_point_lat != 0.0, "Software error in PatCC_Delaunay_Voronoi::load_polars_info");
		shifted_polar_lat = nearest_point_lat + 0.5 * (polar_lat - nearest_point_lat);
		for(unsigned i = 0; i < polars_local_index.size(); i++)
			y_ref[polars_local_index[i]] = shifted_polar_lat;
		if (polar_mode) {
			x_ref[num_points] = 0;
			y_ref[num_points] = polar_lat;
			global_index[num_points++] = -3;
		}
	}
}

#define PAT_COORDHASH_FACTOR (2)
bool PatCC_Delaunay_Voronoi::try_fast_triangulate(double min_x, double max_x, double min_y, double max_y)
{
	double area = (max_x - min_x) * (max_y - min_y);
	double density = std::sqrt(num_points / area);

	int lon_expected_points = std::max(1., density * (max_x - min_x));
	int lat_expected_points = std::max(1., density * (max_y - min_y));

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lon_expected_points > 0, "Software error in try_fast_triangulate, lon_expected_points is %d", lon_expected_points);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lat_expected_points > 0, "Software error in try_fast_triangulate, lat_expected_points is %d", lat_expected_points);

	Coord_Hash lon_table(lon_expected_points * PAT_COORDHASH_FACTOR);
	Coord_Hash lat_table(lat_expected_points * PAT_COORDHASH_FACTOR);

	lon_table.set_hashing_params(min_x, max_x);
	lat_table.set_hashing_params(min_y, max_y);

	bool have_vpolar = vpolar_local_index != -1;
	int vpolar_index = vpolar_local_index; 

	for (int i = 0; i < num_points; i++) {
		if (vpolar_index != -1 && vpolar_index == i)
			continue;

		lon_table.put(x_ref[i]);
		lat_table.put(y_ref[i]);
		if (lon_table.get_num_unique_values() * lat_table.get_num_unique_values() > num_points)
			return false;
	}

	int num_lon = lon_table.get_num_unique_values();
	int num_lat = lat_table.get_num_unique_values();

	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "have vpolar %d the num_lon and num_lat is %d, %d, num_points is %d", have_vpolar, num_lon, num_lat, num_points);
	if (!have_vpolar && num_lon * num_lat != num_points)
		return false;
	if (have_vpolar && num_lon * num_lat != num_points-1)
		return false;

	fast_mode = true;

	lon_table.make_sorted_index();
	lat_table.make_sorted_index();
	
	if (!have_vpolar)
		all_points = (PatCC_Point *)::operator new(num_lon * num_lat * sizeof(PatCC_Point));
	else
		all_points = (PatCC_Point *)::operator new((num_lon * num_lat + 1) * sizeof(PatCC_Point));

	for (int i = 0; i < num_points; i++) {
		if (vpolar_index != -1 && vpolar_index == i)
			continue;

		double lon = x_ref[i];
		double lat = y_ref[i];
		int buf_idx = lon_table.get_index(lon) + lat_table.get_index(lat) * num_lon;
		new(&all_points[buf_idx]) PatCC_Point(lon, lat, i, false, -1, -1);
	}

	/* put vpolar at last space */
	if (have_vpolar)
		new(&all_points[num_lon * num_lat]) PatCC_Point(x_ref[vpolar_index], y_ref[vpolar_index], num_lon * num_lat, false, -1, -1);

	fast_triangulate(num_lon, num_lat, have_vpolar);

	return true; 
}


void PatCC_Delaunay_Voronoi::fast_triangulate(int num_lon, int num_lat, bool have_vpolar)
{
	int last_row_begin_pos = -(2*num_lon-2)*(num_lat-1), current_row_begin_pos;

	if (have_vpolar) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_lon * num_lat + 1 == num_points);
	}
	else
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_lon * num_lat == num_points);

	for (int j = 0; j < num_lat - 1; j++) {
		current_row_begin_pos = all_leaf_triangles.size();
		for (int i = 0; i < num_lon - 1; i++) {
			all_leaf_triangles.push_back(allocate_triangle(j*num_lon+i,	j*num_lon+i+1, (j+1)*num_lon+i));
			if (i > 0)
				build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2], all_leaf_triangles[all_leaf_triangles.size()-1], false);
			if (last_row_begin_pos + i*2 + 1 >= 0)
				build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[last_row_begin_pos+i*2+1], all_leaf_triangles[all_leaf_triangles.size()-1], false);
			all_leaf_triangles.push_back(allocate_triangle((j+1)*num_lon+i, j*num_lon+i+1, (j+1)*num_lon+i+1));
			build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2], all_leaf_triangles[all_leaf_triangles.size()-1], false);
		}
		int i = num_lon - 1;
		if (is_spheric_grid && original_min_lon == 0.0 && original_max_lon == 360.0) {
			all_leaf_triangles.push_back(allocate_triangle(j*num_lon+i, j*num_lon, (j+1)*num_lon+i, true));
			build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2], all_leaf_triangles[all_leaf_triangles.size()-1], false);
			if (last_row_begin_pos + i*2 + 1 >= 0)
				build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[last_row_begin_pos+i*2+1], all_leaf_triangles[all_leaf_triangles.size()-1], false);
			all_leaf_triangles.push_back(allocate_triangle((j+1)*num_lon+i, j*num_lon, (j+1)*num_lon,	true));
			build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2], all_leaf_triangles[all_leaf_triangles.size()-1], false);
			build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2*num_lon], all_leaf_triangles[all_leaf_triangles.size()-1], false);
		}
		last_row_begin_pos = current_row_begin_pos;
	}

	if (have_vpolar) {
		int vpole_idx = num_points-1;
		current_row_begin_pos = all_leaf_triangles.size();
		if (relative_eq_int(all_points[vpole_idx].y, -90, 0.0000002)) {
			/* south pole */
			for (int i = 0; i < num_lon; i++) {
				all_leaf_triangles.push_back(allocate_triangle(vpole_idx, (i+1)%num_lon, i, i == num_lon-1));
				if (i > 0)
					build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2], all_leaf_triangles[all_leaf_triangles.size()-1], false);
				if (current_row_begin_pos > 0)
					build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-1], all_leaf_triangles[i*2], false);
			}
		} else if (relative_eq_int(all_points[vpole_idx].y, 90, 0.0000002)) {
			/* north pole */
			for (int i = 0; i < num_lon; i++) {
				all_leaf_triangles.push_back(allocate_triangle(vpole_idx, (num_lat-1)*num_lon+i, (num_lat-1)*num_lon+(i+1)%num_lon));
				if (i > 0)
					build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-2], all_leaf_triangles[all_leaf_triangles.size()-1], false);
				if (current_row_begin_pos >= 0)
					build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-1], all_leaf_triangles[current_row_begin_pos-2*num_lon+2*i], false);
			}
		} else {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in fast_triangulate");
		}
		build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[all_leaf_triangles.size()-1], all_leaf_triangles[current_row_begin_pos], false);
	}

//	if (report_error_enabled) {
//		for (int i = 0; i < all_leaf_triangles.size(); i ++)
//			for (int j = i+1; j < all_leaf_triangles.size(); j ++)
//				build_check_twin_edge_relationship_between_triangles(all_leaf_triangles[i], all_leaf_triangles[j], true);
//	}
}


void PatCC_Delaunay_Voronoi::triangulate()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_points > 0);

	for (unsigned i = 1; i <= dirty_triangles_count; i++)
		if (triangle_stack[i] && triangle_stack[i]->is_leaf)
			triangulating_process(triangle_stack[i], dirty_triangles_count);

	map_buffer_index_to_point_index();

	make_final_triangle();

	mark_special_triangles();

	delete_irregular_triangles();
}

void PatCC_Delaunay_Voronoi::convert_cyclic_triangle_point(double *converted_lon, int num_lon_points) {
	double max_lon = -99999.0;

	if (is_spheric_grid) {
		for(int i = 0; i < num_lon_points; i ++)
			max_lon = max_lon < converted_lon[i] ? converted_lon[i] : max_lon;
		for(int i = 0; i < num_lon_points; i ++)
			converted_lon[i] = max_lon - converted_lon[i] > PAT_CYCLIC_EDGE_THRESHOLD ? converted_lon[i] + 360.0 : converted_lon[i];
	}
}

void PatCC_Delaunay_Voronoi::delete_irregular_triangles_recursively(PatCC_Triangle* t)
{
	PatCC_Point pt[3], circum_center;
	int i, j, num_null_twin_edge, longest_edge_index;
	PatCC_Triangle *twin_t, temp_t;
	double edge_lengths[3], max_edge_length;


	if (!t->is_leaf)
		return;

	for (i = 0; i < 3; i ++)
		if (all_points[t->v[i]].id >= 0 && global_index[all_points[t->v[i]].id] == -2)
			return;

	if (!t->is_virtual && !on_a_line_lon_lat(&all_points[t->v[0]], &all_points[t->v[1]], &all_points[t->v[2]]) && !on_a_line_projected(&all_points[t->v[0]], &all_points[t->v[1]], &all_points[t->v[2]])) {
		for (i = 0, num_null_twin_edge = 0; i < 3; i ++)
			if (t->edge[i]->twin_edge == NULL)
				num_null_twin_edge ++;
		if (num_null_twin_edge == 0)
			return;

// use coordinate value after projection
		for (i = 0; i < 3; i ++) {
			pt[i] = all_points[t->v[i]];
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, pt[i].x == all_points[t->v[i]].x && pt[i].y == all_points[t->v[i]].y);
			pt[i].x = all_points[t->v[i]].x;
			pt[i].y = all_points[t->v[i]].y;
		}
		temp_t.calulate_circum_circle(&pt[0], &pt[1], &pt[2]);
		circum_center.x = temp_t.circum_center[0];
		circum_center.y = temp_t.circum_center[1];
		circum_center.xx = temp_t.circum_center[0];
		circum_center.yy = temp_t.circum_center[1];		
		circum_center.zz = temp_t.circum_center[2];
		if (circum_center.position_to_triangle(&pt[0], &pt[1], &pt[2], true) != -1) // in this triangle
			return;

// find the longest edge and check its twin_edge
		max_edge_length = -1.0;
		for (i = 0; i < 3; i ++) {
			edge_lengths[i] = head(t->edge[i])->calculate_distance_MD(tail(t->edge[i])->x, tail(t->edge[i])->y);
			if (max_edge_length < edge_lengths[i]) {
				max_edge_length = edge_lengths[i];
				longest_edge_index = i;
			}
		}
		if (t->edge[longest_edge_index]->twin_edge != NULL)
			return;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !relative_eq(edge_lengths[longest_edge_index], edge_lengths[(longest_edge_index+1)%3]) && !relative_eq(edge_lengths[longest_edge_index], edge_lengths[(longest_edge_index+2)%3]), "Software error in delete_irregular_triangles_recursively, equal edge length: %lf, %lf, %lf", edge_lengths[longest_edge_index], edge_lengths[(longest_edge_index+1)%3], edge_lengths[(longest_edge_index+2)%3]);
	}

	t->is_leaf = false;
	for (i = 0; i < 3; i ++) {
		if (t->edge[i]->twin_edge != NULL) {
			twin_t = t->edge[i]->twin_edge->triangle;
			t->edge[i]->twin_edge->twin_edge = NULL;
			delete_irregular_triangles_recursively(twin_t);
		}
	}
}

void PatCC_Delaunay_Voronoi::delete_irregular_triangles()
{
	for (int i = 0; i < all_leaf_triangles.size(); i ++)
		delete_irregular_triangles_recursively(all_leaf_triangles[i]);
}


void PatCC_Delaunay_Voronoi::make_final_triangle()
{
	all_leaf_triangles.clear();
	triangle_allocator->get_all_leaf_triangle(all_leaf_triangles);
}


inline PatCC_Point* PatCC_Delaunay_Voronoi::vertex(const PatCC_Triangle* t, int i) 
{ 
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t->v[i] >= 0 && t->v[i] < num_points, "Software error in PatCC_Delaunay_Voronoi::vertex: %d vs %d", t->v[i], num_points);
	return &all_points[t->v[i]]; 
}


inline PatCC_Point* PatCC_Delaunay_Voronoi::head(const PatCC_Edge* e) 
{ 

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, e->head >= 0 && e->head < num_points, "Software error in PatCC_Delaunay_Voronoi::head: %d vs %d", e->head, num_points);
	return &all_points[e->head]; 
}


inline PatCC_Point* PatCC_Delaunay_Voronoi::tail(const PatCC_Edge* e) 
{ 
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, e->tail >= 0 && e->tail < num_points, "Software error in PatCC_Delaunay_Voronoi::tail: %d vs %d", e->tail, num_points);
	return &all_points[e->tail]; 
}


inline void PatCC_Delaunay_Voronoi::diag_triangle(PatCC_Triangle *triangle)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, triangle->v[0] >= 0 && triangle->v[0] < num_points && triangle->v[1] >= 0 && triangle->v[1] < num_points && triangle->v[2] >= 0 && triangle->v[2] < num_points, "Software error in PatCC_Delaunay_Voronoi::diag_triangle: %d %d %d vs %d", triangle->v[0], triangle->v[1], triangle->v[2], num_points);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, head(triangle->edge[0]) != NULL);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, head(triangle->edge[1]) != NULL);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, head(triangle->edge[2]) != NULL);
}


bool PatCC_Delaunay_Voronoi::are_two_points_the_same(int indx1, int indx2)
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, indx1 != indx2, "Software error in PatCC_Delaunay_Voronoi::are_two_points_the_same");
	return relative_eq_int(all_points[indx1].x,all_points[indx2].x,tolerance) && relative_eq_int(all_points[indx1].y,all_points[indx2].y,tolerance);
}


bool PatCC_Point::is_the_same_with_another(const PatCC_Point *v0) const
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, v0 != this, "Software error in PatCC_Delaunay_Voronoi::are_two_points_the_same");
	return float_eq(v0->x,this->x) && float_eq(v0->y,this->y) && float_eq(v0->xx,this->xx) && float_eq(v0->yy,this->yy) && float_eq(v0->zz,this->zz);
}


void PatCC_Point::transform_coord_from_sphere_to_3D()
{
	if (!is_spheric_grid)
		return;

	double cos_lat = cosl(DEGREE_TO_RADIAN(y));

    xx = cos_lat * cosl(DEGREE_TO_RADIAN(x));
    yy = cos_lat * sinl(DEGREE_TO_RADIAN(x));
    zz = sinl(DEGREE_TO_RADIAN(y)); 
}


static inline unsigned long long hash_two_double(double a, double b)
{
	unsigned long long *aa = reinterpret_cast<unsigned long long*>(&a);
	unsigned long long *bb = reinterpret_cast<unsigned long long*>(&b);
	return (*aa > 1) ^ *bb;
}


PatCC_Edge* PatCC_Delaunay_Voronoi::allocate_edge(int head, int tail)
{
	PatCC_Edge *new_edge = edge_allocator->newElement();
	initialize_edge(new_edge, head, tail);
	num_allocate_edges ++;

	return new_edge;
}


PatCC_Triangle* PatCC_Delaunay_Voronoi::allocate_triangle(PatCC_Edge *edge1, PatCC_Edge *edge2, PatCC_Edge *edge3, bool force)
{
	PatCC_Triangle *new_triangle = triangle_allocator->newElement();
	initialize_triangle_with_edges(new_triangle, edge1, edge2, edge3, force);

	diag_triangle(new_triangle);
	num_allocate_triangles ++;

	return new_triangle;
}


PatCC_Triangle* PatCC_Delaunay_Voronoi::allocate_triangle(int idx1, int idx2, int idx3, bool force)
{
	PatCC_Edge* e1 = allocate_edge(idx1, idx2);
	PatCC_Edge* e2 = allocate_edge(idx2, idx3);
	PatCC_Edge* e3 = allocate_edge(idx3, idx1);
	return allocate_triangle(e1, e2, e3, force);
}


inline double calculate_distence_square(double x1, double y1, double x2, double y2)
{
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
}


unsigned long hash_triangle_by_id(Triangle_inline triangle)
{
	return ((unsigned long)triangle.v[0].id) ^ (((unsigned long)triangle.v[1].id) << 21) ^ (((unsigned long)triangle.v[2].id) << 42) +
		   (unsigned long)triangle.v[0].id + (unsigned long)triangle.v[1].id + (unsigned long)triangle.v[2].id + 
		   (unsigned long)triangle.v[0].id * (unsigned long)triangle.v[1].id * (unsigned long)triangle.v[2].id;
}

/* including triangles intersecting with boundary */
void PatCC_Delaunay_Voronoi::get_triangles_in_region(double min_x, double max_x, double min_y, double max_y, 
											   Triangle_inline *output_triangles, int *output_num_triangles, int buf_len)
{
	int current = 0;
			
	for (unsigned i = 0; i < all_leaf_triangles.size(); i++) {

		if (!all_leaf_triangles[i]->is_leaf || all_leaf_triangles[i]->is_virtual)
			continue;

		bool in = true;
		int* v_idx = all_leaf_triangles[i]->v;
		for (int j = 0; j < 3; j++)
			if (!(x_ref[all_points[v_idx[j]].id] <= max_x && x_ref[all_points[v_idx[j]].id] >= min_x &&
				 y_ref[all_points[v_idx[j]].id] <= max_y && y_ref[all_points[v_idx[j]].id] >= min_y)) {
				in = false;
				break;
			}

		if (in) {
			output_triangles[current] = Triangle_inline(PatCC_Point(origin_lon_values[all_points[v_idx[0]].id], origin_lat_values[all_points[v_idx[0]].id], global_index[all_points[v_idx[0]].id], false),PatCC_Point(origin_lon_values[all_points[v_idx[1]].id], origin_lat_values[all_points[v_idx[1]].id], global_index[all_points[v_idx[1]].id], false), PatCC_Point(origin_lon_values[all_points[v_idx[2]].id], origin_lat_values[all_points[v_idx[2]].id], global_index[all_points[v_idx[2]].id], false));
			for (int j = 0; j < 3; j++) {
				if (is_spheric_grid && calculate_distance(origin_lon_values[all_points[v_idx[j]].id], origin_lat_values[all_points[v_idx[j]].id], origin_lon_values[all_points[v_idx[(j+1)%3]].id], origin_lat_values[all_points[v_idx[(j+1)%3]].id]) > PAT_CYCLIC_EDGE_THRESHOLD) {
					output_triangles[current].is_cyclic = true;
					break;
				}
			}
			current ++;
		}
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, current <= buf_len);
	*output_num_triangles = current;
}


void PatCC_Delaunay_Voronoi::set_regional(bool mode)
{
	is_regional = mode;
}


#ifdef OPENCV
void PatCC_Delaunay_Voronoi::plot_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
	unsigned num_edges;
	double *head_coord[2], *tail_coord[2];

	num_edges = 3 * all_leaf_triangles.size();
	head_coord[0] = new double[num_edges];
	head_coord[1] = new double[num_edges];
	tail_coord[0] = new double[num_edges];
	tail_coord[1] = new double[num_edges];

	num_edges = 0;
	if (x_ref) {
		for(unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
			if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual) {
				for(unsigned j = 0; j < 3; j++) {
					head_coord[0][num_edges]   = x_ref[head(all_leaf_triangles[i]->edge[j])->id];
					head_coord[1][num_edges]   = y_ref[head(all_leaf_triangles[i]->edge[j])->id];
					tail_coord[0][num_edges]   = x_ref[tail(all_leaf_triangles[i]->edge[j])->id];
					tail_coord[1][num_edges++] = y_ref[tail(all_leaf_triangles[i]->edge[j])->id];
				}
				if(all_leaf_triangles[i]->is_cyclic)
					for(unsigned j = num_edges-1; j > num_edges-4; j--) {
						if(head_coord[0][j] > 180) head_coord[0][j] -= 360;
						if(tail_coord[0][j] > 180) tail_coord[0][j] -= 360;
					}
			}
		}
	} else {
		for(unsigned i = 0; i < all_leaf_triangles.size(); i ++) {
			if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual) {
				for(unsigned j = 0; j < 3; j++) {
					head_coord[0][num_edges]   = head(all_leaf_triangles[i]->edge[j])->x;
					head_coord[1][num_edges]   = head(all_leaf_triangles[i]->edge[j])->y;
					tail_coord[0][num_edges]   = tail(all_leaf_triangles[i]->edge[j])->x;
					tail_coord[1][num_edges++] = tail(all_leaf_triangles[i]->edge[j])->y;
				}
				if(all_leaf_triangles[i]->is_cyclic)
					for(unsigned j = num_edges-1; j > num_edges-4; j--) {
						if(head_coord[0][j] > 180) head_coord[0][j] -= 360;
						if(tail_coord[0][j] > 180) tail_coord[0][j] -= 360;
					}
			}
		}
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges%3 == 0);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges <= 3 * all_leaf_triangles.size());
	plot_edge_into_file(filename, head_coord, tail_coord, num_edges, min_x, max_x, min_y, max_y);

	delete head_coord[0];
	delete head_coord[1];
	delete tail_coord[0];
	delete tail_coord[1];
}


void PatCC_Delaunay_Voronoi::plot_current_step_into_file(const char *filename)
{
	unsigned num_edges;
	double *head_coord[2], *tail_coord[2];

	num_edges = 3 * all_leaf_triangles.size();
	head_coord[0] = new double[num_edges];
	head_coord[1] = new double[num_edges];
	tail_coord[0] = new double[num_edges];
	tail_coord[1] = new double[num_edges];

	num_edges = 0;
	for(unsigned i = 0; i < all_leaf_triangles.size(); i ++)
		if(all_leaf_triangles[i] && all_leaf_triangles[i]->is_leaf)
			for(unsigned j = 0; j < 3; j++) {
				head_coord[0][num_edges] = head(all_leaf_triangles[i]->edge[j])->x;
				head_coord[1][num_edges] = head(all_leaf_triangles[i]->edge[j])->y;
				tail_coord[0][num_edges] = tail(all_leaf_triangles[i]->edge[j])->x;
				tail_coord[1][num_edges] = tail(all_leaf_triangles[i]->edge[j])->y;
				num_edges++;
			}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges%3 == 0);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges <= 3 * all_leaf_triangles.size());
	plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_DEFAULT, PDLN_PLOT_FILEMODE_NEW);

	delete head_coord[0];
	delete head_coord[1];
	delete tail_coord[0];
	delete tail_coord[1];
}

void PatCC_Delaunay_Voronoi::plot_projection_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
	unsigned num_edges;
	double *head_coord[2], *tail_coord[2];

	num_edges = 3 * all_leaf_triangles.size();
	head_coord[0] = new double[num_edges];
	head_coord[1] = new double[num_edges];
	tail_coord[0] = new double[num_edges];
	tail_coord[1] = new double[num_edges];

	num_edges = 0;
	for(unsigned i = 0; i < all_leaf_triangles.size(); i ++)
		if(all_leaf_triangles[i]->is_leaf && !all_leaf_triangles[i]->is_virtual)
			for(unsigned j = 0; j < 3; j++) {
				head_coord[0][num_edges] = head(all_leaf_triangles[i]->edge[j])->x;
				head_coord[1][num_edges] = head(all_leaf_triangles[i]->edge[j])->y;
				tail_coord[0][num_edges] = tail(all_leaf_triangles[i]->edge[j])->x;
				tail_coord[1][num_edges] = tail(all_leaf_triangles[i]->edge[j])->y;
				num_edges++;
			}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges <= 3 * all_leaf_triangles.size());
	plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_DEFAULT, PDLN_PLOT_FILEMODE_NEW);

	delete head_coord[0];
	delete head_coord[1];
	delete tail_coord[0];
	delete tail_coord[1];
}


void PatCC_Delaunay_Voronoi::plot_original_points_into_file(const char *filename, double min_x, double max_x, double min_y, double max_y)
{
	double *coord[2];

	coord[0] = new double[num_points - PAT_NUM_LOCAL_VPOINTS];
	coord[1] = new double[num_points - PAT_NUM_LOCAL_VPOINTS];

	int num = 0;
	for(int i = PAT_NUM_LOCAL_VPOINTS; i < num_points; i ++) {
		coord[0][num] = all_points[i].x;
		coord[1][num++] = all_points[i].y;
	}

	plot_points_into_file(filename, coord[0], coord[1], NULL, num_points - PAT_NUM_LOCAL_VPOINTS, min_x, max_x, min_y, max_y);

	delete coord[0];
	delete coord[1];
}


void plot_triangles_into_file(const char *prefix, Triangle_inline *t, int num, bool plot_cyclic_triangles)
{
	int num_edges;
	double *head_coord[2], *tail_coord[2];
	char filename[128];

	num_edges = 3 * num;
	head_coord[0] = new double[num_edges];
	head_coord[1] = new double[num_edges];
	tail_coord[0] = new double[num_edges];
	tail_coord[1] = new double[num_edges];

	num_edges = 0;
	for(int i = 0; i < num; i ++) {
		if (plot_cyclic_triangles && (t[i].v[0].x > 360 || t[i].v[1].x > 360 || t[i].v[2].x > 360)) {
			for(int j = 0; j < 3; j++)
				t[i].v[j].x -= 360;
		} else if (plot_cyclic_triangles &&
				   (t[i].v[0].calculate_distance(&t[i].v[1]) >= PAT_CYCLIC_EDGE_THRESHOLD ||
				   t[i].v[1].calculate_distance(&t[i].v[2]) >= PAT_CYCLIC_EDGE_THRESHOLD ||
				   t[i].v[2].calculate_distance(&t[i].v[0]) >= PAT_CYCLIC_EDGE_THRESHOLD )) {
			for(int j = 0; j < 3; j++)
				if(t[i].v[j].x > 180) t[i].v[j].x -= 360;
		}

		for(int j = 0; j < 3; j++) {
			head_coord[0][num_edges] = t[i].v[j].x;
			head_coord[1][num_edges] = t[i].v[j].y;
			tail_coord[0][num_edges] = t[i].v[(j+1)%3].x;
			tail_coord[1][num_edges++] = t[i].v[(j+1)%3].y;
		}
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges%3 == 0);
	snprintf(filename, 128, "%s.png", prefix);
	plot_edge_into_file(filename, head_coord, tail_coord, num_edges);

	delete head_coord[0];
	delete head_coord[1];
	delete tail_coord[0];
	delete tail_coord[1];

}


void plot_triangles_into_file(const char *prefix, std::vector<PatCC_Triangle*> t, PatCC_Point* whole_points)
{
	unsigned num = t.size();
	int num_edges;
	double *head_coord[2], *tail_coord[2];
	char filename[128];

	num_edges = 3 * num;
	head_coord[0] = new double[num_edges];
	head_coord[1] = new double[num_edges];
	tail_coord[0] = new double[num_edges];
	tail_coord[1] = new double[num_edges];

	num_edges = 0;
	for(unsigned i = 0; i < num; i ++)
		for(int j = 0; j < 3; j++) {
			head_coord[0][num_edges] = whole_points[t[i]->v[j]].x;
			head_coord[1][num_edges] = whole_points[t[i]->v[j]].y;
			tail_coord[0][num_edges] = whole_points[t[i]->v[(j+1)%3]].x;
			tail_coord[1][num_edges] = whole_points[t[i]->v[(j+1)%3]].y;
			num_edges++;
		}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_edges%3 == 0, "Software error in plot_triangles_into_file");
	snprintf(filename, 128, "%s.png", prefix);
	plot_projected_edge_into_file(filename, head_coord, tail_coord, num_edges, PDLN_PLOT_COLOR_WHITE, PDLN_PLOT_FILEMODE_NEW);

	delete head_coord[0];
	delete head_coord[1];
	delete tail_coord[0];
	delete tail_coord[1];

}
#endif


Triangle_inline::Triangle_inline(PatCC_Point p0, PatCC_Point p1, PatCC_Point p2, bool cyclic)
{
	v[0] = p0;
	v[1] = p1;
	v[2] = p2;
	is_cyclic = cyclic;
}


void Triangle_inline::check_cyclic()
{
	if (v[0].calculate_distance(&v[1]) > PAT_CYCLIC_EDGE_THRESHOLD ||
		v[1].calculate_distance(&v[2]) > PAT_CYCLIC_EDGE_THRESHOLD ||
		v[2].calculate_distance(&v[0]) > PAT_CYCLIC_EDGE_THRESHOLD )
		is_cyclic = true;
}


bool operator == (Triangle_inline t1, Triangle_inline t2)
{

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t1.v[0] != t1.v[1] && t1.v[1] != t1.v[2] && t1.v[2] != t1.v[0]);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, t2.v[0] != t2.v[1] && t2.v[1] != t2.v[2] && t2.v[2] != t2.v[0]);

	if(t2.v[0] != t1.v[0] && t2.v[0] != t1.v[1] && t2.v[0] != t1.v[2])
		return false;
	if(t2.v[1] != t1.v[0] && t2.v[1] != t1.v[1] && t2.v[1] != t1.v[2])
		return false;
	if(t2.v[2] != t1.v[0] && t2.v[2] != t1.v[1] && t2.v[2] != t1.v[2])
		return false;
	return true;
}


void save_triangles_info_file(const char *prefix, Triangle_inline *t, int num)
{
	char filename[128];
	FILE* fp;

	snprintf(filename, 128, "%s.txt", prefix);
	fp = fopen(filename, "w");
	for(int i = 0; i < num; i++)
		fprintf(fp, "[%d] (%.15lf, %.15lf), (%.15lf, %.15lf), (%.15lf, %.15lf)\n", i, t[i].v[0].x, t[i].v[0].y, t[i].v[1].x, t[i].v[1].y, t[i].v[2].x, t[i].v[2].y);
	fclose(fp);
}
