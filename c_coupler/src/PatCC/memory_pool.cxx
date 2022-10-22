/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "memory_pool.h"
#include <cstring>

#define PDLN_MEMORY_POOL_SIZE (0x100000)
Triangle_pool::Triangle_pool()
    : pagesize(PDLN_MEMORY_POOL_SIZE)
    , cur_page(NULL)
    , top_chunk(NULL)
    , end_chunk(NULL)
    , bins(NULL)
{
}


Triangle_pool::~Triangle_pool()
{
    char* cur = cur_page;
    while(cur) {
        char* next = (char*)((Bin*)cur)->next;
        operator delete(cur);
        cur = next;
    }
}


size_t padding(char* p, size_t align)
{
    size_t result = (size_t)p;
    return ((align - result) % align);
}


PatCC_Triangle* Triangle_pool::get_page_body(char* page)
{
    char* body = page + sizeof(Bin);
    body += padding(body, sizeof(PatCC_Triangle));
    return (PatCC_Triangle*)body;
}


PatCC_Triangle* Triangle_pool::get_page_end(char* page)
{
    return (PatCC_Triangle*)(page + pagesize - sizeof(PatCC_Triangle));
}


void Triangle_pool::allocNewPage()
{
    char* new_page = (char*)operator new(pagesize);
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli allocate new page: %d", pagesize);
    memset(new_page, 0, pagesize);
    ((Bin*)new_page)->next = (Bin*)cur_page;
    cur_page = new_page;

    char* page_body = new_page + sizeof(Bin);
    page_body += padding(page_body, sizeof(PatCC_Triangle));
    top_chunk = (PatCC_Triangle*)page_body;
    end_chunk = (PatCC_Triangle*)(new_page + pagesize - sizeof(PatCC_Triangle));
}


PatCC_Triangle* Triangle_pool::alloc()
{
    if (bins) {
        PatCC_Triangle* c = (PatCC_Triangle*)bins;
        bins = bins->next;
        return c;
    } else {
        if (top_chunk >= end_chunk)
            allocNewPage();
        return top_chunk++;
    }
}


void Triangle_pool::free(PatCC_Triangle* chunk)
{
    Bin* bin = (Bin*)chunk;
    bin->next = bins;
    bins = bin;
}


PatCC_Triangle* Triangle_pool::newElement()
{
    PatCC_Triangle* result = alloc();
    new (result) PatCC_Triangle();
    return result;
}

void Triangle_pool::deleteElement(PatCC_Triangle* c)
{
    if(c) {
        c->~PatCC_Triangle();
        free(c);
    }
}


void Triangle_pool::get_all_leaf_triangle(std::vector<PatCC_Triangle*>& all/*, std::vector<PatCC_Triangle*>& boundary*/)
{
    char* curpage = cur_page;
    while (curpage) {
        PatCC_Triangle* begin = get_page_body(curpage);
        PatCC_Triangle* end   = get_page_end(curpage);
        for (;begin < end; begin++)
            if (begin->is_leaf) {
                all.push_back(begin);

                //Triangle_inline tp = pack_triangle(PatCC_Triangle);
                //if (is_triangle_intersecting_with_segment(&tp, bound_vertexes[0], bound_vertexes[1], checking_threshold) ||
                //    is_triangle_intersecting_with_segment(&tp, bound_vertexes[1], bound_vertexes[2], checking_threshold) ||
                //    is_triangle_intersecting_with_segment(&tp, bound_vertexes[2], bound_vertexes[3], checking_threshold) ||
                //    is_triangle_intersecting_with_segment(&tp, bound_vertexes[3], bound_vertexes[0], checking_threshold) )
                //    all_leaf_triangles_on_boundary.push_back(PatCC_Triangle);

            }
        curpage = (char*)((Bin*)curpage)->next;
    }
}


Edge_pool::Edge_pool()
    : pagesize(PDLN_MEMORY_POOL_SIZE)
    , cur_page(NULL)
    , top_chunk(NULL)
    , end_chunk(NULL)
    , bins(NULL)
{
}


Edge_pool::~Edge_pool()
{
    char* cur = cur_page;
    while(cur) {
        char* next = (char*)((Bin*)cur)->next;
        operator delete(cur);
        cur = next;
    }
}


void Edge_pool::allocNewPage()
{
    char* new_page = (char*)operator new(pagesize);
    memset(new_page, 0, pagesize);
    ((Bin*)new_page)->next = (Bin*)cur_page;
    cur_page = new_page;

    char* page_body = new_page + sizeof(Bin);
    page_body += padding(page_body, sizeof(PatCC_Edge));
    top_chunk = (PatCC_Edge*)page_body;
    end_chunk = (PatCC_Edge*)(new_page + pagesize - sizeof(PatCC_Edge));
}


PatCC_Edge* Edge_pool::alloc()
{
    if (bins) {
        PatCC_Edge* c = (PatCC_Edge*)bins;
        bins = bins->next;
        return c;
    } else {
        if (top_chunk >= end_chunk)
            allocNewPage();
        return top_chunk++;
    }
}


void Edge_pool::free(PatCC_Edge* chunk)
{
    Bin* bin = (Bin*)chunk;
    bin->next = bins;
    bins = bin;
}

PatCC_Edge* Edge_pool::newElement()
{
    PatCC_Edge* result = alloc();
    new (result) PatCC_Edge();
    return result;
}

void Edge_pool::deleteElement(PatCC_Edge* c)
{
    if(c) {
        c->~PatCC_Edge();
        free(c);
    }
}
