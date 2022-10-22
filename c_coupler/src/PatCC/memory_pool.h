/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef PDLN_MEMORY_POOL_H
#define PDLN_MEMORY_POOL_H

#include <cstddef>
#include <vector>
#include "triangle.h"

/* freed chunk */
struct Bin {
    Bin* next;
};


class Triangle_pool {
    public:
        Triangle_pool();
        ~Triangle_pool();

        PatCC_Triangle* alloc();
        void free(PatCC_Triangle*);
        PatCC_Triangle* newElement();
        void deleteElement(PatCC_Triangle*);

        void get_all_leaf_triangle(std::vector<PatCC_Triangle*>&);
    private:

        PatCC_Triangle* get_page_body(char*);
        PatCC_Triangle* get_page_end(char*);
        void allocNewPage();

        size_t    pagesize;
        char*     cur_page;     // current operating page
        PatCC_Triangle* top_chunk;    // first unallocated chunk
        PatCC_Triangle* end_chunk;
        Bin*      bins;         // freed chunk list
};


class Edge_pool {
    public:
        Edge_pool();
        ~Edge_pool();

        PatCC_Edge* alloc();
        void  free(PatCC_Edge*);
        PatCC_Edge* newElement();
        void  deleteElement(PatCC_Edge*);

    private:

        void allocNewPage();

        size_t pagesize;
        char*  cur_page;     // current operating page
        PatCC_Edge*  top_chunk;    // first unallocated chunk
        PatCC_Edge*  end_chunk;
        Bin*   bins;         // freed chunk list
};
#endif
