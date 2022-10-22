/***************************************************************
  *  Copyright (c) 2019, Tsinghua University.
  *  This is a source file of PatCC.
  *  This file was initially finished by Dr. Li Liu and
  *  Haoyu Yang. If you have any problem,
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef PD_ASSERT_H
#define PD_ASSERT_H

#include <mpi.h>
#include "common_utils.h"


#define PDLN_ABS_TOLERANCE_LOW ((double) 1e-8)
#define PDLN_ABS_TOLERANCE     ((double) 1e-10) // if less than 1e-10, will three point in a line, if more than 1e-15, will not pass check
#define PDLN_ABS_TOLERANCE_HI  ((double) 1e-11) // normal grid less than 1e-11

#define PDLN_RELATIVE_TOLERANCE_LOW ((double) 1e-3)
#define PDLN_RELATIVE_TOLERANCE     ((double) 1e-4)
#define PDLN_RELATIVE_TOLERANCE_HI  ((double) 1e-5)


typedef long double PAT_REAL;


#define PI ((double) 3.1415926535897932384626433)


class PDLN_Timer
{
    public:
        double time;

        PDLN_Timer() : time(0){};
        double tick() {
            double t = MPI_Wtime();
            double result = t - time;
            time = t;
            return result;
        }
};

#endif
