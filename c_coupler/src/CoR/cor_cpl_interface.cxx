/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "cor_cpl_interface.h"
#include "cor_global_data.h"
#include "execution_report.h"


long cpl_get_grid_size(const char *grid_name)
{
    Remap_grid_class *model_grid;


    model_grid = remap_grid_manager->search_remap_grid_with_grid_name(grid_name);    
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid != NULL, "grid %s is undefined\n", grid_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, model_grid->get_grid_size() > 0, "the grid size %%s is not known\n", grid_name);

    model_grid->end_grid_definition_stage(NULL);
    return model_grid->get_grid_size();
}


void cpl_check_remap_weights_format(Remap_weight_of_strategy_class *remap_weights)
{
    remap_weights->check_remap_weights_format();
}

