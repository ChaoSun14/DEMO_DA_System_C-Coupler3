/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef GLOBAL_DATA
#define GLOBAL_DATA

#include "compset_communicators_info_mgt.h"
#include "decomp_info_mgt.h"
#include "field_info_mgt.h"
#include "memory_mgt.h"
#include "timer_mgt.h"
#include "routing_info_mgt.h"
#include "cor_cpl_interface.h"
#include "cor_global_data.h"
#include "fields_gather_scatter_mgt.h"
#include "decomp_grid_mgt.h"
#include "common_utils.h"
#include "execution_report.h"
#include "performance_timing_mgt.h"
#include "ensemble_mgt.h"
#include "object_type_prefix.h"
#include "CCPL_api_mgt.h"
#include "original_grid_mgt.h"
#include "annotation_mgt.h"
#include "inout_interface_mgt.h"
#include "coupling_generator.h"
#include "runtime_trans_algorithm.h"
#include "datamodel_mgt.h"
#include "remapping_configuration_mgt.h"
#include "runtime_remapping_weights_mgt.h"
#include "external_procedures_mgt.h"
#include "distributed_H2D_grid.h"
#include "triangulation.h"
#include "distributed_H2D_wgts_gen.h"
#include "ensemble_procedures_mgt.h"


extern char software_name[];


extern Comp_comm_group_mgt_mgr *comp_comm_group_mgt_mgr;
extern Original_grid_mgt *original_grid_mgr;
extern Routing_info_mgt *routing_info_mgr;
extern Timer_mgt *timer_mgr;
extern Time_mgt *restart_read_timer_mgr;
extern Decomp_info_mgt *decomps_info_mgr;
extern Field_info_mgt *fields_info;
extern Memory_mgt *memory_manager;
extern Remap_mgt *grid_remap_mgr;
extern Fields_gather_scatter_mgt *fields_gather_scatter_mgr;
extern Decomp_grid_mgt *decomp_grids_mgr;
extern Performance_timing_mgt *performance_timing_mgr;
extern Ensemble_mgt *ensemble_mgr;
extern Annotation_mgt *annotation_mgr;
extern Components_time_mgt *components_time_mgrs;
extern Inout_interface_mgt *inout_interface_mgr;
extern Remapping_configuration_mgt *remapping_configuration_mgr;
extern Coupling_generator *coupling_generator;
extern Runtime_remapping_weights_mgt *runtime_remapping_weights_mgr;
extern H2D_remapping_wgt_file_container *all_H2D_remapping_wgt_files_info;
extern External_procedures_mgt *external_procedures_mgr;
extern Distributed_H2D_grid_mgt *distributed_H2D_grid_mgr;
extern PatCC_Delaunay_Voronoi *current_triangulation;
extern Distributed_H2D_weights_generator *current_distributed_H2D_weights_generator;
extern Ensemble_procedures_mgt *ensemble_procedures_mgr;
extern Datamodel_mgt *datamodel_mgr;
extern int *current_remap_local_cell_global_indexes;
extern int num_current_remap_local_cell_global_indexes;
extern Union_comm_mgt *union_comm_mgr;


#endif
