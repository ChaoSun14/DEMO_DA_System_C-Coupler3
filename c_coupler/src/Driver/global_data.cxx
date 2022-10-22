/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"


char software_name[16] = "CoR";


Comp_comm_group_mgt_mgr *comp_comm_group_mgt_mgr = NULL;
Original_grid_mgt *original_grid_mgr = NULL;
Routing_info_mgt *routing_info_mgr = NULL;
Timer_mgt *timer_mgr = NULL;
Time_mgt *restart_read_timer_mgr = NULL;
Decomp_info_mgt *decomps_info_mgr = NULL;
Field_info_mgt *fields_info = NULL;
Memory_mgt *memory_manager = NULL;
Remap_mgt *grid_remap_mgr = NULL;
Fields_gather_scatter_mgt *fields_gather_scatter_mgr = NULL;
Decomp_grid_mgt *decomp_grids_mgr = NULL;
Performance_timing_mgt *performance_timing_mgr = NULL;
Ensemble_mgt *ensemble_mgr = NULL;
Annotation_mgt *annotation_mgr = NULL;
Components_time_mgt *components_time_mgrs = NULL;
Inout_interface_mgt *inout_interface_mgr = NULL;
Remapping_configuration_mgt *remapping_configuration_mgr = NULL;
Coupling_generator *coupling_generator = NULL;
Runtime_remapping_weights_mgt *runtime_remapping_weights_mgr = NULL;
H2D_remapping_wgt_file_container *all_H2D_remapping_wgt_files_info = NULL;
External_procedures_mgt *external_procedures_mgr = NULL;
Distributed_H2D_grid_mgt *distributed_H2D_grid_mgr = NULL;
PatCC_Delaunay_Voronoi *current_triangulation = NULL;
Distributed_H2D_weights_generator *current_distributed_H2D_weights_generator = NULL;
Ensemble_procedures_mgt *ensemble_procedures_mgr = NULL;
Datamodel_mgt *datamodel_mgr = NULL;
int *current_remap_local_cell_global_indexes = NULL;
int num_current_remap_local_cell_global_indexes;
Union_comm_mgt *union_comm_mgr = NULL;



