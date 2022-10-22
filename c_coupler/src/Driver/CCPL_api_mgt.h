/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef CCPL_API_MGT
#define CCPL_API_MGT


#include <mpi.h>
#include "tinyxml.h"


enum 
{
    API_ID_FINALIZE,
    API_ID_COMP_MGT_REG_COMP,
    API_ID_COMP_MGT_RECORD_ENSEMBLE_RELATION,
    API_ID_COMP_MGT_GET_COMP_LOG_FILE_NAME,
    API_ID_COMP_MGT_GET_COMP_LOG_FILE_DEVICE,
    API_ID_COMP_MGT_END_COMP_REG,
    API_ID_COMP_MGT_IS_CURRENT_PROC_IN_COMP,
    API_ID_COMP_MGT_GET_CURRENT_PROC_ID_IN_COMP,
    API_ID_COMP_MGT_GET_NUM_PROC_IN_COMP,
    API_ID_COMP_MGT_GET_COMP_PROC_GLOBAL_ID,
    API_ID_COMP_MGT_GET_COMP_ID,
    API_ID_COMP_MGT_GET_PARENT_COMP_ID,
    API_ID_COMP_MGT_GET_COMP_COMM,
    API_ID_COMP_MGT_IS_COMP_TYPE_COUPLED,
    API_ID_GRID_MGT_REG_H2D_GRID_WITHOUT_DATA,
    API_ID_COMP_MGT_GET_LOCAL_COMP_FULL_NAME,
    API_ID_GRID_MGT_REG_H2D_GRID_VIA_LOCAL_DATA,
    API_ID_GRID_MGT_REG_H2D_GRID_VIA_GLOBAL_DATA,
    API_ID_GRID_MGT_REG_H2D_GRID_VIA_FILE,
    API_ID_GRID_MGT_REG_MD_GRID_VIA_MULTI_GRIDS,
    API_ID_GRID_MGT_REG_GRID_VIA_COR,
    API_ID_GRID_MGT_REG_GRID_VIA_LOCAL,
    API_ID_GRID_MGT_REG_H2D_GRID_VIA_COMP,
    API_ID_GRID_MGT_CMP_GRID_VIA_REMOTE,
    API_ID_GRID_MGT_GET_GRID_ID,
    API_ID_GRID_MGT_SET_GRID_DATA,
    API_ID_GRID_MGT_SET_3D_GRID_3D_VERT_FLD,
    API_ID_GRID_MGT_SET_3D_GRID_DYN_BOT_FLD,
    API_ID_GRID_MGT_SET_3D_GRID_STATIC_BOT_FLD,
    API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD,
    API_ID_GRID_MGT_GET_H2D_GRID_DATA,
    API_ID_GRID_MGT_GET_H2D_GRID_AREA_FROM_WGTS,
    API_ID_GRID_MGT_REG_MID_POINT_GRID,
    API_ID_GRID_MGT_GET_GRID_SIZE,
    API_ID_GRID_MGT_REG_NORMAL_1D_GRID_NO_DATA,
    API_ID_GRID_MGT_REG_V1D_GRID_NO_DATA,
    API_ID_GRID_MGT_REG_V1D_Z_GRID_VIA_MODEL,
    API_ID_GRID_MGT_REG_V1D_SIGMA_GRID_VIA_MODEL,
    API_ID_GRID_MGT_REG_V1D_HYBRID_GRID_VIA_MODEL,
    API_ID_DECOMP_MGT_REG_NORMAL_DECOMP,
    API_ID_DECOMP_MGT_REG_CHUNK_DECOMP,
    API_ID_DECOMP_MGT_GET_NUM_CHUNKS,
    API_ID_DECOMP_MGT_GET_CHUNKS_SIZE,
    API_ID_DECOMP_MGT_GET_GRID_SIZE,
    API_ID_DECOMP_MGT_REG_HALO,
    API_ID_FIELD_MGT_REG_FIELD_INST,
    API_ID_FIELD_MGT_GET_FIELD_INST_GRID_ID,
    API_ID_GRID_MGT_GET_SUB_H2D_GRID_ID,
    API_ID_FIELD_MGT_COPY_FIELD_INST_DATA,
    API_ID_FIELD_MGT_START_CHUNK_FIELD_INST_REG,
    API_ID_FIELD_MGT_REG_FIELD_INST_CHUNK,
    API_ID_FIELD_MGT_FINISH_CHUNK_FIELD_INST_REG,
    API_ID_FIELD_MGT_REPORT_FIELD_INST_CKSUM,
    API_ID_FIELD_MGT_REPORT_ALL_FIELD_INST_CKSUM,
    API_ID_FIELD_MGT_REG_IO_FIELD_from_INST,
    API_ID_FIELD_MGT_REG_IO_FIELDs_from_INSTs,
    API_ID_FIELD_MGT_REG_IO_FIELD_from_BUFFER,
    API_ID_TIME_MGT_SET_NORMAL_TIME_STEP,
    API_ID_TIME_MGT_ADVANCE_TIME,
    API_ID_TIME_MGT_RESET_TIME_TO_START,
    API_ID_TIME_MGT_DEFINE_SINGLE_TIMER,
    API_ID_TIME_MGT_DEFINE_COMPLEX_TIMER,
    API_ID_TIME_MGT_GET_CURRENT_NUM_DAYS_IN_YEAR,
    API_ID_TIME_MGT_GET_CURRENT_YEAR,
    API_ID_TIME_MGT_GET_CURRENT_DATE,
    API_ID_TIME_MGT_GET_CURRENT_SECOND,
    API_ID_TIME_MGT_GET_START_TIME,
    API_ID_TIME_MGT_GET_STOP_TIME,
    API_ID_TIME_MGT_GET_PREVIOUS_TIME,
    API_ID_TIME_MGT_GET_CURRENT_TIME,
    API_ID_TIME_MGT_GET_ELAPSED_DAYS_FROM_REF,
    API_ID_TIME_MGT_GET_ELAPSED_DAYS_FROM_START,
    API_ID_TIME_MGT_IS_END_CURRENT_DAY,
    API_ID_TIME_MGT_IS_END_CURRENT_MONTH,
    API_ID_TIME_MGT_GET_CURRENT_CAL_TIME,
    API_ID_TIME_MGT_IS_FIRST_STEP,
    API_ID_TIME_MGT_IS_FIRST_RESTART_STEP,
    API_ID_TIME_MGT_GET_NUM_CURRENT_STEP,
    API_ID_TIME_MGT_GET_NUM_TOTAL_STEPS,
    API_ID_TIME_MGT_GET_NORMAL_TIME_STEP,
    API_ID_TIME_MGT_CHECK_CURRENT_TIME,
    API_ID_TIME_MGT_IS_TIMER_ON,
    API_ID_TIME_MGT_IS_MODEL_LAST_STEP,
    API_ID_TIME_MGT_IS_MODEL_RUN_ENDED,
    API_ID_INTERFACE_REG_IMPORT,
    API_ID_INTERFACE_REG_EXPORT,
    API_ID_INTERFACE_REG_NORMAL_REMAP,
    API_ID_INTERFACE_REG_FRAC_REMAP,
    API_ID_INTERFACE_REG_HALO_EXCHANGE,
    API_ID_INTERFACE_EXECUTE_WITH_ID,
    API_ID_INTERFACE_EXECUTE_WITH_NAME,
    API_ID_INTERFACE_CHECK_DATA_NAN_INF,
    API_ID_INTERFACE_CHECK_IMPORT_FIELD_CONNECTED,
    API_ID_INTERFACE_DO_HALO_WITH_ID,
    API_ID_INTERFACE_DO_HALO_WITH_NAME,
    API_ID_INTERFACE_FIN_HALO_WITH_ID,
    API_ID_INTERFACE_FIN_HALO_WITH_NAME,
    API_ID_INTERFACE_GET_SENDER_TIME,
    API_ID_EXTERNAL_PROC_INST_INIT,
    API_ID_EXTERNAL_PROC_INST_RUN,
    API_ID_EXTERNAL_PROC_INST_GET_COMM,
    API_ID_EXTERNAL_PROC_INST_GET_COMP,
    API_ID_EXTERNAL_PROC_INST_GET_PROC_ACTIVE,
    API_ID_EXTERNAL_PROC_INST_GET_NUM_CONTROLS,
    API_ID_EXTERNAL_PROC_INST_GET_NUM_TIMERS,
    API_ID_EXTERNAL_PROC_INST_GET_NUM_FIELDS,
    API_ID_EXTERNAL_PROC_INST_GET_CONTROL_VAR,
    API_ID_EXTERNAL_PROC_INST_GET_GRID_ID,
    API_ID_EXTERNAL_PROC_INST_GET_DECOMP_ID,
    API_ID_EXTERNAL_PROC_INST_GET_TIMER_ID,
    API_ID_EXTERNAL_PROC_INST_GET_FIELD_ID,
    API_ID_EXTERNAL_PROC_INST_DECLARE_FIELD,
    API_ID_EXTERNAL_PROC_INST_GET_FIELD_PTR,
    API_ID_EXTERNAL_PROC_INST_FINALIZE,
    API_ID_EXTERNAL_SCRIPT_EXE_KEYWORD,
    API_ID_EXTERNAL_SCRIPT_EXE_FILENAME, 
    API_ID_ENSEMBLE_PROC_INST_INIT,
    API_ID_ENSEMBLE_PROC_INST_RUN,
    API_ID_ENSEMBLE_PROC_INST_GET_ENSEMBLE_NUM,
    API_ID_ENSEMBLE_PROC_INST_GET_ENSEMBLE_ID,
    API_ID_ENSEMBLE_PROC_INST_GET_COMM,
    API_ID_ENSEMBLE_PROC_INST_GET_COMP,
    API_ID_ENSEMBLE_PROC_INST_GET_PROC_ACTIVE,
    API_ID_ENSEMBLE_PROC_INST_GET_NUM_GRID_DECOMPS,
    API_ID_ENSEMBLE_PROC_INST_GET_NUM_CONTROLS,
    API_ID_ENSEMBLE_PROC_INST_GET_NUM_TIMERS,
    API_ID_ENSEMBLE_PROC_INST_GET_NUM_FIELDS,
    API_ID_ENSEMBLE_PROC_INST_GET_CONTROL_VAR,
    API_ID_ENSEMBLE_PROC_INST_GET_GRID_ID,
    API_ID_ENSEMBLE_PROC_INST_GET_DECOMP_ID,
    API_ID_ENSEMBLE_PROC_INST_GET_TIMER_ID,
    API_ID_ENSEMBLE_PROC_INST_GET_FIELD_ID,
    API_ID_ENSEMBLE_PROC_INST_DECLARE_FIELD,
    API_ID_ENSEMBLE_PROC_INST_GET_FIELD_PTR,
    API_ID_REPORT_LOG,
    API_ID_REPORT_ERROR,
    API_ID_REPORT_PROGRESS,
    API_ID_RESTART_MGT_WRITE_IO,
    API_ID_RESTART_MGT_START_READ_IO,
    API_ID_RESTART_MGT_IS_TIMER_ON,
    API_ID_RESTART_MGT_READ_ALL,
    API_ID_RESTART_MGT_READ_INTERFACE,
    API_ID_RESTART_MGT_GET_SETTING,
    API_ID_COUPLING_GEN_FAMILY,
    API_ID_COUPLING_GEN_EXTERNAL,
    API_ID_COUPLING_GEN_INDIVIDUAL,
    API_ID_COUPLING_GEN_GET_COMPS,
    API_ID_HANDLER_DATAMODEL_OUTPUT,
    API_ID_HANDLE_NORMAL_EXPLICIT_OUTPUT,
    API_ID_HANDLE_NORMAL_IMPLICIT_OUTPUT,
    API_ID_HANDLER_FIELD_INSTANCES_OUTPUT,
    API_ID_NORMAL_TIMING_START,
    API_ID_NORMAL_TIMING_STOP,
    API_ID_HANDLER_DATA_FILE_INPUT,
    API_ID_HANDLER_DATAMODEL_INPUT,
    API_ID_HANDLE_NORMAL_EXPLICIT_INPUT
};


extern void synchronize_comp_processes_for_API(int, int, MPI_Comm, const char *, const char *);
extern void check_API_parameter_string(int, int, MPI_Comm, const char*, const char*, const char*, const char*);
extern void check_API_parameter_int(int, int, MPI_Comm, const char*, int, const char*, const char*);
extern void check_API_parameter_float(int, int, MPI_Comm, const char*, float, const char*, const char*);
extern void check_API_parameter_double(int, int, MPI_Comm, const char*, double, const char*, const char*);
extern void check_API_parameter_long(int, int, MPI_Comm, const char*, long, const char*, const char*);
extern void check_API_parameter_bool(int, int, MPI_Comm, const char *, bool, const char *, const char *);
extern void check_API_parameter_data_array(int, int, MPI_Comm, const char *, int, int, const char *, const char *, const char *);
extern void check_API_parameter_timer(int, int, MPI_Comm, const char*, int, const char*, const char*);
extern void check_API_parameter_field_instance(int, int, MPI_Comm, const char*, int, const char*, const char*);
extern void check_API_parameter_decomp(int, int, MPI_Comm, const char *, int, const char *, const char *);
extern void check_API_parameter_comp_or_grid(int, int, MPI_Comm, const char *, const char *, const char *);
extern void get_API_hint(int, int, char*);
extern void check_and_verify_name_format_of_string_for_API(int, const char*, int, const char*, const char*);
extern void check_and_verify_name_format_of_string_for_XML(int, const char*, const char*, const char*, int);
extern const char *get_XML_attribute(int, int, TiXmlElement*, const char*, const char*, int&, const char*, const char*, bool);
extern bool is_XML_setting_on(int, TiXmlElement*, const char*, const char*, const char*);
extern void transfer_array_from_one_comp_to_another(int, int, int, int, MPI_Comm, char **, long &);
extern void gather_array_in_one_comp(int, int, void *, int, int, int *, void **, long &, MPI_Comm);
extern void bcast_array_in_one_comp(int, char **, long &, MPI_Comm);
extern char *check_and_aggregate_local_grid_data(int, int, MPI_Comm, const char *, int, int, int, char *, const char *, int, const int *, int &, const char *);
extern bool does_file_exist(const char *);
extern TiXmlDocument *open_XML_file_to_read(int, const char *, MPI_Comm, bool);
extern TiXmlNode *get_XML_first_child_of_unique_root(int, const char *, TiXmlDocument *);
extern void get_required_XML_element(int, const char *, const char *, TiXmlNode **, const char *, bool, bool, bool);

#endif

