/***************************************************************
 *  Copyright (c) 2017, Tsinghua University.
 *  This is a source file of C-Coupler.
 *  This file was initially finished by Dr. Chao Sun and then
 *  modified by Dr. Chao Sun and Dr. Li Liu.
 *  If you have any problem,
 *  please contact Dr. Chao Sun via sunchaoo@tsinghua.edu.cn
 *  or Dr. Li Liu via liuli-cess@tsinghua.edu.cn
 ***************************************************************/


#ifndef ENSEMBLE_PROCEDURES_MGT
#define ENSEMBLE_PROCEDURES_MGT


#include <mpi.h>
#include <vector>
#include "memory_mgt.h"
#include "inout_interface_mgt.h"
#include "original_grid_mgt.h"
#include "ensemble_field_operation.h"

#define STATISTICAL_METHOD_DEFAULT    "default"
#define STATISTICAL_METHOD_INST       "inst"
#define STATISTICAL_METHOD_AVER       "aver"
#define STATISTICAL_METHOD_ACCUM      "accum"
#define STATISTICAL_METHOD_MAX        "max"
#define STATISTICAL_METHOD_MIN        "min"

#define ENSEMBLE_OP_DEFAULT           "default"
#define ENSEMBLE_OP_NONE              "none"
#define ENSEMBLE_OP_GATHER            "gather"
#define ENSEMBLE_OP_AVER              "aver"
#define ENSEMBLE_OP_ANOM              "anom"
#define ENSEMBLE_OP_MAX               "max"
#define ENSEMBLE_OP_MIN               "min"
#define ENSEMBLE_OP_MEM               "mem"

#define GET_ENS_PROCEDURES_INST_INDEX(ID)       (((ID & 0x00FF0000) >> 16) & 0x00000FFF)
#define GET_ENS_INST_PROCEDURE_INDEX(ID)        ((ID & 0x00FFF000) >> 12)

class Ensemble_procedures_inst;

struct field_op
{
	char field_name[NAME_STR_SIZE];
	char field_statistical_method[NAME_STR_SIZE];          // inst/aver/accum/max/min
	char field_ensemble_op[NAME_STR_SIZE];                 // none/gather/aver/anom/max/min/mem_%d
	int member_id;                                         // if field_ensemble_op="mem_%d", member_id=%d; else member_id=-1;
	int statistical_line_number;
	int op_line_number;
};


class Field_instances_operation
{
private:
	const char *field_instances_statistical_method;  // default statistical method of all field instances
	const char *field_instances_ensemble_op;         // default ensemble_op of all field instances
	int field_instances_member_id;                                  // default member_id
	int statistical_periodic_timer_id;
	bool do_field_instances_operation;
	bool use_statistical_method;
	bool do_ensemble_OP;
	bool use_set_grids;
	std::vector<field_op> fields_op;
	Ensemble_procedures_inst *ensemble_procedures_inst;
public:
	Field_instances_operation(Ensemble_procedures_inst*, TiXmlElement*);
	~Field_instances_operation() {}
	bool if_do_field_instances_operation() { return do_field_instances_operation; }
	bool if_use_statistical_method() { return use_statistical_method; }
	bool if_do_ensemble_OP() { return do_ensemble_OP; }
	bool if_use_set_grids() { return use_set_grids; }
	field_op *search_field_op(const char*);
	void add_field_op(const char*, const char*, int, int, const char*, int);
	void check_statistical_method_format(const char*, const char*, const char*, int, const char*);
	int check_ensemble_op_format(const char*, const char*, const char*, int, const char*);
	const char *get_field_op_statistical_method(const char *);
	const char *get_field_op_ensemble_op(const char *);
	int get_field_op_member_id(const char *);
};


class Ensemble_procedures_inst
{
private:
	friend class Ensemble_procedures_operation;
	friend class Field_instances_operation;
	int instance_id;
	int external_instance_id;
	int member_comp_id;
	int ensemble_components_id;
	int num_ens_members;
	int proc_member_id;
	const char *instance_name;
	const char *procedures_name;
	const char *procedures_type;
	const char *procedures_dl_name;
	const char *ensemble_components_name;
	int ensemble_components_num;
	std::vector<const char*> ensemble_components_full_name;
	const char *XML_file_name;
	TiXmlDocument *XML_file;
	Field_instances_operation *field_instances_op;
	const char *default_working_directory;
	const char *working_directory;
	const char *pre_instance_script;
	const char *post_instance_script;
	int periodic_timer_id;
	bool use_ensemble_components;
	bool use_periodic_timer;
	MPI_Comm local_comm;
	MPI_Comm host_comm;
	std::vector<int> control_vars;

	std::vector<Field_mem_info*> API_specified_field_insts;
	std::vector<Field_mem_info*> mirror_API_specified_field_insts;
	std::vector<Field_mem_info*> mirror_procedures_field_insts;
	std::vector<Field_mem_info*> mirror_procedures_import_field_insts;
	std::vector<Field_mem_info*> mirror_procedures_export_field_insts;

	std::vector<int> copy_in_index;
	std::vector<int> copy_out_index;

	int ens_member_virtual_grid_id;
	int ensemble_set_decomp_id;
	int ensemble_set_timer_id;
	int **ensemble_member_field_insts_id;
	std::vector<Field_mem_info*> ensemble_set_import_field_insts_before_ens_op;
	std::vector<Field_mem_info*> ensemble_set_export_field_insts_before_ens_op;
	std::vector<Field_mem_info*> ensemble_set_import_field_insts_one_member;
	std::vector<Field_mem_info*> ensemble_set_export_field_insts_one_member;
	std::vector<Field_mem_info*> ensemble_set_import_field_insts_after_ens_op;
	std::vector<Field_mem_info*> ensemble_set_export_field_insts_after_ens_op;
	int *model_import_field_insts_id;
	int *model_export_field_insts_id;
	int **ensemble_member_import_field_insts_id;
	int **ensemble_member_export_field_insts_id;

	int *ensemble_member_export_interface_id;
	int *ensemble_member_import_interface_id;
	int model_export_interface_id;
	int model_import_interface_id;
	int *model_export_field_update_status;
	int *model_import_field_update_status;
	int **ensemble_member_export_field_update_status;
	int **ensemble_member_import_field_update_status;
	std::vector<Member_to_set_operation*> fields_member_to_set_operation;


public:
	Ensemble_procedures_inst(int, const char *, int, int, int, int, int, const int *, const int *, const char *);
	int get_ensemble_components_id() { return ensemble_components_id; }
	int get_host_comp_id() { return ensemble_components_id; }
	int get_member_comp_id() { return member_comp_id; }
	const char *get_instance_name() { return instance_name; }
	const char *get_procedures_name() { return procedures_name; }
	const char *get_procedures_type() { return procedures_type; }
	const char *get_procedures_dl_name() { return procedures_dl_name; }
	const char *get_XML_file_name() { return XML_file_name; }
	const char *get_ensemble_components_name() { return ensemble_components_name; }
	~Ensemble_procedures_inst();
	int get_ens_member_virtual_grid_id() { return ens_member_virtual_grid_id; }
	int get_num_control_vars() { return control_vars.size(); }
	int get_ensemble_size() { return num_ens_members; }
	int get_ensemble_member_index() { return proc_member_id; }
	int get_instance_id() { return instance_id; }
	int get_external_instance_id() { return external_instance_id; }
	MPI_Comm get_local_comm() { return local_comm; }
	MPI_Comm get_host_comm() { return host_comm; }
	bool if_use_ensemble_components() { return use_ensemble_components; }
	int get_control_var(int control_var_index) { return control_vars[control_var_index - 1]; }
	int get_specified_field_instance_id(int field_index) { return API_specified_field_insts[field_index - 1]->get_field_instance_id(); }
	TiXmlElement *get_XML_file_with_configuration();
	Field_instances_operation *get_field_instances_op() { return field_instances_op; }
	int get_num_specified_field_instances() { return API_specified_field_insts.size(); }
	Field_mem_info *get_specified_field_instances(int specified_field_instances_index) { return API_specified_field_insts[specified_field_instances_index - 1]; }

	void do_copy_in();
	void do_copy_out();
	void do_none_ensemble_op_initialize();
	void do_ensemble_OP(std::vector<Field_mem_info*> &, std::vector<Field_mem_info*> &, std::vector<Field_mem_info*> &, int **, std::vector<Field_mem_info*> &);
	void Ensemble_procedures_inst::alloc_ensemble_set_mem(std::vector<Field_mem_info*> &, std::vector<Field_mem_info*> &, std::vector<Field_mem_info*> &, int **, bool );
	void do_ensemble_OP_initialize();
	void run(const char*);
	void execute_config_script(int, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*);
	void change_to_default_work_dir(int);
	void change_to_work_dir(int, const char*);
	//void finalize(const char*);

};

class Ensemble_procedures_mgt
{
private:
	std::vector<Ensemble_procedures_inst *> registered_ensemble_procedures_insts;

public:
	Ensemble_procedures_mgt();
	~Ensemble_procedures_mgt();
	int initialize_ensemble_procedures_inst(const char *, int , int, int, int, int, const int *, const int *, const char *);
	Ensemble_procedures_inst *get_procedures_inst(int, int, const char*);
	int get_ensemble_size(int, const char *);
	int get_ensemble_member_index(int, const char *);
	void add_ensemble_procedures_inst(Ensemble_procedures_inst *instance) { registered_ensemble_procedures_insts.push_back(instance); }
	int get_num_registered_ensemble_procedures_insts() { return registered_ensemble_procedures_insts.size(); }
	Ensemble_procedures_inst *get_registered_ensemble_procedures_inst(int registered_ensemble_procedures_inst_index) { return registered_ensemble_procedures_insts[registered_ensemble_procedures_inst_index - 1]; }
};

#endif