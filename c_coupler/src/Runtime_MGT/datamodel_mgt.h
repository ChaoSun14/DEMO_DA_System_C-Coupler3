/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Miss Xinzhu Yu and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef DATAMODEL_MGT_H
#define DATAMODEL_MGT_H

#include "common_utils.h"
#include <netcdf.h>
#include "coupling_generator.h"
#include "cor_global_data.h"
#include "timer_mgt.h"
#include "memory_mgt.h"
#include "tinyxml.h"
#include <string>
#include <vector>
#include "compset_communicators_info_mgt.h"
#include "io_pnetcdf.h"
#include "object_type_prefix.h"
#include "runtime_cumulate_average_algorithm.h"


#define INPUT_DATAMODEL            ((int)0)
#define OUTPUT_DATAMODEL           ((int)1)
#define RESTART_HANDLER            ((int)0)
#define OUTPUT_HANDLER             ((int)1)
#define INPUT_HANDLER              ((int)2)


struct Time_field_param
{
	const char *varname;
	const char *time_format_in_data_file;
	int id_time_format_in_data_file;
	Time_field_param() {
		varname = NULL;
		time_format_in_data_file = NULL;
	}
	~Time_field_param() {
		delete varname;
		delete time_format_in_data_file;
	}
};


struct Time_field_value
{
	std::vector<long> time_fields;//all time_field times
	std::vector<int> last_time_global_ind;//the size of this vector is the same with the num of files
};

struct Input_time_field_setting
{
	const char *specification;//file_name, file_field
	int time_point_type;
	char time_field_full_time_format[32];
	std::vector<Time_field_param*> time_fields;
	int id_full_time_format;
	Input_time_field_setting() {
		specification = NULL;
		time_field_full_time_format[0] = '\0';
	}
	~Input_time_field_setting() {
		if (specification != NULL) delete specification;
		for (int i = 0; i < time_fields.size(); i++) delete time_fields[i];
	}
};


class Input_file_time_info
{
public:
	char full_file_name[NAME_STR_SIZE];
#ifdef USE_PARALLEL_IO
	IO_pnetcdf *netcdf_file_object;
#else
	IO_netcdf *netcdf_file_object;
#endif
	int time_ind_in_dir;
	int time_position;//local time position for time_field
	int global_ind;//for time_field
	bool read_flag;
	bool use_flag;

	Input_file_time_info() {
		full_file_name[0] = '\0';
		netcdf_file_object = NULL;
		time_position = -1;
		time_ind_in_dir = -1;
		global_ind = -1;
		read_flag = false;
	}
	~Input_file_time_info() {
		if (netcdf_file_object != NULL)
			delete netcdf_file_object;
	}
	void copy_info(Input_file_time_info*);
	void set_value(void*, int, int, int, bool);
};


struct Datamodel_file_info
{
	char *datamodel_files_dir_name;
	char *file_dir;
	char *file_name_prefix;
	char *file_name_suffix;
	char *time_format_in_file_name;
	char *file_type;
	int id_time_format_in_file_name;
	
	bool move_flag;
	bool perfect_match;
	Input_time_field_setting *time_field_setting;
	std::vector<long> file_name_times;//time in file_name
	Time_field_value *time_fields_info;
	std::vector<Input_file_time_info*> next_input_file_times;//at most two

	~Datamodel_file_info() {
		delete datamodel_files_dir_name;
		delete file_dir;
		delete file_name_prefix;
		delete file_name_suffix;
		delete time_format_in_file_name;
		delete file_type;
		delete time_field_setting;
		for (int i = 0; i < next_input_file_times.size(); i++)
			delete next_input_file_times[i];
	}
};

struct Field_config_info
{
	bool datamodel_type;//input: 0, output: 1
	char *name_in_model;
	char *name_in_file;
	char *grid_name;
	char *float_datatype;
	char *integer_datatype;
	char *operation;
	char *unit;
	double add_offset;
	double scale_factor;
	char *default_output_grid_name;

	double min_bound_value;
	double max_bound_value;
	bool have_min_bound;
	bool have_max_bound;

	Datamodel_file_info *file_set_for_use;

	~Field_config_info() {
		delete name_in_model;
		delete name_in_file;
		delete grid_name;
		delete float_datatype;
		delete integer_datatype;
		delete operation;
	}
};


struct Fields_output_common_setting
{
	char *default_operation;
	char *default_float_type;
	char *default_integer_type;
	int default_output_grid_id;
	char *default_output_grid_name;
	int id_time_format_in_data_file;
	char *file_mid_name;
	char *field_specification;
	void init() {
		this->default_operation = strdup("inst");
		this->default_float_type = strdup(DATA_TYPE_FLOAT);
		this->default_integer_type = strdup(DATA_TYPE_INT);
		this->default_output_grid_id = -1;
		this->default_output_grid_name = strdup("");
		this->file_mid_name = strdup("");
		this->field_specification = strdup("fields");
	}
	~Fields_output_common_setting() {
		delete default_operation;
		delete default_integer_type;
		delete default_output_grid_name;
		delete file_mid_name;
		delete field_specification;
	}
};


struct V3d_grid_info {
	char *grid_name;
	char *h2d_subgrid_name;
	char *v1d_subgrid_name;
	char *surface_field_name;
	char *dimension_order;
	char *mid_point_grid_name;
	int line_number;
	~V3d_grid_info() {
		delete grid_name;
		delete h2d_subgrid_name;
		delete v1d_subgrid_name;
		delete surface_field_name;
		delete dimension_order;
		delete mid_point_grid_name;
	}
};


struct Field_grid_info
{
	int grid_id;
	bool is_sigma;
	char *grid_name;
	char *surface_field_name;
	~Field_grid_info() {
		delete grid_name;
		delete surface_field_name;
	}
};


struct Import_field_info
{
	int io_field_instance_id;
	int model_field_instance_id;
	char *import_field_operation;
	char *field_name_in_file;
	Datamodel_file_info *file_set_for_use;
	std::vector<Field_mem_info*> averaging_field_mems;//size=two on model decomposition, returned by scatter
	Import_field_info() {
		file_set_for_use = NULL;
	}

	~Import_field_info() {
		if (import_field_operation != NULL)
			delete import_field_operation;
		if (field_name_in_file != NULL)
			delete field_name_in_file;
	}
};

struct Period_setting
{
	long period_data_start_time;//long?
	const char *period_time_format;
	int id_period_time_format;
	const char *period_unit;
	int id_period_unit;
	int period_count;
	~Period_setting() {
		delete period_time_format;
		delete period_unit;
	}
};

struct Offset_setting
{
	const char *offset_unit;
	int offset_count;
	int id_offset_unit;
	~Offset_setting() {
		delete offset_unit;
	}
};

struct Field_rename_map
{
	const char *name_in_model;
	const char *name_in_file;
	~Field_rename_map() {
		delete name_in_model;
		delete name_in_file;
	}
};

struct V3d_grid_surface_map
{
	int v3d_grid_id;
	char *surface_field_name_in_file;
	Datamodel_file_info *file_set_for_use;
};

class Inout_datamodel
{
private:
	friend class Output_handler;
	friend class Input_handler_controller;
	bool datamodel_type;//input: 0, output: 1
	char datamodel_config_dir[NAME_STR_SIZE];//config dir
	char datamodel_output_data_dir[NAME_STR_SIZE];//output data dir
	char datamodel_input_data_dir[NAME_STR_SIZE];//input data dir
	char XML_file_name[NAME_STR_SIZE];//dir+filename
	char *annotation;
	char *datamodel_name;
	int host_comp_id;

	std::vector<int> io_timer_ids;
	std::vector<int> file_timer_ids;

	std::vector<Fields_output_common_setting*> default_settings;
	std::vector<Field_grid_info*> fields_grid_info;
	std::vector<V3d_grid_info*> special_v3ds;
	std::vector<std::vector<Field_config_info*> > fields_config_info;
	std::vector<V3d_grid_surface_map*> v3d_grid_surfaces;

public:
	std::vector<Datamodel_file_info*> file_sets_info;
	std::vector<Field_mem_info*> input_surface_mems;

	Inout_datamodel(int, const char*, bool, const char*);
	Inout_datamodel(int, const char*, const char*, const char*, int, int, int, int, const char*, const char*, bool, const char*);
	Inout_datamodel(int, const char*, const char*);
	~Inout_datamodel();
	const char *get_datamodel_name() { return datamodel_name; }
	int get_host_comp_id() { return host_comp_id; }
	char *get_datamodel_output_data_dir() { return datamodel_output_data_dir; }
	char *get_datamodel_input_data_dir() { return datamodel_input_data_dir; }
	char *get_XML_file_name() { return XML_file_name; }
	int get_datamodel_grid_id(const char*, int&, bool, bool);
	void generate_datamodel_files_info(const char*, const char*, int&, Datamodel_file_info*);
	void config_data_files_for_datamodel(TiXmlNode*);
	void config_horizontal_grids_for_datamodel(TiXmlNode*);
	void config_vertical_grids_for_datamodel(TiXmlNode*);
	void config_v3d_grids_for_datamodel(TiXmlNode*);
	int config_horizontal_grid_via_CCPL_grid_file(TiXmlNode*, const char*, const char*);
	int config_horizontal_grid_via_grid_data_file_field(TiXmlNode*, const char*, const char*);
	int config_horizontal_grid_via_uniform_lonlat_grid(TiXmlNode*, const char*, const char*);
	void config_vertical_z_grid(MPI_Comm, TiXmlNode*, bool, IO_netcdf*, char*, void **coord_values, int &, char*, bool);
	void config_vertical_sigma_grid(MPI_Comm, TiXmlNode*, bool, IO_netcdf*, char*, char*, void**, int&, char*, bool);
	void config_vertical_hybrid_grid(MPI_Comm, TiXmlNode*, bool, const char*, IO_netcdf*, char*, char*, void**, void**, int&, char*, bool);
	void config_field_output_settings_for_datamodel(TiXmlNode*);
	void config_field_info(TiXmlNode*, int);
	void config_output_time_setting(TiXmlNode*);
	void config_fields_output_common_setting(TiXmlNode*);
	int register_common_H2D_grid_for_datamodel(const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char*, const char *, const char*, const char*, int, const char*);
	void register_datamodel_mid_point_v3d_grid(const char*, int);
	int register_datamodel_v3d_grid(const char*, char*, int, int, int);
	void push_special_v3d_grid_info(const char*, const char*, const char*, char*, char*, char*, int);
	bool is_special_v3d_grid(const char*, V3d_grid_info*&);
	Field_grid_info *search_datamodel_field_grid_info(int);
	Coupling_timer *config_period_slots_points_node(TiXmlNode*, const char*, int);
	void get_time_slots_points_from_str(std::vector<std::pair<long, long> >&, const char*, TiXmlElement*, bool);
	void config_data_file_node_info(TiXmlNode*, Datamodel_file_info*);
	void config_time_fields_node_info(TiXmlNode*, Datamodel_file_info*);
	void config_datamodel_configuration_file_info(const char*, bool);
	void find_target_root_node_of_datamodel(int, const char*, TiXmlNode*&, const char*);
	Datamodel_file_info *config_input_datamodel_data_time_series(TiXmlNode*);
	void config_input_fields_for_datamodel(TiXmlNode*);
	const char *randomly_match_a_data_file_in_datamodel(Datamodel_file_info*);
	void initialize_input_datamodel_file_sets_info();
	void read_time_fields_from_a_data_file(Datamodel_file_info*, char*, int);
	void allocate_surface_fields_for_input_datamodel();
};


class Input_instance
{
private:
	int host_comp_id;
	char *XML_file_name;
	const char *instance_name;
	const char *datamodel_name;
	std::vector<Field_rename_map*> field_renamings;
	const char *annotation;
public:
	Inout_datamodel *input_datamodel;
	Period_setting *period_info;
	Offset_setting *offset_info;

	Input_instance(int, const char*, const char*);
	~Input_instance();
	const char *get_datamodel_name() { return datamodel_name; }
	const char *get_input_instance_name() { return instance_name; }
	void config_input_instance_from_xml_file(char*, MPI_Comm, const char*);
	void config_time_mapping_info_from_xml(TiXmlNode*);
	void config_offset_setting(TiXmlElement*);
	void config_period_setting(TiXmlElement*);
	void config_input_instance_field_renamings(TiXmlNode*);
	Inout_datamodel *get_input_datamodel() { return input_datamodel; }
	const char *get_field_name_in_file(const char*);
};


struct Output_file_procedure
{
#ifdef USE_PARALLEL_IO
	IO_pnetcdf *netcdf_file_object;
#else
	IO_netcdf *netcdf_file_object;
#endif
	char *file_mid_name;
	Inout_interface *interface;
	int interface_index;
	int file_timer_id;
	int interface_timer_id;
	int inst_or_aver;
	std::vector<char*> fields_name;
	std::vector<int> dst_field_instance_ids;
	std::vector<int> src_field_instance_ids;
	std::vector<int> dst_grid_id;
	std::vector<Field_mem_info*> averaging_field_mems;
	std::vector<Runtime_cumulate_average_algorithm*> outer_level_averaging_algorithm;

	~Output_file_procedure() {
		if (netcdf_file_object != NULL)
			delete netcdf_file_object;
		delete [] file_mid_name;
		for (int i = 0; i < fields_name.size(); i++)
			delete fields_name[i];
		for (int i = 0; i < averaging_field_mems.size(); i++) {
			if (averaging_field_mems[i] != NULL) {
				delete averaging_field_mems[i];
				delete outer_level_averaging_algorithm[i];
			}
		}
	}
};

class Inout_handler
{
public:
	Inout_datamodel *inout_datamodel;//common
	char *handler_name;
	int host_comp_id;//common
	int handler_id;
	int implicit_or_explicit;//1 for explicit, 0 for implicit
	int num_fields;
	char *datamodel_name;
	int pio_proc_num;//common
	int io_proc_mark;//common
	MPI_Comm io_comm;//common
	int total_num_procs;//common
	int io_proc_stride;//common
	char *annotation;
	int handler_type;//common
	int level;//1 or 2: 1 for outer level, 2 for inner level

	Time_mgt *time_mgr;
	std::vector<std::vector<Field_config_info*> > fields_config_info;

	Inout_handler() {}
	~Inout_handler() {}
	Field_mem_info *find_handler_field_instance(char*, int*, int, bool);
	bool find_field_info_from_datamodel_config(const char*, std::vector<Field_config_info*>, int&);
};

class Input_handler_controller: public Inout_handler
{
private:
	int input_timer_id;
	int calculate_unit;
	long period_time_value_add_count;
	long file_time_value;

	Input_instance *input_instance;
	std::vector<Import_field_info*> input_fields_info;
	std::vector<Input_file_time_info*> next_input_file_times;
public:
	Input_handler_controller(const char*, const char*, int, int*, int, int, int, const char*, int*, int*, const char*, const char*, int, char*, char*, bool, bool, const char*, const char*);
	~Input_handler_controller();
	void common_checking_for_read_in_field_handler(const char*, int, int*, int, int, const char*, const char*, int*, int, const char*);
	void config_input_field_info(char*, char*, bool, bool, const char*);
	void config_next_input_file_time_info(const char*);
	void execute_input_handler_controller(const char*);//execute
	const char *get_handler_name() { return handler_name; }
	const char *get_handler_annotation() { return this->annotation; }
	int get_input_handler_controller_id() { return this->handler_id; }
	int get_host_comp_id() { return this->host_comp_id; }
	void initialize_input_file_time_info_structure(const char*);
	void config_read_in_single_file_time(const char*);
	void config_readin_multiple_file_time_infos();
	void find_readin_files_for_target_time(int, int, int, int);
	void config_input_fields_model_io_info(void*, int, const char*, const char*, int, const char*);

	void get_full_file_name_from_time_value_on_target_format(Datamodel_file_info*, int, char*);
	void config_file_name_time_left_right_match(int, int, Datamodel_file_info*);
	void config_file_name_time_perfect_match(int, Datamodel_file_info*);
	void config_time_field_perfect_match(int, Datamodel_file_info*);
	void config_time_field_left_right_match(int, int, Datamodel_file_info*);
	void interpolate_read_in_value_on_target_time(Field_mem_info*, Import_field_info*);
	template<typename T> void linear_time_interpolation(Field_mem_info*, std::vector<float>, std::vector<Field_mem_info*>);
};

class Input_handler_operator: public Inout_handler
{
private:
	int input_model_timer_id;
	int input_io_timer_id;
	Inout_interface *handler_import_interface;
	Inout_interface *handler_export_interface;

	Import_field_info* export_interface_field_info;
public:
	Input_handler_operator(const char*, int, int, int, const char*, const char*, const char*);
	~Input_handler_operator();
	void config_input_handler_operator_export_interface(int, int, const char*, const char*);
	Field_mem_info *get_unique_IO_field_mem();
	int determine_field_io_grid(Field_mem_info*, int);
	void execute_input_handler_operator(const char*, Input_file_time_info*, const char*);
	void execute_read_in_field(const char*, Input_file_time_info*, const char*);
	void execute_input_handler_operator_interface(bool, const char*);
};


class Output_handler: public Inout_handler
{
private:
	int sampling_timer_id;
	int output_grid_id;
	int handler_output_H2D_grid_id;
	int handler_output_V1D_grid_id;
	int export_timer_id;

	int *expanded_export_field_ids;

	Inout_interface *handler_export_interface;
	std::vector<Import_field_info*> import_fields_info;
	std::vector<Output_file_procedure*> output_procedures;
	std::map<const char*, Field_mem_info*> surface_field_src_dst_map;
public:
	Output_handler(const char*, const char*, Inout_datamodel*, int, int, int*, int, int, int, int, int, int, int, const char *, bool, bool, int);
	~Output_handler();
	void initiate_handler_info(const char*, const char*, int, int, int, int, int, int, int, int*, int, const char*);
	int get_handler_type() { return handler_type; }
	bool check_handler_field_instance_match(Field_mem_info*, Field_mem_info*&);
	void config_datamodel_import_interface_parameters(int, int*, bool);
	int get_handler_id() { return handler_id; }
	int get_host_comp_id() { return host_comp_id; }
	bool get_handler_implicit_or_explicit() { return implicit_or_explicit; }
	const char *get_handler_annotation() { return this->annotation; }
	const char *get_handler_name() { return handler_name; }
	Inout_interface *get_handler_interface(int);
	int get_output_procedures_size() { return output_procedures.size(); }
	Output_file_procedure *get_handler_output_procedure(int procedure_index) { return output_procedures[procedure_index]; }
	int *add_surface_fields_to_export_interface(int&, const int*, int&);
	void determine_field_output_grid(Field_mem_info*, const char*, int&);
	void generate_a_coupling_procedure(int, int);
	void register_operation_import_interface(int, int*, int*, int, int, int, const char*, const char*, int);
	void check_grid_dim(int &, Field_mem_info*);
	void get_or_generate_field_real_output_grid(const char*, Field_mem_info*, int&);
	char *modify_special_v3d_grid_name(char*, const char*);
	void execute_handler(bool, int, const char*, void*, const char*);

	MPI_Comm get_io_comm() { return io_comm; }
	int get_io_proc_mark() { return io_proc_mark; }
	int get_pio_proc_num() { return pio_proc_num; }
	bool search_import_field_interface_surface_field(const char*);
	bool search_export_field_id(int*, int, int);
	void execute_handler_coupling_procedure(bool, const char*);
	void execute_handler_average_procedure(bool);
	Field_mem_info *get_unique_IO_field_mem();
};


class Datamodel_mgt
{
private:
	friend class Output_handler;
	friend class Input_handler_controller;
	std::vector<Output_handler*> output_handlers;
	std::vector<Inout_datamodel*> output_datamodels;
	std::vector<Inout_datamodel*> input_datamodels;
	std::vector<Input_handler_controller*> input_handler_controllers;
	std::vector<Input_handler_operator*> input_handler_operators;
	std::vector<int> implicit_output_handler_ids;
	std::vector<Input_instance*> input_instances;
	int handlers_field_mem_buf_mark;
public:
	Datamodel_mgt() { handlers_field_mem_buf_mark = 0; }
	~Datamodel_mgt();
	int get_output_handlers_size() { return output_handlers.size(); }
	void add_handlers_field_mem_buf_mark() { handlers_field_mem_buf_mark ++; }
	int get_handlers_field_mem_buf_mark() { return handlers_field_mem_buf_mark; }
	int register_datamodel_output_handler(const char*, int, int *, const char*, int, int, int, int, int, int, int, const char*);
	int register_field_instances_output_handler(int, int*, const char*, const char*, int, int, int, int, int, int, const char*, const char*, int, bool, int, const char*);
	int get_next_output_handler_id() { return TYPE_OUTPUT_HANDLER_ID_PREFIX|output_handlers.size(); }
	int get_next_input_handler_controller_id() { return TYPE_INPUT_HANDLER_ID_PREFIX|input_handler_controllers.size(); }
	int get_next_input_handler_operator_id() { return TYPE_INPUT_HANDLER_ID_PREFIX|input_handler_operators.size(); }
	void common_checking_for_datamodel_handler_registration(int, int*, int, int, int, int, int, int, const char*, const char*, bool);
	void handle_normal_output(int, bool, int, const char*);
	void handle_explicit_output(int, bool, int, const char*);
	void handle_implicit_output(const char*);
	bool is_legal_output_handler_id(int);
	void get_all_implicit_handler_ids();
	Output_handler *get_output_handler(int);
	Inout_datamodel *search_output_datamodel(int, char*);
	Output_handler *search_output_handler(int handler_id) { return output_handlers[handler_id&TYPE_ID_SUFFIX_MASK]; }
	int get_comp_PIO_proc_setting(int, int&, int&, MPI_Comm&);
	Output_handler *get_output_handler_by_index(int index) { return output_handlers[index]; }
	void read_in_field_from_datafile(int, const char*, const char*, const char*, int, int, char*, char*, bool, bool, const char*, int&, const char*);
	int get_num_input_handler_operators() { return input_handler_operators.size(); }
	Input_handler_operator *search_input_handler_operator(int handler_id) { return input_handler_operators[handler_id&TYPE_ID_SUFFIX_MASK]; }
	int register_input_handler_operator(int, int, const char*, const char*);
	int register_datamodel_input_handler(const char*, int*, int, int, const char*, int, int*, int*, const char*);
	Input_instance *search_input_instance(const char*);
	Inout_datamodel *search_input_datamodel(int, const char*);
	void handle_explicit_input(int, const char*);
	Input_handler_controller *get_input_handler_controller(int);
	bool is_legal_input_handler_controller_id(int);
};


bool at_most_one_node_of(const char*, TiXmlNode*, TiXmlNode*&, int&);
void get_time_string(char*, const char*, bool, const char*, int, int, int, int);
void get_num_elapsed_time(int, int&, int&, int&, int&);
void add_time_offset(int, int&, int&, int&, int&, int, int);

#endif
