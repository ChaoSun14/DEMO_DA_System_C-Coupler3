/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include "runtime_cumulate_average_algorithm.h"
#include "cor_global_data.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



template<typename T> void template_cumulate_or_average(T* dst, const T* src, const int length, 
        const int computing_count, const bool do_average)
{
    if (computing_count == 1) {
        for (int i = 0; i < length; i++)
            dst[i] = src[i];
    } 
    else {
        for (int i = 0; i < length; i++)
            dst[i] += src[i];
    }
    if (do_average && computing_count != 1) {
        /// a trick
        T frac = 1 / ((T)computing_count);
        if (frac == 0) {
            /// not a float number
            for (int i = 0; i < length; i++)
                dst[i] = dst[i] / computing_count;

        } else {
            /// float number
            for (int i = 0; i < length; i++)
                dst[i] = dst[i] * frac;
        }
    }
}


Runtime_cumulate_average_algorithm::Runtime_cumulate_average_algorithm(Connection_coupling_procedure *coupling_procedure, Field_mem_info *field_src, Field_mem_info *field_dst)
{
    cumulate_average_field_info *cumulate_average_field  = new cumulate_average_field_info;
    this->coupling_procedure = coupling_procedure;
	cumulate_average_field->original_mem_info_src = NULL;
    cumulate_average_field->mem_info_src = field_src;
    cumulate_average_field->mem_info_dst = field_dst;
    cumulate_average_field->timer = NULL;
    cumulate_average_field->num_elements_in_field = field_src->get_size_of_field();
    cumulate_average_field->field_data_type = field_src->get_data_type();
    cumulate_average_field->current_computing_count = 0;
	if (field_src->get_field_data()->get_coord_value_grid() != field_dst->get_field_data()->get_coord_value_grid()) {
		cumulate_average_field->original_mem_info_src = field_src;
		if (field_src->get_num_chunks() > 0)
			cumulate_average_field->mem_info_src = memory_manager->alloc_mem(field_dst, BUF_MARK_AVERAGED_BACKUP, 0, NULL, false, false);
	}
    cumulate_average_fields.push_back(cumulate_average_field);

    comp_id = field_src->get_comp_id();
}


void Runtime_cumulate_average_algorithm::handler_cumulate_or_average(const char *data_type, void *data_dst, const void *data_src, const int length, const int computing_count, const bool do_average)
{
	if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
		template_cumulate_or_average<float>((float *) data_dst, (float *) data_src, length, computing_count, do_average);
	else if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
		template_cumulate_or_average<double>((double *) data_dst, (double *) data_src, length, computing_count, do_average);
	else if (words_are_the_same(data_type, DATA_TYPE_INT))
		template_cumulate_or_average<int>((int *) data_dst, (int *) data_src, length, computing_count, do_average);
    else if (words_are_the_same(data_type, DATA_TYPE_SHORT))
        template_cumulate_or_average<short>((short *) data_dst, (short *) data_src, length, computing_count, do_average);
	else EXECUTION_REPORT(REPORT_ERROR, -1, false, "error data type in cumulate_average algorithm\n"); 
}


void Runtime_cumulate_average_algorithm::cumulate_or_average(bool do_average)
{
	void *data_src, *data_dst;

	
    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "before cumulate or average");
    for (int i = 0; i < cumulate_average_fields.size(); i ++) {
		if (cumulate_average_fields[i]->original_mem_info_src != NULL) {
			if (cumulate_average_fields[i]->original_mem_info_src == cumulate_average_fields[i]->mem_info_src) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cumulate_average_fields[i]->original_mem_info_src->get_num_chunks() == 0, "Software error in Runtime_cumulate_average_algorithm::cumulate_or_average");
				cumulate_average_fields[i]->mem_info_src->get_field_data()->interchange_grid_data(cumulate_average_fields[i]->mem_info_dst->get_field_data()->get_coord_value_grid());
			}
			else {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cumulate_average_fields[i]->original_mem_info_src->get_num_chunks() > 0, "Software error in Runtime_cumulate_average_algorithm::cumulate_or_average");
				cumulate_average_fields[i]->original_mem_info_src->transformation_between_chunks_array(true);
				cumulate_average_fields[i]->original_mem_info_src->get_field_data()->interchange_grid_data(cumulate_average_fields[i]->mem_info_dst->get_field_data()->get_coord_value_grid());
				cumulate_average_fields[i]->mem_info_src->confirm_overall_data_buf_for_chunks();
				memcpy(cumulate_average_fields[i]->mem_info_src->get_data_buf(), cumulate_average_fields[i]->original_mem_info_src->get_data_buf(), cumulate_average_fields[i]->mem_info_src->get_size_of_field()*get_data_type_size(cumulate_average_fields[i]->mem_info_src->get_data_type()));
				cumulate_average_fields[i]->mem_info_src->transformation_between_chunks_array(false);
			}
		}
		if (cumulate_average_fields[i]->mem_info_dst->get_field_data()->get_coord_value_grid() != cumulate_average_fields[i]->mem_info_src->get_field_data()->get_coord_value_grid()) {
			cumulate_average_fields[i]->mem_info_src->transformation_between_chunks_array(true);
			cumulate_average_fields[i]->mem_info_src->get_field_data()->interchange_grid_data(cumulate_average_fields[i]->mem_info_dst->get_field_data()->get_coord_value_grid());
			cumulate_average_fields[i]->mem_info_src->transformation_between_chunks_array(false);
		}
        cumulate_average_fields[i]->mem_info_src->check_field_sum(report_internal_log_enabled, true, "(src value) before cumulate or average");
		if (cumulate_average_fields[i]->current_computing_count > 0)
        cumulate_average_fields[i]->mem_info_dst->check_field_sum(report_internal_log_enabled, true, "(dst value) before cumulate or average");
    }
    
    for (int i = 0; i < cumulate_average_fields.size(); i ++) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cumulate_average_fields[i]->mem_info_src->get_size_of_field() == cumulate_average_fields[i]->mem_info_dst->get_size_of_field(), "Software error in Runtime_cumulate_average_algorithm::Runtime_cumulate_average_algorithm");
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(cumulate_average_fields[i]->mem_info_src->get_data_type(), cumulate_average_fields[i]->mem_info_dst->get_data_type()), "Software error in Runtime_cumulate_average_algorithm::Runtime_cumulate_average_algorithm");
        cumulate_average_fields[i]->current_computing_count ++;
		int num_chunks = cumulate_average_fields[i]->mem_info_src->get_num_chunks();
		if (num_chunks == 0)
			handler_cumulate_or_average(cumulate_average_fields[i]->field_data_type, cumulate_average_fields[i]->mem_info_dst->get_data_buf(), cumulate_average_fields[i]->mem_info_src->get_data_buf(), cumulate_average_fields[i]->num_elements_in_field, cumulate_average_fields[i]->current_computing_count, do_average);
		else for (int j = 0; j < num_chunks; j ++) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, cumulate_average_fields[i]->mem_info_src->get_chunk_data_buf_size(j) == cumulate_average_fields[i]->mem_info_dst->get_chunk_data_buf_size(j), "Software error Runtime_cumulate_average_algorithm::cumulate_or_average");
			handler_cumulate_or_average(cumulate_average_fields[i]->field_data_type, cumulate_average_fields[i]->mem_info_dst->get_chunk_buf(j), cumulate_average_fields[i]->mem_info_src->get_chunk_buf(j), cumulate_average_fields[i]->mem_info_src->get_chunk_data_buf_size(j), cumulate_average_fields[i]->current_computing_count, do_average);
		}
        if (do_average) {
            EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "do average at computing count is %d", cumulate_average_fields[i]->current_computing_count);
            cumulate_average_fields[i]->current_computing_count = 0;            
        }
		if (cumulate_average_fields[i]->original_mem_info_src != NULL)
			cumulate_average_fields[i]->original_mem_info_src->use_field_values(NULL);
        else cumulate_average_fields[i]->mem_info_src->use_field_values(NULL);
        cumulate_average_fields[i]->mem_info_dst->define_field_values(false);
    }

    EXECUTION_REPORT_LOG(REPORT_LOG, comp_id, true, "after cumulate or average");
    for (int i = 0; i < cumulate_average_fields.size(); i ++) {
		if (cumulate_average_fields[i]->original_mem_info_src != NULL)
			cumulate_average_fields[i]->original_mem_info_src->get_field_data()->generate_grid_info(cumulate_average_fields[i]->original_mem_info_src->get_field_data()->get_coord_value_grid());
        cumulate_average_fields[i]->mem_info_dst->check_field_sum(report_internal_log_enabled, true, "(dst value) after cumulate or average");
    }
}


bool Runtime_cumulate_average_algorithm::run(bool do_average)
{
    cumulate_or_average(do_average);
    
    return true;
}


Runtime_cumulate_average_algorithm::~Runtime_cumulate_average_algorithm()
{
    for (int i = 0; i < cumulate_average_fields.size(); i ++) {
        delete cumulate_average_fields[i]->timer;
        delete cumulate_average_fields[i];
    }
}


void Runtime_cumulate_average_algorithm::restart_write(Restart_buffer_container *restart_buffer, const char *label)
{
    for (int i = cumulate_average_fields.size()-1; i >= 0; i --) {
        restart_buffer->dump_in_data(&(cumulate_average_fields[i]->current_computing_count), sizeof(int));
        restart_buffer->dump_in_data(&(cumulate_average_fields[i]->num_elements_in_field), sizeof(int));
        if (cumulate_average_fields[i]->current_computing_count != 0) {
            EXECUTION_REPORT_LOG(REPORT_LOG, coupling_procedure->get_inout_interface()->get_comp_id(), true, "Detect the requirements for writing intermediate field data (%s) into restart data file for a Runtime_cumulate_average_algorithm of the export interface \"%s\" to the component model \"%s\" (%d)", label, coupling_procedure->get_inout_interface()->get_interface_name(), coupling_procedure->get_coupling_connection()->get_dst_comp_full_name(), cumulate_average_fields[i]->current_computing_count);
            char temp_label[NAME_STR_SIZE*2];
            sprintf(temp_label, "%s_%s_%s", label, coupling_procedure->get_coupling_connection()->get_dst_comp_full_name(), coupling_procedure->get_coupling_connection()->get_dst_interface_name());
            restart_buffer->get_restart_mgr()->write_restart_field_data(cumulate_average_fields[i]->mem_info_dst, coupling_procedure->get_inout_interface()->get_interface_name(), temp_label, false);
            sprintf(temp_label, "restart write of field \"%s\" average algorithm of export interface \"%s\" to the import interface \"%s\" of the component model \"%s\"", cumulate_average_fields[i]->mem_info_dst->get_field_name(), coupling_procedure->get_inout_interface()->get_interface_name(), coupling_procedure->get_coupling_connection()->get_dst_interface_name(), coupling_procedure->get_coupling_connection()->get_dst_comp_full_name());
            cumulate_average_fields[i]->mem_info_dst->check_field_sum(report_internal_log_enabled, true, temp_label);
        }
    }
    int temp_int = cumulate_average_fields.size();
    restart_buffer->dump_in_data(&temp_int, sizeof(int));
}


void Runtime_cumulate_average_algorithm::restart_read(Restart_buffer_container *restart_buffer, const char *label)
{
    int num_fields, temp_current_computing_count, temp_num_elements_in_field;
    Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(comp_id);


    restart_buffer->load_restart_data(&num_fields, sizeof(int));
    if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH || time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE)
        EXECUTION_REPORT(REPORT_ERROR, comp_id, num_fields == cumulate_average_fields.size(), "Fail to restart the simulation in a branch run or continue run: the coupling configuration has been changed from the original simulation run. Please verify.");
    for (int i = 0; i < num_fields; i ++) {
        restart_buffer->load_restart_data(&temp_num_elements_in_field, sizeof(int));
        restart_buffer->load_restart_data(&temp_current_computing_count, sizeof(int));
        if (time_mgr->get_runtype_mark() == RUNTYPE_MARK_BRANCH || time_mgr->get_runtype_mark() == RUNTYPE_MARK_CONTINUE) {
            EXECUTION_REPORT(REPORT_ERROR, comp_id, temp_num_elements_in_field == cumulate_average_fields[i]->num_elements_in_field, "Fail to restart the simulation in a branch run or continue run: the coupling configuration has been changed from the original simulation run. Please verify.");
            cumulate_average_fields[i]->current_computing_count = temp_current_computing_count;
            if (temp_current_computing_count != 0) {
                char temp_label[NAME_STR_SIZE*2];
                sprintf(temp_label, "%s_%s_%s", label, coupling_procedure->get_coupling_connection()->get_dst_comp_full_name(), coupling_procedure->get_coupling_connection()->get_dst_interface_name());
                restart_buffer->get_restart_mgr()->read_restart_field_data(cumulate_average_fields[i]->mem_info_dst, coupling_procedure->get_inout_interface()->get_interface_name(), temp_label, false, NULL, false, "");
                sprintf(temp_label, "restart read of field \"%s\" average algorithm of export interface \"%s\" to the import interface \"%s\" of the component model \"%s\"", cumulate_average_fields[i]->mem_info_dst->get_field_name(), coupling_procedure->get_inout_interface()->get_interface_name(), coupling_procedure->get_coupling_connection()->get_dst_interface_name(), coupling_procedure->get_coupling_connection()->get_dst_comp_full_name());
                cumulate_average_fields[i]->mem_info_dst->check_field_sum(report_internal_log_enabled, true, temp_label);
            }
        }
    }
}

