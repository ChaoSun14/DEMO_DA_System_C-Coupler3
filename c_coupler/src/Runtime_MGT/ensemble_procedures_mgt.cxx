/***************************************************************
 *  Copyright (c) 2017, Tsinghua University.
 *  This is a source file of C-Coupler.
 *  This file was initially finished by Dr. Chao Sun and then
 *  modified by Dr. Chao Sun and Dr. Li Liu.
 *  If you have any problem,
 *  please contact Dr. Chao Sun via sunchaoo@tsinghua.edu.cn
 *  or Dr. Li Liu via liuli-cess@tsinghua.edu.cn
 ***************************************************************/


#include <mpi.h>
#include <dlfcn.h>
#include "global_data.h"
#include "ensemble_procedures_mgt.h"
#include "ensemble_field_operation.h"
#include <algorithm>
#include <unistd.h>
#include <cstring>


Field_instances_operation::Field_instances_operation(Ensemble_procedures_inst *ensemble_procedures_inst, TiXmlElement *field_instances_XML_element)
{
    int line_number;
    int member_id;
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start Field_instances_operation::Field_instances_operation");
    this->ensemble_procedures_inst = ensemble_procedures_inst;
    int member_comp_id = ensemble_procedures_inst->get_member_comp_id();
    const char *XML_file_name = ensemble_procedures_inst->get_XML_file_name();
    const char *inst_name = ensemble_procedures_inst->get_instance_name();
    this->do_field_instances_operation = true;
    this->use_statistical_method = false;
    this->statistical_periodic_timer_id = -1;
    this->do_ensemble_OP = false;
    this->field_instances_ensemble_op = strdup(ENSEMBLE_OP_DEFAULT);
    this->field_instances_statistical_method = strdup(STATISTICAL_METHOD_DEFAULT);
    char tmp_annotation[NAME_STR_SIZE];
    sprintf(tmp_annotation, "getting DA configuration of %s", inst_name);
    TiXmlNode *time_processing_XML_node = field_instances_XML_element->FirstChildElement();
    TiXmlNode *temp_node = time_processing_XML_node;
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "time_processing", true, true, false);
    if (temp_node != NULL) {
        TiXmlElement *time_processing_XML_element = temp_node->ToElement();
        const char *time_processing_type = get_XML_attribute(-1, CCPL_NAME_STR_LEN, time_processing_XML_element, "type", XML_file_name, line_number, "the the statistical method of current ensemble procedure", "configuration of the attributes of the statistical method for ensemble procedures", true);
        check_statistical_method_format(time_processing_type, NULL, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of default statistical method of all field instances");
        this->field_instances_statistical_method = strdup(time_processing_type);
        if (this->use_statistical_method) {
            int period_count, local_lag_count;
            const char *period_unit = get_XML_attribute(-1, CCPL_NAME_STR_LEN, time_processing_XML_element, "period_unit", XML_file_name, line_number, "period unit of periodic timer", tmp_annotation, true);
            const char *period_count_string = get_XML_attribute(-1, CCPL_NAME_STR_LEN, time_processing_XML_element, "period_count", XML_file_name, line_number, "period count of periodic timer", tmp_annotation, true);
            EXECUTION_REPORT(REPORT_ERROR, -1, is_string_decimal_number(period_count_string), "Error happens when using the XML configuration file \"%s\": the value (\"%s\") of the attribute \"%s\" is not a decimal integer. Please verify the XML file around the line %d.", XML_file_name, period_count_string, "period_count", line_number);
            EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(period_count_string, "%d", &period_count) == 1, "Software error in Field_instances_operation::Field_instances_operation");
            const char *local_lag_count_string = get_XML_attribute(-1, CCPL_NAME_STR_LEN, time_processing_XML_element, "local_lag_count", XML_file_name, line_number, "local lag count of periodic timer", tmp_annotation, true);
            EXECUTION_REPORT(REPORT_ERROR, -1, is_string_decimal_number(local_lag_count_string), "Error happens when using the XML configuration file \"%s\": the value (\"%s\") of the attribute \"%s\" is not a decimal integer. Please verify the XML file around the line %d.", XML_file_name, local_lag_count_string, "local_lag_count", line_number);
            EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(local_lag_count_string, "%d", &local_lag_count) == 1, "Software error in Field_instances_operation::Field_instances_operation");
            components_time_mgrs->get_time_mgr(ensemble_procedures_inst->get_member_comp_id())->check_timer_format_XML(period_unit, period_count, local_lag_count, 0, true, inst_name, XML_file_name, time_processing_XML_element->Row());
            this->statistical_periodic_timer_id = timer_mgr->define_timer(ensemble_procedures_inst->get_member_comp_id(), period_unit, period_count, local_lag_count, 0, "define statistical periodic timer of ensemble procedures instance");

            for (TiXmlNode *time_processing_field_XML_node = time_processing_XML_element->FirstChildElement(); time_processing_field_XML_node != NULL; time_processing_field_XML_node = time_processing_field_XML_node->NextSibling()) {
                if (time_processing_field_XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                    continue;
                TiXmlElement *time_processing_field_XML_element = time_processing_field_XML_node->ToElement();
                EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(time_processing_field_XML_element->Value(), "field"), "Error happens when using the XML configuration file \"%s\": the name of the attribute (currently is \"%s\") should be \"field\". Please verify the XML file around the line %d.", XML_file_name, time_processing_field_XML_element->Value(), time_processing_field_XML_element->Row());
                const char *field_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, time_processing_field_XML_element, "name", ensemble_procedures_inst->get_XML_file_name(), line_number, "name of a field", "configuration of the attributes of shared fields for ensemble procedures", true);
                bool field_existed = true;
                for (int i = 0; i < ensemble_procedures_inst->get_num_specified_field_instances(); i ++) {
                    if (words_are_the_same(ensemble_procedures_inst->API_specified_field_insts[i]->get_field_name(), field_name)) {
                        field_existed = true;
                        break;
                    }
                    else field_existed = false;
                }
                EXECUTION_REPORT(REPORT_ERROR, -1, field_existed, "Error happens when using the XML configuration file \"%s\": the field specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" has not been registered in C-Coupler. Please verify the XML file around the line %d.", XML_file_name, field_name, "name", "field", line_number);
                const field_op *existing_field = search_field_op(field_name);
                EXECUTION_REPORT(REPORT_ERROR, -1, existing_field == NULL, "Error happens when using the XML configuration file \"%s\": Can not specify the attributes of field \"%s\" again as it has already been specified around the line %d. Please verify the XML file around the line %d.", XML_file_name, field_name, existing_field->statistical_line_number, line_number);
                const char *field_statistical_method = get_XML_attribute(-1, CCPL_NAME_STR_LEN, time_processing_field_XML_element, "type", ensemble_procedures_inst->get_XML_file_name(), line_number, "statistical method of a field", "configuration of the attributes of shared fields for ensemble procedures", true);
                check_statistical_method_format(field_statistical_method, field_name, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of statistical method of the field instance");
                add_field_op(field_name, this->field_instances_ensemble_op, -1, -1, field_statistical_method, line_number);
            }
        }
    }
    else {
        this->field_instances_statistical_method = strdup(STATISTICAL_METHOD_DEFAULT);
        this->use_statistical_method = false;
        this->statistical_periodic_timer_id = -1;
        EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"time_processing\" in the XML configuration file \"%s\" is \"off\", which means the time processing is not used in ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    }
    TiXmlNode *ensemble_operation_XML_node = time_processing_XML_node->NextSibling();
    temp_node = ensemble_operation_XML_node;
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "ensemble_operation", true, true, false);
    if (temp_node != NULL) {
        TiXmlElement *ensemble_operation_XML_element = temp_node->ToElement();
        const char *ensemble_operation_type = get_XML_attribute(-1, CCPL_NAME_STR_LEN, ensemble_operation_XML_element, "type", XML_file_name, line_number, "the ensemble operation of current ensemble procedure", "configuration of the attributes of the ensemble operation for ensemble procedures", true);
        this->field_instances_member_id = check_ensemble_op_format(ensemble_operation_type, NULL, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of default ensemble operation of all field instances");
        EXECUTION_REPORT(REPORT_ERROR, -1, this->do_ensemble_OP == ensemble_procedures_inst->if_use_ensemble_components(), "Error happens when using the XML configuration file \"%s\": the ensemble operation specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" conflicts with the attribute of ensemble components. Please verify the XML file around the line %d.", XML_file_name, ensemble_operation_type, "type", "ensemble_operation", line_number);
        this->field_instances_ensemble_op = strdup(ensemble_operation_type);
        for (TiXmlNode *ensemble_operation_field_XML_node = ensemble_operation_XML_element->FirstChildElement(); ensemble_operation_field_XML_node != NULL; ensemble_operation_field_XML_node = ensemble_operation_field_XML_node->NextSibling()) {
            if (ensemble_operation_field_XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                continue;
            TiXmlElement *ensemble_operation_field_XML_element = ensemble_operation_field_XML_node->ToElement();
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(ensemble_operation_field_XML_element->Value(), "field"), "Error happens when using the XML configuration file \"%s\": the name of the attribute (currently is \"%s\") should be \"field\". Please verify the XML file around the line %d.", XML_file_name, ensemble_operation_field_XML_element->Value(), ensemble_operation_field_XML_element->Row());
            const char *field_name_op = get_XML_attribute(-1, CCPL_NAME_STR_LEN, ensemble_operation_field_XML_element, "name", ensemble_procedures_inst->get_XML_file_name(), line_number, "name of a field", "configuration of the attributes of shared fields for ensemble procedures", true);
            bool field_existed = true;
            for (int i = 0; i < ensemble_procedures_inst->get_num_specified_field_instances(); i ++) {
                if (words_are_the_same(ensemble_procedures_inst->API_specified_field_insts[i]->get_field_name(), field_name_op)) {
                    field_existed = true;
                    break;
                }
                else field_existed = false;
            }
            EXECUTION_REPORT(REPORT_ERROR, -1, field_existed, "Error happens when using the XML configuration file \"%s\": the field specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" has not been registered in C-Coupler. Please verify the XML file around the line %d.", XML_file_name, field_name_op, "name", "field", line_number);
            const field_op *existing_field_1 = search_field_op(field_name_op);
            if (existing_field_1 != NULL && existing_field_1->op_line_number != -1)
                EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when using the XML configuration file \"%s\": Can not specify the attributes of field \"%s\" again as it has already been specified around the line %d. Please verify the XML file around the line %d.", XML_file_name, field_name_op, existing_field_1->op_line_number, line_number);
            const char *field_ensemble_op = get_XML_attribute(-1, CCPL_NAME_STR_LEN, ensemble_operation_field_XML_element, "type", ensemble_procedures_inst->get_XML_file_name(), line_number, "ensemble operation of a field", "configuration of the attributes of shared fields for ensemble procedures", true);
            member_id = check_ensemble_op_format(field_ensemble_op, field_name_op, ensemble_procedures_inst->get_XML_file_name(), line_number, "configuration of the attributes of ensemble operation of the field instance");
            EXECUTION_REPORT(REPORT_ERROR, -1, this->do_ensemble_OP == ensemble_procedures_inst->if_use_ensemble_components(), "Error happens when using the XML configuration file \"%s\": the ensemble operation specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" conflicts with the attribute of ensemble components. Please verify the XML file around the line %d.", XML_file_name, field_ensemble_op, "type", "field", line_number);
            add_field_op(field_name_op, field_ensemble_op, member_id, line_number, this->field_instances_statistical_method, -1);
        }
    }
    else {
        this->field_instances_ensemble_op = strdup(ENSEMBLE_OP_NONE);
        this->do_ensemble_OP = false;
        this->field_instances_member_id = -1;
        EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"ensemble_operation\" in the XML configuration file \"%s\" is \"off\", which means the ensemble operation is not used in ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    }
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish Field_instances_operation::Field_instances_operation");
}


field_op *Field_instances_operation::search_field_op(const char *field_name)
{
    for (int i = 0; i < fields_op.size(); i ++)
        if (words_are_the_same(field_name, fields_op[i].field_name))
            return &fields_op[i];
    return NULL;
}


void Field_instances_operation::add_field_op(const char *field_name, const char *field_ensemble_op, int member_id, int line_number, const char *field_statistical_method, int statistical_line_number)
{
    field_op new_local_op;
    if (search_field_op(field_name) == NULL) {
        strcpy(new_local_op.field_name, field_name);
        strcpy(new_local_op.field_statistical_method, field_statistical_method);
        strcpy(new_local_op.field_ensemble_op, field_ensemble_op);
        new_local_op.member_id = member_id;
        new_local_op.statistical_line_number = statistical_line_number;
        new_local_op.op_line_number = line_number;
        fields_op.push_back(new_local_op);
    }
    else {
        for (int i = 0; i < fields_op.size(); i ++) {
            if (words_are_the_same(field_name, fields_op[i].field_name)) {
                if (words_are_the_same(fields_op[i].field_ensemble_op, ENSEMBLE_OP_DEFAULT)) {
                    strcpy(fields_op[i].field_ensemble_op, field_ensemble_op);
                    fields_op[i].member_id = member_id;
                    fields_op[i].op_line_number = line_number;
                }
                if (words_are_the_same(fields_op[i].field_statistical_method, STATISTICAL_METHOD_DEFAULT)) {
                    strcpy(fields_op[i].field_statistical_method, field_statistical_method);
                    fields_op[i].statistical_line_number = statistical_line_number;
                }
            }
        }
    }
}


void Field_instances_operation::check_statistical_method_format(const char *statistical_method, const char *field_name, const char *XML_file_name, int line_number, const char *statistical_method_annotation)
{
    if (field_name == NULL) {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(statistical_method, STATISTICAL_METHOD_INST) || words_are_the_same(statistical_method, STATISTICAL_METHOD_AVER) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MAX) || words_are_the_same(statistical_method, STATISTICAL_METHOD_ACCUM) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MIN), "Error happens when using the XML configuration file \"%s\": the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" must be one of \"inst/aver/accum/max/min\". Please verify the XML file around the line %d.", XML_file_name, statistical_method, "type", "time_processing", line_number);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(statistical_method, "inst"), "Error happens when using the XML configuration file \"%s\": the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" must be \"inst\", while orther statistical operations are not currently well supported. Please verify the XML file around the line %d and contact liuli-cess@tsinghua.edu.cn for more information", XML_file_name, statistical_method, "type", "time_processing", line_number);
        if (words_are_the_same(statistical_method, STATISTICAL_METHOD_INST))
            this->use_statistical_method = false;
        else this->use_statistical_method = true;
    }
    else {
        EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(statistical_method, STATISTICAL_METHOD_INST) || words_are_the_same(statistical_method, STATISTICAL_METHOD_AVER) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MAX) || words_are_the_same(statistical_method, STATISTICAL_METHOD_ACCUM) || words_are_the_same(statistical_method, STATISTICAL_METHOD_MIN), "Error happens when using the XML configuration file \"%s\": the statistical operation for filed \"%s\" specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" must be one of \"inst/aver/accum/max/min\". Please verify the XML file around the line %d.", XML_file_name, field_name, statistical_method, "type", "field", line_number);
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(statistical_method, "inst"), "Error happens when using the XML configuration file \"%s\": the statistical operation for filed \"%s\" specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" must be \"inst\", while orther statistical operations are not currently well supported. Please verify the XML file around the line %d and contact liuli-cess@tsinghua.edu.cn for more information", XML_file_name, field_name, statistical_method, "type", "field", line_number);
    }
}


int Field_instances_operation::check_ensemble_op_format(const char *ensemble_op, const char *field_name, const char *XML_file_name, int line_number, const char *ensemble_op_annotation)
{
    int member_id = -1;
    char tmp_annotation[NAME_STR_SIZE], tmp1_annotation[NAME_STR_SIZE], tmp2_annotation[NAME_STR_SIZE];
    if (field_name == NULL) {
        sprintf(tmp_annotation, "Error happens when using the XML configuration file \"%s\": the ensemble operation specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" must be one of \"none/gather/aver/anom/max/min/mem_[id]\". Please verify the XML file around the line %d.", XML_file_name, ensemble_op, "type", "ensemble_operation", line_number);
        sprintf(tmp1_annotation, "Error happens when using the XML configuration file \"%s\": the ensemble operation specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" is wrong, the ensemble member ID should >= 1 and <= %d. Please verify the XML file around the line %d.", XML_file_name, ensemble_op, "type", "ensemble_operation", ensemble_procedures_inst->num_ens_members, line_number);
        if (words_are_the_same(ensemble_op, ENSEMBLE_OP_NONE))
            this->do_ensemble_OP = false;
        else if (words_are_the_same(ensemble_op, ENSEMBLE_OP_GATHER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_ANOM)) {
            this->do_ensemble_OP = true;
            this->use_set_grids = true;
        }
        else if (words_are_the_same(ensemble_op, ENSEMBLE_OP_AVER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MAX) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MIN)) {
            this->do_ensemble_OP = true;
            this->use_set_grids = false;
        }
        else if (strstr(ensemble_op, "mem_") != NULL) {
            EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(ensemble_op, "mem_%d", &member_id) == 1, tmp_annotation);
            EXECUTION_REPORT(REPORT_ERROR, -1, member_id > 0 && member_id <= ensemble_procedures_inst->num_ens_members, tmp1_annotation);
            this->do_ensemble_OP = true;
            this->use_set_grids = false;
        }
        else EXECUTION_REPORT(REPORT_ERROR, -1, false, tmp_annotation);
    }
    else {
        sprintf(tmp_annotation, "Error happens when using the XML configuration file \"%s\": the ensemble operation for field \"%s\" specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" conflicts with the common value (\"%s\") specified by the attribute \"%s\". Please verify the XML file around the line %d.", XML_file_name, field_name, ensemble_op, "type", "field", field_instances_ensemble_op, "ensemble_operation", line_number);
        sprintf(tmp1_annotation, "Error happens when using the XML configuration file \"%s\": the ensemble operation for field \"%s\" specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" must be one of \"none/gather/aver/anom/max/min/mem_[id]\". Please verify the XML file around the line %d.", XML_file_name, field_name, ensemble_op, "type", "ensemble_operation", line_number);
        sprintf(tmp2_annotation, "Error happens when using the XML configuration file \"%s\": the ensemble operation for field \"%s\" specified by the value (currently is \"%s\") of the attribute \"%s\" in \"%s\" is wrong, the ensemble member ID should >= 1 and <= %d. Please verify the XML file around the line %d.", XML_file_name, field_name, ensemble_op, "type", "field", ensemble_procedures_inst->num_ens_members, line_number);
        if (words_are_the_same(ensemble_op, ENSEMBLE_OP_NONE))
            EXECUTION_REPORT(REPORT_ERROR, -1, !this->do_ensemble_OP, tmp_annotation);
        else {
            if (strstr(ensemble_op, "mem_") != NULL) {
                EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(ensemble_op, "mem_%d", &member_id) == 1, tmp1_annotation);
                if (member_id == 0 || member_id > ensemble_procedures_inst->num_ens_members)
                    EXECUTION_REPORT(REPORT_ERROR, -1, false, tmp2_annotation);
            }
            else EXECUTION_REPORT(REPORT_ERROR, -1, (words_are_the_same(ensemble_op, ENSEMBLE_OP_GATHER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_ANOM) || words_are_the_same(ensemble_op, ENSEMBLE_OP_AVER) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MAX) || words_are_the_same(ensemble_op, ENSEMBLE_OP_MIN)), tmp1_annotation);
            EXECUTION_REPORT(REPORT_ERROR, -1, this->do_ensemble_OP, tmp_annotation);
        }
    }
    return member_id;
}


const char *Field_instances_operation::get_field_op_statistical_method(const char *field_name)
{
    if (search_field_op(field_name) == NULL)
        return field_instances_statistical_method;
    else if (words_are_the_same(search_field_op(field_name)->field_statistical_method, STATISTICAL_METHOD_DEFAULT))
        return field_instances_statistical_method;
    else return search_field_op(field_name)->field_statistical_method;
}


const char *Field_instances_operation::get_field_op_ensemble_op(const char *field_name)
{
    if (search_field_op(field_name) == NULL)
        return field_instances_ensemble_op;
    else if (words_are_the_same(search_field_op(field_name)->field_ensemble_op, ENSEMBLE_OP_DEFAULT))
        return field_instances_statistical_method;
    else return search_field_op(field_name)->field_ensemble_op;
}


int Field_instances_operation::get_field_op_member_id(const char *field_name)
{
    if (search_field_op(field_name) == NULL)
        return field_instances_member_id;
    return search_field_op(field_name)->member_id;
}


Ensemble_procedures_inst::Ensemble_procedures_inst(int instance_id, const char *inst_name, int member_comp_id, int ensemble_nums, int ensemble_id, int size_field_inst, int size_controls, const int *field_inst_ids, const int *control_vars, const char *annotation)
{
    int line_number;
    int tmp_int;
    double time1, time2, time3;
    wtime(&time1);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start Ensemble_procedures_inst");
    check_and_verify_name_format_of_string_for_API(member_comp_id, inst_name, API_ID_ENSEMBLE_PROC_INST_INIT, "the name of the ensemble procedure(s) instance", annotation);
    MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(member_comp_id, "Ensemble_procedures_inst::Ensemble_procedures_inst");
    synchronize_comp_processes_for_API(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", annotation);
    check_API_parameter_string(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", inst_name, "inst_name", annotation);
    check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_controls, "array size of \"control_vars\" (array size)", annotation);
    check_API_parameter_int(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_field_inst, "array size of \"field_inst_ids\" (array size)", annotation);
    check_API_parameter_data_array(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", size_controls, sizeof(int), (const char*) control_vars, "control_vars", annotation);
    for (int i = 0; i < size_controls; i ++)
        this->control_vars.push_back(control_vars[i]);
    for (int i = 0; i < size_field_inst; i ++) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, memory_manager->get_field_instance(field_inst_ids[i])->get_comp_id() == member_comp_id, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the number %d field instance specified by the input parameter \"field_ids\"\" does not belong to the component model corresponding to the input parameter of \"member_comp_id\". Please verify at the model code with the annotation \"%s\".", inst_name, i + 1, annotation);
        check_API_parameter_field_instance(member_comp_id, API_ID_ENSEMBLE_PROC_INST_INIT, comm, "registering an instance of ensemble procedures", field_inst_ids[i], "field_inst_ids", annotation);
        this->API_specified_field_insts.push_back(memory_manager->get_field_instance(field_inst_ids[i]));
    }
    EXECUTION_REPORT(REPORT_ERROR, member_comp_id, components_time_mgrs->get_time_mgr(member_comp_id)->get_time_step_in_second() > 0, "Error happers when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the time step of the corresponding component model \"%s\" has not been set yet. Please specify the time step before the model code with the annotation \"%s\".", inst_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(member_comp_id, true, "Ensemble_procedures_inst::Ensemble_procedures_inst")->get_comp_full_name(), annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish API parameters check in Ensemble_procedures_inst");
    this->member_comp_id = member_comp_id;
    this->instance_name = strdup(inst_name);
    this->instance_id = instance_id;
    this->local_comm = comm;
    int member_id;
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "%s", comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name());
    if (strstr(comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name(), "_member") != NULL) {
        if (sscanf(strstr(comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name(), "_member"), "_member%d", &member_id) != 1 )
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "ERROR happens to get the ensemble id of current process in Ensemble_procedures_inst::Ensemble_procedures_inst.");
    }
    else EXECUTION_REPORT(REPORT_ERROR, -1, false, "ERROR happens to get the ensemble id of current process in Ensemble_procedures_inst::Ensemble_procedures_inst: the componant name should be extented with \"_member[id]\", while the current name is \"%s\".", comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_full_name());
    EXECUTION_REPORT(REPORT_ERROR, -1, member_id == ensemble_id, "ERROR happens to get the ensemble id of current process in Ensemble_procedures_inst::Ensemble_procedures_inst: the member ID from registered component name (%d) is different that from MPI run parameter (%d).", member_id, ensemble_id);
    this->proc_member_id = member_id;
    this->num_ens_members = ensemble_nums;
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "The number of ensemble members is %d, and the ensemble id of current process is %d", this->num_ens_members, this->proc_member_id);
    annotation_mgr->add_annotation(instance_id, "registering an instance of ensemble procedures", annotation);
    for (int i = 0; i < API_specified_field_insts.size(); i ++)
        for (int j = i + 1; j < API_specified_field_insts.size(); j ++)
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, !words_are_the_same(API_specified_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()), "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\": the input parameter \"field_inst_ids\" contains multiple elements corresponding to the same field \"%s\" (at number %d and %d elements). Please verify at the model code with the annotation \"%s\" (please specify at most one element for this field).", inst_name, API_specified_field_insts[i]->get_field_name(), i + 1, j + 1, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations in Ensemble_procedures_inst");
    this->XML_file_name = NULL;
    TiXmlElement *config_XML_element = get_XML_file_with_configuration();
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, config_XML_element != NULL, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance \"%s\" of ensemble procedures at the model code with the annotation \"%s\": no corresponding configuration information can be found from the corresponding XML files named \"%s\".", inst_name, annotation, this->XML_file_name);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about external_procedures in Ensemble_procedures_inst");
    char tmp_annotation[NAME_STR_SIZE], dll_name_path[NAME_STR_SIZE];
    sprintf(tmp_annotation, "getting DA configuration of %s", this->instance_name);
    TiXmlNode *external_procedures_node = config_XML_element->FirstChildElement();
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &external_procedures_node, "external_procedures", true, true, true);
    TiXmlElement *external_procedures_XML_element = external_procedures_node->ToElement();
    const char *p_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, external_procedures_XML_element, "procedures_name", XML_file_name, line_number, "name of external procedures", tmp_annotation, true);
    const char *p_dl_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, external_procedures_XML_element, "dll_name", XML_file_name, line_number, "dynamic library name of external procedures", tmp_annotation, true);
    const char *p_type = get_XML_attribute(-1, CCPL_NAME_STR_LEN, external_procedures_XML_element, "type", XML_file_name, line_number, "type of external procedures", tmp_annotation, true);
    sprintf(dll_name_path, "%s/%s", comp_comm_group_mgt_mgr->get_external_procedure_lib_dir(), p_dl_name);
    EXECUTION_REPORT(REPORT_ERROR, -1, access(dll_name_path, F_OK) == 0, "Error happens when using the XML configuration file \"%s\": the dynamic library specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" dose not exist. Please verify the XML file around the line %d.", XML_file_name, dll_name_path, "dll_name", "external_procedures", line_number);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(p_type, "individual") || words_are_the_same(p_type, "package"), "Error happens when using the XML configuration file \"%s\":  the dynamic library type specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" is not \"individual\" and \"package\". Please verify the XML file around the line %d.", XML_file_name, p_type, "type", "external_procedures", line_number);
    //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "procedures_name, procedures_type, procedures_dl_name, %s, %s, %s", p_name, p_type, p_dl_name);
    this->procedures_name = strdup(p_name);
    this->procedures_type = strdup(p_type);
    this->procedures_dl_name = strdup(p_dl_name);

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about ensemble_components in Ensemble_procedures_inst");
    this->ensemble_components_name = NULL;
    this->ensemble_components_num = 0;
    this->use_ensemble_components = false;
    char tmp_ensemble_components_name[NAME_STR_SIZE], XML_file_common_name_command[NAME_STR_SIZE];
    int status;
    int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(this->member_comp_id, "Ensemble_procedures_inst");
    sprintf(tmp_ensemble_components_name, "ENSEMBLE_COMP");
    this->ensemble_components_name = strdup(tmp_ensemble_components_name);
    TiXmlNode *ensemble_components_node = external_procedures_node->NextSibling();
    TiXmlNode *temp_node = ensemble_components_node;
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "ensemble_components", true, true, false);
    if (temp_node != NULL) {
        TiXmlElement *ensemble_components_XML_element = temp_node->ToElement();
        for (TiXmlNode *ensemble_components_name_XML_node = ensemble_components_XML_element->FirstChildElement(); ensemble_components_name_XML_node != NULL; ensemble_components_name_XML_node = ensemble_components_name_XML_node->NextSibling()) {
            if (ensemble_components_name_XML_node->Type() != TiXmlNode::TINYXML_ELEMENT)
                continue;
            TiXmlElement *ensemble_components_name_XML_element = ensemble_components_name_XML_node->ToElement();
            EXECUTION_REPORT(REPORT_ERROR, -1, words_are_the_same(ensemble_components_name_XML_element->Value(), "comp_ensemble_full_name"), "Error happens when using the XML configuration file \"%s\": the name of the attribute (currently is \"%s\") should be \"comp_ensemble_full_name\". Please verify the XML file arround the line number %d.", XML_file_name, ensemble_components_name_XML_element->Value(), ensemble_components_name_XML_element->Row());
            const char *comp_ensemble_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, ensemble_components_name_XML_element, "name", XML_file_name, line_number, "name of an ensemble component", "configuration of the attributes of the ensemble component name for ensemble procedures", true);
            if (local_proc_id == 0) {
                sprintf(XML_file_common_name_command, "%s %s/%s.basic_info.xml", "ls 1>/dev/null", comp_comm_group_mgt_mgr->get_components_processes_dir(), comp_ensemble_name);
                status = system(XML_file_common_name_command);
                if (status == -1 || WEXITSTATUS(status) != 0)
                    EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when using the XML configuration file \"%s\": there is no member information of component model that matches the common component name of the ensemble (\"%s\"). Please verify the XML file around the line %d.", XML_file_name, comp_ensemble_name, line_number);
            }
            EXECUTION_REPORT_LOG(REPORT_LOG, member_comp_id, true, "Dead lock may happen when tring to load the members information in the ensemble of the component model with wrong common component name of the ensemble (\"%s\"). Please verify the XML file \"%s\" around the line %d.", comp_ensemble_name, XML_file_name, line_number);
            this->ensemble_components_num = this->ensemble_components_num + 1;
            sprintf(tmp_ensemble_components_name, "%s_%s", this->ensemble_components_name, comp_ensemble_name);
            this->ensemble_components_name = strdup(tmp_ensemble_components_name);
            this->ensemble_components_full_name.push_back(strdup(comp_ensemble_name));

        }
        EXECUTION_REPORT(REPORT_ERROR, -1, this->ensemble_components_num == 1, "Error happens when using the XML configuration file \"%s\": the number of the attributes \"comp_ensemble_full_name\" (currently is %d) in \"ensemble_components\" must be 1 currently. Please verify the XML file around the line %d and contact Dr. Li Liu for special support.", XML_file_name, this->ensemble_components_num, line_number);
        this->use_ensemble_components = true;
    }
    else {
        this->use_ensemble_components = false;
        EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "The status of \"ensemble_components\" in the XML configuration file \"%s\" is \"off\", which means the ensemble host is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    }
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about periodic_timer in Ensemble_procedures_inst");
    TiXmlNode *periodic_timer_node = ensemble_components_node->NextSibling();
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &periodic_timer_node, "periodic_timer", true, true, true);
    TiXmlElement *periodic_timer_XML_element = periodic_timer_node->ToElement();
    int period_count, local_lag_count;
    const char *period_unit = get_XML_attribute(-1, CCPL_NAME_STR_LEN, periodic_timer_XML_element, "period_unit", XML_file_name, line_number, "period unit of periodic timer", tmp_annotation, true);
    const char *period_count_string = get_XML_attribute(-1, CCPL_NAME_STR_LEN, periodic_timer_XML_element, "period_count", XML_file_name, line_number, "period count of periodic timer", tmp_annotation, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, is_string_decimal_number(period_count_string), "Error happens when using the XML configuration file \"%s\": the value (\"%s\") of the attribute \"%s\" is not a decimal integer. Please verify the XML file around the line %d.", XML_file_name, period_count_string, "period_count", line_number);
    EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(period_count_string, "%d", &period_count) == 1, "Software error in Ensemble_procedures_inst::Ensemble_procedures_inst");
    const char *local_lag_count_string = get_XML_attribute(-1, CCPL_NAME_STR_LEN, periodic_timer_XML_element, "local_lag_count", XML_file_name, line_number, "local lag count of periodic timer", tmp_annotation, true);
    EXECUTION_REPORT(REPORT_ERROR, -1, is_string_decimal_number(local_lag_count_string), "Error happens when using the XML configuration file \"%s\": the value (\"%s\") of the attribute \"%s\" is not a decimal integer. Please verify the XML file around the line %d.", XML_file_name, local_lag_count_string, "local_lag_count", line_number);
    EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(local_lag_count_string, "%d", &local_lag_count) == 1, "Software error in Ensemble_procedures_inst::Ensemble_procedures_inst");
    components_time_mgrs->get_time_mgr(member_comp_id)->check_timer_format_XML(period_unit, period_count, local_lag_count, 0, true, inst_name, XML_file_name, periodic_timer_XML_element->Row());
    this->periodic_timer_id = timer_mgr->define_timer(member_comp_id, period_unit, period_count, local_lag_count, 0, "define periodic timer of ensemble procedures instance");

    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about field_instances in Ensemble_procedures_inst");
    TiXmlNode *field_instances_node = periodic_timer_node->NextSibling();
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &field_instances_node, "field_instances", true, true, true);
    TiXmlElement *field_instances_XML_element = field_instances_node->ToElement();
    this->field_instances_op = new Field_instances_operation(this, field_instances_XML_element);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to get XML configurations about processing_control in Ensemble_procedures_inst");
    this->default_working_directory = getcwd(NULL, 0);
    this->working_directory = NULL;
    this->pre_instance_script = NULL;
    this->post_instance_script = NULL;
    char full_dir[NAME_STR_SIZE * 32], full_script_name[NAME_STR_SIZE * 32 + NAME_STR_SIZE];
    TiXmlNode *processing_control_node = field_instances_node->NextSibling();
    temp_node = processing_control_node;
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "processing_control", true, true, false);
    if (temp_node != NULL) {
        TiXmlElement *processing_control_XML_element = temp_node->ToElement();
        TiXmlNode *working_directory_XML_node = temp_node->FirstChildElement();
        temp_node = working_directory_XML_node;
        get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "working_directory", true, true, false);
        if (temp_node != NULL) {
            TiXmlElement *working_directory_XML_element = temp_node->ToElement();
            const char *working_directory_path = get_XML_attribute(-1, CCPL_NAME_STR_LEN, working_directory_XML_element, "path", XML_file_name, line_number, "the working directory of current ensemble procedure", "configuration of the attributes of working directory for ensemble procedures", true);
            if (working_directory_path[0] == '/')
                strcpy(full_dir, working_directory_path);
            else sprintf(full_dir, "%s/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), working_directory_path);
            EXECUTION_REPORT(REPORT_ERROR, -1, access(full_dir, F_OK) == 0, "Error happens when using the XML configuration file \"%s\": the path specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" dose not exist. Please verify the XML file around the line %d.", XML_file_name, full_dir, "path", "working_directory", line_number);
            this->working_directory = strdup(full_dir);
        }
        else {
            this->working_directory = NULL;
            EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "This XML element \"working_directory\" in the XML configuration file \"%s\" dose not exist or its status is \"off\", which means the working directory of component model named \"%s\" is default used as the the working directory of ensemble procedures instance \"%s\".", XML_file_name, comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id)->get_comp_name(), inst_name);
        }
        TiXmlNode *config_scripts_XML_node = working_directory_XML_node->NextSibling();
        temp_node = config_scripts_XML_node;
        get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "config_scripts", true, true, false);
        if (temp_node != NULL ) {
            TiXmlElement *config_scripts_XML_element = temp_node->ToElement();
            TiXmlNode *pre_instance_script_XML_node = config_scripts_XML_element->FirstChildElement();
            temp_node = pre_instance_script_XML_node;
            get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "pre_instance_script", true, true, false);
            if (temp_node != NULL) {
                TiXmlElement *pre_instance_script_XML_element = temp_node->ToElement();
                const char *pre_instance_script_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, pre_instance_script_XML_element, "name", XML_file_name, line_number, "name of configuration script before instance run", "configuration of the attributes of config scripts for ensemble procedures", true);
                if (pre_instance_script_name[0] == '/')
                    strcpy(full_script_name, pre_instance_script_name);
                else if (this->working_directory != NULL)
                    sprintf(full_script_name, "%s/%s", this->working_directory, pre_instance_script_name);
                else sprintf(full_script_name, "%s/%s", this->default_working_directory, pre_instance_script_name);
                EXECUTION_REPORT(REPORT_ERROR, -1, access(full_script_name, F_OK) == 0 && access(full_script_name, X_OK) == 0, "Error happens when using the XML configuration file \"%s\": the file specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" dose not exist or is not executale. Please verify the XML file around the line %d.", XML_file_name, full_script_name, "name", "pre_instance_script", line_number);
                this->pre_instance_script = strdup(full_script_name);
            }
            else {
                this->pre_instance_script = NULL;
                EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "This XML element \"pre_instance_script\" in the XML configuration file \"%s\" dose not exist or its status is \"off\", which means the configuration script before instance run is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
            }
            TiXmlNode *post_instance_script_XML_node = pre_instance_script_XML_node->NextSibling();
            temp_node = post_instance_script_XML_node;
            get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &temp_node, "post_instance_script", true, true, false);
            if (temp_node != NULL) {
                TiXmlElement *post_instance_script_XML_element = temp_node->ToElement();
                const char *post_instance_script_name = get_XML_attribute(-1, CCPL_NAME_STR_LEN, post_instance_script_XML_element, "name", XML_file_name, line_number, "name of configuration script after instance run", "configuration of the attributes of shared fields for ensemble procedures", true);
                if (post_instance_script_name[0] == '/')
                    strcpy(full_script_name, post_instance_script_name);
                else if (this->working_directory != NULL)
                    sprintf(full_script_name, "%s/%s", this->working_directory, post_instance_script_name);
                else sprintf(full_script_name, "%s/%s", this->default_working_directory, post_instance_script_name);
                EXECUTION_REPORT(REPORT_ERROR, -1, access(full_script_name, F_OK) == 0 && access(full_script_name, X_OK) == 0, "Error happens when using the XML configuration file \"%s\": the file specified by the value (\"%s\") of the attribute \"%s\" in \"%s\" dose not exist or is not executale. Please verify the XML file around the line %d.", XML_file_name, full_script_name, "name", "post_instance_script", line_number);
                this->post_instance_script = strdup(full_script_name);
            }
            else {
                this->post_instance_script = NULL;
                EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "This XML element \"post_instance_script\" in the XML configuration file \"%s\" dose not exist or its status is \"off\", which means the configuration script after instance run is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
            }
        }
        else EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "This XML element \"config_scripts\" in the XML configuration file \"%s\" dose not exist or its status is \"off\", which means the configuration scripts before and after instance run are not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    }
    else EXECUTION_REPORT(REPORT_LOG, member_comp_id, true, "This XML element \"processing_control\" in the XML configuration file \"%s\" dose not exist or its status is \"off\", which means the processing control is not used in the ensemble procedures instance \"%s\".", XML_file_name, inst_name);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish getting XML configurations in Ensemble_procedures_inst");
    ensemble_procedures_mgr->add_ensemble_procedures_inst(this);
    wtime(&time2);
    if (!this->get_field_instances_op()->if_do_ensemble_OP())
        do_none_ensemble_op_initialize();
    if (this->get_field_instances_op()->if_do_ensemble_OP())
        do_ensemble_OP_initialize();
    wtime(&time3);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in Ensemble_procedures_inst: all (%lf),  ensemble_OP_initialize (%lf)", time3 - time1, time3 - time2);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish Ensemble_procedures_inst");
}


Ensemble_procedures_inst::~Ensemble_procedures_inst()
{
    delete [] XML_file_name;
    delete [] procedures_name;
    delete [] procedures_type;
    delete [] procedures_dl_name;
    delete [] ensemble_components_name;
    for (int i = 0; i < ensemble_components_full_name.size(); i ++)
        delete ensemble_components_full_name[i];
}


TiXmlElement *Ensemble_procedures_inst::get_XML_file_with_configuration()
{
    char current_XML_file_name[NAME_STR_SIZE];
    Comp_comm_group_mgt_node *member_comp_node = comp_comm_group_mgt_mgr->search_global_node(member_comp_id);
    TiXmlDocument *current_XML_file;
    bool da_instance_status = false;
    if (words_are_the_same(member_comp_node->get_comp_type(), COMP_TYPE_PSEUDO_COUPLED))
        return NULL;
    sprintf(current_XML_file_name, "%s/%s_DA_config.xml", comp_comm_group_mgt_mgr->get_ensemble_procedure_config_dir(), this->instance_name);
    XML_file_name = strdup(current_XML_file_name);
    current_XML_file = open_XML_file_to_read(member_comp_id, current_XML_file_name, member_comp_node->get_comm_group(), false);
    if (current_XML_file == NULL)
        return NULL;
    TiXmlNode *da_instance_node = get_XML_first_child_of_unique_root(member_comp_id, XML_file_name, current_XML_file);
    char tmp_annotation[NAME_STR_SIZE];
    sprintf(tmp_annotation, "getting DA configuration of %s", this->instance_name);
    get_required_XML_element(member_comp_id, XML_file_name, tmp_annotation, &da_instance_node, "da_instance", true, true, true);
    if (da_instance_node != NULL) {
        XML_file = current_XML_file;
        TiXmlElement *XML_element = da_instance_node->ToElement();
        da_instance_status = true;
        EXECUTION_REPORT_LOG(REPORT_LOG, member_comp_id, true, "In the process of calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an data assimilation instance \"%s\": find the XML file \"%s\" with the configuration information", this->instance_name, XML_file_name);
        return XML_element;
    }
    else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Ensemble_procedures_inst::get_XML_file_with_configuration");
}


void Ensemble_procedures_inst::change_to_default_work_dir(int comp_id)
{
    change_to_work_dir(comp_id, this->default_working_directory);
}


void Ensemble_procedures_inst::change_to_work_dir(int comp_id, const char *work_dir)
{
    char dir_after[NAME_STR_SIZE * 32], *tmp_dir;
    int ierr;
    MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::change_to_work_dir"));
    tmp_dir = getcwd(NULL, 0);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "before directory is \"%s\"", tmp_dir);
    free(tmp_dir);
    if (work_dir[0] == '/')
        strcpy(dir_after, work_dir);
    else sprintf(dir_after, "%s/%s", comp_comm_group_mgt_mgr->get_root_working_dir(), work_dir);
    ierr = chdir(dir_after);
    EXECUTION_REPORT(REPORT_ERROR, -1, ierr != -1, "Error happens when changing the work directory to \"%s\". Please make sure this directory existing or verify the corresponding configuration file", dir_after);
    MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::change_to_work_dir"));
    tmp_dir = getcwd(NULL, 0);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Change the ensemble work directory to \"%s\"", tmp_dir);
    free(tmp_dir);
}


void Ensemble_procedures_inst::do_none_ensemble_op_initialize()
{
    double time1, time2, time3, time4;
    wtime(&time1);
    int *fields_id = new int [API_specified_field_insts.size()];
    for (int i = 0; i < API_specified_field_insts.size(); i ++) {
        const char *data_type = API_specified_field_insts[i]->get_data_type();
        mirror_API_specified_field_insts.push_back(memory_manager->alloc_mem(API_specified_field_insts[i], BUF_MARK_ENS_DATA_TRANSFER, this->instance_id, data_type, false, false));
        memory_manager->copy_field_data_values(mirror_API_specified_field_insts[i], API_specified_field_insts[i]);
        mirror_API_specified_field_insts[i]->define_field_values(false);
        fields_id[i] = mirror_API_specified_field_insts[i]->get_field_instance_id();
    }
    int tmp_control_vars[control_vars.size() + 1];
    for (int i = 0; i < control_vars.size(); i++)
        tmp_control_vars[i] = control_vars[i];
    tmp_control_vars[control_vars.size()] = this->instance_id;
    if (this->working_directory != NULL)
        change_to_work_dir(this->member_comp_id, this->working_directory);
    wtime(&time2);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run initialize_external_procedures_inst: %s@%s", this->instance_name, this->get_procedures_name());
    this->external_instance_id = external_procedures_mgr->initialize_external_procedures_inst(instance_name, procedures_name, procedures_type, member_comp_id, procedures_dl_name, 1, control_vars.size() + 1, API_specified_field_insts.size(), 1, tmp_control_vars, fields_id, &(this->periodic_timer_id), "Initialize external data assimilation procedures");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running initialize_external_procedures_inst: %s@%s, instance_id: %x, external_instance_id: %x", this->instance_name, this->get_procedures_name(), this->instance_id, this->external_instance_id);
    wtime(&time3);
    External_procedures_inst *external_procedures_inst = external_procedures_mgr->get_procedures_inst(this->external_instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMM, " Obtain registered external procedures insts");
    external_procedures_inst->get_procedures_import_field_insts(mirror_procedures_import_field_insts);
    external_procedures_inst->get_procedures_export_field_insts(mirror_procedures_export_field_insts);

    for (int i = 0; i < mirror_procedures_import_field_insts.size(); i ++)
        for (int j = 0; j < API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_import_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()))
                this->copy_in_index.push_back(j);
    for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++)
        for (int j = 0; j < API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_export_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()))
                this->copy_out_index.push_back(j);
    if (this->working_directory != NULL)
        change_to_default_work_dir(this->member_comp_id);
    wtime(&time4);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_none_ensemble_op_initialize: all (%lf),  before_initialize_external (%lf), initialize_external (%lf), after_initialize_external (%lf)", time4 - time1, time2 - time1, time3 - time2, time4 - time3);
}


void Ensemble_procedures_inst::do_copy_in()
{
    for (int i = 0; i < copy_in_index.size(); i ++)
        memory_manager->copy_field_data_values(this->mirror_API_specified_field_insts[copy_in_index[i]], this->API_specified_field_insts[copy_in_index[i]]);
}


void Ensemble_procedures_inst::do_copy_out()
{
    for (int i = 0; i < copy_out_index.size(); i ++)
        memory_manager->copy_field_data_values(this->API_specified_field_insts[copy_out_index[i]], this->mirror_API_specified_field_insts[copy_out_index[i]]);
}


void Ensemble_procedures_inst::do_ensemble_OP(std::vector<Field_mem_info*> &field_insts_list, std::vector<Field_mem_info*> &tmp_ensemble_set_field_insts_before_ens_op, std::vector<Field_mem_info*> &tmp_ensemble_set_field_insts_one_member, int **tmp_ensemble_member_field_insts_id, std::vector<Field_mem_info*> &tmp_ensemble_set_import_field_insts_after_ens_op)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to ensemble_procedures_inst::do_ensemble_OP");
    double time1, time2, time3, time4;
    int ensemble_operation, mem_id;
    char op_mem[NAME_STR_SIZE];
    bool do_member_to_set_operation = false;
    std::vector<Field_mem_info*> tmp_member_field_insts;
    wtime(&time1);
    if (fields_member_to_set_operation.empty())
        do_member_to_set_operation = true;
    for (int i = 0; i < field_insts_list.size(); i ++) {
        if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_GATHER)) {
            tmp_ensemble_set_import_field_insts_after_ens_op.push_back(tmp_ensemble_set_field_insts_before_ens_op[i]);
            tmp_ensemble_set_import_field_insts_after_ens_op.back()->define_field_values(false);
            fields_member_to_set_operation.push_back(NULL);
        }
        else {
            mem_id = this->field_instances_op->get_field_op_member_id(field_insts_list[i]->get_field_name());
            sprintf(op_mem, "mem_%d", mem_id);
            if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_AVER))
                ensemble_operation = ENSEMBLE_OP_TYPE_MEAN;
            else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_MAX))
                ensemble_operation = ENSEMBLE_OP_TYPE_MAX;
            else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_MIN))
                ensemble_operation = ENSEMBLE_OP_TYPE_MIN;
            //else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), ENSEMBLE_OP_SUM))
            //ensemble_operation = ENSEMBLE_OP_TYPE_SUM;
            else if (words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(field_insts_list[i]->get_field_name()), op_mem))
                ensemble_operation = ENSEMBLE_OP_TYPE_ANY;
            else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Software error in Ensemble_procedures_inst::do_ensemble_OP");
            tmp_member_field_insts.clear();
            for (int n = 0; n < this->num_ens_members; n ++)
                tmp_member_field_insts.push_back(memory_manager->get_field_instance(tmp_ensemble_member_field_insts_id[n][i]));
            if (do_member_to_set_operation)
                fields_member_to_set_operation.push_back(new Member_to_set_operation(tmp_member_field_insts, tmp_ensemble_set_field_insts_one_member[i], ensemble_operation, mem_id));
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, fields_member_to_set_operation[i] != NULL, "Software error in Ensemble_procedures_inst::do_ensemble_OP");
            fields_member_to_set_operation[i]->execute();
            tmp_ensemble_set_import_field_insts_after_ens_op.push_back(tmp_ensemble_set_field_insts_one_member[i]);
        }
    }
    wtime(&time2);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP: all (%lf)", time2 - time1);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish ensemble_procedures_inst::do_ensemble_OP");
}


void Ensemble_procedures_inst::do_ensemble_OP_initialize()
{
    double time1, time2, time3, time4, time5, time6, time7, time8;
    wtime(&time1);
    if (comp_comm_group_mgt_mgr->search_comp_with_comp_name(this->ensemble_components_name) == NULL) {
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register the ensemble host component named \"%s\"", this->ensemble_components_name);
        comp_comm_group_mgt_mgr->clear_temp_comps_name();
        for (int i = 0; i < this->ensemble_components_full_name.size(); i ++)
            comp_comm_group_mgt_mgr->add_temp_comp_name(this->ensemble_components_full_name[i]);
        this->ensemble_components_id = comp_comm_group_mgt_mgr->register_ensemble_set_of_comps(this->ensemble_components_name, this->member_comp_id, this->num_ens_members, "Register the ensemble host component");
        components_time_mgrs->clone_parent_comp_time_mgr(this->ensemble_components_id, this->member_comp_id, "Register the ensemble host component");
        components_time_mgrs->set_component_time_step(this->ensemble_components_id, components_time_mgrs->get_time_mgr(this->member_comp_id)->get_time_step_in_second(), "Set host component time step");
        int num_ensemble = comp_comm_group_mgt_mgr->get_num_members_in_ensemble_set(comp_comm_group_mgt_mgr->search_global_node(this->ensemble_components_id), comp_comm_group_mgt_mgr->search_global_node(this->member_comp_id));
        EXECUTION_REPORT(REPORT_ERROR, -1, this->num_ens_members == num_ensemble, "ERROR happens in Ensemble_procedures_inst::do_ensemble_OP_initialize: the number of ensemble members get from the global communication (\"%d\") is different from that read from MPI run parameter (\"%d\").", num_ensemble, this->num_ens_members);
        this->ens_member_virtual_grid_id = original_grid_mgr->register_V1D_grid_via_data(API_ID_GRID_MGT_REG_NORMAL_1D_GRID_NO_DATA, this->ensemble_components_id, "ens_member_virtual_grid", 0, NULL, this->num_ens_members, 0.0, NULL, NULL, false, "register the virtual grid with the number of ensemble members for ensemble set");
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish registering the ensemble host component named \"%s\"", this->ensemble_components_name);
    }
    else {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, comp_comm_group_mgt_mgr->search_comp_with_comp_name(this->ensemble_components_name) != NULL, "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
        this->ensemble_components_id = comp_comm_group_mgt_mgr->search_comp_with_comp_name(this->ensemble_components_name)->get_comp_id();
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->search_grid_info("ens_member_virtual_grid", this->ensemble_components_id) != NULL, "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
        this->ens_member_virtual_grid_id = original_grid_mgr->search_grid_info("ens_member_virtual_grid", this->ensemble_components_id)->get_grid_id();
    }
    wtime(&time2);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP_initialize: register the ensemble host component (%lf)", time2 - time1);
    std::vector<Field_mem_info*> tmp_ensemble_set_field_insts_before_ens_op;
    std::vector<Field_mem_info*> tmp_ensemble_set_field_insts_one_member;
    int **tmp_ensemble_member_field_insts_id;
    tmp_ensemble_member_field_insts_id = new int *[this->num_ens_members];
    for (int i = 0; i < this->num_ens_members; i++)
        tmp_ensemble_member_field_insts_id[i] = new int [API_specified_field_insts.size()];
    alloc_ensemble_set_mem(API_specified_field_insts, tmp_ensemble_set_field_insts_before_ens_op, tmp_ensemble_set_field_insts_one_member, tmp_ensemble_member_field_insts_id, true);
    wtime(&time3);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP_initialize: alloc_ensemble_set_mem 1 (%lf)", time3 - time2);
    //Generate and execute the model temporary output interface and the temporary ensemble member input interface of Ensemble set
    int tmp_field_update_status_size = API_specified_field_insts.size() + API_specified_field_insts.size() + 1;
    int tmp_field_update_status[tmp_field_update_status_size];
    int tmp_ensemble_member_import_interface_id[this->num_ens_members];
    int tmp_model_src_field_ids[API_specified_field_insts.size()];
    int tmp_ensemble_member_dst_field_ids[API_specified_field_insts.size()];
    for (int i = 0; i < API_specified_field_insts.size(); i ++)
        tmp_model_src_field_ids[i] = API_specified_field_insts[i]->get_field_instance_id();
    for (int i = 0; i < API_specified_field_insts.size(); i ++)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "tmp_model_src_field_ids: %s", memory_manager->get_field_instance(tmp_model_src_field_ids[i])->get_field_name());
    char interface_name[NAME_STR_SIZE], tmp_annotation[NAME_STR_SIZE];
    sprintf(interface_name, "%s_tmp_model_export_interface_for_member_%d", this->instance_name, this->proc_member_id);
    sprintf(tmp_annotation, "%s register model export interface for ensemble member_%d", this->instance_name, this->proc_member_id);
    int tmp_model_export_interface_id = inout_interface_mgr->register_inout_interface(interface_name, 1, API_specified_field_insts.size(), tmp_model_src_field_ids, API_specified_field_insts.size(), this->periodic_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    this->ensemble_set_timer_id = timer_mgr->define_timer(this->ensemble_components_id, timer_mgr->get_timer(this->periodic_timer_id)->get_frequency_unit(), timer_mgr->get_timer(this->periodic_timer_id)->get_frequency_count(), timer_mgr->get_timer(this->periodic_timer_id)->get_local_lag_count(), timer_mgr->get_timer(this->periodic_timer_id)->get_remote_lag_count(), "define timer of ensemble set");
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "if this->ensemble_set_timer_id: %d", this->ensemble_set_timer_id);
    for (int i = 0; i < this->num_ens_members; i ++) {
        for (int j = 0; j < API_specified_field_insts.size(); j ++) {
            tmp_ensemble_member_dst_field_ids[j] = tmp_ensemble_member_field_insts_id[i][j];
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "tmp_ensemble_member_dst_field_ids: %s", memory_manager->get_field_instance(tmp_ensemble_member_dst_field_ids[j])->get_field_name());
        }
        sprintf(interface_name, "%s_tmp_ensemble_set_import_interface_for_member_%d", this->instance_name, i + 1);
        sprintf(tmp_annotation, "Ensemble set: %s register import interface for ensemble member_%d", this->instance_name, i + 1);

        tmp_ensemble_member_import_interface_id[i] = inout_interface_mgr->register_inout_interface(interface_name, 0, API_specified_field_insts.size(), tmp_ensemble_member_dst_field_ids, API_specified_field_insts.size(), ensemble_set_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    }
    int mem_id = this->proc_member_id;
    for (int i = 0; i < this->num_ens_members; i ++) {
        if (mem_id - 1 == i)
            EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(tmp_ensemble_member_import_interface_id[i]), inout_interface_mgr->get_interface(tmp_model_export_interface_id), -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
        else EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(tmp_ensemble_member_import_interface_id[i]), NULL, -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
    }
    int num_dst_fields_model, num_dst_fields_ensemble_member;
    sprintf(tmp_annotation, "%s execute model export interface for ensemble member_%d", this->instance_name, mem_id);
    inout_interface_mgr->execute_interface(tmp_model_export_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, tmp_field_update_status, tmp_field_update_status_size, &num_dst_fields_model, tmp_annotation);
    for (int i = 0; i < this->num_ens_members; i ++) {
        sprintf(tmp_annotation, "Ensemble set: %s execute import interface for ensemble member_%d", this->instance_name, i);
        inout_interface_mgr->execute_interface(tmp_ensemble_member_import_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, tmp_field_update_status, tmp_field_update_status_size, &num_dst_fields_ensemble_member, tmp_annotation);
    }
    //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special1");
    inout_interface_mgr->get_interface(tmp_model_export_interface_id)->get_unique_data_send_algorithm()->wait_sending_data();

    for (int i = 0; i < tmp_ensemble_set_field_insts_before_ens_op.size(); ++i)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "tmp_ensemble_set_field_insts_before_ens_op: \"%s\"", tmp_ensemble_set_field_insts_before_ens_op[i]->get_field_name());

    wtime(&time4);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP_initialize: Generate and execute temporary interface 1 (%lf)", time4 - time3);
    //Do ensemble operations based on the XML configuration
    fields_member_to_set_operation.clear();
    do_ensemble_OP(API_specified_field_insts, tmp_ensemble_set_field_insts_before_ens_op, tmp_ensemble_set_field_insts_one_member, tmp_ensemble_member_field_insts_id, ensemble_set_import_field_insts_after_ens_op);
    for (int i = 0; i < this->num_ens_members; i++)
        delete [] tmp_ensemble_member_field_insts_id[i];
    delete [] tmp_ensemble_member_field_insts_id;

    for (int i = 0; i < ensemble_set_import_field_insts_after_ens_op.size(); ++i)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "ensemble_set_import_field_insts_after_ens_op: \"%s\"", ensemble_set_import_field_insts_after_ens_op[i]->get_field_name());

    //Call external procedures initialization API by ensemble set
    int tmp_control_vars[control_vars.size() + 1];
    for (int i = 0; i < control_vars.size(); ++i)
        tmp_control_vars[i] = control_vars[i];
    tmp_control_vars[control_vars.size()] = this->instance_id;
    int tmp_field_ids[ensemble_set_import_field_insts_after_ens_op.size()];
    for (int i = 0; i < ensemble_set_import_field_insts_after_ens_op.size(); i++)
        tmp_field_ids[i] = ensemble_set_import_field_insts_after_ens_op[i]->get_field_instance_id();

    wtime(&time5);
    if (this->working_directory != NULL)
        change_to_work_dir(this->ensemble_components_id, this->working_directory);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run initialize_external_procedures_inst: %s@%s", this->instance_name, this->get_procedures_name());
    this->external_instance_id = external_procedures_mgr->initialize_external_procedures_inst(instance_name, procedures_name, procedures_type, ensemble_components_id, procedures_dl_name, 1, control_vars.size() + 1, ensemble_set_import_field_insts_after_ens_op.size(), 1, tmp_control_vars, tmp_field_ids, &(this->ensemble_set_timer_id), "Initialize external data assimilation procedures");
    external_procedures_mgr->get_procedures_inst(this->external_instance_id, API_ID_ENSEMBLE_PROC_INST_INIT, "C-Coupler internal")->set_are_ensemble_procedures();
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running initialize_external_procedures_inst: %s@%s, instance_id: %x, external_instance_id: %x", this->instance_name, this->get_procedures_name(), this->instance_id, this->external_instance_id);
    if (this->working_directory != NULL)
        change_to_default_work_dir(this->ensemble_components_id);
    wtime(&time6);
    //Get the list of input and output field instances for the external process and generate the list of input and output field instances for model and ensmble member accordingly
    External_procedures_inst *external_procedures_inst = external_procedures_mgr->get_procedures_inst(this->external_instance_id, API_ID_ENSEMBLE_PROC_INST_GET_COMM, "Obtain registered external procedures insts");
    external_procedures_inst->get_procedures_import_field_insts(mirror_procedures_import_field_insts);
    external_procedures_inst->get_procedures_export_field_insts(mirror_procedures_export_field_insts);
    for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++)
        if (!words_are_the_same(this->field_instances_op->get_field_op_ensemble_op(mirror_procedures_export_field_insts[i]->get_field_name()), ENSEMBLE_OP_GATHER))
            EXECUTION_REPORT(REPORT_ERROR, -1, false, "ERROR happens in Ensemble_procedures_inst::do_ensemble_OP_initialize: the ensemble operation \"%s\" for the field \"%s\" in DA instance \"%s\" is not supported to be an \"CCPL_PARA_TYPE_INOUT\" or \"CCPL_PARA_TYPE_OUT\" type. Please change the inout_type when declaring this field in \"%s\" or change the ensemble operation in the XML file \"%s\".", this->field_instances_op->get_field_op_ensemble_op(mirror_procedures_export_field_insts[i]->get_field_name()), mirror_procedures_export_field_insts[i]->get_field_name(), this->instance_name, this->procedures_dl_name, this->XML_file_name);

    model_export_field_insts_id = new int [mirror_procedures_import_field_insts.size()];
    ensemble_member_import_field_insts_id = new int *[this->num_ens_members];
    for (int i = 0; i < this->num_ens_members; i++)
        ensemble_member_import_field_insts_id[i] = new int [mirror_procedures_import_field_insts.size()];
    alloc_ensemble_set_mem(mirror_procedures_import_field_insts, this->ensemble_set_import_field_insts_before_ens_op, this->ensemble_set_import_field_insts_one_member, this->ensemble_member_import_field_insts_id, false);

    for (int i = 0; i < mirror_procedures_import_field_insts.size(); i ++)
        for (int j = 0; j < API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_import_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()))
                model_export_field_insts_id[i] = API_specified_field_insts[j]->get_field_instance_id();

    model_import_field_insts_id = new int [mirror_procedures_export_field_insts.size()];
    ensemble_member_export_field_insts_id = new int *[this->num_ens_members];
    for (int i = 0; i < this->num_ens_members; i++)
        ensemble_member_export_field_insts_id[i] = new int [mirror_procedures_export_field_insts.size()];
    alloc_ensemble_set_mem(mirror_procedures_export_field_insts, this->ensemble_set_export_field_insts_before_ens_op, this->ensemble_set_export_field_insts_one_member, this->ensemble_member_export_field_insts_id, false);
    wtime(&time7);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP_initialize: alloc_ensemble_set_mem 2 (%lf)", time7 - time6);
    for (int i = 0; i < mirror_procedures_export_field_insts.size(); i ++)
        for (int j = 0; j < API_specified_field_insts.size(); j ++)
            if (words_are_the_same(mirror_procedures_export_field_insts[i]->get_field_name(), API_specified_field_insts[j]->get_field_name()))
                model_import_field_insts_id[i] = API_specified_field_insts[j]->get_field_instance_id();

    //Generate the normal model output interface and the normal ensemble member input interface of Ensemble set
    sprintf(interface_name, "%s_model_export_interface_for_member_%d", this->instance_name, this->proc_member_id);
    sprintf(tmp_annotation, "%s register normal model export interface for ensemble member_%d", this->instance_name, this->proc_member_id);
    this->model_export_interface_id = inout_interface_mgr->register_inout_interface(interface_name, 1, mirror_procedures_import_field_insts.size(), model_export_field_insts_id, mirror_procedures_import_field_insts.size(), this->periodic_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    if (!mirror_procedures_export_field_insts.empty()) {
        sprintf(interface_name, "%s_model_import_interface_for_member_%d", this->instance_name, this->proc_member_id);
        sprintf(tmp_annotation, "%s register normal model import interface for ensemble member_%d", this->instance_name, this->proc_member_id);
        this->model_import_interface_id = inout_interface_mgr->register_inout_interface(interface_name, 0, mirror_procedures_export_field_insts.size(), model_import_field_insts_id, mirror_procedures_export_field_insts.size(), this->periodic_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
    }
    int tmp_ensemble_member_import_field_ids[mirror_procedures_import_field_insts.size()];
    ensemble_member_import_interface_id = new int [this->num_ens_members];
    int tmp_ensemble_member_export_field_ids[mirror_procedures_export_field_insts.size()];
    ensemble_member_export_interface_id = new int [this->num_ens_members];
    for (int i = 0; i < this->num_ens_members; i ++) {
        for (int j = 0; j < mirror_procedures_import_field_insts.size(); j ++)
            tmp_ensemble_member_import_field_ids[j] = ensemble_member_import_field_insts_id[i][j];
        sprintf(interface_name, "%s_ensemble_set_import_interface_for_member_%d", this->instance_name, i + 1);
        sprintf(tmp_annotation, "Ensemble set: %s register normal import interface for ensemble member_%d", this->instance_name, i + 1);
        ensemble_member_import_interface_id[i] = inout_interface_mgr->register_inout_interface(interface_name, 0, mirror_procedures_import_field_insts.size(), tmp_ensemble_member_import_field_ids, mirror_procedures_import_field_insts.size(), ensemble_set_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
        if (!mirror_procedures_export_field_insts.empty()) {
            for (int j = 0; j < mirror_procedures_export_field_insts.size(); j ++) tmp_ensemble_member_export_field_ids[j] = ensemble_member_export_field_insts_id[i][j];
            sprintf(interface_name, "%s_ensemble_set_export_interface_for_member_%d", this->instance_name, i + 1);
            sprintf(tmp_annotation, "Ensemble set: %s register normal export interface for ensemble member_%d", this->instance_name, i + 1);
            ensemble_member_export_interface_id[i] = inout_interface_mgr->register_inout_interface(interface_name, 1, mirror_procedures_export_field_insts.size(), tmp_ensemble_member_export_field_ids, mirror_procedures_export_field_insts.size(), ensemble_set_timer_id, 0, tmp_annotation, INTERFACE_SOURCE_REGISTER);
        }
    }
    //Do coupling connection between the model normal input and output interface and the normal ensemble member input and output interfaces
    for (int i = 0; i < this->num_ens_members; i ++) {
        if (mem_id - 1 == i) {
            EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_import_interface_id[i]), inout_interface_mgr->get_interface(this->model_export_interface_id), -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
            if (!mirror_procedures_export_field_insts.empty())
                EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_export_interface_id[i]), inout_interface_mgr->get_interface(this->model_import_interface_id), -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
        }
        else {
            EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_import_interface_id[i]), NULL, -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
            if (!mirror_procedures_export_field_insts.empty())
                EXECUTION_REPORT(REPORT_ERROR, -1, inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(inout_interface_mgr->get_interface(ensemble_member_export_interface_id[i]), NULL, -1, false), "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
        }
    }
    int model_export_field_update_status[mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1];
    sprintf(tmp_annotation, "%s execute model export interface for ensemble member_%d", this->instance_name, this->proc_member_id);
    inout_interface_mgr->execute_interface(this->model_export_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_export_field_update_status, mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1, &num_dst_fields_model, tmp_annotation);
    for (int i = 0; i < this->num_ens_members; i ++) {
        sprintf(tmp_annotation, "Ensemble set: %s execute import interface for ensemble member_%d", this->instance_name, i);
        inout_interface_mgr->execute_interface(ensemble_member_import_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_export_field_update_status, mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1, &num_dst_fields_ensemble_member, tmp_annotation);
    }
    ensemble_set_import_field_insts_after_ens_op.clear();
    fields_member_to_set_operation.clear();
    do_ensemble_OP(mirror_procedures_import_field_insts, this->ensemble_set_import_field_insts_before_ens_op, this->ensemble_set_import_field_insts_one_member, this->ensemble_member_import_field_insts_id, ensemble_set_import_field_insts_after_ens_op);
    //Reset external procedures based on new ensemble set field instances
    if (this->working_directory != NULL)
        change_to_work_dir(this->ensemble_components_id, this->working_directory);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Reset external procedures based on new ensemble set field instances: %s@%s", this->instance_name, this->get_procedures_name());
    external_procedures_mgr->get_procedures_inst(this->external_instance_id, API_ID_ENSEMBLE_PROC_INST_INIT, "Reset external procedures based on new ensemble set field instances")->reset_procedures_inst(ensemble_set_import_field_insts_after_ens_op);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish reset of external procedures based on new ensemble set field instances: %s@%s, instance_id: %x, external_instance_id: %x", this->instance_name, this->get_procedures_name(), this->instance_id, this->external_instance_id);
    if (this->working_directory != NULL)
        change_to_default_work_dir(this->ensemble_components_id);
    wtime(&time8);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP_initialize: Generate and execute normal interface 2 (%lf)", time8 - time7);
    EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in do_ensemble_OP_initialize: all (%lf), before_initialize_external (%lf), initialize_external (%lf), after_initialize_external (%lf)", time8 - time1, time5 - time1, time6 - time5, time8 - time6);
}


void Ensemble_procedures_inst::alloc_ensemble_set_mem(std::vector<Field_mem_info*> &reference_field_insts, std::vector<Field_mem_info*> &tmp_ensemble_set_field_insts_all_member, std::vector<Field_mem_info*> &tmp_ensemble_set_field_insts_one_member, int **tmp_ensemble_member_field_insts_id, bool use_default_parallel_decomp)
{
    int tmp_member_comp_or_grid_id;
    char tmp_annotation[NAME_STR_SIZE];
    std::vector<int> ensemble_set_grid_ids;
    std::vector<int> ensemble_member_grid_ids;
    std::vector<int> ensemble_set_decomp_ids;
    std::vector<Field_mem_info*> tmp_field_insts;
    ensemble_set_grid_ids.clear();
    ensemble_member_grid_ids.clear();
    ensemble_set_decomp_ids.clear();
    tmp_field_insts.clear();
    char set_md_grid_name[NAME_STR_SIZE];
    for (int i = 0; i < reference_field_insts.size(); i ++) {
        ensemble_set_grid_ids.push_back(this->ens_member_virtual_grid_id);
        ensemble_set_decomp_ids.push_back(-1);
        ensemble_member_grid_ids.push_back(-1);
        if (reference_field_insts[i]->get_grid_id() != -1) {
            if (use_default_parallel_decomp) {
                ensemble_member_grid_ids[i] = original_grid_mgr->promote_ensemble_member_grid_to_set(this->ensemble_components_id, original_grid_mgr->search_grid_info(reference_field_insts[i]->get_grid_id()), true)->get_grid_id();
                sprintf(set_md_grid_name, "set_md_grid_%s", original_grid_mgr->search_grid_info(reference_field_insts[i]->get_grid_id())->get_grid_name());
                Original_grid_info *set_md_grid = original_grid_mgr->search_grid_info(set_md_grid_name, this->ensemble_components_id);
                if (set_md_grid == NULL)
                    ensemble_set_grid_ids[i] = original_grid_mgr->register_md_grid_via_multi_grids(this->ensemble_components_id, set_md_grid_name, ensemble_member_grid_ids[i], this->ens_member_virtual_grid_id, -1, -1, NULL, false, "register grid with the number of ensemble members for ensemble set");
                else ensemble_set_grid_ids[i] = set_md_grid->get_grid_id();
                if (original_grid_mgr->search_grid_info(ensemble_set_grid_ids[i])->get_H2D_sub_CoR_grid() != NULL) {
                    Original_grid_info *set_H2D_grid = original_grid_mgr->search_grid_info(ensemble_set_grid_ids[i])->get_H2D_sub_grid();
                    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, set_H2D_grid != NULL, "Software error in Ensemble_procedures_inst::do_ensemble_OP_initialize.");
                    ensemble_set_decomp_ids[i] = decomps_info_mgr->generate_default_parallel_decomp(set_H2D_grid, -1)->get_decomp_id();
                }
            }
            else {
                Original_grid_info *field_grid = original_grid_mgr->search_grid_info(reference_field_insts[i]->get_grid_id());
                if (field_grid->get_Tracer1D_sub_grid() == NULL)  {
                    ensemble_member_grid_ids[i] = reference_field_insts[i]->get_grid_id();
                    sprintf(set_md_grid_name, "set_md_grid_%s", original_grid_mgr->search_grid_info(reference_field_insts[i]->get_grid_id())->get_grid_name());
                    Original_grid_info *set_md_grid = original_grid_mgr->search_grid_info(set_md_grid_name, this->ensemble_components_id);
                    if (set_md_grid == NULL)
                        ensemble_set_grid_ids[i] = original_grid_mgr->register_md_grid_via_multi_grids(this->ensemble_components_id, set_md_grid_name, ensemble_member_grid_ids[i], this->ens_member_virtual_grid_id, -1, -1, NULL, false, "register grid with the number of ensemble members for ensemble set");
                    else ensemble_set_grid_ids[i] = set_md_grid->get_grid_id();
                }
                else {
                    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, field_grid->get_Tracer1D_sub_grid()->get_original_CoR_grid()->get_grid_size() == this->num_ens_members, "Software error in Ensemble_procedures_inst::alloc_ensemble_set_mem.");
                    ensemble_set_grid_ids[i] = reference_field_insts[i]->get_grid_id();
                    if (field_grid->get_original_CoR_grid()->get_num_dimensions() > 1)
                        ensemble_member_grid_ids[i] = field_grid->get_max_sub_grid_under_V3D()->get_grid_id();
                    else ensemble_member_grid_ids[i] = -1;
                }
                if (original_grid_mgr->search_grid_info(ensemble_set_grid_ids[i])->get_H2D_sub_CoR_grid() != NULL)
                    ensemble_set_decomp_ids[i] = reference_field_insts[i]->get_decomp_id();
                else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, reference_field_insts[i]->get_decomp_id() == -1, "Software error in Ensemble_procedures_inst::alloc_ensemble_set_mem.");
            }
        }
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to register ensemble set field instance for %s", reference_field_insts[i]->get_field_name());
        tmp_member_comp_or_grid_id = ensemble_member_grid_ids[i];
        if (tmp_member_comp_or_grid_id == -1)
            tmp_member_comp_or_grid_id = this->ensemble_components_id;
        sprintf(tmp_annotation, "Ensemble set: register field instance of %s", reference_field_insts[i]->get_field_name());
        tmp_ensemble_set_field_insts_all_member.push_back(memory_manager->alloc_mem(reference_field_insts[i]->get_field_name(), ensemble_set_decomp_ids[i], ensemble_set_grid_ids[i], -(this->instance_id), reference_field_insts[i]->get_data_type(), reference_field_insts[i]->get_unit(), tmp_annotation, false, false));
        sprintf(tmp_annotation, "Ensemble set: register field instance of %s for ensemble member", reference_field_insts[i]->get_field_name());
        tmp_ensemble_set_field_insts_one_member.push_back(memory_manager->alloc_mem(reference_field_insts[i]->get_field_name(), ensemble_set_decomp_ids[i], tmp_member_comp_or_grid_id, -(this->instance_id), reference_field_insts[i]->get_data_type(), reference_field_insts[i]->get_unit(), tmp_annotation, false, false));
        tmp_ensemble_set_field_insts_one_member.back()->define_field_values(false);
        Field_mem_info *existing_field_instance_instance;
        for (int j = 0; j < this->num_ens_members; j ++) {
            tmp_field_insts.clear();
            sprintf(tmp_annotation, "Ensemble set: register field instance of %s for ensemble member_%d", reference_field_insts[i]->get_field_name(), j + 1);
            existing_field_instance_instance = memory_manager->search_field_instance(reference_field_insts[i]->get_field_name(), ensemble_set_decomp_ids[i], tmp_member_comp_or_grid_id, -(this->instance_id + 16 + j));
            if (existing_field_instance_instance == NULL) {
                tmp_field_insts.push_back(memory_manager->alloc_mem(reference_field_insts[i]->get_field_name(), ensemble_set_decomp_ids[i], tmp_member_comp_or_grid_id, -(this->instance_id + 16 + j), reference_field_insts[i]->get_data_type(), reference_field_insts[i]->get_unit(), tmp_annotation, false, false));
                tmp_field_insts.back()->reset_mem_buf((char *)(tmp_ensemble_set_field_insts_all_member.back()->get_data_buf()) + j * (tmp_field_insts.back()->get_size_of_field()) * get_data_type_size(tmp_field_insts.back()->get_data_type()), true, -1);
            }
            else tmp_field_insts.push_back(existing_field_instance_instance);
            EXECUTION_REPORT(REPORT_ERROR, -1, memory_manager->check_is_legal_field_instance_id(tmp_field_insts.back()->get_field_instance_id()), "Error happens when %s." , tmp_annotation);
            tmp_ensemble_member_field_insts_id[j][i] = tmp_field_insts.back()->get_field_instance_id();
        }
    }
}


void Ensemble_procedures_inst::execute_config_script(int comp_id, const char *file_name, const char *str_para0, const char *str_para1, const char *str_para2, const char *str_para3, const char *str_para4, const char *str_para5, const char *annotation)
{
    int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(comp_id, "Ensemble_procedures_inst::execute_config_script");
    MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::execute_config_script"));
    if (local_proc_id == 0) {
        char full_file_name[NAME_STR_SIZE * 16], working_directory[NAME_STR_SIZE * 16], full_command[NAME_STR_SIZE * 32], tmp_full_command[NAME_STR_SIZE * 32];
        const char *str_paras[6];
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "execute_config_script for \"%s\" begins at the model code with the annotation \"%s\"", file_name, annotation);
        realpath(file_name, full_file_name);
        getcwd(working_directory, NAME_STR_SIZE * 16);
        EXECUTION_REPORT(REPORT_ERROR, -1, does_file_exist(file_name), "ERROR happens when executing the execute_config_script at the model code with the annotation \"%s\" for the script \"%s\": this script (the absolute file name is \"%s\") does not exist. Please verify (please note that the current working directory is \"%s\").", annotation, file_name, full_file_name,    working_directory);
        str_paras[0] = str_para0;
        str_paras[1] = str_para1;
        str_paras[2] = str_para2;
        str_paras[3] = str_para3;
        str_paras[4] = str_para4;
        str_paras[5] = str_para5;
        sprintf(full_command, "\"%s\"", full_file_name);
        for (int i = 0; i < 5; i ++)
            if (strlen(str_paras[i]) > 0) {
                strcpy(tmp_full_command, full_command);
                sprintf(full_command, "%s \"%s\"", tmp_full_command, str_paras[i]);
            }
        EXECUTION_REPORT(REPORT_ERROR, -1, system(full_command) != -1, "ERROR happens when executing the execute_config_script at the model code with the annotation \"%s\" for the script \"%s\"    (\"%s\"): fail to execute this script. Please verify.", annotation, file_name, full_file_name);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "execute_config_script for \"%s\" ends", file_name);
    }
    MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "Ensemble_procedures_inst::execute_config_script"));
}


void Ensemble_procedures_inst::run(const char *annotation)
{
    double time1, time2, time3, time4, sum_time, mean_time;
    check_for_ccpl_managers_allocated(API_ID_TIME_MGT_IS_TIMER_ON, annotation);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to Ensemble_procedures_inst::run: %s@%s", this->instance_name, this->get_procedures_name());
    if (ensemble_components_id != member_comp_id) {
        Time_mgt *ensemble_time_mgr = components_time_mgrs->get_time_mgr(ensemble_components_id);
        Time_mgt *member_time_mgr = components_time_mgrs->get_time_mgr(member_comp_id);
        member_time_mgr->clone_time_mgr(ensemble_time_mgr, ensemble_components_id);
    }
    if (timer_mgr->is_timer_on(this->periodic_timer_id, annotation)) {
        char full_date[15], date[9], hours[3], minutes[3], seconds[3];
        sprintf(date, "%08d", components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_date());
        std::strncpy(date, date, 8);
        sprintf(hours, "%02d", (components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second()) / 3600);
        std::strncpy(hours, hours, 2);
        sprintf(minutes, "%02d", ((components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second()) % 3600) / 60);
        std::strncpy(minutes, minutes, 2);
        sprintf(seconds, "%02d", ((components_time_mgrs->get_time_mgr(this->member_comp_id)->get_current_second()) % 3600) % 60);
        std::strncpy(seconds, seconds, 2);
        sprintf(full_date, "%08s%02s%02s%02s", date, hours, minutes, seconds);
        std::strncpy(full_date, full_date, 14);
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "current full date: %14s", full_date);
        if (this->get_field_instances_op()->if_do_field_instances_operation() && !this->get_field_instances_op()->if_do_ensemble_OP()) {
            MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->member_comp_id, "Ensemble_procedures_inst::run"));
            wtime(&time1);
            this->do_copy_in();
            wtime(&time2);
            mean_time = time2 - time1;
            MPI_Reduce(&mean_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->member_comp_id, "Ensemble_procedures_inst::run"));
            EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in Ensemble_procedures_inst::run: individual run copy in (%lf)", sum_time/comp_comm_group_mgt_mgr->get_num_proc_in_comp(this->member_comp_id, "Ensemble_procedures_inst::run"));
            if (this->working_directory != NULL)
                change_to_work_dir(this->member_comp_id, this->working_directory);
            if (this->pre_instance_script != NULL )
                execute_config_script(this->member_comp_id, this->pre_instance_script, full_date, "", "", "", "", "", "before instance configuration script");
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Start to run external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            check_for_ccpl_managers_allocated(API_ID_EXTERNAL_PROC_INST_RUN, annotation);
            external_procedures_mgr->get_procedures_inst(this->get_external_instance_id(), API_ID_EXTERNAL_PROC_INST_RUN, annotation)->run(-1, annotation);
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            if (this->post_instance_script != NULL )
                execute_config_script(this->member_comp_id, this->post_instance_script, full_date, "", "", "", "", "", "after instance configuration script");
            if (this->working_directory != NULL)
                change_to_default_work_dir(this->member_comp_id);
            MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->member_comp_id, "Ensemble_procedures_inst::run"));
            wtime(&time3);
            this->do_copy_out();
            wtime(&time4);
            mean_time = time4 - time3;
            MPI_Reduce(&mean_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->member_comp_id, "Ensemble_procedures_inst::run"));
            EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in Ensemble_procedures_inst::run: individual run copy out (%lf)", sum_time/comp_comm_group_mgt_mgr->get_num_proc_in_comp(this->member_comp_id, "Ensemble_procedures_inst::run"));
        }
        if (this->get_field_instances_op()->if_do_field_instances_operation() && this->get_field_instances_op()->if_do_ensemble_OP()) {
            char tmp_annotation[NAME_STR_SIZE];
            int num_dst_fields_model, num_dst_fields_ensemble_member;
            int model_export_field_update_status[mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1];
            int mem_id = this->proc_member_id;
            MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->ensemble_components_id, "Ensemble_procedures_inst::run"));
            wtime(&time1);
            sprintf(tmp_annotation, "%s execute model export interface for ensemble member_%d", this->instance_name, this->proc_member_id);
            inout_interface_mgr->execute_interface(this->model_export_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_export_field_update_status, mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1, &num_dst_fields_model, tmp_annotation);
            for (int i = 0; i < this->num_ens_members; i ++) {
                sprintf(tmp_annotation, "Ensemble set: %s execute import interface for ensemble member_%d", this->instance_name, i);
                inout_interface_mgr->execute_interface(ensemble_member_import_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_export_field_update_status, mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1, &num_dst_fields_ensemble_member, tmp_annotation);
            }
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special2");
            inout_interface_mgr->get_interface(this->model_export_interface_id)->get_unique_data_send_algorithm()->wait_sending_data();
            wtime(&time2);
            mean_time = time2 - time1;
            MPI_Reduce(&mean_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->ensemble_components_id, "Ensemble_procedures_inst::run"));
            EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in Ensemble_procedures_inst::run: ensemble run before run %s (%lf)", annotation, sum_time/comp_comm_group_mgt_mgr->get_num_proc_in_comp(this->ensemble_components_id, "Ensemble_procedures_inst::run"));
            ensemble_set_import_field_insts_after_ens_op.clear();
            do_ensemble_OP(mirror_procedures_import_field_insts, this->ensemble_set_import_field_insts_before_ens_op, this->ensemble_set_import_field_insts_one_member, this->ensemble_member_import_field_insts_id, ensemble_set_import_field_insts_after_ens_op);
            if (this->working_directory != NULL)
                change_to_work_dir(this->ensemble_components_id, this->working_directory);
            if (this->pre_instance_script != NULL )
                execute_config_script(this->ensemble_components_id, this->pre_instance_script, full_date, "", "", "", "", "", "before instance configuration script");
            EXECUTION_REPORT_LOG(REPORT_LOG, this->member_comp_id, true, "Start to run external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            check_for_ccpl_managers_allocated(API_ID_EXTERNAL_PROC_INST_RUN, annotation);
            external_procedures_mgr->get_procedures_inst(this->get_external_instance_id(), API_ID_EXTERNAL_PROC_INST_RUN, annotation)->run(-1, annotation);
            EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Finish running external procedures: %s@%s", this->instance_name, this->get_procedures_name());
            if (this->post_instance_script != NULL )
                execute_config_script(this->ensemble_components_id, this->post_instance_script, full_date, "", "", "", "", "", "after instance configuration script");
            if (this->working_directory != NULL)
                change_to_default_work_dir(this->ensemble_components_id);
            MPI_Barrier(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->ensemble_components_id, "Ensemble_procedures_inst::run"));
            wtime(&time3);
            if (!mirror_procedures_export_field_insts.empty()) {
                int model_import_field_update_status[mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1];
                for (int i = 0; i < this->num_ens_members; i ++) {
                    sprintf(tmp_annotation, "Ensemble set: execute export interface for ensemble member_%d", i);
                    inout_interface_mgr->execute_interface(ensemble_member_export_interface_id[i], API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_import_field_update_status, mirror_procedures_export_field_insts.size() + mirror_procedures_export_field_insts.size() + 1, &num_dst_fields_ensemble_member, tmp_annotation);
                    if (mem_id - 1 == i) {
                        sprintf(tmp_annotation, "execute model import interface for ensemble member_%d", this->proc_member_id);
                        inout_interface_mgr->execute_interface(this->model_import_interface_id, API_ID_INTERFACE_EXECUTE_WITH_ID, true, model_import_field_update_status, mirror_procedures_import_field_insts.size() + mirror_procedures_export_field_insts.size() + 1, &num_dst_fields_model, tmp_annotation);
                    }
                    //EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special3");
                    //inout_interface_mgr->get_interface(ensemble_member_export_interface_id[i])->get_unique_data_send_algorithm()->wait_sending_data();
                }
            }
            if (!mirror_procedures_export_field_insts.empty()) {
                for (int i = 0; i < this->num_ens_members; i ++) {
                    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "liuli special3");
                    inout_interface_mgr->get_interface(ensemble_member_export_interface_id[i])->get_unique_data_send_algorithm()->wait_sending_data();
                }
            }
            wtime(&time4);
            mean_time = time4 - time3;
            MPI_Reduce(&mean_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(this->ensemble_components_id, "Ensemble_procedures_inst::run"));
            EXECUTION_REPORT(REPORT_LOG, -1, true, "TIME in Ensemble_procedures_inst::run: ensemble run after run %s (%lf)", annotation, sum_time/comp_comm_group_mgt_mgr->get_num_proc_in_comp(this->ensemble_components_id, "Ensemble_procedures_inst::run"));
        }
    }
}


Ensemble_procedures_mgt::Ensemble_procedures_mgt()
{
}


Ensemble_procedures_mgt::~Ensemble_procedures_mgt()
{
    for (int i = 0; i < registered_ensemble_procedures_insts.size(); i ++)
        delete registered_ensemble_procedures_insts[i];
}


int Ensemble_procedures_mgt::initialize_ensemble_procedures_inst(const char *inst_name, int member_comp_id, int ensemble_nums, int ensemble_id, int size_field_inst, int size_controls,
        const int *field_inst_ids, const int *control_vars, const char *annotation)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start Ensemble_procedures_mgt::initialize_ensemble_procedures_inst");
    for (int i = 0; i < registered_ensemble_procedures_insts.size(); i ++)
        if (registered_ensemble_procedures_insts[i]->get_member_comp_id() == member_comp_id && words_are_the_same(registered_ensemble_procedures_insts[i]->get_instance_name(), inst_name))
            EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, member_comp_id, false, "ERROR happens when calling the API \"CCPL_ensemble_procedures_inst_init\" to initialize an instance named \"%s\" at the model code with the annotation \"%s\": another instance with the same name has been initialized before (at the model code with the annotation \"%s\"). Please verify.", inst_name, annotation, annotation_mgr->get_annotation(registered_ensemble_procedures_insts[i]->get_instance_id(), "registering an instance of ensemble procedures"));
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "registered_ensemble_procedures_insts.size is \"%d\"", registered_ensemble_procedures_insts.size());
    int instance_id = TYPE_ENS_PROCEDURE_PREFIX | (registered_ensemble_procedures_insts.size() << 16);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to Ensemble_procedures_inst (instance_id = %x)", instance_id);
    Ensemble_procedures_inst *procedures_inst = new Ensemble_procedures_inst(instance_id, inst_name, member_comp_id, ensemble_nums, ensemble_id, size_field_inst, size_controls,
            field_inst_ids, control_vars, annotation);
    return instance_id;
}


Ensemble_procedures_inst *Ensemble_procedures_mgt::get_procedures_inst(int instance_id, int API_id, const char *annotation)
{
    char API_label[NAME_STR_SIZE];
    int instance_index;
    get_API_hint(-1, API_id, API_label);
    instance_index = GET_ENS_PROCEDURES_INST_INDEX(instance_id);
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "Ensemble_procedures_inst::get_procedures_inst (instance_id = %x, instance_index = %d)", instance_id, instance_index);
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (instance_id & TYPE_ID_PREFIX_MASK) == TYPE_ENS_PROCEDURE_PREFIX && instance_index < registered_ensemble_procedures_insts.size(), "ERROR happens when calling the API \"%s\": this API is not called in the native code of an ensemble procedure. Please verify the model code with the anotation \"%s\"", API_label, annotation);
    return registered_ensemble_procedures_insts[instance_index];
}


int Ensemble_procedures_mgt::get_ensemble_size(int external_inst_id, const char *annotation)
{
    int num_control_vars, ens_inst_id;
    num_control_vars = external_procedures_mgr->get_instance_num_control_vars(external_inst_id, "API \"CCPL_ensemble_procedures_external_inst_get_ensemble_num\"");
    ens_inst_id = external_procedures_mgr->get_instance_control_var(external_inst_id, num_control_vars, "API \"CCPL_ensemble_procedures_external_inst_get_ensemble_num\"");
    return get_procedures_inst(ens_inst_id, API_ID_ENSEMBLE_PROC_INST_GET_ENSEMBLE_NUM, annotation)->get_ensemble_size();
}


int Ensemble_procedures_mgt::get_ensemble_member_index(int external_inst_id, const char *annotation)
{
    int num_control_vars, ens_inst_id;
    num_control_vars = external_procedures_mgr->get_instance_num_control_vars(external_inst_id, "API \"CCPL_ensemble_procedures_external_inst_get_ensemble_num\"");
    ens_inst_id = external_procedures_mgr->get_instance_control_var(external_inst_id, num_control_vars, "API \"CCPL_ensemble_procedures_external_inst_get_ensemble_num\"");
    return get_procedures_inst(ens_inst_id, API_ID_ENSEMBLE_PROC_INST_GET_ENSEMBLE_ID, annotation)->get_ensemble_member_index();
}
