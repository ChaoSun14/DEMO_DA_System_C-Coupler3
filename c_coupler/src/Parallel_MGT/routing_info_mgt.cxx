/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu and then
  *  modified by Dr. Cheng Zhang and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn or
  *  Dr. Cheng Zhang via zhangc-cess@tsinghua.edu.cn
  ***************************************************************/


#include <mpi.h>
#include "routing_info_mgt.h"
#include "global_data.h"
#include "cor_global_data.h"
#include "CCPL_api_mgt.h"
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <time.h>


Routing_info *Routing_info_mgt::search_or_add_router(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    Routing_info *router;

    router = search_router(src_comp_id, dst_comp_id, src_decomp_name, dst_decomp_name);

    if (router != NULL) {
		EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_id, true, "The router from (%s %s) to (%s %s) already exists", comp_comm_group_mgt_mgr->search_global_node(src_comp_id)->get_full_name(), src_decomp_name, comp_comm_group_mgt_mgr->search_global_node(dst_comp_id)->get_full_name(), dst_decomp_name);
        return router;
    }

    router = new Routing_info(src_comp_id, dst_comp_id, src_decomp_name, dst_decomp_name);
    routers.push_back(router);

    return router;
}


Routing_info_mgt::~Routing_info_mgt()
{
    for (int i = 0; i < routers.size(); i ++)
        delete routers[i];
}


Routing_info *Routing_info_mgt::search_router(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    for (int i = 0; i < routers.size(); i ++)
        if (routers[i]->match_router(src_comp_id, dst_comp_id, src_decomp_name, dst_decomp_name))
            return routers[i];

    return NULL;
}


Routing_info::Routing_info(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
	double time1, time2;

	wtime(&time1);
    src_decomp_info = decomps_info_mgr->search_decomp_info(src_decomp_name, src_comp_id);
    dst_decomp_info = decomps_info_mgr->search_decomp_info(dst_decomp_name, dst_comp_id);
    this->src_comp_id = src_comp_id;
    this->dst_comp_id = dst_comp_id;
    src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_comp_id);
    dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_comp_id);
    strcpy(this->index_dst_comp_full_name, dst_comp_node->get_comp_full_name());
	if (src_decomp_info != NULL)
		src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_decomp_info->get_host_comp_id());
    if (dst_decomp_info != NULL)
        dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(dst_decomp_info->get_host_comp_id());
    strcpy(this->src_comp_full_name, src_comp_node->get_comp_full_name());
    strcpy(this->true_dst_comp_full_name, dst_comp_node->get_comp_full_name());
    src_comp_node_id = src_comp_node->get_comp_id();
    dst_comp_node_id = dst_comp_node->get_comp_id();
    strcpy(this->src_decomp_name, src_decomp_name);
    strcpy(this->dst_decomp_name, dst_decomp_name);
    src_decomp_size = 0;
    dst_decomp_size = 0;
    current_proc_id_src_comp = src_comp_node->get_current_proc_local_id();
    current_proc_id_dst_comp = dst_comp_node->get_current_proc_local_id();

    src_local_routing_mapping_table_entries = NULL;
    dst_local_routing_mapping_table_entries = NULL;
    src_remote_routing_mapping_table_entries = NULL;
    dst_remote_routing_mapping_table_entries = NULL;
    num_src_local_routing_mapping_table_entries = 0;
    num_dst_local_routing_mapping_table_entries = 0;
    num_src_remote_routing_mapping_table_entries = 0;
    num_dst_remote_routing_mapping_table_entries = 0;

    distribute_sorting = new Distribute_merge_sort<routing_mapping_table_entry>(current_proc_id_src_comp, current_proc_id_dst_comp, src_comp_node, dst_comp_node);

    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_id, true, "Start to generate router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_id, true, "Start to generate router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);

    if (words_are_the_same(src_decomp_name, "NULL")) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(dst_decomp_name, "NULL"), "for router of scalar variables, the local and remote decompositions must be \"NULL\"\n");
        num_dimensions = 0;
        if (current_proc_id_src_comp != -1) 
            src_decomp_size = 1;
        if (current_proc_id_dst_comp != -1) 
            dst_decomp_size = 1;
    }
    else {
        if (current_proc_id_src_comp != -1)
            src_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information new");
        else dst_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information new");
        calculate_routing_mapping_tables();
		if (!report_error_enabled)
			build_router_based_on_routing_mapping_tables();
        if (current_proc_id_src_comp != -1)
            src_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information new");
        else dst_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information new");        
        num_dimensions = 2;
        if (current_proc_id_src_comp != -1)
            src_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information old");
        else dst_comp_node->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information old");
		if (report_error_enabled)
	        build_2D_router();
        if (current_proc_id_src_comp != -1)
            src_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information old");
        else dst_comp_node->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_COMPUTATION, -1, -1, "generating routing information old");
        if (current_proc_id_src_comp != -1) 
            src_decomp_size = src_decomp_info->get_num_local_cells();
        if (current_proc_id_dst_comp != -1) 
            dst_decomp_size = dst_decomp_info->get_num_local_cells();
    }

    if (current_proc_id_src_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, src_comp_id, true, "Finish generating router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_LOG(REPORT_LOG, dst_comp_id, true, "Finish generating router from (%s %s) to (%s %s)", src_comp_full_name, src_decomp_name, index_dst_comp_full_name, dst_decomp_name);

    delete distribute_sorting;
	distribute_sorting = NULL;

	wtime(&time2);
	EXECUTION_REPORT(REPORT_LOG, -1, true, "use %lf seconds for routing information generation", time2 - time1);
}


void Routing_info::output_routing_mapping_table(common_sort_struct<routing_mapping_table_entry> *routing_mapping_table_entries, char *hint, int num_local_routing_mapping_table_entries, int current_proc_id)
{
    EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "output routing_mapping_table with the hint \"%s\"", hint);
    for (int i = 0; i < num_local_routing_mapping_table_entries; i ++)
        EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "table entry (%d): %d %d %d %d %d\n", i, routing_mapping_table_entries[i].key, routing_mapping_table_entries[i].content.global_index, routing_mapping_table_entries[i].content.local_process_id, routing_mapping_table_entries[i].content.local_index, routing_mapping_table_entries[i].content.remote_local_process_id, routing_mapping_table_entries[i].content.remote_local_index);
}

void Routing_info::exchange_routing_mapping_tables_between_components(int num_src_procs_adjust, int num_dst_procs_adjust, Decomp_info *decomp_info, bool exchange_src)
{
    int comp_max_num_proc_adjust, src_proc_scale_factor, dst_proc_scale_factor;
    int remote_proc_local_id_recv_from, remote_proc_local_id_send_to, *set_remote_proc_global_id_send_to, *set_remote_proc_global_id_recv_from, j, k;
    int send_addr_disp = 0, *set_send_entry_count;
    int min_grid_index, max_grid_index;
    common_sort_struct<routing_mapping_table_entry> **set_send_routing_mapping_table_entries_buffer, **set_recv_routing_mapping_table_entries_buffer, *send_routing_mapping_table_entries, **recv_routing_mapping_table_entries;
    int *set_num_recv_routing_mapping_table_entries, num_send_routing_mapping_table_entries, *num_recv_routing_mapping_table_entries;
    int *set_send_recv_mark, send_proc_scale_factor, max_proc_scale_factor;
    
    if (!(current_proc_id_src_comp != -1 && current_proc_id_src_comp < num_src_procs_adjust || current_proc_id_dst_comp != -1 && current_proc_id_dst_comp < num_dst_procs_adjust))
        return;

    comp_max_num_proc_adjust = std::max(num_src_procs_adjust, num_dst_procs_adjust);
    src_proc_scale_factor = comp_max_num_proc_adjust / num_src_procs_adjust;
    dst_proc_scale_factor = comp_max_num_proc_adjust / num_dst_procs_adjust;
    max_proc_scale_factor = std::max(src_proc_scale_factor, dst_proc_scale_factor);
    set_send_recv_mark = new int [max_proc_scale_factor];
    set_send_routing_mapping_table_entries_buffer = new common_sort_struct<routing_mapping_table_entry> *[max_proc_scale_factor];
    set_recv_routing_mapping_table_entries_buffer = new common_sort_struct<routing_mapping_table_entry> *[max_proc_scale_factor];
	for (int i = 0; i < max_proc_scale_factor; i ++)
		set_recv_routing_mapping_table_entries_buffer[i] = NULL;
    set_num_recv_routing_mapping_table_entries = new int [max_proc_scale_factor];
    set_send_entry_count = new int [max_proc_scale_factor];
    set_remote_proc_global_id_send_to = new int [max_proc_scale_factor];
    set_remote_proc_global_id_recv_from = new int [max_proc_scale_factor];
    send_routing_mapping_table_entries = NULL;
    recv_routing_mapping_table_entries = NULL;
    num_recv_routing_mapping_table_entries = NULL;
    
    if (exchange_src) {
        if (current_proc_id_src_comp != -1 && current_proc_id_src_comp < num_src_procs_adjust) {
            send_routing_mapping_table_entries = src_local_routing_mapping_table_entries;
            num_send_routing_mapping_table_entries = num_src_local_routing_mapping_table_entries;
            send_proc_scale_factor = src_proc_scale_factor;
        }
        if (current_proc_id_dst_comp != -1 && current_proc_id_dst_comp < num_dst_procs_adjust) {
            recv_routing_mapping_table_entries = &src_remote_routing_mapping_table_entries;
            num_recv_routing_mapping_table_entries = &num_src_remote_routing_mapping_table_entries;
        }
    }
    else {
        if (current_proc_id_src_comp != -1 && current_proc_id_src_comp < num_src_procs_adjust) {
            recv_routing_mapping_table_entries = &dst_remote_routing_mapping_table_entries;
            num_recv_routing_mapping_table_entries = &num_dst_remote_routing_mapping_table_entries;
        }
        if (current_proc_id_dst_comp != -1 && current_proc_id_dst_comp < num_dst_procs_adjust) {
            send_routing_mapping_table_entries = dst_local_routing_mapping_table_entries;
            num_send_routing_mapping_table_entries = num_dst_local_routing_mapping_table_entries;
            send_proc_scale_factor = dst_proc_scale_factor;
        }
    }

    for (int i = 0; i < max_proc_scale_factor; i ++) {
        set_send_recv_mark[i] = 0;
        set_remote_proc_global_id_send_to[i] = -1;
        set_remote_proc_global_id_recv_from[i] = -1; 
        if (exchange_src) {
            if (current_proc_id_src_comp != -1 && current_proc_id_src_comp < num_src_procs_adjust && (i%dst_proc_scale_factor) == 0) {
                set_send_recv_mark[i] = set_send_recv_mark[i] | ROUTER_SEND;
                remote_proc_local_id_send_to = current_proc_id_src_comp * num_dst_procs_adjust / num_src_procs_adjust + i;
                set_remote_proc_global_id_send_to[i] = dst_comp_node->get_local_proc_global_id(remote_proc_local_id_send_to);
            }
            if (current_proc_id_dst_comp != -1 && current_proc_id_dst_comp < num_dst_procs_adjust && (i%src_proc_scale_factor) == 0) {
                set_send_recv_mark[i] = set_send_recv_mark[i] | ROUTER_RECV;
                remote_proc_local_id_recv_from = current_proc_id_dst_comp * num_src_procs_adjust / num_dst_procs_adjust + i;
                set_remote_proc_global_id_recv_from[i] = src_comp_node->get_local_proc_global_id(remote_proc_local_id_recv_from);
            }
        }
        else {
            if (current_proc_id_src_comp != -1 && current_proc_id_src_comp < num_src_procs_adjust && (i%dst_proc_scale_factor) == 0) {
                set_send_recv_mark[i] = set_send_recv_mark[i] | ROUTER_RECV;
                remote_proc_local_id_recv_from = current_proc_id_src_comp * num_dst_procs_adjust / num_src_procs_adjust + i;
                set_remote_proc_global_id_recv_from[i] = dst_comp_node->get_local_proc_global_id(remote_proc_local_id_recv_from);
            }
            if (current_proc_id_dst_comp != -1 && current_proc_id_dst_comp < num_dst_procs_adjust && (i%src_proc_scale_factor) == 0) {
                set_send_recv_mark[i] = set_send_recv_mark[i] | ROUTER_SEND;
                remote_proc_local_id_send_to = current_proc_id_dst_comp * num_src_procs_adjust / num_dst_procs_adjust + i;
                set_remote_proc_global_id_send_to[i] = src_comp_node->get_local_proc_global_id(remote_proc_local_id_send_to);
            }
        }
        if (set_send_recv_mark[i] == 0)
            continue;
        if ((set_send_recv_mark[i] & ROUTER_SEND) > 0) {
            if (send_proc_scale_factor == 1) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, send_addr_disp == 0, "Software error in Routing_info::exchange_routing_mapping_tables_between_components");
                set_send_entry_count[i] = num_send_routing_mapping_table_entries;
            } else {
                distribute_sorting->calculate_min_max_grid_index(decomp_info->get_num_global_cells(), remote_proc_local_id_send_to, remote_proc_local_id_send_to, comp_max_num_proc_adjust, 1, &min_grid_index, &max_grid_index);
                for (; send_addr_disp < num_send_routing_mapping_table_entries; send_addr_disp ++) {
                    if(send_routing_mapping_table_entries[send_addr_disp].content.global_index >= min_grid_index)
                        break;
                }
                set_send_entry_count[i] = 0;
                for (j = send_addr_disp; j < num_send_routing_mapping_table_entries; j ++, set_send_entry_count[i] ++) {
                    if((send_routing_mapping_table_entries)[j].content.global_index >= max_grid_index)
                        break;
                }
				if (i == send_proc_scale_factor - 1)
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, j == num_send_routing_mapping_table_entries, "Software error in Routing_info::exchange_routing_mapping_tables_between_components");					
            }
			set_send_routing_mapping_table_entries_buffer[i] = send_routing_mapping_table_entries + send_addr_disp;
        }
    }

    distribute_sorting->exchange_sorted_data_between_processes(max_proc_scale_factor, set_send_routing_mapping_table_entries_buffer, set_recv_routing_mapping_table_entries_buffer, set_send_entry_count, set_num_recv_routing_mapping_table_entries, set_remote_proc_global_id_send_to, set_remote_proc_global_id_recv_from, MPI_COMM_WORLD, set_send_recv_mark);

    if (num_recv_routing_mapping_table_entries != NULL) {
        *num_recv_routing_mapping_table_entries = 0;
        for (int i = 0; i < max_proc_scale_factor; i ++) 
            if ((set_send_recv_mark[i] & ROUTER_RECV) > 0)
                *num_recv_routing_mapping_table_entries += set_num_recv_routing_mapping_table_entries[i];
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, *recv_routing_mapping_table_entries == NULL, "Software error in Routing_info::exchange_routing_mapping_tables_between_components");
        *recv_routing_mapping_table_entries = new common_sort_struct<routing_mapping_table_entry> [*num_recv_routing_mapping_table_entries];
        *num_recv_routing_mapping_table_entries = 0;
        for (int i = 0; i < max_proc_scale_factor; i ++) 
            if ((set_send_recv_mark[i] & ROUTER_RECV) > 0 && set_recv_routing_mapping_table_entries_buffer[i] != NULL) {
                memcpy((*recv_routing_mapping_table_entries)+*num_recv_routing_mapping_table_entries, set_recv_routing_mapping_table_entries_buffer[i], sizeof(struct common_sort_struct<routing_mapping_table_entry>)*set_num_recv_routing_mapping_table_entries[i]);
                delete [] set_recv_routing_mapping_table_entries_buffer[i];
                *num_recv_routing_mapping_table_entries += set_num_recv_routing_mapping_table_entries[i];
            }
    }

    delete [] set_send_recv_mark;
    delete [] set_send_routing_mapping_table_entries_buffer;
    delete [] set_recv_routing_mapping_table_entries_buffer;
    delete [] set_num_recv_routing_mapping_table_entries;
    delete [] set_send_entry_count;
    delete [] set_remote_proc_global_id_send_to;
    delete [] set_remote_proc_global_id_recv_from;
}


void Routing_info::intersect_routing_mapping_tables_between_components(common_sort_struct<routing_mapping_table_entry> **local_routing_mapping_table_entries, common_sort_struct<routing_mapping_table_entry> *remote_routing_mapping_table_entries, int *num_local_routing_mapping_table_entries, int num_remote_routing_mapping_table_entries, bool is_src)
{
    std::vector<common_sort_struct<routing_mapping_table_entry> > mapping_routing_mapping_table_entries_pairs;
    int local_pointer = 0, remote_pointer = 0;

    if (*num_local_routing_mapping_table_entries == 0 || num_remote_routing_mapping_table_entries == 0) {
        if ((*local_routing_mapping_table_entries) != NULL)
            delete [] (*local_routing_mapping_table_entries);
        (*local_routing_mapping_table_entries) = NULL;
        *num_local_routing_mapping_table_entries = 0;
        return;
    }

    while (local_pointer < *num_local_routing_mapping_table_entries && remote_pointer < num_remote_routing_mapping_table_entries) {
        while((local_pointer < *num_local_routing_mapping_table_entries) && ((*local_routing_mapping_table_entries)[local_pointer].content.global_index < remote_routing_mapping_table_entries[remote_pointer].content.global_index))
            local_pointer ++;
        while ((remote_pointer < num_remote_routing_mapping_table_entries) && ((*local_routing_mapping_table_entries)[local_pointer].content.global_index > remote_routing_mapping_table_entries[remote_pointer].content.global_index))
            remote_pointer ++;
        
        if ((*local_routing_mapping_table_entries)[local_pointer].content.global_index == remote_routing_mapping_table_entries[remote_pointer].content.global_index) {
            int num_local_routing_mapping_table_entry_same = 1, num_remote_routing_mapping_table_entry_same = 1;
            while ((local_pointer+num_local_routing_mapping_table_entry_same < *num_local_routing_mapping_table_entries) && ((*local_routing_mapping_table_entries)[local_pointer].content.global_index == (*local_routing_mapping_table_entries)[local_pointer+num_local_routing_mapping_table_entry_same].content.global_index))
                num_local_routing_mapping_table_entry_same ++;
            while ((remote_pointer+num_remote_routing_mapping_table_entry_same < num_remote_routing_mapping_table_entries) && (remote_routing_mapping_table_entries[remote_pointer].content.global_index == remote_routing_mapping_table_entries[remote_pointer+num_remote_routing_mapping_table_entry_same].content.global_index))
                num_remote_routing_mapping_table_entry_same ++;
            
            if (is_src) {
                for (int j = 0; j < num_remote_routing_mapping_table_entry_same; j ++) {
                    common_sort_struct<routing_mapping_table_entry> temp_routing_mapping_table_entry;
					temp_routing_mapping_table_entry.key = -1;
					temp_routing_mapping_table_entry.content.global_index = (*local_routing_mapping_table_entries)[local_pointer + j%num_local_routing_mapping_table_entry_same].content.global_index;
					temp_routing_mapping_table_entry.content.local_index = (*local_routing_mapping_table_entries)[local_pointer + j%num_local_routing_mapping_table_entry_same].content.local_index;
					temp_routing_mapping_table_entry.content.local_process_id = (*local_routing_mapping_table_entries)[local_pointer + j%num_local_routing_mapping_table_entry_same].content.local_process_id;
					temp_routing_mapping_table_entry.content.remote_local_index = remote_routing_mapping_table_entries[remote_pointer + j].content.local_index;
					temp_routing_mapping_table_entry.content.remote_local_process_id = remote_routing_mapping_table_entries[remote_pointer + j].content.local_process_id;
                    mapping_routing_mapping_table_entries_pairs.push_back(temp_routing_mapping_table_entry);
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (*local_routing_mapping_table_entries)[local_pointer + j%num_local_routing_mapping_table_entry_same].content.global_index == remote_routing_mapping_table_entries[remote_pointer + j].content.global_index, "Software error in Routing_info::intersect_routing_mapping_tables_between_components");
					if (num_local_routing_mapping_table_entry_same > 1)
						EXECUTION_REPORT_LOG(REPORT_LOG, true, src_comp_node->get_comp_id(), "Router generation detects multiple input of the same cell (%d), and then generates a connection from src (%d %d) to dst (%d %d)", temp_routing_mapping_table_entry.content.global_index, temp_routing_mapping_table_entry.content.local_process_id, temp_routing_mapping_table_entry.content.local_index, temp_routing_mapping_table_entry.content.remote_local_process_id, temp_routing_mapping_table_entry.content.remote_local_index);
                }
            }
            else { 
                for (int j = 0; j < num_local_routing_mapping_table_entry_same; j ++) {
                    common_sort_struct<routing_mapping_table_entry> temp_routing_mapping_table_entry;
					temp_routing_mapping_table_entry.key = -1;
					temp_routing_mapping_table_entry.content.global_index = (*local_routing_mapping_table_entries)[local_pointer + j].content.global_index;
					temp_routing_mapping_table_entry.content.local_index = (*local_routing_mapping_table_entries)[local_pointer + j].content.local_index;
					temp_routing_mapping_table_entry.content.local_process_id = (*local_routing_mapping_table_entries)[local_pointer + j].content.local_process_id;
					temp_routing_mapping_table_entry.content.remote_local_index = remote_routing_mapping_table_entries[remote_pointer + j%num_remote_routing_mapping_table_entry_same].content.local_index;
					temp_routing_mapping_table_entry.content.remote_local_process_id = remote_routing_mapping_table_entries[remote_pointer + j%num_remote_routing_mapping_table_entry_same].content.local_process_id;
                    mapping_routing_mapping_table_entries_pairs.push_back(temp_routing_mapping_table_entry);
					EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (*local_routing_mapping_table_entries)[local_pointer + j].content.global_index == remote_routing_mapping_table_entries[remote_pointer + j%num_remote_routing_mapping_table_entry_same].content.global_index, "Software error in Routing_info::intersect_routing_mapping_tables_between_components");
                }
            }
            local_pointer += num_local_routing_mapping_table_entry_same;
            remote_pointer += num_remote_routing_mapping_table_entry_same;
        }
    }

    *num_local_routing_mapping_table_entries = mapping_routing_mapping_table_entries_pairs.size();
    if ((*local_routing_mapping_table_entries) != NULL)
        delete [] (*local_routing_mapping_table_entries);
    (*local_routing_mapping_table_entries) = NULL;
    if (*num_local_routing_mapping_table_entries > 0 ) {
        (*local_routing_mapping_table_entries) = new common_sort_struct<routing_mapping_table_entry> [*num_local_routing_mapping_table_entries];
        for (int i = 0; i < *num_local_routing_mapping_table_entries; i ++) {
            (*local_routing_mapping_table_entries)[i] = mapping_routing_mapping_table_entries_pairs[i];
        }
    }
}

void Routing_info::initialize_routing_mapping_table(common_sort_struct<routing_mapping_table_entry> **routing_mapping_table_entries, Decomp_info *decomp_info, Comp_comm_group_mgt_node *comp_node, int *num_local_routing_mapping_table_entries, int current_proc_id, int num_local_procs, int num_local_procs_adjust)
{
    int i = 0, j = 0;
    const int *local_cells_gobal_index = decomp_info->get_local_cell_global_indx();

    if (*num_local_routing_mapping_table_entries > 0) {
        (*routing_mapping_table_entries) = new common_sort_struct<routing_mapping_table_entry> [*num_local_routing_mapping_table_entries];
        for (int i = 0; i < *num_local_routing_mapping_table_entries; i ++) {
            if (local_cells_gobal_index[i] == CCPL_NULL_INT)
                continue;
            (*routing_mapping_table_entries)[j].content.global_index = local_cells_gobal_index[i];    
            (*routing_mapping_table_entries)[j].key = (*routing_mapping_table_entries)[j].content.global_index;
            (*routing_mapping_table_entries)[j].content.local_index = i;
            (*routing_mapping_table_entries)[j].content.local_process_id = current_proc_id;
            (*routing_mapping_table_entries)[j].content.remote_local_index = -1;
            (*routing_mapping_table_entries)[j].content.remote_local_process_id = -1;
            j++;
        }
        *num_local_routing_mapping_table_entries = j;
        do_quick_sort(*routing_mapping_table_entries, (int*)NULL, 0, *num_local_routing_mapping_table_entries-1);
    }
    //output_routing_mapping_table(*routing_mapping_table_entries, "origin", *num_local_routing_mapping_table_entries, current_proc_id);
    distribute_sorting->adjust_sorted_data_among_processes(routing_mapping_table_entries, decomp_info, comp_node, num_local_routing_mapping_table_entries, current_proc_id, num_local_procs, num_local_procs_adjust);
}


// compare_option is 1: based on global_index, 2: based on local_process_id
void Routing_info::sort_routing_mapping_table_entries(common_sort_struct<routing_mapping_table_entry>** routing_mapping_table_entries, Decomp_info* decomp_info, Comp_comm_group_mgt_node* comp_node, int* num_local_routing_mapping_table_entries, int current_proc_id, int num_local_procs, int num_local_procs_adjust, int num_remote_procs_adjust, int compare_option)
{    
    int temp_num_routing_mapping_table_entries;
    common_sort_struct<routing_mapping_table_entry> *temp_routing_mapping_table_entries = NULL;
    
    distribute_sorting->do_data_sorting(routing_mapping_table_entries, decomp_info, comp_node, num_local_routing_mapping_table_entries, current_proc_id, num_local_procs, num_local_procs_adjust, num_remote_procs_adjust, compare_option);

    if(compare_option == 2 && num_local_procs != num_local_procs_adjust) {
        
        int proc_id_send_to = current_proc_id + num_local_procs_adjust;
        int proc_id_recv_from = current_proc_id - num_local_procs_adjust;
        int send_recv_mark;
        if (current_proc_id < num_local_procs-num_local_procs_adjust) {
            send_recv_mark = ROUTER_SEND;
            distribute_sorting->exchange_sorted_data_between_processes(1, routing_mapping_table_entries, &temp_routing_mapping_table_entries, num_local_routing_mapping_table_entries, NULL, &proc_id_send_to, &proc_id_recv_from, comp_node->get_comm_group(), &send_recv_mark);
        } 
        else if (current_proc_id >= num_local_procs_adjust) {
            send_recv_mark = ROUTER_RECV;
            distribute_sorting->exchange_sorted_data_between_processes(1, &temp_routing_mapping_table_entries, routing_mapping_table_entries, NULL, num_local_routing_mapping_table_entries, &proc_id_send_to, &proc_id_recv_from, comp_node->get_comm_group(), &send_recv_mark);
        }       

        for (int i = 0; i < *num_local_routing_mapping_table_entries; i ++)
            (*routing_mapping_table_entries)[i].key = (*routing_mapping_table_entries)[i].content.local_process_id;
        distribute_sorting->merge_sorted_data(routing_mapping_table_entries, (common_sort_struct<routing_mapping_table_entry>*)NULL, *num_local_routing_mapping_table_entries, 0, current_proc_id, current_proc_id+1, num_local_routing_mapping_table_entries);
        
        if (temp_routing_mapping_table_entries != NULL)
            delete [] temp_routing_mapping_table_entries;
    }
}


void Routing_info::calculate_routing_mapping_tables()
{
    int num_src_procs = src_comp_node->get_num_procs();
    int num_src_procs_adjust = distribute_sorting->calculate_max_power2(num_src_procs);
    int num_dst_procs = dst_comp_node->get_num_procs();
    int num_dst_procs_adjust = distribute_sorting->calculate_max_power2(num_dst_procs);
    
    int num_src_global_routing_mapping_table_entries, num_dst_global_routing_mapping_table_entries;
 
    if (current_proc_id_src_comp != -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_decomp_info != NULL, "Software error in Routing_info::calculate_routing_mapping_tables: NULL src decomp info");
        num_src_local_routing_mapping_table_entries = src_decomp_info->get_num_local_cells();
        initialize_routing_mapping_table(&src_local_routing_mapping_table_entries, src_decomp_info, src_comp_node, &num_src_local_routing_mapping_table_entries, current_proc_id_src_comp, num_src_procs, num_src_procs_adjust);
        sort_routing_mapping_table_entries(&src_local_routing_mapping_table_entries, src_decomp_info, src_comp_node, &num_src_local_routing_mapping_table_entries, current_proc_id_src_comp, num_src_procs, num_src_procs_adjust, num_dst_procs_adjust, 1);
    }

    if (current_proc_id_dst_comp != -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_decomp_info != NULL, "Software error in Routing_info::calculate_routing_mapping_tables: NULL dst decomp info");
        num_dst_local_routing_mapping_table_entries = dst_decomp_info->get_num_local_cells();
        initialize_routing_mapping_table(&dst_local_routing_mapping_table_entries, dst_decomp_info, dst_comp_node, &num_dst_local_routing_mapping_table_entries, current_proc_id_dst_comp, num_dst_procs, num_dst_procs_adjust);        
        sort_routing_mapping_table_entries(&dst_local_routing_mapping_table_entries, dst_decomp_info, dst_comp_node, &num_dst_local_routing_mapping_table_entries, current_proc_id_dst_comp, num_dst_procs, num_dst_procs_adjust, num_src_procs_adjust, 1);
    }

    exchange_routing_mapping_tables_between_components(num_src_procs_adjust, num_dst_procs_adjust, src_decomp_info, true);        
    exchange_routing_mapping_tables_between_components(num_src_procs_adjust, num_dst_procs_adjust, dst_decomp_info, false);
    
    if (current_proc_id_src_comp != -1) 
        intersect_routing_mapping_tables_between_components(&src_local_routing_mapping_table_entries, dst_remote_routing_mapping_table_entries, &num_src_local_routing_mapping_table_entries, num_dst_remote_routing_mapping_table_entries, true);
    if (current_proc_id_dst_comp != -1)
        intersect_routing_mapping_tables_between_components(&dst_local_routing_mapping_table_entries, src_remote_routing_mapping_table_entries, &num_dst_local_routing_mapping_table_entries, num_src_remote_routing_mapping_table_entries, false);
    if (src_remote_routing_mapping_table_entries != NULL)
        delete [] src_remote_routing_mapping_table_entries;
    if (dst_remote_routing_mapping_table_entries != NULL)
        delete [] dst_remote_routing_mapping_table_entries;

    if (current_proc_id_src_comp != -1) {
        for (int i = 0; i < num_src_local_routing_mapping_table_entries; i ++)
            src_local_routing_mapping_table_entries[i].key = src_local_routing_mapping_table_entries[i].content.local_process_id % num_src_procs_adjust;
        do_quick_sort(src_local_routing_mapping_table_entries, (int*)NULL, 0, num_src_local_routing_mapping_table_entries-1);
        sort_routing_mapping_table_entries(&src_local_routing_mapping_table_entries, src_decomp_info, src_comp_node, &num_src_local_routing_mapping_table_entries, current_proc_id_src_comp, num_src_procs, num_src_procs_adjust, num_src_procs_adjust, 2);        
        for (int i = 0; i < num_src_local_routing_mapping_table_entries; i ++)
            src_local_routing_mapping_table_entries[i].key = src_local_routing_mapping_table_entries[i].content.remote_local_process_id;
        do_quick_sort(src_local_routing_mapping_table_entries, (int*)NULL, 0, num_src_local_routing_mapping_table_entries-1);
	}

    if (current_proc_id_dst_comp != -1) {
        for (int i = 0; i < num_dst_local_routing_mapping_table_entries; i ++)
            dst_local_routing_mapping_table_entries[i].key = dst_local_routing_mapping_table_entries[i].content.local_process_id % num_dst_procs_adjust;
        do_quick_sort(dst_local_routing_mapping_table_entries, (int*)NULL, 0, num_dst_local_routing_mapping_table_entries-1);
        sort_routing_mapping_table_entries(&dst_local_routing_mapping_table_entries, dst_decomp_info, dst_comp_node, &num_dst_local_routing_mapping_table_entries, current_proc_id_dst_comp, num_dst_procs, num_dst_procs_adjust, num_dst_procs_adjust, 2);
        for (int i = 0; i < num_dst_local_routing_mapping_table_entries; i ++)
            dst_local_routing_mapping_table_entries[i].key = dst_local_routing_mapping_table_entries[i].content.remote_local_process_id;
        do_quick_sort(dst_local_routing_mapping_table_entries, (int*)NULL, 0, num_dst_local_routing_mapping_table_entries-1);
    }
}


Routing_info::~Routing_info()
{
    for (int i = 0; i < recv_from_remote_procs_routing_info.size(); i ++) {
        if (recv_from_remote_procs_routing_info[i]->num_elements_transferred > 0) {
            delete [] recv_from_remote_procs_routing_info[i]->local_indx_segment_starts;
            delete [] recv_from_remote_procs_routing_info[i]->local_indx_segment_lengths;
        }
		delete recv_from_remote_procs_routing_info[i];
    }
    for (int i = 0; i < send_to_remote_procs_routing_info.size(); i ++) {
        if (send_to_remote_procs_routing_info[i]->num_elements_transferred > 0) {
            delete [] send_to_remote_procs_routing_info[i]->local_indx_segment_starts;
            delete [] send_to_remote_procs_routing_info[i]->local_indx_segment_lengths;
        }
		delete send_to_remote_procs_routing_info[i];
    }
}


bool Routing_info::match_router(const int src_comp_id, const int dst_comp_id, const char *src_decomp_name, const char *dst_decomp_name)
{
    Comp_comm_group_mgt_node *src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_comp_id);
	Decomp_info *src_decomp_info = decomps_info_mgr->search_decomp_info(src_decomp_name, src_comp_id);
	if (src_decomp_info != NULL)
		src_comp_node = comp_comm_group_mgt_mgr->search_global_node(src_decomp_info->get_host_comp_id());
    return (words_are_the_same(src_comp_full_name, src_comp_node->get_full_name()) && words_are_the_same(index_dst_comp_full_name, comp_comm_group_mgt_mgr->search_global_node(dst_comp_id)->get_full_name()) &&
            words_are_the_same(this->src_decomp_name, src_decomp_name) && words_are_the_same(this->dst_decomp_name, dst_decomp_name));
}


Routing_info_with_one_process *Routing_info::generate_routing_info_between_procs(common_sort_struct<routing_mapping_table_entry>* local_routing_mapping_table_entries, int num_local_routing_mapping_table_entries, int remote_proc_global_id, const int *local_cells_chunk_id, bool is_src)
{
    Routing_info_with_one_process *routing_info;
	int last_reference_local_index;

    routing_info = new Routing_info_with_one_process;
    routing_info->num_elements_transferred = num_local_routing_mapping_table_entries;
    routing_info->num_local_indx_segments = 0;
    routing_info->remote_proc_global_id = remote_proc_global_id;

	if (num_local_routing_mapping_table_entries == 0)
		return routing_info;

	if (report_error_enabled) 
		for (int i = 1; i < num_local_routing_mapping_table_entries; i ++)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, local_routing_mapping_table_entries[0].content.remote_local_process_id == local_routing_mapping_table_entries[i].content.remote_local_process_id, "Software error in Routing_info::generate_routing_info_between_procs");

    if (is_src)
        for (int i = 0; i < num_local_routing_mapping_table_entries; i ++)
            local_routing_mapping_table_entries[i].key = local_routing_mapping_table_entries[i].content.remote_local_index;
    else
        for (int i = 0; i < num_local_routing_mapping_table_entries; i ++)
            local_routing_mapping_table_entries[i].key = local_routing_mapping_table_entries[i].content.local_index;
    do_quick_sort(local_routing_mapping_table_entries, (int*)NULL, 0, num_local_routing_mapping_table_entries-1);

	last_reference_local_index = -100;
	for (int i = 0; i < num_local_routing_mapping_table_entries; i ++) {
		if (last_reference_local_index + 1 != local_routing_mapping_table_entries[i].content.local_index || local_cells_chunk_id != NULL && local_cells_chunk_id[last_reference_local_index] != local_cells_chunk_id[local_routing_mapping_table_entries[i].content.local_index])
            routing_info->num_local_indx_segments ++;
        last_reference_local_index = local_routing_mapping_table_entries[i].content.local_index;
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, routing_info->num_local_indx_segments > 0, "Software error in Routing_info::generate_routing_info_between_procs");
	routing_info->local_indx_segment_starts = new int [routing_info->num_local_indx_segments];
	routing_info->local_indx_segment_lengths = new int [routing_info->num_local_indx_segments];
	routing_info->num_local_indx_segments = 0;
	last_reference_local_index = -100;
	for (int i = 0; i < num_local_routing_mapping_table_entries; i ++) {
		if (last_reference_local_index + 1 != local_routing_mapping_table_entries[i].content.local_index || local_cells_chunk_id != NULL && local_cells_chunk_id[last_reference_local_index] != local_cells_chunk_id[local_routing_mapping_table_entries[i].content.local_index]) {
			routing_info->local_indx_segment_starts[routing_info->num_local_indx_segments] = local_routing_mapping_table_entries[i].content.local_index;
            routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments] = 1;
            routing_info->num_local_indx_segments ++;
		}
		else routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments - 1] ++;
        last_reference_local_index = local_routing_mapping_table_entries[i].content.local_index;
	}

    return routing_info;
}


void Routing_info::build_router_based_on_routing_mapping_tables()
{
    Routing_info_with_one_process *routing_info;

    
    if (current_proc_id_src_comp != -1) {
		int *remote_proc_num_routing_mapping_table_entries = new int [dst_comp_node->get_num_procs()];
		for (int i = 0; i < dst_comp_node->get_num_procs(); i ++)
			remote_proc_num_routing_mapping_table_entries[i] = 0;
        for (int i = 0; i < num_src_local_routing_mapping_table_entries; i ++)
			remote_proc_num_routing_mapping_table_entries[src_local_routing_mapping_table_entries[i].content.remote_local_process_id] ++;
		int offset = 0;
        for (int i = 0; i < dst_comp_node->get_num_procs(); i ++) {
            routing_info = generate_routing_info_between_procs(src_local_routing_mapping_table_entries+offset, remote_proc_num_routing_mapping_table_entries[i], dst_comp_node->get_local_proc_global_id(i), src_decomp_info->get_local_cell_chunk_id(), true);
            send_to_remote_procs_routing_info.push_back(routing_info);
			offset += remote_proc_num_routing_mapping_table_entries[i];
        }
		delete [] remote_proc_num_routing_mapping_table_entries;
    }

    if (current_proc_id_dst_comp != -1) {
		int *remote_proc_num_routing_mapping_table_entries = new int [src_comp_node->get_num_procs()];
		for (int i = 0; i < src_comp_node->get_num_procs(); i ++)
			remote_proc_num_routing_mapping_table_entries[i] = 0;
        for (int i = 0; i < num_dst_local_routing_mapping_table_entries; i ++)
			remote_proc_num_routing_mapping_table_entries[dst_local_routing_mapping_table_entries[i].content.remote_local_process_id] ++;
		int offset = 0;
        for (int i = 0; i < src_comp_node->get_num_procs(); i ++) {
            routing_info = generate_routing_info_between_procs(dst_local_routing_mapping_table_entries+offset, remote_proc_num_routing_mapping_table_entries[i], src_comp_node->get_local_proc_global_id(i), dst_decomp_info->get_local_cell_chunk_id(), false);
            recv_from_remote_procs_routing_info.push_back(routing_info);
			offset += remote_proc_num_routing_mapping_table_entries[i];
        }
		delete [] remote_proc_num_routing_mapping_table_entries;
    }

    if (src_local_routing_mapping_table_entries != NULL)
        delete [] src_local_routing_mapping_table_entries;
    if (dst_local_routing_mapping_table_entries != NULL)
        delete [] dst_local_routing_mapping_table_entries;
}


void Routing_info::build_2D_router()
{
    int num_src_procs = src_comp_node->get_num_procs();
    int *num_cells_each_src_proc = new int [num_src_procs];
    int num_dst_procs = dst_comp_node->get_num_procs();
    int * num_cells_each_dst_proc = new int [num_dst_procs];
    int num_local_src_cells, num_local_dst_cells, *num_global_src_cells = new int [1], *num_global_dst_cells = new int [1];
    int *cells_indx_each_src_proc = NULL;
    int *cells_indx_each_dst_proc = NULL;
    int src_comp_root_proc_global_id = src_comp_node->get_root_proc_global_id();
    int dst_comp_root_proc_global_id = dst_comp_node->get_root_proc_global_id();
    Routing_info_with_one_process *routing_info;
    long total_src_cells, total_dst_cells;

    if (current_proc_id_src_comp != -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_decomp_info != NULL, "Software error in Routing_info::build_2D_router: NULL src decomp info");
        num_local_src_cells = src_decomp_info->get_num_local_cells();
        *num_global_src_cells = src_decomp_info->get_num_global_cells();
        gather_array_in_one_comp(num_src_procs, current_proc_id_src_comp, (void*)src_decomp_info->get_local_cell_global_indx(), num_local_src_cells, 
                                 sizeof(int), num_cells_each_src_proc, (void**)(&cells_indx_each_src_proc), total_src_cells, src_comp_node->get_comm_group());
    }
    if (current_proc_id_dst_comp != -1) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, dst_decomp_info != NULL, "Software error in Routing_info::build_2D_router: NULL dst decomp info");
        num_local_dst_cells = dst_decomp_info->get_num_local_cells();
        *num_global_dst_cells = dst_decomp_info->get_num_global_cells();
        gather_array_in_one_comp(num_dst_procs, current_proc_id_dst_comp, (void*)dst_decomp_info->get_local_cell_global_indx(), num_local_dst_cells, 
                                 sizeof(int), num_cells_each_dst_proc, (void**)(&cells_indx_each_dst_proc), total_dst_cells, dst_comp_node->get_comm_group());
    }

    long temp_size = num_src_procs*sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), (char**)(&num_cells_each_src_proc), temp_size);
    temp_size = num_dst_procs*sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), (char**)(&num_cells_each_dst_proc), temp_size);
    temp_size = sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), (char**)(&num_global_src_cells), temp_size);    
    if (current_proc_id_dst_comp != -1)
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, *num_global_src_cells == *num_global_dst_cells, "Software error in Routing_info::build_2D_router: different global decomp grid size: %d vs %d", *num_global_src_cells, *num_global_dst_cells);
    total_src_cells = 0;
    for (int i = 0; i < num_src_procs; i ++) 
        total_src_cells += num_cells_each_src_proc[i] * sizeof(int);
    total_dst_cells = 0;
    for (int i = 0; i < num_dst_procs; i ++) 
        total_dst_cells += num_cells_each_dst_proc[i] * sizeof(int);
    transfer_array_from_one_comp_to_another(current_proc_id_src_comp, src_comp_root_proc_global_id, current_proc_id_dst_comp, dst_comp_root_proc_global_id, dst_comp_node->get_comm_group(), (char**)(&cells_indx_each_src_proc), total_src_cells);
    transfer_array_from_one_comp_to_another(current_proc_id_dst_comp, dst_comp_root_proc_global_id, current_proc_id_src_comp, src_comp_root_proc_global_id, src_comp_node->get_comm_group(), (char**)(&cells_indx_each_dst_proc), total_dst_cells);

	logical_indx_lookup_table_remote = new int [*num_global_src_cells];
	logical_indx_lookup_table_local = new int [*num_global_src_cells];
	for (int j = 0; j < *num_global_src_cells; j ++) {
			logical_indx_lookup_table_local[j] = -1;
			logical_indx_lookup_table_remote[j] = -1;
	}
	
    if (current_proc_id_src_comp != -1) {
        int *temp_num_cells_each_dst_proc_disp = new int [num_dst_procs];
        int *temp_num_cells_each_dst_proc_count = new int [num_dst_procs];
        for (int i = 0; i < num_dst_procs; i ++) {
            temp_num_cells_each_dst_proc_disp[i] = 0;
            temp_num_cells_each_dst_proc_count[i] = 0;
        }
        int last_dst_proc = -1;
        for (int i = 0; i < num_src_local_routing_mapping_table_entries; i ++) {
            if (src_local_routing_mapping_table_entries[i].content.remote_local_process_id != last_dst_proc) {
                last_dst_proc = src_local_routing_mapping_table_entries[i].content.remote_local_process_id;
                temp_num_cells_each_dst_proc_disp[last_dst_proc] = i; 
            }
            temp_num_cells_each_dst_proc_count[last_dst_proc] ++;
        }

        int tmp_displs = 0;
        for (int i = 0; i < num_dst_procs; i ++) {
            routing_info = compute_routing_info_between_decomps(num_local_src_cells, src_decomp_info->get_local_cell_global_indx(), src_decomp_info->get_local_cell_chunk_id(), num_cells_each_dst_proc[i], cells_indx_each_dst_proc+tmp_displs, 
                                                                src_decomp_info->get_num_global_cells(), comp_comm_group_mgt_mgr->get_current_proc_global_id(), dst_comp_node->get_local_proc_global_id(i), true,
                                                                src_local_routing_mapping_table_entries+temp_num_cells_each_dst_proc_disp[i], temp_num_cells_each_dst_proc_count[i]);
            tmp_displs += num_cells_each_dst_proc[i];
			if (routing_info != NULL)
	            send_to_remote_procs_routing_info.push_back(routing_info);
        }
		delete [] temp_num_cells_each_dst_proc_disp;
		delete [] temp_num_cells_each_dst_proc_count;
    }

    if (current_proc_id_dst_comp != -1) {
        int *temp_num_cells_each_src_proc_disp = new int [num_src_procs];
        int *temp_num_cells_each_src_proc_count = new int [num_src_procs];
        for (int i = 0; i < num_src_procs; i ++) {
            temp_num_cells_each_src_proc_disp[i] = 0;
            temp_num_cells_each_src_proc_count[i] = 0;
        }
        int last_src_proc = -1;
        for (int i = 0; i < num_dst_local_routing_mapping_table_entries; i ++) {
            if (dst_local_routing_mapping_table_entries[i].content.remote_local_process_id != last_src_proc) {
                last_src_proc = dst_local_routing_mapping_table_entries[i].content.remote_local_process_id;
                temp_num_cells_each_src_proc_disp[last_src_proc] = i; 
            }
            temp_num_cells_each_src_proc_count[last_src_proc] ++;
        }

        int tmp_displs = 0;
        for (int i = 0; i < num_src_procs; i ++) {
            routing_info = compute_routing_info_between_decomps(num_local_dst_cells, dst_decomp_info->get_local_cell_global_indx(), dst_decomp_info->get_local_cell_chunk_id(), num_cells_each_src_proc[i], cells_indx_each_src_proc+tmp_displs, 
                                                                dst_decomp_info->get_num_global_cells(), comp_comm_group_mgt_mgr->get_current_proc_global_id(), src_comp_node->get_local_proc_global_id(i), false,
                                                                dst_local_routing_mapping_table_entries+temp_num_cells_each_src_proc_disp[i], temp_num_cells_each_src_proc_count[i]);
            tmp_displs += num_cells_each_src_proc[i];
			if (routing_info != NULL)
	            recv_from_remote_procs_routing_info.push_back(routing_info);
        }
		delete [] temp_num_cells_each_src_proc_disp;
		delete [] temp_num_cells_each_src_proc_count;
    }

    delete [] logical_indx_lookup_table_remote;
    delete [] logical_indx_lookup_table_local;

    if (src_local_routing_mapping_table_entries != NULL)
        delete [] src_local_routing_mapping_table_entries;
    if (dst_local_routing_mapping_table_entries != NULL)
        delete [] dst_local_routing_mapping_table_entries;
    if (cells_indx_each_src_proc != NULL) 
        delete [] cells_indx_each_src_proc;
    if (cells_indx_each_dst_proc != NULL) 
        delete [] cells_indx_each_dst_proc;
    delete [] num_cells_each_src_proc; 
    delete [] num_cells_each_dst_proc;
    delete [] num_global_src_cells;
    delete [] num_global_dst_cells;
}


Routing_info_with_one_process *Routing_info::compute_routing_info_between_decomps(int num_local_cells_local, const int *local_cells_global_indexes_local, const int *local_cells_chunk_id_local,
                                                  int num_local_cells_remote, const int *local_cells_global_indexes_remote, 
                                                  int num_global_cells, int local_proc_id, int remote_proc_id, bool is_src,
                                                  common_sort_struct<routing_mapping_table_entry>* local_routing_mapping_table_entries, int routing_mapping_table_entries_count)
{
    Routing_info_with_one_process *routing_info;
    const int *reference_cell_indx;
    int num_reference_cells;
    int last_local_logical_indx;
    int j;

    routing_info = new Routing_info_with_one_process;
    routing_info->num_elements_transferred = 0;
    routing_info->num_local_indx_segments = 0;
    routing_info->remote_proc_global_id = remote_proc_id;

    if (num_local_cells_local == 0)
        return routing_info;
    
    /* Determine the reference cell index table according to the table size */
    if (is_src) {
        reference_cell_indx = local_cells_global_indexes_remote;
        num_reference_cells = num_local_cells_remote;  
    }
    else {
        reference_cell_indx = local_cells_global_indexes_local;
        num_reference_cells = num_local_cells_local; 
    }
    
    for (j = 0; j < num_local_cells_local; j ++)
        if (local_cells_global_indexes_local[j] >= 0)
            if (local_cells_global_indexes_local[j] != CCPL_NULL_INT)
                logical_indx_lookup_table_local[local_cells_global_indexes_local[j]] = j;
    for (j = 0; j < num_local_cells_remote; j ++)
        if (local_cells_global_indexes_remote[j] >= 0)
            if (local_cells_global_indexes_remote[j] != CCPL_NULL_INT)
                logical_indx_lookup_table_remote[local_cells_global_indexes_remote[j]] = j;

    if (is_src)
        for (int i = 0; i < routing_mapping_table_entries_count; i ++)
            local_routing_mapping_table_entries[i].key = local_routing_mapping_table_entries[i].content.remote_local_index;
    else
        for (int i = 0; i < routing_mapping_table_entries_count; i ++)
            local_routing_mapping_table_entries[i].key = local_routing_mapping_table_entries[i].content.local_index;
    do_quick_sort(local_routing_mapping_table_entries, (int*)NULL, 0, routing_mapping_table_entries_count-1);

    /* Compute the number of common cells and the number of segments of common cells */
    last_local_logical_indx = -100;
    for (j = 0; j < num_reference_cells; j ++) 
        if (reference_cell_indx[j] != CCPL_NULL_INT && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, local_routing_mapping_table_entries[routing_info->num_elements_transferred].content.global_index == reference_cell_indx[j] && j == local_routing_mapping_table_entries[routing_info->num_elements_transferred].key, "ERROR happens in Routing_info::compute_routing_info_between_decomps: different index between router info: %d %d, %d %d\n", local_routing_mapping_table_entries[routing_info->num_elements_transferred].content.global_index, local_routing_mapping_table_entries[routing_info->num_elements_transferred].key, reference_cell_indx[j], j);
			if (reference_cell_indx == local_cells_global_indexes_local) {
	            if (last_local_logical_indx + 1 != j || local_cells_chunk_id_local != NULL && local_cells_chunk_id_local[last_local_logical_indx] != local_cells_chunk_id_local[j]) 
	                routing_info->num_local_indx_segments ++;
	            last_local_logical_indx = j;
	            routing_info->num_elements_transferred ++;				
			}
			else {
	            if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]] || local_cells_chunk_id_local != NULL && local_cells_chunk_id_local[last_local_logical_indx] != local_cells_chunk_id_local[logical_indx_lookup_table_local[reference_cell_indx[j]]]) 
	                routing_info->num_local_indx_segments ++;
	            last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
	            routing_info->num_elements_transferred ++;
			}
        }
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, routing_info->num_elements_transferred == routing_mapping_table_entries_count, "ERROR happens in Routing_info::compute_routing_info_between_decomps: router_info->num_elements_transferred(%d) is different from cells_count(%d)\n", routing_info->num_elements_transferred, routing_mapping_table_entries_count);

    /* Compute the info of segments when there are common cells */
    last_local_logical_indx = -100;
    if (routing_info->num_elements_transferred > 0) {
        routing_info->local_indx_segment_starts = new int [routing_info->num_local_indx_segments];
        routing_info->local_indx_segment_lengths = new int [routing_info->num_local_indx_segments];
        routing_info->num_local_indx_segments = 0;
        for (j = 0; j < num_reference_cells; j ++) 
            if (reference_cell_indx[j] != CCPL_NULL_INT && logical_indx_lookup_table_local[reference_cell_indx[j]] != -1 && logical_indx_lookup_table_remote[reference_cell_indx[j]] != -1) {
				if (reference_cell_indx == local_cells_global_indexes_local) {
					if (last_local_logical_indx + 1 != j || local_cells_chunk_id_local != NULL && local_cells_chunk_id_local[last_local_logical_indx] != local_cells_chunk_id_local[j]) {
	                    routing_info->local_indx_segment_starts[routing_info->num_local_indx_segments] = j;
	                    routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments] = 1;
	                    routing_info->num_local_indx_segments ++;
	                }
	                else routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments - 1] ++;
	                last_local_logical_indx = j;
				}
				else {
					if (last_local_logical_indx + 1 != logical_indx_lookup_table_local[reference_cell_indx[j]] || local_cells_chunk_id_local != NULL && local_cells_chunk_id_local[last_local_logical_indx] != local_cells_chunk_id_local[logical_indx_lookup_table_local[reference_cell_indx[j]]]) {
	                    routing_info->local_indx_segment_starts[routing_info->num_local_indx_segments] = logical_indx_lookup_table_local[reference_cell_indx[j]];
	                    routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments] = 1;
	                    routing_info->num_local_indx_segments ++;
	                }
	                else routing_info->local_indx_segment_lengths[routing_info->num_local_indx_segments - 1] ++;
	                last_local_logical_indx = logical_indx_lookup_table_local[reference_cell_indx[j]];
				}
            }
    }


    for (j = 0; j < num_local_cells_local; j ++)
        if (local_cells_global_indexes_local[j] >= 0 && local_cells_global_indexes_local[j] != CCPL_NULL_INT)
                        logical_indx_lookup_table_local[local_cells_global_indexes_local[j]] = -1;
    for (j = 0; j < num_local_cells_remote; j ++)
        if (local_cells_global_indexes_remote[j] >= 0 && local_cells_global_indexes_remote[j] != CCPL_NULL_INT)
            logical_indx_lookup_table_remote[local_cells_global_indexes_remote[j]] = -1;

    return routing_info;
}


Routing_info_with_one_process *Routing_info::get_routing_info(bool is_send, int i)
{
    if (is_send) {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, i >= 0 && i < send_to_remote_procs_routing_info.size(), "Software error in Routing_info::get_num_elements_transferred_with_remote_proc: wrong i at sender");
        return send_to_remote_procs_routing_info[i];
    }
    else {
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, i >= 0 && i < recv_from_remote_procs_routing_info.size(), "Software error in Routing_info::get_num_elements_transferred_with_remote_proc: wrong i at receiver");
        return recv_from_remote_procs_routing_info[i];
    }
}


Comp_comm_group_mgt_node *Routing_info::get_src_comp_node() 
{ 
    return comp_comm_group_mgt_mgr->search_global_node(src_comp_full_name); 
}


Comp_comm_group_mgt_node *Routing_info::get_dst_comp_node() 
{ 
    return comp_comm_group_mgt_mgr->search_global_node(true_dst_comp_full_name); 
}
