/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Mr. Hao Yu and Dr. Li 
  *  Liu and then modified by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn 
  ***************************************************************/

#ifndef DISTRIBUTED_MERGE_SORT
#define DISTRIBUTED_MERGE_SORT


#include "compset_communicators_info_mgt.h"
#include "execution_report.h"
#include <mpi.h>
#include <math.h>


#define ROUTER_SEND 1
#define ROUTER_RECV 2

template <typename T>
struct common_sort_struct
{
    long key;
	int target_proc_id;
    T content;
    bool operator<(const common_sort_struct<T> b) const { return this->key < b.key; }
    bool operator>(const common_sort_struct<T> b) const { return this->key > b.key; }   
    bool operator<=(const common_sort_struct<T> b) const { return this->key <= b.key; }
    bool operator>=(const common_sort_struct<T> b) const { return this->key >= b.key; }   
    bool operator==(const common_sort_struct<T> b) const { return this->key == b.key; }
};


template<class T>
class Distribute_merge_sort 
{
    private:
        int current_proc_id_src_comp;
        int current_proc_id_dst_comp;
        Comp_comm_group_mgt_node * src_comp_node;
        Comp_comm_group_mgt_node * dst_comp_node;        
    public:
        Distribute_merge_sort(int, int, Comp_comm_group_mgt_node *, Comp_comm_group_mgt_node *);
        void exchange_sorted_data_between_processes(int, common_sort_struct<T> **, common_sort_struct<T> **, int *, int *, int *, int *, MPI_Comm comm, int *);
        void merge_sorted_data(common_sort_struct<T>**, common_sort_struct<T>*, int, int, int, int, int*);        
        void check_data_sorting_order(common_sort_struct<T>*, int, int, int, char*);
        int calculate_max_power2(int);
        void calculate_min_max_grid_index(int, int, int, int, int, int*, int*); 

        // re-write
        void do_data_sorting(common_sort_struct<T>**, Decomp_info*, Comp_comm_group_mgt_node*, int*, int, int, int, int, int);
        void adjust_sorted_data_among_processes(common_sort_struct<T>**, Decomp_info*, Comp_comm_group_mgt_node*, int*, int, int, int);
		void do_data_sorting_with_target_process_id(common_sort_struct<T> **, Decomp_info *, Comp_comm_group_mgt_node *, int *);
};


template<typename T>
Distribute_merge_sort<T>::Distribute_merge_sort(int current_proc_id_src_comp, int current_proc_id_dst_comp, Comp_comm_group_mgt_node * src_comp_node, Comp_comm_group_mgt_node * dst_comp_node)
{
    this->current_proc_id_src_comp = current_proc_id_src_comp;
    this->current_proc_id_dst_comp = current_proc_id_dst_comp;
    this->src_comp_node = src_comp_node;
    this->dst_comp_node = dst_comp_node;
}


template<typename T>
int Distribute_merge_sort<T>::calculate_max_power2(int origin_num_procs)
{
    int new_num_procs = origin_num_procs;
    if (origin_num_procs > 1 && (origin_num_procs & (origin_num_procs-1)) != 0)
        new_num_procs = pow(2, (int)log2(origin_num_procs));
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, new_num_procs > 0 && new_num_procs <= origin_num_procs && 2*new_num_procs > origin_num_procs, "Software error in Distribute_merge_sort::calculate_max_power2");
        
    return new_num_procs;
}


template<typename T>
void Distribute_merge_sort<T>::merge_sorted_data(common_sort_struct<T> **local_sorted_data, common_sort_struct<T> *remote_sorted_data, int num_local_sorted_data, int num_remote_sorted_data, int min_key_value, int max_key_value, int *num_return_sorted_data)
{
    std::vector<common_sort_struct<T> > temp_sorted_data;
    int local_data_pointer = 0;
    int remote_data_pointer = 0;

    while (local_data_pointer < num_local_sorted_data && remote_data_pointer < num_remote_sorted_data) {
        if ((*local_sorted_data)[local_data_pointer].key <= remote_sorted_data[remote_data_pointer].key) {
            if ((*local_sorted_data)[local_data_pointer].key >= min_key_value && (*local_sorted_data)[local_data_pointer].key < max_key_value)
                temp_sorted_data.push_back((*local_sorted_data)[local_data_pointer]);
            local_data_pointer ++;
        }
        else {
            if (remote_sorted_data[remote_data_pointer].key >= min_key_value && remote_sorted_data[remote_data_pointer].key < max_key_value)
                temp_sorted_data.push_back(remote_sorted_data[remote_data_pointer]);
            remote_data_pointer++;
        }
    }
    for (; local_data_pointer < num_local_sorted_data; local_data_pointer ++)
        if ( (*local_sorted_data)[local_data_pointer].key >= min_key_value && (*local_sorted_data)[local_data_pointer].key < max_key_value)
            temp_sorted_data.push_back((*local_sorted_data)[local_data_pointer]);
    for (; remote_data_pointer < num_remote_sorted_data; remote_data_pointer ++)
        if (remote_sorted_data[remote_data_pointer].key >= min_key_value && remote_sorted_data[remote_data_pointer].key < max_key_value)
            temp_sorted_data.push_back(remote_sorted_data[remote_data_pointer]);
    
    *num_return_sorted_data = temp_sorted_data.size();
    if ((*local_sorted_data) != NULL)
        delete [] (*local_sorted_data);
    (*local_sorted_data) = NULL;
    if (*num_return_sorted_data > 0) {
        (*local_sorted_data) = new common_sort_struct<T> [*num_return_sorted_data];
        for (int i = 0; i < *num_return_sorted_data; i ++)
            (*local_sorted_data)[i] = temp_sorted_data[i];
    }
}


template<typename T>
void Distribute_merge_sort<T>::exchange_sorted_data_between_processes(int set_size, common_sort_struct<T> **set_send_sorted_data, common_sort_struct<T> **set_recv_sorted_data, int *set_num_send_sorted_data, int *set_num_recv_sorted_data, int *set_proc_id_send_to, int *set_proc_id_recv_from, MPI_Comm comm, int *set_send_recv_mark)
{
    MPI_Status status;
    MPI_Request *mpi_requests;
    int num_mpi_requests = 0;

    mpi_requests = new MPI_Request [2*set_size];
    for (int i = 0; i < set_size; i ++) {
        if ((set_send_recv_mark[i] & ROUTER_SEND) > 0)
            MPI_Isend(set_num_send_sorted_data+i, 1, MPI_INT, set_proc_id_send_to[i], 0, comm, &mpi_requests[num_mpi_requests++]);
        if ((set_send_recv_mark[i] & ROUTER_RECV) > 0)
            MPI_Irecv(set_num_recv_sorted_data+i, 1, MPI_INT, set_proc_id_recv_from[i], 0, comm, &mpi_requests[num_mpi_requests++]);
    }
    for (int i = 0; i < num_mpi_requests; i ++)
        MPI_Wait(&mpi_requests[i], &status);

    num_mpi_requests = 0;
    for (int i = 0; i < set_size; i ++) {
        if ((set_send_recv_mark[i] & ROUTER_SEND) > 0 && set_num_send_sorted_data[i] > 0)
            MPI_Isend((char*)(set_send_sorted_data[i]), sizeof(struct common_sort_struct<T>)*set_num_send_sorted_data[i], MPI_CHAR, set_proc_id_send_to[i], 0, comm, &mpi_requests[num_mpi_requests++]);
		if ((set_send_recv_mark[i] & ROUTER_RECV) > 0) {
			if (set_recv_sorted_data[i] != NULL)
				delete [] set_recv_sorted_data[i];
			set_recv_sorted_data[i] = NULL;
		}
        if ((set_send_recv_mark[i] & ROUTER_RECV) > 0 && set_num_recv_sorted_data[i] > 0) {
            set_recv_sorted_data[i] = new common_sort_struct<T> [set_num_recv_sorted_data[i]];
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, sizeof(struct common_sort_struct<T>) == (char*)(set_recv_sorted_data[i]+1)-(char*)(set_recv_sorted_data[i]), "Software error in Distribute_merge_sort<T>::exchange_sorted_data_between_processes");
            MPI_Irecv((char*)(set_recv_sorted_data[i]), sizeof(struct common_sort_struct<T>)*(set_num_recv_sorted_data[i]), MPI_CHAR, set_proc_id_recv_from[i], 0, comm, &mpi_requests[num_mpi_requests++]);
        }
    }
    for (int i = 0; i < num_mpi_requests; i ++)
        MPI_Wait(&mpi_requests[i], &status);

    delete [] mpi_requests;
}


template<typename T>
void Distribute_merge_sort<T>::calculate_min_max_grid_index(int num_grid_cells, int min_proc_id, int max_proc_id, int ref_proc_num, int proc_scale_factor, int* min_grid_index, int* max_grid_index)
{
    *min_grid_index = (num_grid_cells / ref_proc_num) * (min_proc_id * proc_scale_factor);
    *min_grid_index += std::min(num_grid_cells % ref_proc_num, min_proc_id * proc_scale_factor);
    *max_grid_index = (num_grid_cells / ref_proc_num) * (max_proc_id + 1) * proc_scale_factor;
    *max_grid_index += std::min(num_grid_cells % ref_proc_num, (max_proc_id + 1) * proc_scale_factor);
}


template<typename T>
void Distribute_merge_sort<T>::check_data_sorting_order(common_sort_struct<T> *sorted_data, int num_local_sorted_data, int min_key_value, int max_key_value, char* hint)
{
    if (!report_error_enabled)
        return;
    
    for (int i = 0; i < num_local_sorted_data - 1; i ++)
        EXECUTION_REPORT(REPORT_ERROR, -1, sorted_data[i].key >= min_key_value && sorted_data[i].key < max_key_value && sorted_data[i+1].key >= min_key_value && sorted_data[i+1].key < max_key_value && sorted_data[i].key <= sorted_data[i+1].key, "ERROR in check_data_sorting_order (%s): %d th is %ld, while the range is %d %d, next value key is %d\n", hint, i, sorted_data[i].key, min_key_value, max_key_value, sorted_data[i+1].key);
}


template<typename T>
void Distribute_merge_sort<T>::adjust_sorted_data_among_processes(common_sort_struct<T> **sorted_data, Decomp_info *decomp_info, Comp_comm_group_mgt_node *comp_node, int *num_local_sorted_data, int current_proc_id, int num_local_procs, int num_local_procs_adjust)
{
    common_sort_struct<T> *temp_sorted_data = NULL;
    int temp_num_sorted_data;
    int send_recv_mark, proc_id_send_to, proc_id_recv_from;
    if (num_local_procs != num_local_procs_adjust) {
        proc_id_send_to = current_proc_id-num_local_procs_adjust;
        proc_id_recv_from = current_proc_id+num_local_procs_adjust;
        if (current_proc_id >= num_local_procs_adjust) {
            send_recv_mark = ROUTER_SEND;
            exchange_sorted_data_between_processes(1, sorted_data, &temp_sorted_data, num_local_sorted_data, &temp_num_sorted_data, &proc_id_send_to, &proc_id_recv_from, comp_node->get_comm_group(), &send_recv_mark);
            *num_local_sorted_data = 0; // proc_id > 2^n, so the sorted_data number is 0
        } 
        else if (current_proc_id <= (num_local_procs-num_local_procs_adjust-1)) {
            send_recv_mark = ROUTER_RECV;
            exchange_sorted_data_between_processes(1, sorted_data, &temp_sorted_data, num_local_sorted_data, &temp_num_sorted_data, &proc_id_send_to, &proc_id_recv_from, comp_node->get_comm_group(), &send_recv_mark);
            merge_sorted_data(sorted_data, temp_sorted_data, *num_local_sorted_data, temp_num_sorted_data, 0, decomp_info->get_num_global_cells(), num_local_sorted_data);
            if (temp_sorted_data != NULL)
                delete [] temp_sorted_data;
        }
    }
}


template<typename T>
void Distribute_merge_sort<T>::do_data_sorting_with_target_process_id(common_sort_struct<T> **sorted_data, Decomp_info *decomp_info, Comp_comm_group_mgt_node *comp_node, int *num_local_sorted_data)
{
	int send_recv_mark, proc_id_send_to, proc_id_recv_from, num_procs_in_sorting = calculate_max_power2(comp_node->get_num_procs());
	int original_num_total_data, new_num_total_data;

	if (report_error_enabled) 
		MPI_Allreduce(num_local_sorted_data, &original_num_total_data, 1, MPI_INT, MPI_SUM, comp_node->get_comm_group());
	for (int i = 0; i < *num_local_sorted_data; i ++)
		(*sorted_data)[i].key = (*sorted_data)[i].target_proc_id % num_procs_in_sorting;
	adjust_sorted_data_among_processes(sorted_data, decomp_info, comp_node, num_local_sorted_data, comp_node->get_current_proc_local_id(), comp_node->get_num_procs(), num_procs_in_sorting);
	do_quick_sort(*sorted_data, (int*)NULL, 0, *num_local_sorted_data-1);
	do_data_sorting(sorted_data, decomp_info, comp_node, num_local_sorted_data, comp_node->get_current_proc_local_id(), comp_node->get_num_procs(), num_procs_in_sorting, num_procs_in_sorting, 2);
	if (comp_node->get_current_proc_local_id() < comp_node->get_num_procs()-num_procs_in_sorting) {
		send_recv_mark = ROUTER_SEND;
		proc_id_send_to = comp_node->get_current_proc_local_id() + num_procs_in_sorting;
		exchange_sorted_data_between_processes(1, sorted_data, sorted_data, num_local_sorted_data, NULL, &proc_id_send_to, &proc_id_recv_from, comp_node->get_comm_group(), &send_recv_mark);
	} 
	else if (comp_node->get_current_proc_local_id() >= num_procs_in_sorting) {
		send_recv_mark = ROUTER_RECV;
		proc_id_recv_from = comp_node->get_current_proc_local_id() - num_procs_in_sorting;
		exchange_sorted_data_between_processes(1, sorted_data, sorted_data, NULL, num_local_sorted_data, &proc_id_send_to, &proc_id_recv_from, comp_node->get_comm_group(), &send_recv_mark);
	}
	for (int i = 0; i < *num_local_sorted_data; i ++)
		(*sorted_data)[i].key = (*sorted_data)[i].target_proc_id;
	merge_sorted_data(sorted_data, (common_sort_struct<T>*)NULL, *num_local_sorted_data, 0, comp_node->get_current_proc_local_id(), comp_node->get_current_proc_local_id()+1, num_local_sorted_data); 

	if (report_error_enabled) {
		MPI_Allreduce(num_local_sorted_data, &new_num_total_data, 1, MPI_INT, MPI_SUM, comp_node->get_comm_group());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_num_total_data == new_num_total_data, "Software error in do_data_sorting_with_target_process_id: %d vs %d", original_num_total_data, new_num_total_data);
	}
}


template<typename T>
void Distribute_merge_sort<T>::do_data_sorting(common_sort_struct<T> **sorted_data, Decomp_info *decomp_info, Comp_comm_group_mgt_node *comp_node, int *num_local_sorted_data, int current_proc_id, int num_local_procs, int num_local_procs_adjust, int num_remote_procs_adjust, int compare_option)
{
    int temp_num_sorted_data;
    common_sort_struct<T> *temp_sorted_data = NULL;

    int ref_proc_num = std::max(num_local_procs_adjust, num_remote_procs_adjust);
    int proc_scale_factor = ref_proc_num / num_local_procs_adjust;
    int min_proc_id = 0, max_proc_id = num_local_procs_adjust - 1, mid_proc_id, remote_proc_id;
    int min_key_value, mid_key_value, max_key_value, send_recv_mark;
    int proc_id_send_to, proc_id_recv_from;

    if (current_proc_id < num_local_procs_adjust) {
        while (min_proc_id != max_proc_id) {
            mid_proc_id = (max_proc_id - min_proc_id + 1) / 2 + min_proc_id - 1;            
            if (compare_option == 1 ) {
                calculate_min_max_grid_index(decomp_info->get_num_global_cells(), min_proc_id, max_proc_id, ref_proc_num, proc_scale_factor, &min_key_value, &max_key_value);
                mid_key_value = (decomp_info->get_num_global_cells() / ref_proc_num) * (mid_proc_id + 1) * proc_scale_factor;
                mid_key_value += std::min(decomp_info->get_num_global_cells() % ref_proc_num, (mid_proc_id + 1) * proc_scale_factor);
            }
            else if (compare_option == 2) {
                min_key_value = min_proc_id;
                mid_key_value = mid_proc_id + 1;
                max_key_value = max_proc_id + 1;
            }

            send_recv_mark = ROUTER_SEND|ROUTER_RECV;
            if (current_proc_id <= mid_proc_id ) {
                remote_proc_id = current_proc_id + (max_proc_id - min_proc_id + 1) / 2;
                exchange_sorted_data_between_processes(1, sorted_data, &temp_sorted_data, num_local_sorted_data, &temp_num_sorted_data, &remote_proc_id, &remote_proc_id, comp_node->get_comm_group(), &send_recv_mark);
                merge_sorted_data(sorted_data, temp_sorted_data, *num_local_sorted_data, temp_num_sorted_data, min_key_value, mid_key_value, num_local_sorted_data); 
                max_proc_id = mid_proc_id;
            } else {
                remote_proc_id = current_proc_id - (max_proc_id - min_proc_id +  1) / 2;
                exchange_sorted_data_between_processes(1, sorted_data, &temp_sorted_data, num_local_sorted_data, &temp_num_sorted_data, &remote_proc_id, &remote_proc_id, comp_node->get_comm_group(), &send_recv_mark);
                merge_sorted_data(sorted_data, temp_sorted_data, *num_local_sorted_data, temp_num_sorted_data, mid_key_value, max_key_value, num_local_sorted_data);
                min_proc_id = mid_proc_id + 1;
            }
            check_data_sorting_order(*sorted_data, *num_local_sorted_data, min_key_value, max_key_value, "in do_data_sorting"); 

            if (temp_sorted_data != NULL) {
                delete [] temp_sorted_data;
				temp_sorted_data = NULL;
            }
        }
        EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, min_proc_id == current_proc_id && current_proc_id == max_proc_id, "Software error in Distribute_merge_sort::do_data_sorting");
    }
}


#endif

