/**********************************************M*****************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Miss Xinzhu Yu and Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#include "global_data.h"
#include <dirent.h>
#include "field_info_mgt.h"
#include <netcdf.h>
#ifdef USE_PARALLEL_IO
#include <pnetcdf.h>
#endif
#include "quick_sort.h"
#include "datamodel_mgt.h"

#define rename_datamodel_grid(grid_name_str2,datamodel_name, grid_name_str) sprintf(grid_name_str2, "datamodel-%s-grid-%s", datamodel_name, grid_name_str)
extern int elapsed_days_on_start_of_month_of_nonleap_year[];
extern int elapsed_days_on_start_of_month_of_leap_year[];
extern int num_days_of_month_of_nonleap_year[];
extern int num_days_of_month_of_leap_year[];

bool is_string_a_value(const char* string) 
{
	int num_dots = 0;

    for (int i = 0; i < strlen(string); i ++) {
        if (!(((string[i] >= '0') && (string[i] <= '9')) || (string[i] == '.' && num_dots == 0) || (string[i] == '-' && i == 0))) 
            return false;//not decimal value
		if (string[i] == '.') 
			num_dots ++;
	}

    return true;//decimal value
}


bool is_string_a_nonnegative_integer(const char* string) 
{
	if (!is_string_a_value(string))
		return false;

	for (int i = 0; i < strlen(string); i ++)
		if (string[i] == '.' || string[0] == '-')
			return false;

	return true;
}


char *characters_to_lower_case(const char *string) 
{
	if (string == NULL) 
		return NULL;

	char *new_str = strdup(string);
	int str_len = strlen(new_str);
	for (int i = 0; i < str_len; i ++)
		if (new_str[i]<='Z' && new_str[i]>='A')
			new_str[i]+=32;

	return new_str;
}


bool words_are_similar(const char *str1, const char *str2) 
{
	char *word1 = characters_to_lower_case(str1);
	char *word2 = characters_to_lower_case(str2);
	bool result;

	if (words_are_the_same(word1, word2))
		result = true;
	else result = false;

	delete [] word1;
	delete [] word2;
	return result;
}


int set_unit(const char* input_unit, const char *report_hint) 
{
    int output_unit;

    if (words_are_the_same(input_unit,"year") || words_are_the_same(input_unit,"years") || words_are_the_same(input_unit,"nyear") || words_are_the_same(input_unit,"nyears"))
        output_unit = 1;
    else if (words_are_the_same(input_unit,"month") || words_are_the_same(input_unit,"months") || words_are_the_same(input_unit,"nmonth") || words_are_the_same(input_unit,"nmonths"))
        output_unit = 2;
    else if (words_are_the_same(input_unit,"day") || words_are_the_same(input_unit,"days") || words_are_the_same(input_unit,"nday") || words_are_the_same(input_unit,"ndays"))
        output_unit = 3;
    else if (words_are_the_same(input_unit,"second") || words_are_the_same(input_unit,"seconds") || words_are_the_same(input_unit,"nsecond") || words_are_the_same(input_unit,"nseconds"))
        output_unit = 4;
    else if (words_are_the_same(input_unit,"step") || words_are_the_same(input_unit,"steps") || words_are_the_same(input_unit,"nstep") || words_are_the_same(input_unit,"nsteps"))
        output_unit = 5;
    else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "the unit \"%s\" for \"%s\" is incorrect, Please check.", input_unit, report_hint);

    return output_unit;
}


int check_time_format(const char *time_format, const char *report_hint) 
{
    int time_format_ID;

    if (words_are_the_same(time_format, "YYYY"))
        time_format_ID = TIME_FORMAT_YYYY;
    else if (words_are_the_same(time_format, "YYYYMM") || words_are_the_same(time_format, "YYYY-MM"))
        time_format_ID = TIME_FORMAT_YYYYMM;
    else if (words_are_the_same(time_format, "YYYYMMDD") || words_are_the_same(time_format, "YYYY-MM-DD"))
        time_format_ID = TIME_FORMAT_YYYYMMDD;
    else if (words_are_the_same(time_format, "YYYYMMDD.SSSSS") || words_are_the_same(time_format, "YYYY-MM-DD.SSSSS") || 
             words_are_the_same(time_format, "YYYYMMDD-SSSSS") || words_are_the_same(time_format, "YYYY-MM-DD-SSSSS") ||
             words_are_the_same(time_format, "YYYYMMDDSSSSS"))
        time_format_ID = TIME_FORMAT_YYYYMMDDSSSSS;
    else if (words_are_the_same(time_format, "YYYYMMDD.HHMMSS") || words_are_the_same(time_format, "YYYY-MM-DD.HH-MM-SS") ||
             words_are_the_same(time_format, "YYYYMMDD.HH-MM-SS") || words_are_the_same(time_format, "YYYY-MM-DD.HHMMSS") ||
             words_are_the_same(time_format, "YYYY-MM-DD-HH-MM-SS") || words_are_the_same(time_format, "YYYYMMDDHHMMSS"))
        time_format_ID = TIME_FORMAT_YYYYMMDDHHMMSS;
    else if (words_are_the_same(time_format, "SSSSS"))
        time_format_ID = TIME_FORMAT_SSSSS;
    else if (words_are_the_same(time_format, "HHMMSS") || words_are_the_same(time_format, "HH-MM-SS") || words_are_the_same(time_format, "HH:MM:SS"))
        time_format_ID = TIME_FORMAT_HHMMSS;
    else if (words_are_the_same(time_format, "MMDD.HHMMSS") || words_are_the_same(time_format, "MMDD-HHMMSS") || words_are_the_same(time_format, "MM-DD.HH-MM-SS") ||
             words_are_the_same(time_format, "MM-DD-HH-MM-SS") || words_are_the_same(time_format, "MMDDHHMMSS"))
        time_format_ID = TIME_FORMAT_MMDDHHMMSS;
    else if (words_are_the_same(time_format, "MMDD.SSSSS") || words_are_the_same(time_format, "MMDD-SSSSS") || words_are_the_same(time_format, "MM-DD.SSSSS") ||
             words_are_the_same(time_format, "MM-DD-SSSSS") || words_are_the_same(time_format, "MMDDSSSSS"))
        time_format_ID = TIME_FORMAT_MMDDSSSSS;
    else if (words_are_the_same(time_format, "DD.HHMMSS") || words_are_the_same(time_format, "DD-HHMMSS") || words_are_the_same(time_format, "DD.HH-MM-SS") ||
             words_are_the_same(time_format, "DD-HH-MM-SS") || words_are_the_same(time_format, "DDHHMMSS"))
        time_format_ID = TIME_FORMAT_DDHHMMSS;
    else if (words_are_the_same(time_format, "DD.SSSSS") || words_are_the_same(time_format, "DD-SSSSS") || words_are_the_same(time_format, "DDSSSSS"))
        time_format_ID = TIME_FORMAT_DDSSSSS;
    else if (words_are_the_same(time_format, "YYYYMMDD.HHMM") || words_are_the_same(time_format, "YYYY-MM-DD.HH-MM") ||
             words_are_the_same(time_format, "YYYYMMDD.HH-MM") || words_are_the_same(time_format, "YYYY-MM-DD.HHMM") ||
             words_are_the_same(time_format, "YYYY-MM-DD-HH-MM") || words_are_the_same(time_format, "YYYYMMDDHHMM"))
        time_format_ID = TIME_FORMAT_YYYYMMDDHHMM;
    else if (words_are_the_same(time_format, "HHMM") || words_are_the_same(time_format, "HH-MM") || words_are_the_same(time_format, "HH:MM"))
        time_format_ID = TIME_FORMAT_HHMM;
    else if (words_are_the_same(time_format, "MMDD.HHMM") || words_are_the_same(time_format, "MMDD-HHMM") || words_are_the_same(time_format, "MM-DD.HH-MM") ||
             words_are_the_same(time_format, "MM-DD-HH-MM") || words_are_the_same(time_format, "MMDDHHMM"))
        time_format_ID = TIME_FORMAT_MMDDHHMM;
    else if (words_are_the_same(time_format, "DD.HHMM") || words_are_the_same(time_format, "DD-HHMM") || words_are_the_same(time_format, "DD.HH-MM") ||
             words_are_the_same(time_format, "DD-HH-MM") || words_are_the_same(time_format, "DDHHMM"))
        time_format_ID = TIME_FORMAT_DDHHMM;
    else if (words_are_the_same(time_format, "YYYYMMDD.HH") || words_are_the_same(time_format, "YYYY-MM-DD.HH") ||
             words_are_the_same(time_format, "YYYYMMDD.HH") || words_are_the_same(time_format, "YYYY-MM-DD.HH") ||
             words_are_the_same(time_format, "YYYY-MM-DD-HH") || words_are_the_same(time_format, "YYYYMMDDHH"))
        time_format_ID = TIME_FORMAT_YYYYMMDDHH;
    else if (words_are_the_same(time_format, "HH"))
        time_format_ID = TIME_FORMAT_HH;
    else if (words_are_the_same(time_format, "MMDD.HH") || words_are_the_same(time_format, "MMDD-HH") || words_are_the_same(time_format, "MM-DD.HH") ||
             words_are_the_same(time_format, "MM-DD-HH") || words_are_the_same(time_format, "MMDDHH"))
        time_format_ID = TIME_FORMAT_MMDDHH;
    else if (words_are_the_same(time_format, "DD.HH") || words_are_the_same(time_format, "DD-HH") || words_are_the_same(time_format, "DDHH"))
        time_format_ID = TIME_FORMAT_DDHH;
    else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "The time format \"%s\" in the XML configuration \"%s\" is incorrect. Please verify", time_format, report_hint);//测一下

    return time_format_ID;
}

int get_time_format_smallest_unit(int time_format_id)
{
	if (time_format_id == TIME_FORMAT_YYYY) return 1;
	if (time_format_id == TIME_FORMAT_YYYYMM) return 2;
	if (time_format_id == TIME_FORMAT_YYYYMMDD || time_format_id == TIME_FORMAT_DD) return 3;
	return 4;
}

bool check_time_format_id_legal(int time_format_id) {
	if (time_format_id&((long)(0x00000FFFFFFF)) != 0) return false;
	return time_format_id == TIME_FORMAT_YYYY || time_format_id == TIME_FORMAT_YYYYMM || time_format_id == TIME_FORMAT_YYYYMMDD || time_format_id == TIME_FORMAT_YYYYMMDDHHMMSS || time_format_id == TIME_FORMAT_YYYYMMDDSSSSS || time_format_id == TIME_FORMAT_MMO || time_format_id == TIME_FORMAT_MMDDSSSSS || time_format_id == TIME_FORMAT_MMDDHHMM || time_format_id == TIME_FORMAT_DD || time_format_id == TIME_FORMAT_DDHHMMSS || time_format_id == TIME_FORMAT_DDSSSSS || time_format_id == TIME_FORMAT_SSSSS || time_format_id == TIME_FORMAT_HHMMSS || time_format_id == TIME_FORMAT_YYYYMMDDHHMM || time_format_id == TIME_FORMAT_MMDD || time_format_id == TIME_FORMAT_MMDDHHMM || time_format_id == TIME_FORMAT_DDHHMM || time_format_id == TIME_FORMAT_HHMM || time_format_id == TIME_FORMAT_YYYYMMDDHH || time_format_id == TIME_FORMAT_MMDDHH || time_format_id == TIME_FORMAT_DDHH || time_format_id == TIME_FORMAT_HH;
}


template <class T> void get_values_from_string(const char *str_buffer, T **coord_data_buffer, int &array_size, const char *report_string, const char *delimiter, const char *data_type) 
{
	std::vector<T> array;
	char *element;
	*coord_data_buffer = NULL;
	array_size = 0;
	if (str_buffer == NULL || strlen(str_buffer) == 0)
		return;

	char *str_buffer2 = strdup(str_buffer);

	element = strtok(str_buffer2, delimiter);
	while (element) {
		T element_temp;
		if (words_are_the_same(data_type, DATA_TYPE_DOUBLE))
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(element, "%lf", &element_temp) == 1, report_string);
		else if (words_are_the_same(data_type, DATA_TYPE_FLOAT))
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(element, "%f", &element_temp) == 1, report_string);
		else if (words_are_the_same(data_type, DATA_TYPE_INT))
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(element, "%d", &element_temp) == 1, report_string);
		else if (words_are_the_same(data_type, DATA_TYPE_LONG))
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(element, "%ld", &element_temp) == 1, report_string);
		array.push_back(element_temp);
		element = strtok(NULL, delimiter);
	}
	delete str_buffer2;

	if (array.size() == 0)
		return;

	*coord_data_buffer = new T [array.size()];
	for (int i = 0; i < array.size(); i ++)
		(*coord_data_buffer)[i] = array[i];
	array_size = array.size();
}


void get_time_from_string(int host_comp_id, long time_string_value, const char *time_format_str, int time_format_id, bool add_tail, int time_point_type, int &current_year, int &current_month, int &current_day, int &current_second)
{
	//time_string_str: 20160101000000, time_format_str: YYYYMMDDHHMMSS
	//default full_time_format: YYYYMMDDSSSSS, default time_point_type: start(0)
	bool set_minute = false;
	int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0, i = strlen(time_format_str)-1;
	Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);
	while (i >= 0) {
		if (time_format_str[i] == '-' || time_format_str[i] == '.' || time_format_str[i] == ':') continue;
		else if (time_format_str[i] == 'S' && time_format_str[i-2] == 'S') {
			second = time_string_value%100000;
			time_string_value /= 100000;
			set_minute = true;
			i -=5;
		}
		else if (time_format_str[i] == 'S') {
			second = time_string_value%100;
			time_string_value /= 100;
			i -=2;
		}
		else if (time_format_str[i] == 'M' && time_format_str[i-2] == 'H') {//MM means minute
			minute = time_string_value%100;
			time_string_value /=100;
			i -= 2;
		}
		else if (time_format_str[i] == 'H') {
			hour = time_string_value%100;
			time_string_value /= 100;
			i -=2;
		}
		else if (time_format_str[i] == 'D') {
			day = time_string_value%100;
			time_string_value /= 100;
			i -=2;
		}
		else if (time_format_str[i] == 'M') {
			month = time_string_value%100;
			time_string_value /= 100;
			i -=2;
		}
		else if (time_format_str[i] == 'Y') {
			year = time_string_value%10000;
			time_string_value /= 10000;
			i -=4;
		}
		else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error, abnormal time_format %s undetected", time_format_str);
	}
	if (add_tail && get_time_format_smallest_unit(time_format_id) != 4) {//only expand at smaller unit
		if (time_format_id == TIME_FORMAT_YYYY) {
	        if (time_point_type == 0) {
	            month = 1;
	            day = 1;
	        }
	        else if (time_point_type == 2) {
	            month = 12;
	            day = 31;
	            hour = 23;
	            minute = 59;
	            second = 59;
	        }
	        else {
	            month = 7;
	            day = 2;
	            if (time_mgr->get_is_leap_year_on() && (((year%4) == 0 && (year%100) != 0) || (year%400) == 0))
	                hour = 0;
	            else hour = 12;
	        }
		}
		else if (time_format_id == TIME_FORMAT_YYYYMM){
            if (time_point_type == 0)
                day = 1;
            else {
                if (time_mgr->get_is_leap_year_on() && (((year%4) == 0 && (year%100) != 0) || (year%400) == 0)) {
                    if (time_point_type == 1) {
                        day = num_days_of_month_of_leap_year[month-1]/2 + 1;
                        if ((num_days_of_month_of_leap_year[month-1]%2) == 1)
                            hour = 12;
                    }
                    else {
                        day = num_days_of_month_of_leap_year[month-1];
                        hour = 23;
                        minute = 59;
                        second = 59;
                    }
                }
                else {
                    if (time_point_type == 1) {
                        day = num_days_of_month_of_nonleap_year[month-1]/2 + 1;
                        if ((num_days_of_month_of_nonleap_year[month-1]%2) == 1)
                            hour = 12;
                    }
                    else {
                        day = num_days_of_month_of_nonleap_year[month-1];
                        hour = 23;
                        minute = 59;
                        second = 59;
                    }
                }
            }
		}
		else if (time_format_id == TIME_FORMAT_YYYYMMDD || time_format_id == TIME_FORMAT_DD) {
	        if (time_point_type == 1)
	            hour = 12;
	        else if (time_point_type == 2) {
	            hour = 23;
	            minute = 59;
	            second = 59;
	        }
		}
	}
	current_second = hour*3600 + minute*60 + second;
	current_year = year;
	current_month = month == 0? 1: month;
	current_day = day == 0? 1: month;
}


bool time_string_and_format_match(int host_comp_id, const char *time_string, const char *time_format, int time_format_id, bool add_tail, int time_point_type, long &time_value) {
	//time_point_type: 0 start, 1 middle, 2 end; default 0
	//add_tail full_format:YYYYMMDDSSSSS
	int time_string_len = strlen(time_format), i = 0, ind = 0;
	char time_str[time_string_len], full_time_str[32];
	int current_year=0, current_month=0, current_day=0, current_second=0;
	for (i = 0; i < time_string_len; i++) {
		if ((time_string[i] < '0' || time_string[i] > '9')) {
			if (time_string[i] != time_format[i]) return false;
			else time_str[ind++] = time_format[i];
		}
		else if (time_format[i] > 'Z' || time_format[i] < 'A') return false;
		else time_str[ind++] = time_string[i];
	}
	time_str[ind] = '\0';
	time_value = atoi(time_str);

	if (add_tail) {
		get_time_from_string(host_comp_id, time_value, time_format, time_format_id, true, time_point_type, current_year, current_month, current_day, current_second);
		EXECUTION_REPORT(REPORT_LOG, host_comp_id ,true, "1:current_year:%d, current_month:%d, current_day:%d, current_second:%d",current_year, current_month, current_day, current_second);
		get_time_string(full_time_str, "YYYYMMDDSSSSS", false, NULL, current_year, current_month, current_day, current_second);
		EXECUTION_REPORT(REPORT_LOG, host_comp_id ,true, "full_time_str:%s, current_year:%d, current_month:%d, current_day:%d, current_second:%d", full_time_str, current_year, current_month, current_day, current_second);
		bool sure_match = time_string_and_format_match(host_comp_id, full_time_str, "YYYYMMDDSSSSS", TIME_FORMAT_YYYYMMDDSSSSS, false, 0, time_value);
	}
	return true;
}

void get_time_string_from_time_value_and_format(long time_value, char *time_value_format, char *time_string)
{//time_string preallocated
	int i = 0;
	char time_value_str[8];
	sprintf(time_value_str, "%ld", time_value);
	int time_ind = strlen(time_value_str)-1;
	for (i = strlen(time_value_format)-1; i >= 0 ; i--) {
		if (time_value_format[i] > 'A' && time_value_format[i] < 'Z') {
			if (time_ind < 0) time_string[i] = '0';
			else time_string[i] = time_value_str[time_ind--];
		}
		else time_string[i] = time_value_format[i];
	}
	time_string[strlen(time_value_format)] = '\0';
}


void get_time_string_from_time_value_and_format(int host_comp_id, long time_value, char *time_value_format, char *target_time_format, char *time_string)//usually target_time_format is smaller than time_value_format
{
	//time_string preallocated
	int current_year=0, current_month=0, current_day=0, current_second=0;
	get_time_from_string(host_comp_id, time_value, time_value_format, -1, false, -1, current_year, current_month, current_day, current_second);
	get_time_string(time_string, target_time_format, false, NULL, current_year, current_month, current_day, current_second);
}


void get_num_elapsed_time(int host_comp_id, int &elapsed_year, int &elapsed_month, int &elapsed_day, int &elapsed_second)
{
	Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);
	int current_year = time_mgr->get_current_year(), current_month = time_mgr->get_current_month(), current_day = time_mgr->get_current_day(), current_second = time_mgr->get_current_second();
	int start_year = time_mgr->get_start_year(), start_month = time_mgr->get_start_month(), start_day = time_mgr->get_start_day(), start_second = time_mgr->get_start_second();
	elapsed_year = current_year - start_year;
	elapsed_month = current_month - start_month;
	elapsed_day = current_day - start_day;
	elapsed_second = current_second - start_second;
	if (elapsed_second < 0) {
		elapsed_second += 86400;
		elapsed_day --;
	}
	if (elapsed_day < 0) {
		if (time_mgr->get_is_leap_year_on() && time_mgr->is_a_leap_year(current_year)) elapsed_day += num_days_of_month_of_leap_year[current_month-1];
		else elapsed_day += num_days_of_month_of_nonleap_year[current_month-1];
		elapsed_month --;
	}
	if (elapsed_month < 0) {
		elapsed_month += 12;
		elapsed_year --;
	}
	if (elapsed_year < 0) 
		EXECUTION_REPORT(REPORT_ERROR, -1, false, "current_time smaller than model start_time: current_year:%d, current_month:%d, current_day:%d, current_second:%d, start_year:%d, start_month:%d, start_day:%d, start_second:%d", current_year, current_month, current_day, current_second, start_year, start_month, start_day, start_second);
}

int generate_time_on_target_unit(int elapsed_year, int elapsed_month, int elapsed_day, int elapsed_second, int target_unit_id)
{
	int unit_count = 0;
	if (target_unit_id == 1) return elapsed_year;
	if (target_unit_id == 2) return 12*elapsed_year+elapsed_month;
}

void get_time_offset_from_unit(int id_offset_unit, int offset_count, int &offset_year, int &offset_month, int &offset_day, int &offset_second)
{
	if (id_offset_unit == 1) offset_year += offset_count;
	else if (id_offset_unit == 2) offset_month += offset_count;
	else if (id_offset_unit == 3) offset_day += offset_count;
	else if (id_offset_unit == 4) offset_second += offset_count;
}

void add_time_offset(int host_comp_id, int &base_year, int &base_month, int &base_day, int &base_second, int add_count, int unit_id)//base_time is a calendar, not a time offset
{
	Time_mgt *time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);
	EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "before add_offset: base_year:%d, base_month:%d, base_day:%d, base_second:%d", base_year, base_month, base_day, base_second);

	if (unit_id == 1) base_year += add_count;
	else if (unit_id == 2) {
		base_month += add_count;
		while (base_month > 12) {
			base_year ++;
			base_month -= 12;
		}
	}
	else if (unit_id == 3 || unit_id == 4) {
		int num_elapsed_days_at_start_year = time_mgr->get_is_leap_year_on() && time_mgr->is_a_leap_year(base_year)? elapsed_days_on_start_of_month_of_leap_year[base_month-1]+base_day: elapsed_days_on_start_of_month_of_nonleap_year[base_month-1]+base_day;
		if (unit_id == 4) {
			base_second += add_count;
			num_elapsed_days_at_start_year += base_second/SECONDS_PER_DAY;
			base_second %= SECONDS_PER_DAY;
		}
		while (num_elapsed_days_at_start_year >= 365) {
			if (time_mgr->get_is_leap_year_on() && time_mgr->is_a_leap_year(base_year))
				{
					if (num_elapsed_days_at_start_year >= 366)
						{ num_elapsed_days_at_start_year -= 366; base_year++; }
				}
			else { num_elapsed_days_at_start_year -= 365; base_year++; }
		}
		base_month = 0;
		if (time_mgr->get_is_leap_year_on() && time_mgr->is_a_leap_year(base_year)) {
			while (elapsed_days_on_start_of_month_of_leap_year[base_month] < num_elapsed_days_at_start_year) base_month ++;
			base_day = num_elapsed_days_at_start_year - elapsed_days_on_start_of_month_of_leap_year[base_month-1];
		}
		else {
			while (elapsed_days_on_start_of_month_of_nonleap_year[base_month] < num_elapsed_days_at_start_year) base_month ++;
			base_day = num_elapsed_days_at_start_year - elapsed_days_on_start_of_month_of_leap_year[base_month-1];
		}
	}
	EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "after add_offset: base_year:%d, base_month:%d, base_day:%d, base_second:%d", base_year, base_month, base_day, base_second);
}


int get_time_format_biggest_unit(const char *time_format)
{
	if (time_format[0] == 'Y') return 1;
	if (time_format[0] == 'M') return 2;
	if (time_format[0] == 'D') return 3;
	return 4;
}

void get_time_string(char *time_string, const char *time_format_str, bool chop_head, const char *ref_time_format, int current_year, int current_month, int current_day, int current_second)
{//ref_time_format is smaller than time_format_str
	int i = 0, using_hour = 0;
	if (chop_head) {
		if (ref_time_format[0] == 'M') current_year = 0;
		else if (ref_time_format[0] == 'D') {
			current_year = 0;
			current_month = 0;
		}
		else if (ref_time_format[0] == 'H') {
			current_year = 0;
			current_month = 0;
			current_day = 0;
		}
	}

	time_string[0] = '\0';
	while (i < strlen(time_format_str)) {
		if (time_format_str[i] == '-' || time_format_str[i] == '.' || time_format_str[i] == ':')
			sprintf(time_string+strlen(time_string), "%c", time_format_str[i++]);
		else if (time_format_str[i] == 'Y') {
			sprintf(time_string+strlen(time_string), "%04d", current_year);
			i = i + 4;
		}
		else if (time_format_str[i] == 'M' && using_hour == 0) {
			sprintf(time_string+strlen(time_string), "%02d", current_month);
			i = i + 2;
		}
		else if (time_format_str[i] == 'M') {
			sprintf(time_string+strlen(time_string), "%02d", (current_second%3600)/60);
			i = i + 2;
		}
		else if (time_format_str[i] == 'D') {
			sprintf(time_string+strlen(time_string), "%02d", current_day);
			i = i + 2;
		}
		else if (time_format_str[i] == 'H') {
			sprintf(time_string+strlen(time_string), "%02d", current_second/3600);
			i = i + 2;
			using_hour = 1;
		}
		else if (time_format_str[i] == 'S' && using_hour == 0) {
			sprintf(time_string+strlen(time_string), "%05d", current_second);
			i = i + 5;
		}
		else if (time_format_str[i] == 'S') {
			sprintf(time_string+strlen(time_string), "%02d", current_second%60);
			i = i + 2;
		}
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "software error in get_time_string");
	}
}


void Input_file_time_info::copy_info(Input_file_time_info *next_input_file_time) {
	strcpy(this->full_file_name, next_input_file_time->full_file_name);
	this->time_ind_in_dir = next_input_file_time->time_ind_in_dir;
	this->time_position = next_input_file_time->time_position;
	this->global_ind = next_input_file_time->global_ind;
	this->read_flag = next_input_file_time->read_flag;
}

void Input_file_time_info::set_value(void *netcdf_file_object, int time_ind_in_dir, int time_position, int global_ind, bool read_flag)
{
	strcpy(this->full_file_name, full_file_name);
#ifdef USE_PARALLEL_IO
	this->netcdf_file_object = (IO_pnetcdf*)netcdf_file_object;
#else
	this->netcdf_file_object = (IO_netcdf*)netcdf_file_object;
#endif
	this->time_ind_in_dir = time_ind_in_dir;
	this->time_position = time_position;
	this->global_ind = global_ind;
	this->read_flag = read_flag;
}


int Datamodel_mgt::register_datamodel_output_handler(const char *handler_name, int num_fields, int *field_ids, const char *output_datamodel_name, int implicit_or_explicit, int output_grid_id, int handler_output_H2D_grid_id, int handler_output_V1D_grid_id, int sampling_timer_id, int field_instance_ids_size, int handler_level, const char *annotation) 
{
	int API_id = API_ID_HANDLER_DATAMODEL_OUTPUT, i, parallel_io_proc_num;;
	Inout_datamodel *output_datamodel;
	char API_label[256];

	get_API_hint(-1, API_ID_HANDLER_DATAMODEL_OUTPUT, API_label);

	common_checking_for_datamodel_handler_registration(num_fields, field_ids, implicit_or_explicit, output_grid_id, handler_output_H2D_grid_id, handler_output_V1D_grid_id, sampling_timer_id, field_instance_ids_size, handler_name, annotation, true);
	int host_comp_id = memory_manager->get_field_instance(field_ids[0])->get_comp_id();

	for (int m = 0; m < output_handlers.size(); m ++)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !words_are_the_same(handler_name, output_handlers[m]->get_handler_name()) && host_comp_id == output_handlers[m]->get_host_comp_id(), "Error happens when calling \"%s\" at the model code with the annotation \"%s\", a handler with the name \"%s\" has already been registered at the model code with the annotation \"%s\", Please Verify.", API_label, annotation, handler_name, output_handlers[m]->get_handler_annotation());

	for (i = 0; i < output_datamodels.size(); i ++)
		if (words_are_the_same(output_datamodel_name, output_datamodels[i]->get_datamodel_name()) && host_comp_id == output_datamodels[i]->get_host_comp_id()) {
			output_datamodel = output_datamodels[i];
			EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "The output_datamodel \"%s\" under the DIR %s, has already been configured.", output_datamodel_name, output_datamodel->get_datamodel_output_data_dir());
			break;
		}

	if (i == output_datamodels.size()) {
		output_datamodel = new Inout_datamodel(host_comp_id, output_datamodel_name,  OUTPUT_DATAMODEL, annotation);
		output_datamodels.push_back(output_datamodel);
	}

	Output_handler *new_output_handler = new Output_handler(handler_name, output_datamodel_name, NULL, get_next_output_handler_id(), num_fields, field_ids, NULL, implicit_or_explicit, output_grid_id, handler_output_H2D_grid_id, handler_output_V1D_grid_id, sampling_timer_id, field_instance_ids_size, annotation, true, true, handler_level);
	output_handlers.push_back(new_output_handler);

	return new_output_handler->get_handler_id();
}


int Datamodel_mgt::register_field_instances_output_handler(int num_fields, int *field_instance_ids, const char *file_name, const char *file_type, int implicit_or_explicit, int sampling_timer_id, int output_timer_id, int file_timer_id, int output_grid_id, int inst_or_aver, const char *float_datatype, const char *integer_datatype, int field_instance_ids_size, bool considering_surface_fields, int handler_level, const char *annotation) 
{
	int API_id = API_ID_HANDLER_FIELD_INSTANCES_OUTPUT, i, parallel_io_proc_num;;
	char handler_name[NAME_STR_SIZE], output_datamodel_name[NAME_STR_SIZE];

	common_checking_for_datamodel_handler_registration(num_fields, field_instance_ids, implicit_or_explicit, output_grid_id, -1, -1, sampling_timer_id, field_instance_ids_size, NULL, annotation, false);
	int host_comp_id = memory_manager->get_field_instance(field_instance_ids[0])->get_comp_id();

	int handler_id = get_next_output_handler_id();
	sprintf(handler_name, "output_handler_to_file_%s", file_name);
	sprintf(output_datamodel_name, "output_datamodel_to_file_%s", file_name);

	Inout_datamodel *output_datamodel = new Inout_datamodel(host_comp_id, output_datamodel_name, file_name, file_type, output_timer_id, file_timer_id, output_grid_id, inst_or_aver, float_datatype, integer_datatype, OUTPUT_DATAMODEL, annotation);
	output_datamodels.push_back(output_datamodel);

	Output_handler *new_output_handler = new Output_handler(handler_name, output_datamodel_name, output_datamodel, handler_id, num_fields, field_instance_ids, inst_or_aver, implicit_or_explicit, output_grid_id, -1, -1, sampling_timer_id, field_instance_ids_size, annotation, false, considering_surface_fields, handler_level);
	output_handlers.push_back(new_output_handler);

	return handler_id;
}


void Output_handler::initiate_handler_info(const char *output_handler_name, const char *datamodel_name_str, int handler_id, int sampling_timer_id, int implicit_or_explicit, int num_fields, int output_grid_id, int handler_output_H2D_grid_id, int handler_output_V1D_grid_id, int *export_field_ids, int level, const char *annotation) 
{
	this->handler_name = strdup(output_handler_name);
	this->datamodel_name =strdup(datamodel_name_str);
	this->handler_id = handler_id;
	this->sampling_timer_id = sampling_timer_id;
	this->implicit_or_explicit = implicit_or_explicit;
	this->num_fields = num_fields;
	this->output_grid_id = output_grid_id;
	this->handler_output_H2D_grid_id = handler_output_H2D_grid_id;
	this->handler_output_V1D_grid_id = handler_output_V1D_grid_id;

	this->host_comp_id = memory_manager->get_field_instance(export_field_ids[0])->get_comp_id();
	this->pio_proc_num = -1;
	this->io_comm = MPI_COMM_NULL;
	this->io_proc_mark = 0;
	this->total_num_procs = comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_num_procs();
	this->level = level;
	this->annotation = strdup(annotation);
}


Output_handler::Output_handler(const char *output_handler_name, const char *datamodel_name_str, Inout_datamodel *existing_datamodel, int handler_id, int num_fields, int *export_field_ids, int inst_or_aver, int implicit_or_explicit, int output_grid_id, int handler_output_H2D_grid_id, int handler_output_V1D_grid_id, int sampling_timer_id, int field_instance_ids_size, const char *annotation, bool is_datamodel_handler, bool considering_surface_fields, int level) 
{
	char handler_export_interface_name[NAME_STR_SIZE];
	int num_expanded_export_fields;

	handler_type = OUTPUT_HANDLER;
	initiate_handler_info(output_handler_name, datamodel_name_str, handler_id, sampling_timer_id, implicit_or_explicit, num_fields, output_grid_id, handler_output_H2D_grid_id, handler_output_V1D_grid_id, export_field_ids, level, annotation);
	this->inout_datamodel = existing_datamodel != NULL? existing_datamodel : datamodel_mgr->search_output_datamodel(host_comp_id, this->datamodel_name);
	time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);
	this->fields_config_info = inout_datamodel->fields_config_info;

#ifdef USE_PARALLEL_IO
	pio_proc_num = datamodel_mgr->get_comp_PIO_proc_setting(host_comp_id, io_proc_stride, io_proc_mark, io_comm);
#endif

	//register export_interface for handler
	this->export_timer_id = (sampling_timer_id != -1? sampling_timer_id: timer_mgr->define_timer(host_comp_id, "steps", 1, 0, 0, this->annotation));
	sprintf(handler_export_interface_name, "%s_export_interface", output_handler_name);
	num_expanded_export_fields = num_fields;
	if (considering_surface_fields) {
		expanded_export_field_ids = add_surface_fields_to_export_interface(num_expanded_export_fields, export_field_ids, field_instance_ids_size);
	}
	else {
		expanded_export_field_ids = new int [num_fields];
		for (int i = 0; i < num_fields; i++)
			expanded_export_field_ids[i] = export_field_ids[i];
			num_expanded_export_fields = num_fields;
	}
	if (level == 2)
		handler_export_interface = new Inout_interface(handler_export_interface_name, TYPE_OUTPUT_HANDLER_ID_PREFIX|output_procedures.size(), 1, num_expanded_export_fields, expanded_export_field_ids, field_instance_ids_size, export_timer_id, 0, "field_instance_IDs", this->annotation, API_ID_HANDLER_DATAMODEL_OUTPUT, INTERFACE_SOURCE_IO_OUTPUT, true);
	else handler_export_interface = NULL;

	config_datamodel_import_interface_parameters(num_fields, expanded_export_field_ids, considering_surface_fields);//or expanded_export_field_ids
	if (considering_surface_fields)
		delete [] expanded_export_field_ids;
}


void Output_handler::config_datamodel_import_interface_parameters(int num_fields, int *export_field_ids, bool considering_surface_fields) 
{
	bool field_output_config_from_common_setting_or_not;
	Import_field_info *import_field_info;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, import_fields_info.size() == 0, "Software error in Output_handler::config_datamodel_import_interface_parameters");
	for (int i = 0; i < inout_datamodel->default_settings.size(); i ++) {
		datamodel_mgr->add_handlers_field_mem_buf_mark();
		for (int j = 0; j < num_fields; j++) {
			const char *dst_field_grid_name = NULL, *config_dst_field_int_type = NULL, *config_dst_field_float_type = NULL, *dst_field_data_type;
			int dst_field_grid_id = -1, dst_field_decomp_id = -1;
			Field_mem_info *handler_field_info = memory_manager->get_field_instance(export_field_ids[j]);
			const char *dst_field_unit = handler_field_info->get_unit();
			int field_index_in_xml = -1, field_instance_id;
			if (words_are_the_same(inout_datamodel->default_settings[i]->field_specification, "fields")) {
				if (find_field_info_from_datamodel_config(handler_field_info->get_field_name(), fields_config_info[i], field_index_in_xml))
					field_output_config_from_common_setting_or_not = false;
				else continue;
			}
			else field_output_config_from_common_setting_or_not = (words_are_the_same(inout_datamodel->default_settings[i]->field_specification, "default") || !find_field_info_from_datamodel_config(handler_field_info->get_field_name(), fields_config_info[i], field_index_in_xml));

			if (field_output_config_from_common_setting_or_not) {
				dst_field_grid_name = inout_datamodel->default_settings[i]->default_output_grid_name;
				config_dst_field_float_type = inout_datamodel->default_settings[i]->default_float_type;
				config_dst_field_int_type = inout_datamodel->default_settings[i]->default_integer_type;
			}
			else {
				dst_field_grid_name = fields_config_info[i][field_index_in_xml]->grid_name;
				config_dst_field_float_type = fields_config_info[i][field_index_in_xml]->float_datatype;
				config_dst_field_int_type = fields_config_info[i][field_index_in_xml]->integer_datatype;
				if (!words_are_the_same(fields_config_info[i][field_index_in_xml]->unit, ""))
					dst_field_unit = fields_config_info[i][field_index_in_xml]->unit;
			}
			determine_field_output_grid(handler_field_info, dst_field_grid_name, dst_field_grid_id);
			if (dst_field_grid_id != -1) {
#ifdef USE_PARALLEL_IO
				Decomp_info *new_decomp = decomps_info_mgr->generate_parallel_decomp_for_parallel_IO(original_grid_mgr->search_grid_info(dst_field_grid_id), this->io_proc_stride, this->pio_proc_num, this->io_proc_mark);
#else
				Decomp_info *new_decomp = decomps_info_mgr->generate_default_parallel_decomp_serial(original_grid_mgr->search_grid_info(dst_field_grid_id));
#endif
				if (new_decomp != NULL)
					dst_field_decomp_id = new_decomp->get_decomp_id();
			}
			else dst_field_grid_id = host_comp_id;

			if (words_are_the_same(handler_field_info->get_data_type(), DATA_TYPE_DOUBLE) || words_are_the_same(handler_field_info->get_data_type(), DATA_TYPE_FLOAT)) {
				if (!words_are_the_same(config_dst_field_float_type, ""))
					dst_field_data_type = config_dst_field_float_type;
				else dst_field_data_type = handler_field_info->get_data_type();
			}
			else if (words_are_the_same(handler_field_info->get_data_type(), DATA_TYPE_INT) || words_are_the_same(handler_field_info->get_data_type(), DATA_TYPE_SHORT) || words_are_the_same(handler_field_info->get_data_type(), DATA_TYPE_LONG)) {
				if (!words_are_the_same(config_dst_field_int_type, ""))
					dst_field_data_type = config_dst_field_int_type;
				else dst_field_data_type = handler_field_info->get_data_type();
			}
			else dst_field_data_type = handler_field_info->get_data_type();

			EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "field_name:%s, data_type:%s, original_data_type:%s", handler_field_info->get_field_name(), dst_field_data_type, handler_field_info->get_data_type());
			import_field_info = new Import_field_info();
			if (level != 2 && dst_field_decomp_id != -1) {
				int new_dst_field_decomp_id = decomps_info_mgr->generate_empty_decomp(dst_field_decomp_id);
				dst_field_decomp_id = new_dst_field_decomp_id;
			}

			Field_mem_info *dst_field_instance = memory_manager->alloc_mem(handler_field_info->get_field_name(), dst_field_decomp_id, dst_field_grid_id, BUF_MARK_IO_FIELD_MIRROR ^ datamodel_mgr->get_handlers_field_mem_buf_mark(), dst_field_data_type, dst_field_unit, this->annotation, false, false);
			dst_field_instance->set_usage_tag(REG_FIELD_TAG_NONE);
			import_field_info->model_field_instance_id = handler_field_info->get_field_instance_id();
			import_field_info->io_field_instance_id = dst_field_instance->get_field_instance_id();
			import_field_info->import_field_operation = strdup(inout_datamodel->default_settings[i]->default_operation);
			import_field_info->field_name_in_file = NULL;
			import_fields_info.push_back(import_field_info);
		}

		if (considering_surface_fields) {
			for (int j = 0; j < import_fields_info.size(); j ++) {
				Field_mem_info *dst_field_instance = memory_manager->get_field_instance(import_fields_info[j]->io_field_instance_id);
				int dst_grid_id = dst_field_instance->get_grid_id();
				if (dst_grid_id != -1) {
					Original_grid_info *dst_original_grid = original_grid_mgr->search_grid_info(dst_grid_id);
					if (dst_original_grid->has_sigma_sub_grid()) {
						Original_grid_info *src_original_grid = NULL;
						for (int k = 0; k < num_fields; k ++)
							if (words_are_the_same(dst_field_instance->get_field_name(), memory_manager->get_field_instance(export_field_ids[k])->get_field_name())) {
								src_original_grid = original_grid_mgr->search_grid_info(memory_manager->get_field_instance(export_field_ids[k])->get_grid_id());
								break;
							}
						EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, src_original_grid->has_sigma_sub_grid(), "Software error in Output_handler::config_datamodel_import_interface_parameters");
						Field_mem_info *src_surface_field = memory_manager->get_field_instance(src_original_grid->get_bottom_field_id());
						if (!search_import_field_interface_surface_field(src_surface_field->get_field_name())) {
							Import_field_info *import_field_info = new Import_field_info();
								Field_mem_info *new_surface_field = memory_manager->alloc_mem(src_surface_field->get_field_name(), dst_field_instance->get_decomp_id(), dst_original_grid->get_H2D_sub_grid()->get_grid_id(), dst_field_instance->get_buf_mark(), src_surface_field->get_data_type(), src_surface_field->get_unit(), this->annotation, false, false);//interpolate to field_instance grid
								new_surface_field->set_usage_tag(REG_FIELD_TAG_NONE);
								import_field_info->model_field_instance_id = src_surface_field->get_field_instance_id();
								import_field_info->io_field_instance_id = new_surface_field->get_field_instance_id();
								import_field_info->import_field_operation = strdup(import_fields_info[j]->import_field_operation);
								import_fields_info.push_back(import_field_info);
						}
					}
				}
			}
		}

		EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "import_fields_info size:%d", import_fields_info.size());
		int *import_inst_field_ids = new int [import_fields_info.size()];
		int *import_aver_field_ids = new int [import_fields_info.size()];
		int *import_inst_src_field_ids = new int [import_fields_info.size()];
		int *import_aver_src_field_ids = new int [import_fields_info.size()];

		int inst_id = 0, aver_id = 0;
		for (int m = 0; m < import_fields_info.size(); m++) {
			if (words_are_the_same(import_fields_info[m]->import_field_operation, "inst")) {
				import_inst_field_ids[inst_id] = import_fields_info[m]->io_field_instance_id;
				import_inst_src_field_ids[inst_id++] = import_fields_info[m]->model_field_instance_id;
			}
			else if (words_are_the_same(import_fields_info[m]->import_field_operation, "aver")) {
				import_aver_field_ids[aver_id] = import_fields_info[m]->io_field_instance_id;
				import_aver_src_field_ids[aver_id++] = import_fields_info[m]->model_field_instance_id;
			}
		}
		//register import interfaces
		int import_timer_id = inout_datamodel->io_timer_ids[i];
		EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "inst_id:%d, aver_id:%d, num_fields:%d, level:%d", inst_id, aver_id, num_fields, level);
		register_operation_import_interface(inst_id, import_inst_field_ids, import_inst_src_field_ids, import_fields_info.size(), import_timer_id, inout_datamodel->file_timer_ids[i], inout_datamodel->default_settings[i]->file_mid_name, "inst", i);
		register_operation_import_interface(aver_id, import_aver_field_ids, import_aver_src_field_ids, import_fields_info.size(), import_timer_id, inout_datamodel->file_timer_ids[i], inout_datamodel->default_settings[i]->file_mid_name, "aver", i);
		import_fields_info.clear();
		delete [] import_inst_field_ids, import_aver_field_ids, import_inst_src_field_ids, import_aver_src_field_ids;
	}
}


int *Output_handler::add_surface_fields_to_export_interface(int &num_fields, const int *field_ids, int &field_instance_ids_size) 
{
	int *new_field_ids = new int [field_instance_ids_size*2];
	int num_new_field_ids = num_fields;
	int new_field_instance_ids_size = field_instance_ids_size * 2;

	for (int i = 0; i < num_fields; i ++)
		new_field_ids[i] = field_ids[i];

	for (int i = 0; i < num_fields; i ++) {
		Field_mem_info *field_inst = memory_manager->get_field_instance(field_ids[i]);
		if (field_inst->get_grid_id() == -1)
			continue;
		int src_surface_field_id = original_grid_mgr->search_grid_info(field_inst->get_grid_id())->get_bottom_field_id();
		if (src_surface_field_id != -1 && !search_export_field_id(new_field_ids, num_new_field_ids, src_surface_field_id))
			new_field_ids[num_new_field_ids++] = src_surface_field_id;
	}

	num_fields = num_new_field_ids;
	field_instance_ids_size = new_field_instance_ids_size;
	return new_field_ids;
}


bool Output_handler::search_export_field_id(int *field_ids, int num_fields, int target_field_id) {
	for (int i = 0; i < num_fields; i ++) {
		if (target_field_id == field_ids[i])
			return true;
	}
	return false;
}


void Output_handler::check_grid_dim(int &dest_grid_id, Field_mem_info *src_field_instance) 
{
	if (src_field_instance->get_grid_id() == -1) {
		dest_grid_id = -1;
		return;
	}

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(dest_grid_id), "Software error in Output_handler::check_grid_dim grid_id not legal:%d", dest_grid_id);

	Original_grid_info *specified_grid = original_grid_mgr->search_grid_info(dest_grid_id);
	Original_grid_info *field_grid = original_grid_mgr->search_grid_info(src_field_instance->get_grid_id());
	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "field_name: %s, dest grid_name: %s, src grid_name: %s", src_field_instance->get_field_name(), specified_grid->get_grid_name(),field_grid->get_grid_name());
	Original_grid_info *output_grid = specified_grid->get_sub_original_grid_corresponding_to_another(field_grid);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, output_grid != NULL, "Error happens when the data model \"%s\" tries to output the field \"%s\": the dimensions of the grid of this field (grid is \"%s\") does not match the grid specified (\"%s\") for the output. Please verify.", handler_name, field_grid->get_grid_name(), specified_grid->get_grid_name());

	dest_grid_id = output_grid->get_grid_id();
}


void Output_handler::register_operation_import_interface(int num_insts, int *field_instance_ids, int *src_field_instance_ids, int field_instance_size, int output_timer_id, int file_timer_id, const char *file_mid_name, const char *name_tag, int interface_index) 
{
	int interface_id;
	Inout_interface *new_interface;
	Output_file_procedure *output_procedure = new Output_file_procedure();

	if (num_insts > 0) {
		int inst_or_aver = words_are_the_same(name_tag, "inst")?0:1;
		if (level == 2) {
			char interface_name[NAME_STR_SIZE];
			int *field_necessity = new int [num_insts];
			sprintf(interface_name, "PIO_output_import_interface_for_%s_%s", handler_name, name_tag);
			for (int i = 0; i < num_insts; i ++)
				field_necessity[i] = 1;
			new_interface = new Inout_interface(interface_name, TYPE_OUTPUT_HANDLER_ID_PREFIX|output_procedures.size(), 0, num_insts, field_instance_ids, field_instance_size, output_timer_id, inst_or_aver, "field_instance_IDs", this->annotation, API_ID_HANDLER_DATAMODEL_OUTPUT, INTERFACE_SOURCE_IO_OUTPUT, true);
			new_interface->set_fields_necessity(field_necessity, num_insts, this->annotation);
			output_procedure->interface = new_interface;
			delete field_necessity;
		}
		else {
			output_procedure->interface = NULL;
			output_procedure->interface_timer_id = output_timer_id;
		}

		output_procedure->inst_or_aver = inst_or_aver;
		for (int i = 0; i < num_insts; i++) {
			output_procedure->fields_name.push_back(strdup(memory_manager->get_field_instance(field_instance_ids[i])->get_field_name()));
			output_procedure->dst_field_instance_ids.push_back(field_instance_ids[i]);//level = 2: parallel decomp, dst attributes, size >= 1; level == 1: null decomp, dst attributes;
			output_procedure->src_field_instance_ids.push_back(src_field_instance_ids[i]);
		}
		if (level == 1) {
			for (int i = 0; i < output_procedure->src_field_instance_ids.size(); i ++) {
				output_procedure->averaging_field_mems.push_back(NULL);
				output_procedure->outer_level_averaging_algorithm.push_back(NULL);
				if (output_procedure->inst_or_aver == 1) {
					datamodel_mgr->add_handlers_field_mem_buf_mark();
					int new_buf_mark = datamodel_mgr->get_handlers_field_mem_buf_mark() | ((((get_handler_id() + 16)&TYPE_ID_SUFFIX_MASK)) << 12);
					Field_mem_info *src_field_mem = memory_manager->get_field_instance(src_field_instance_ids[i]);
					output_procedure->averaging_field_mems[i] = memory_manager->alloc_mem(src_field_mem, BUF_MARK_AVERAGED_INNER, new_buf_mark, NULL, false, false);
					EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "output_handler_field_mem_buf_mark:%x, TYPE_ID_SUFFIX_MASK:%x,handler_id:%x, and suffix:%x, base:%x, buf_mark:%x, real aver_buf_mark:%x", datamodel_mgr->get_handlers_field_mem_buf_mark(), TYPE_ID_SUFFIX_MASK, get_handler_id(), get_handler_id()&TYPE_ID_SUFFIX_MASK, ((((get_handler_id() +16)&TYPE_ID_SUFFIX_MASK)) << 12), new_buf_mark, BUF_MARK_AVERAGED_INNER ^ new_buf_mark);
					output_procedure->outer_level_averaging_algorithm[i] = new Runtime_cumulate_average_algorithm(NULL, src_field_mem, output_procedure->averaging_field_mems[i]);
				}
			}
		}

		output_procedure->interface_index = interface_index;
		output_procedure->netcdf_file_object = NULL;
		output_procedure->file_timer_id = file_timer_id;
		output_procedure->file_mid_name = strdup(file_mid_name);
		output_procedures.push_back(output_procedure);
		if (level == 2)
			inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(new_interface, handler_export_interface, -1, false);//import, export
	}
}


Inout_interface *Output_handler::get_handler_interface(int handler_interface_id) 
{
	bool is_interface_legal = true;
    if ((handler_interface_id & TYPE_ID_PREFIX_MASK) != TYPE_OUTPUT_HANDLER_ID_PREFIX)
        is_interface_legal = false;
    else
    	is_interface_legal = (handler_interface_id&TYPE_ID_SUFFIX_MASK) < output_procedures.size();
    EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_interface_legal, "Software error in Output_handler::get_handler_interface: NULL interface");
    return output_procedures[handler_interface_id&TYPE_ID_SUFFIX_MASK]->interface;
}


void Output_handler::generate_a_coupling_procedure(int export_interface_id, int import_interface_id) 
{
	coupling_generator->synchronize_latest_connection_id(comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, ""));
	Coupling_connection *coupling_connection = new Coupling_connection(coupling_generator->apply_connection_id());
	coupling_connection->export_interface = get_handler_interface(export_interface_id);
	coupling_connection->import_interface = get_handler_interface(import_interface_id);
	strcpy(coupling_connection->dst_comp_full_name, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(host_comp_id, false, "in Output_handler::generate_a_coupling_procedure")->get_full_name());
	strcpy(coupling_connection->dst_interface_name, coupling_connection->import_interface->get_interface_name());

	std::vector<const char*> import_fields_name;
	coupling_connection->import_interface->get_fields_name(&import_fields_name);
	for (int k = 0; k < import_fields_name.size(); k ++) {
		coupling_connection->fields_name.push_back(strdup(import_fields_name[k]));
	}

	std::pair<const char*, const char*> src_comp_interface;
	src_comp_interface.first = strdup(coupling_connection->dst_comp_full_name);
	src_comp_interface.second = strdup(coupling_connection->export_interface->get_interface_name());
	coupling_connection->src_comp_interfaces.push_back(src_comp_interface);

	coupling_connection->src_comp_node = comp_comm_group_mgt_mgr->search_global_node(coupling_connection->src_comp_interfaces[0].first);
	coupling_connection->dst_comp_node = comp_comm_group_mgt_mgr->search_global_node(coupling_connection->dst_comp_full_name);
	coupling_connection->current_proc_id_src_comp = coupling_connection->src_comp_node->get_current_proc_local_id();
	coupling_connection->current_proc_id_dst_comp = coupling_connection->dst_comp_node->get_current_proc_local_id();
	coupling_connection->create_union_comm();
	coupling_connection->exchange_connection_fields_info();
	coupling_connection->generate_interpolation(false);
	coupling_connection->export_procedure = new Connection_coupling_procedure(coupling_connection->export_interface, coupling_connection);
	coupling_connection->export_interface->add_coupling_procedure(coupling_connection->export_procedure);
	coupling_connection->import_procedure = new Connection_coupling_procedure(coupling_connection->import_interface, coupling_connection);
	coupling_connection->import_interface->add_coupling_procedure(coupling_connection->import_procedure);
	coupling_connection->generate_data_transfer();
}


void Output_handler::determine_field_output_grid(Field_mem_info *source_field_instance, const char *dst_grid_name, int &grid_id) 
{	
	if (source_field_instance->get_grid_id() == -1) {
		grid_id = -1;
		return;
	}
	if (!words_are_the_same(dst_grid_name, "") && !words_are_the_same(dst_grid_name, "original_grid"))
		get_or_generate_field_real_output_grid(dst_grid_name, source_field_instance, grid_id);		
	else if (!words_are_the_same(dst_grid_name, "original_grid") && output_grid_id != -1)
		grid_id = output_grid_id;
	else grid_id = source_field_instance->get_grid_id();

	check_grid_dim(grid_id, source_field_instance);

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, original_grid_mgr->is_grid_id_legal(grid_id), "Software error in Output_handler::determine_field_output_grid grid_id:%d", grid_id);
}


char *Output_handler::modify_special_v3d_grid_name(char *field_instance_grid_name, const char *field_name) {
	//char *new_field_instance_grid_name = new char [64+strlen(field_instance_grid_name)+strlen(field_name)];
	char *new_field_instance_grid_name = new char [NAME_STR_SIZE];
	sprintf(new_field_instance_grid_name, "%s_%s", field_instance_grid_name, field_name);
	return new_field_instance_grid_name;
}


void Output_handler::get_or_generate_field_real_output_grid(const char *field_instance_grid_name, Field_mem_info *source_field_instance, int &grid_id) 
{
	V3d_grid_info *v3d_grid_info;
	char *new_field_instance_grid_name, *current_field_instance_grid_name = NULL;
	int h2d_subgrid_id, v1d_subgrid_id, special_v3d_grid_id, grid_index = -1;

	if (inout_datamodel->is_special_v3d_grid(field_instance_grid_name, v3d_grid_info)) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, source_field_instance->get_grid_id() != -1 && original_grid_mgr->search_grid_info(source_field_instance->get_grid_id())->get_V3D_sub_CoR_grid() != NULL, "Error happens when registering the output handler \"%s\" corresponding to the data model XML configuration file \"%s\" at the model code with the annotation \"%s\": the grid of the model field \"%s\" is not a V3D grid while the configuration file has been specified to output it to a V3D grid. Please verify.", handler_name, inout_datamodel->get_XML_file_name(), annotation, source_field_instance->get_field_name());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !(words_are_the_same(v3d_grid_info->h2d_subgrid_name, "handler_output_H2D_grid") && words_are_the_same(v3d_grid_info->v1d_subgrid_name, "handler_output_V1D_grid")), "Error happens when registering the output handler \"%s\" corresponding to the data model XML configuration file \"%s\" at the model code with the annotation \"%s\": both \"handler_output_H2D_grid\" and \"handler_output_V1D_grid\" have been specified to the 3D grid for outputting the model field \"%s\" in the configuration file. That is not allowed currently. Please verify.", handler_name, inout_datamodel->get_XML_file_name(), annotation, source_field_instance->get_field_name());
		char *new_field_instance_grid_name = new char [NAME_STR_SIZE];
		sprintf(new_field_instance_grid_name, "%s_%s_%s", inout_datamodel->datamodel_name, source_field_instance->get_grid_name(), v3d_grid_info->grid_name);
		current_field_instance_grid_name = strdup(new_field_instance_grid_name);
		if (original_grid_mgr->search_grid_info(current_field_instance_grid_name, host_comp_id) == NULL) {
			if (words_are_the_same(v3d_grid_info->h2d_subgrid_name, "handler_output_H2D_grid"))
				h2d_subgrid_id = handler_output_H2D_grid_id != -1? handler_output_H2D_grid_id: original_grid_mgr->search_grid_info(source_field_instance->get_grid_id())->get_H2D_sub_grid()->get_grid_id();
			else {
				char *hg_subgrid_name = new char [NAME_STR_SIZE];
				rename_datamodel_grid(hg_subgrid_name, inout_datamodel->datamodel_name, v3d_grid_info->h2d_subgrid_name);
				h2d_subgrid_id = original_grid_mgr->search_grid_info(hg_subgrid_name, host_comp_id)->get_grid_id();
				delete [] hg_subgrid_name;
			}
			if (words_are_the_same(v3d_grid_info->v1d_subgrid_name, "handler_output_V1D_grid"))
				v1d_subgrid_id = handler_output_V1D_grid_id != -1? handler_output_V1D_grid_id: original_grid_mgr->search_grid_info(original_grid_mgr->search_grid_info(source_field_instance->get_grid_id())->get_V1D_sub_CoR_grid()->get_grid_name(), host_comp_id)->get_grid_id();
			else {
				char *vg_subgrid_name = new char [NAME_STR_SIZE];
				rename_datamodel_grid(vg_subgrid_name, inout_datamodel->datamodel_name, v3d_grid_info->v1d_subgrid_name);
				v1d_subgrid_id = original_grid_mgr->search_grid_info(vg_subgrid_name, host_comp_id)->get_grid_id();
				delete [] vg_subgrid_name;
			}

			new_field_instance_grid_name = new char [NAME_STR_SIZE];
			sprintf(new_field_instance_grid_name, "%s_%s_%s", inout_datamodel->datamodel_name, source_field_instance->get_grid_name(), v3d_grid_info->grid_name);
			special_v3d_grid_id = inout_datamodel->register_datamodel_v3d_grid(v3d_grid_info->dimension_order, new_field_instance_grid_name, h2d_subgrid_id , v1d_subgrid_id, v3d_grid_info->line_number);
			delete [] new_field_instance_grid_name;
			if (!words_are_the_same(v3d_grid_info->mid_point_grid_name, "")) {
				new_field_instance_grid_name = new char [NAME_STR_SIZE];
				sprintf(new_field_instance_grid_name, "%s_%s", source_field_instance->get_grid_name(), v3d_grid_info->mid_point_grid_name);
				inout_datamodel->register_datamodel_mid_point_v3d_grid(new_field_instance_grid_name, special_v3d_grid_id);
				delete [] new_field_instance_grid_name;
			}
			if (original_grid_mgr->search_grid_info(v1d_subgrid_id)->get_original_CoR_grid()->is_sigma_grid()) {
				char API_label[NAME_STR_SIZE];
				get_API_hint(-1, API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD, API_label);
				original_grid_mgr->set_3d_grid_bottom_field(host_comp_id, special_v3d_grid_id, -1, BOTTOM_FIELD_VARIATION_EXTERNAL, API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD, API_label, this->annotation);
			}
			if (!words_are_the_same(v3d_grid_info->mid_point_grid_name, "")) {
				inout_datamodel->fields_grid_info[inout_datamodel->fields_grid_info.size()-1]->surface_field_name = strdup(v3d_grid_info->surface_field_name);
				inout_datamodel->fields_grid_info[inout_datamodel->fields_grid_info.size()-2]->surface_field_name = strdup(v3d_grid_info->surface_field_name);
			}
			else inout_datamodel->fields_grid_info[inout_datamodel->fields_grid_info.size()-1]->surface_field_name = strdup(v3d_grid_info->surface_field_name);
		}
	}

	if (current_field_instance_grid_name == NULL)
		current_field_instance_grid_name = strdup(field_instance_grid_name);

	grid_id = inout_datamodel->get_datamodel_grid_id(current_field_instance_grid_name, grid_index, true, false);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, grid_id != -1, "Software error at Output_handler::get_or_generate_field_real_output_grid \"%s\"", current_field_instance_grid_name);
	delete [] current_field_instance_grid_name;
}


bool Output_handler::search_import_field_interface_surface_field(const char *surface_field_name) 
{
	for (int i = 0; i < import_fields_info.size(); i ++) {
		if (words_are_the_same(surface_field_name, memory_manager->get_field_instance(import_fields_info[i]->io_field_instance_id)->get_field_name()))
			return true;
	}
	return false;
}

Field_mem_info *Inout_handler::find_handler_field_instance(char *field_name, int *export_field_ids, int num_fields, bool report_error) 
{
	for (int i = 0; i < num_fields; i ++) {
		Field_mem_info *handler_field_info = memory_manager->get_field_instance(export_field_ids[i]);
		if (words_are_the_same(handler_field_info->get_field_name(), field_name)) {
			return handler_field_info;
		}
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !report_error, "Error happens when registering Inout_handler \"%s\" at the model code with the annotation \"%s\": the corresponding data model \"%s\" (the XML configuration file is \"%s\") wants to output the field \"%s\" while the registration of the handler does not specify the corresponding field instance. Please verify.", handler_name, annotation, datamodel_name, inout_datamodel->get_XML_file_name(), field_name);
	return NULL;
}


bool Inout_handler::find_field_info_from_datamodel_config(const char *handler_field_name, std::vector<Field_config_info*> field_info_list, int &index) {
	for (int i = 0; i < field_info_list.size(); i ++) {
		if (words_are_the_same(field_info_list[i]->name_in_model, handler_field_name)) {
			index = i;
			return true;
		}
	}
	return false;
}


Field_mem_info *Output_handler::get_unique_IO_field_mem()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, output_procedures.size() == 1, "Software error in Output_handler::get_unique_IO_field_mem size:%d", output_procedures.size());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, output_procedures[0]->interface->get_num_fields_mem_registered() == 1, "Software error in Output_handler::get_unique_IO_field_mem");

	return output_procedures[0]->interface->get_field_mem(0);
}

Inout_datamodel *Datamodel_mgt::search_output_datamodel(int host_comp_id, char *datamodel_name) 
{
	for (int i = 0; i < output_datamodels.size(); i ++) {
		if (output_datamodels[i]->get_host_comp_id() == host_comp_id && words_are_the_same(output_datamodels[i]->get_datamodel_name(), datamodel_name))
			return output_datamodels[i];
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Error happens when searching for datamodel named \"%s\", Please verify the corresponding XML configuration file.", datamodel_name);
	return NULL;
}


void Datamodel_mgt::common_checking_for_datamodel_handler_registration(int num_fields, int *field_ids, int implicit_or_explicit, int output_grid_id, int handler_output_H2D_grid_id, int handler_output_V1D_grid_id, int sampling_timer_id, int field_instance_ids_size, const char *handler_name, const char *annotation, bool is_handler_for_datamodel)
{
	int comp_id = -1;
	char str[NAME_STR_SIZE], API_label[256];

	get_API_hint(-1, API_ID_HANDLER_DATAMODEL_OUTPUT, API_label);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_fields > 0, "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\", the parameter \"num_field_instances\" (currently is %d) cannot be smaller than 1. Please verify the model code.", API_label, annotation, num_fields);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_fields <= field_instance_ids_size, "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the array size (currently is %d) of parameter \"%s\" cannot be smaller than the parameter \"num_field_instances\" (currently is %d). Please verify the corresponding model code.", API_label, annotation, field_instance_ids_size, "field_instance_ids", num_fields);
	comp_id = memory_manager->get_field_instance(field_ids[0])->get_comp_id();

	for (int i = 0; i < num_fields; i ++) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, memory_manager->check_is_legal_field_instance_id(field_ids[i]) , "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"%s\" contains wrong field instance id (the %dth element of the array is wrong). Please verify the corresponding model code.", API_label, annotation, "field_instance_ids", i+1);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, comp_id == memory_manager->get_field_instance(field_ids[i])->get_comp_id(), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the field instances specified via the parameter \"%s\" must correspond to the same component model while the first field instance corresponds to the component model \"%s\" and the %dth field instance corresponds to the component model \"%s\". Please verify the corresponding model code.", API_label, annotation, "field_instance_ids", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_comp_full_name(), i+1, comp_comm_group_mgt_mgr->get_global_node_of_local_comp(memory_manager->get_field_instance(field_ids[i])->get_comp_id(), false, "")->get_comp_full_name());
	}
	sprintf(str, "registering an Inout_handler with the annotation \"%s\"", annotation);
	synchronize_comp_processes_for_API(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), str, annotation);
	int temp_int = (sampling_timer_id == -1)? 0 : 1;
	check_API_parameter_int(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), NULL, temp_int, "sampling_timer_id", annotation);
	if (sampling_timer_id != -1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, timer_mgr->check_is_legal_timer_id(sampling_timer_id), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"sampling_timer_id\" (currently is 0x%x) is not the legal id of a timer. Please verify the corresponding model code.", API_label, annotation, sampling_timer_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, comp_id == timer_mgr->get_timer(sampling_timer_id)->get_comp_id(), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"sampling_timer_id\" and the parameter \"%s\" do not correspond to the same component model (the parameter \"sampling_timer_id\" corresponds to the component model \"%s\" while \"%s\" corresponds to the component model \"%s\"). Please verify.", API_label, annotation, "field_instance_ids", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(timer_mgr->get_timer(sampling_timer_id)->get_comp_id(), false, "")->get_comp_full_name(), "field_instance_ids", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_comp_full_name());//timer和field_instance必须对应相同的comp_id
		check_API_parameter_timer(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", sampling_timer_id, "sampling_timer_id (the information of the timer)", annotation);
	}
	temp_int = (output_grid_id == -1)? 0 : 1;
	check_API_parameter_int(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", temp_int, "output_grid_id", annotation);
	if (output_grid_id != -1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(output_grid_id), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"output_grid_id\" (currently is 0x%x) is not the legal id of a grid. Please verify the corresponding model code.", API_label, annotation, output_grid_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, comp_id == original_grid_mgr->search_grid_info(output_grid_id)->get_comp_id(), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"output_grid_id\" and teh parameter \"%s\" do not correspond to the same component model (the parameter \"%s\" corresponds to the component model \"%s\"). Please verify.", API_label, annotation, "field_instance_ids", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(original_grid_mgr->search_grid_info(output_grid_id)->get_comp_id(), false, "")->get_comp_full_name(), "field_instance_dis", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_comp_full_name());
		check_API_parameter_comp_or_grid(output_grid_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", "output_grid_id", annotation);
	}
	temp_int = (handler_output_H2D_grid_id == -1)? 0 : 1;
	check_API_parameter_int(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", temp_int, "handler_output_H2D_grid_id", annotation);
	if (handler_output_H2D_grid_id != -1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(handler_output_H2D_grid_id), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"handler_output_H2D_grid_id\" (currently is 0x%x) is not the legal id of a grid. Please verify the corresponding model code.", API_label, annotation, handler_output_H2D_grid_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, comp_id == original_grid_mgr->search_grid_info(handler_output_H2D_grid_id)->get_comp_id(), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"handler_output_H2D_grid_id\" and the parameter \"%s\" do not correspond to the same component model (the parameter \"%s\" corresponds to the component model \"%s\"). Please verify.", API_label, annotation, "field_instance_ids", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(original_grid_mgr->search_grid_info(handler_output_H2D_grid_id)->get_comp_id(), false, "")->get_comp_full_name(), "field_instance_dis", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_comp_full_name());
		check_API_parameter_comp_or_grid(handler_output_H2D_grid_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Output_handler::Output_handler"), "registering an output datamodel handler", "handler_output_H2D_grid_id", annotation);
	}
	temp_int = (handler_output_V1D_grid_id == -1)? 0 : 1;
	check_API_parameter_int(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", temp_int, "handler_output_V1D_grid_id", annotation);
	if (handler_output_V1D_grid_id != -1) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, original_grid_mgr->is_grid_id_legal(handler_output_V1D_grid_id), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"handler_output_V1D_grid_id\" (currently is 0x%x) is not the legal id of a grid. Please verify the corresponding model code.", API_label, annotation, handler_output_V1D_grid_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, comp_id, comp_id == original_grid_mgr->search_grid_info(handler_output_V1D_grid_id)->get_comp_id(), "Error happens when calling the API \"%s\" to register a datamodel handler with the annotation \"%s\": the parameter \"handler_output_V1D_grid_id\" and teh parameter \"%s\" do not correspond to the same component model (the parameter \"%s\" corresponds to the component model \"%s\"). Please verify.", API_label, annotation, "field_instance_ids", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(original_grid_mgr->search_grid_info(handler_output_V1D_grid_id)->get_comp_id(), false, "")->get_comp_full_name(), "field_instance_dis", comp_comm_group_mgt_mgr->get_global_node_of_local_comp(comp_id, false, "")->get_comp_full_name());
		check_API_parameter_comp_or_grid(handler_output_V1D_grid_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Output_handler::Output_handler"), "registering an output datamodel handler", "handler_output_V1D_grid_id", annotation);
	}

	check_API_parameter_int(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Output_handler:Output_handler"), NULL, num_fields, "num_field_instances", annotation);
	sprintf(str, "\"%s\" (the information of the field instances)", "field_instance_ids");
	for (int i = 0; i < num_fields; i ++)
		check_API_parameter_field_instance(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", field_ids[i], str, annotation);
	check_API_parameter_int(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", implicit_or_explicit, "implicit_or_explicit (the tag of executing the handler implicitly or explicitly)", annotation);
	if (is_handler_for_datamodel)
		check_API_parameter_string(comp_id, API_ID_HANDLER_DATAMODEL_OUTPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Inout_handler:Inout_handler"), "registering an output datamodel handler", handler_name, "handler_name", annotation);
}


Inout_datamodel::Inout_datamodel(int host_comp_id, const char *output_datamodel_name, const char *file_name, const char *file_type, int output_timer_id, int file_timer_id, int output_grid_id, int inst_or_aver, const char *float_datatype, const char *integer_datatype, bool datamodel_type_label, const char *annotation) //for output API
{
	int num_stars = 0;
	this->host_comp_id = host_comp_id;
	this->datamodel_type = datamodel_type_label;
	this->datamodel_name = strdup(output_datamodel_name);
	this->annotation = strdup(annotation);
	Datamodel_file_info *current_file_info = new Datamodel_file_info();
	sprintf(this->datamodel_config_dir, "%s/CCPL_dir/datamodel/config",comp_comm_group_mgt_mgr->get_root_working_dir());
	sprintf(this->datamodel_input_data_dir, "%s/CCPL_dir/datamodel/data", comp_comm_group_mgt_mgr->get_root_working_dir());
	sprintf(this->datamodel_output_data_dir, "%s/%s", comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_working_dir(), "output_datamodel");
	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in register_output_datamodel");
	generate_datamodel_files_info(file_name, file_type, num_stars, current_file_info);
	current_file_info->time_format_in_file_name = strdup("YYYY-MM-DD-SSSSS");
	this->file_sets_info.push_back(current_file_info);

	if (output_timer_id == -1)
		io_timer_ids.push_back(timer_mgr->define_timer(host_comp_id, "steps", 1, 0, 0, this->annotation));
	else io_timer_ids.push_back(output_timer_id);
	if (file_timer_id == -1)
		file_timer_ids.push_back(timer_mgr->define_timer(host_comp_id, "steps", 1, 0, 0, this->annotation));
	else file_timer_ids.push_back(file_timer_id);
	//config default_setting
	Fields_output_common_setting *default_setting = new Fields_output_common_setting();
	default_setting->default_operation = inst_or_aver == 0? strdup("inst"): strdup("aver");
	if (words_are_the_same(float_datatype, DATA_TYPE_DOUBLE))
		default_setting->default_float_type = strdup(DATA_TYPE_DOUBLE);
	else if (words_are_the_same(float_datatype, DATA_TYPE_FLOAT))
		default_setting->default_float_type = strdup(DATA_TYPE_FLOAT);
	else if (words_are_the_same(float_datatype, ""))
		default_setting->default_float_type = strdup("");
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "the float_datatype for a datamodel should be \"%s\" or \"%s\", Please check the registration code for datamodel \"%s\"", DATA_TYPE_DOUBLE, DATA_TYPE_FLOAT, output_datamodel_name);
	
	if (words_are_the_same(integer_datatype, DATA_TYPE_SHORT))
		default_setting->default_integer_type = strdup(DATA_TYPE_SHORT);
	else if (words_are_the_same(integer_datatype, DATA_TYPE_LONG))
		default_setting->default_integer_type = strdup(DATA_TYPE_LONG);
	else if (words_are_the_same(integer_datatype, DATA_TYPE_INT))
		default_setting->default_integer_type = strdup(DATA_TYPE_INT);
	else if (words_are_the_same(integer_datatype, ""))
		default_setting->default_integer_type = strdup("");
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "the integer_datatype for a datamodel should be \"%s\", \"%s\" or \"%s\", Please check the registration code for datamodel \"%s\"", DATA_TYPE_INT, DATA_TYPE_SHORT, DATA_TYPE_LONG, output_datamodel_name);

	default_setting->default_output_grid_id = output_grid_id;
	//default_setting->default_output_grid_name = output_grid_id != -1? strdup(original_grid_mgr->search_grid_info(output_grid_id)->get_grid_name()): strdup("");
	if (output_grid_id != -1) {
		default_setting->default_output_grid_name = strdup(original_grid_mgr->search_grid_info(output_grid_id)->get_grid_name());
		Field_grid_info *grid_info = new Field_grid_info();
		grid_info->grid_name = strdup(default_setting->default_output_grid_name);
		grid_info->grid_id = output_grid_id;
		fields_grid_info.push_back(grid_info);
	}
	else default_setting->default_output_grid_name = strdup("");
	default_setting->file_mid_name = strdup("");
	default_setting->field_specification = strdup("default");
	default_settings.push_back(default_setting);
}

Inout_datamodel::Inout_datamodel(int host_comp_id, const char *input_datamodel_name, const char *annotation)//for input datamodel
{
	int line_number, num_nodes = 0;
	TiXmlNode *input_datamodel_node, *sub_node;

	this->host_comp_id = host_comp_id;
	this->datamodel_type = INPUT_DATAMODEL;
	this->datamodel_name = strdup(input_datamodel_name);
	this->annotation = strdup(annotation);

	config_datamodel_configuration_file_info(input_datamodel_name, INPUT_DATAMODEL);
	find_target_root_node_of_datamodel(host_comp_id, input_datamodel_name, input_datamodel_node, "input_datamodel");

	at_most_one_node_of("data_time_series", input_datamodel_node, sub_node, num_nodes);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_nodes == 1, "Error happens when reading the XML configuration file \"%s\" for the datamodel \"%s\": there should be one and only one XML node named \"data_time_series\" under XML node \"input_datamodel\" while there are %d XML nodes detected. Please verify.", XML_file_name, datamodel_name, num_nodes);
	Datamodel_file_info *current_file_info = config_input_datamodel_data_time_series(sub_node);
	//config input_datamodel grids
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("horizontal_grids", input_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"horizontal_grids\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL)//if only scalar fields are needed to be read in
		config_horizontal_grids_for_datamodel(sub_node);
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("vertical_grids", input_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"vertical_grids\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL)
		config_vertical_grids_for_datamodel(sub_node);
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("V3D_grids", input_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"V3D_grids\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL) 
		config_v3d_grids_for_datamodel(sub_node);

	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("input_fields", input_datamodel_node, sub_node, num_nodes) && sub_node != NULL, "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"input_fields\" are found under node \"input_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	config_input_fields_for_datamodel(sub_node);
	initialize_input_datamodel_file_sets_info();
	allocate_surface_fields_for_input_datamodel();
}

void Inout_datamodel::initialize_input_datamodel_file_sets_info()
{
	long times_in_filename[2048], time_value;
	char full_file_name[NAME_STR_SIZE], time_string[32];


	for (int i = 0; i < file_sets_info.size(); i++) {
		int size_times_in_filename = 0;
		Datamodel_file_info *current_file_info = file_sets_info[i];
		int size_time_format = strlen(current_file_info->time_format_in_file_name);
		DIR *dir = opendir(current_file_info->file_dir);
		struct dirent *ent;
		if (size_time_format != 0) {
			while((ent=readdir(dir)) != NULL) {
				int size_prefix = strlen(current_file_info->file_name_prefix);
				if (strncmp(ent->d_name, current_file_info->file_name_prefix, size_prefix) != 0) continue;
				if (strlen(current_file_info->file_name_suffix) != 0 && strncmp(current_file_info->file_name_suffix, ent->d_name+size_prefix+size_time_format, strlen(current_file_info->file_name_suffix)) != 0) continue;
				if (time_string_and_format_match(host_comp_id, ent->d_name+size_prefix, current_file_info->time_format_in_file_name, current_file_info->id_time_format_in_file_name, true, 0, time_value)) {//already expanded time_format in file_name
					times_in_filename[size_times_in_filename] = time_value;
					size_times_in_filename ++;
				}
			}
			do_quick_sort(times_in_filename, (long*)NULL, 0, size_times_in_filename-1);
			for (int j = 0; j < size_times_in_filename; j++)
				current_file_info->file_name_times.push_back(times_in_filename[j]);
			if (words_are_the_same(current_file_info->time_field_setting->specification, "file_field")) {
				for (int j = 0; j < size_times_in_filename; j++) {
					get_time_string_from_time_value_and_format(host_comp_id, current_file_info->file_name_times[j], "YYYYMMDDSSSSS", current_file_info->time_format_in_file_name, time_string);
					sprintf(full_file_name, "%s/%s%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, time_string, current_file_info->file_name_suffix);
					read_time_fields_from_a_data_file(current_file_info, full_file_name, j);
				}
			}
		}
		else {
			sprintf(full_file_name, "%s/%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, current_file_info->file_name_suffix);
			read_time_fields_from_a_data_file(current_file_info, full_file_name, 0);
		}
	}
}


void Inout_datamodel::read_time_fields_from_a_data_file(Datamodel_file_info *current_file_info, char *full_file_name, int file_ind)
{
	int size_time_from_file_field;
	std::vector<int*> time_field_bufs;
	char time_field_data_type[16], time_string[32];
	IO_netcdf *netcdf_file_object = new IO_netcdf(full_file_name, full_file_name, "r", false);
	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in read_time_fields_from_a_data_file");
	bool is_root_proc = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "in read_time_fields_from_a_data_file") == 0;

	for (int i = 0; i < current_file_info->time_field_setting->time_fields.size(); i++) {
		char *time_field_buf;
		Time_field_param *current_time_field = current_file_info->time_field_setting->time_fields[i];
		netcdf_file_object->read_file_field(current_time_field->varname, (void**)(&time_field_buf), &size_time_from_file_field, time_field_data_type, comm, is_root_proc);
		int *time_field_buf2 = new int[size_time_from_file_field];
		transform_datatype_of_arrays(time_field_buf, (char*)time_field_buf2, time_field_data_type, DATA_TYPE_INT, size_time_from_file_field);
		time_field_bufs.push_back(time_field_buf2);
		delete [] time_field_buf;
	}

	for (int i = 0; i < size_time_from_file_field; i ++) {
		long time_value = 0;
		for (int j = current_file_info->time_field_setting->time_fields.size()-1; j >= 0; j --) {
			int time_format_len = strlen(current_file_info->time_field_setting->time_fields[j]->time_format_in_data_file);
			time_value = time_value * (10^time_format_len) + time_field_bufs[j][i];
			if (current_file_info->time_field_setting->id_full_time_format != TIME_FORMAT_YYYYMMDDSSSSS) {
				get_time_string_from_time_value_and_format(time_value, current_file_info->time_field_setting->time_field_full_time_format, time_string);
				time_string_and_format_match(host_comp_id, time_string, current_file_info->time_field_setting->time_field_full_time_format, current_file_info->time_field_setting->id_full_time_format, true, 0, time_value);
			}
		}
		current_file_info->time_fields_info->time_fields.push_back(time_value);
	}
	if (current_file_info->time_fields_info->last_time_global_ind.size() == 0)
		current_file_info->time_fields_info->last_time_global_ind.push_back(size_time_from_file_field-1);
	else current_file_info->time_fields_info->last_time_global_ind.push_back(size_time_from_file_field + current_file_info->time_fields_info->last_time_global_ind[current_file_info->time_fields_info->last_time_global_ind.size()-1]);
	delete netcdf_file_object;
	for (int i = 0; i < current_file_info->time_field_setting->time_fields.size(); i ++) {
		delete [] time_field_bufs[i];
	}
}


void Inout_datamodel::config_input_fields_for_datamodel(TiXmlNode *input_fields_node)
{
	int temp_index;
	std::vector<Field_config_info*> current_set_of_fields;

	for (TiXmlNode *input_field_node = input_fields_node->FirstChild(); input_field_node != NULL; input_field_node = input_field_node->NextSibling()) {
		TiXmlElement *input_field_element = input_field_node->ToElement();
		Field_config_info *input_field_info = new Field_config_info();

		const char *name_in_model_str = strdup(get_XML_attribute(host_comp_id, 80, input_field_element, "name_in_model", XML_file_name, line_number, "The \"name_in_model\" of the output_datamodel field","output datamodel xml file",true));
		const char *grid_name_str = get_XML_attribute(host_comp_id, 80, input_field_element, "grid_name", XML_file_name, line_number, "The \"grid_name\" of the output_datamodel field","output datamodel xml file",true);

		input_field_info->name_in_model = strdup(name_in_model_str);

		char grid_name_str2[NAME_STR_SIZE];
		rename_datamodel_grid(grid_name_str2, datamodel_name, grid_name_str);
		input_field_info->grid_name = strdup(grid_name_str2);
		int grid_id = get_datamodel_grid_id(grid_name_str2, temp_index, true, true);

		TiXmlNode *data_time_series_node = input_field_node->FirstChild();
		if (data_time_series_node != NULL) {
			if (words_are_the_same(data_time_series_node->ToElement()->Value(), "data_time_series"))
				input_field_info->file_set_for_use = config_input_datamodel_data_time_series(data_time_series_node);
			else EXECUTION_REPORT(REPORT_ERROR, -1, false, "Error happens when configuring fields information for datamodel \"%s\": the node under \"input_field\" node shoule be named \"data_time_series\" or there should be no other inner node. Please check the XML configuration file \"%s\" around line_number %d", datamodel_name, XML_file_name, data_time_series_node->Row());
		}
		else input_field_info->file_set_for_use = file_sets_info[0];
		//check if fields with the same name_in_model exists
		for (int i = 0; i < current_set_of_fields.size(); i ++) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !words_are_the_same(name_in_model_str, current_set_of_fields[i]->name_in_model), "Error happens when configuring fields information for datamodel \"%s\": the field \"%s\" has been specified more than once under the same XML node \"input_fields\". Please check the XML configuration file \"%s\" around line_number %d", datamodel_name, name_in_model_str, XML_file_name, input_field_element->Row());
		}
		if (fields_config_info.size() == 0)
			current_set_of_fields.push_back(input_field_info);
		else fields_config_info[0].push_back(input_field_info);
	}
	if (fields_config_info.size() == 0)
		fields_config_info.push_back(current_set_of_fields);
}

Datamodel_file_info* Inout_datamodel::config_input_datamodel_data_time_series(TiXmlNode *data_time_series_node)
{
	bool data_files_configured = false, time_fields_configured = false;
	Datamodel_file_info *current_file_info = new Datamodel_file_info();

	TiXmlElement *data_time_series_element = data_time_series_node->ToElement();
	for (TiXmlNode *sub_node = data_time_series_node->FirstChild(); sub_node != NULL; sub_node = sub_node->NextSibling()) {
		TiXmlElement *sub_element = sub_node->ToElement();
		if (is_XML_setting_on(host_comp_id, sub_element, XML_file_name,"the status of a data_time_series node of an input_datamodel", "input_datamodel XML file")) {
			if (words_are_the_same(sub_element->Value(), "data_files")) {
				EXECUTION_REPORT(REPORT_ERROR, host_comp_id, !data_files_configured, "Error happens when reading the XML configuration file \"%s\" for the datamodel \"%s\": there should be one and only one XML node named \"data_files\" under XML node \"data_time_series\". Please verify.", XML_file_name, datamodel_name);
				config_data_file_node_info(sub_node, current_file_info);
				data_files_configured = true;
			}
			else if (words_are_the_same(sub_element->Value(), "time_fields")) {
				EXECUTION_REPORT(REPORT_ERROR, host_comp_id, !time_fields_configured, "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\": there should be one and only one XML node named \"time_fields\" under XML node \"data_time_series\". Please verify.", XML_file_name, datamodel_name);
				config_time_fields_node_info(sub_node, current_file_info);
				time_fields_configured = true;
			}
			else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\": the node under \"data_time_series\" node should only be named \"data_files\" or \"time_fields\". Please verify the XML file around line %d", XML_file_name, datamodel_name, sub_element->Row());
		}
	}
	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, time_fields_configured && data_files_configured, "Error happens when reading the XML configuration file for datamodel \"%s\": both node \"data_files\" and \"time_fields\" must be configured. Please verify the XML file \"%s\" around line %d", datamodel_name, XML_file_name, data_time_series_element->Row());

	file_sets_info.push_back(current_file_info);
	return current_file_info;
}


void Inout_datamodel::config_time_fields_node_info(TiXmlNode *time_fields_node, Datamodel_file_info *current_file_info)
{
	int line_number;
	TiXmlElement *time_fields_element = time_fields_node->ToElement();
	bool time_fields_configured = false;
	current_file_info->time_field_setting = new Input_time_field_setting();
	current_file_info->time_field_setting->specification = get_XML_attribute(host_comp_id, NAME_STR_SIZE, time_fields_element, "specification", XML_file_name, line_number, "The \"time_filed\" of the datamodel","input datamodel XML configuration file", true);
	//time_point_type
	const char *time_point_type_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, time_fields_element, "time_point_type", XML_file_name, line_number, "The \"time_filed\" of the datamodel","input datamodel XML configuration file", false);
	if (time_point_type_str == NULL) current_file_info->time_field_setting->time_point_type = 0;
	else if (words_are_the_same(time_point_type_str, "start")) current_file_info->time_field_setting->time_point_type = 0;
	else if (words_are_the_same(time_point_type_str, "middle")) current_file_info->time_field_setting->time_point_type = 1;
	else if (words_are_the_same(time_point_type_str, "end")) current_file_info->time_field_setting->time_point_type = 2;
	else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "Error happens when configuring \"time_fields\" for datamodel \"%s\", \"time_point_type\" should be \"start\", \"middle\" or \"end\"", datamodel_name);

	if (words_are_the_same(current_file_info->time_field_setting->specification, "file_field")) {
		int time_format_id_sum = 0, time_format_id_or = 0;
		for (TiXmlNode *time_field_node = time_fields_node->FirstChild(); time_field_node != NULL; time_field_node = time_field_node->NextSibling()) {
			TiXmlElement *time_field_element = time_field_node->ToElement();

			Time_field_param *time_field_info = new Time_field_param();
			time_field_info->varname = get_XML_attribute(host_comp_id, NAME_STR_SIZE, time_fields_element, "variable", XML_file_name, line_number, "The \"time_filed\" of the datamodel","input datamodel XML configuration file", false);
			const char *time_format_in_datafile_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, time_fields_element, "time_format_in_datafile", XML_file_name, line_number, "The \"time_filed\" of the datamodel","input datamodel XML configuration file", false);
			time_field_info->id_time_format_in_data_file = check_time_format(time_format_in_datafile_str, "time_format_in_datafile");
			time_field_info->time_format_in_data_file = strdup(time_format_in_datafile_str);
			current_file_info->time_field_setting->time_fields.push_back(time_field_info);
			time_format_id_sum += time_field_info->id_time_format_in_data_file;
			time_format_id_or |= time_field_info->id_time_format_in_data_file;
		}
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, time_format_id_sum == time_format_id_or, "Error happens in \"time_field\" node configuration in XML file \"%s\", the \"time_format\" of the time_fields overlaps. Please verify.", XML_file_name);
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, check_time_format_id_legal(time_format_id_sum), "Error happens in \"time_field\" node configuration in XML file \"%s\", the \"time_format\" of the time_fields combined is not a legal time_format. Please verify", XML_file_name);
		current_file_info->time_field_setting->id_full_time_format = time_format_id_sum;
		//sort by time_field time_formats from small to large
		for (int i = 0; i < current_file_info->time_field_setting->time_fields.size(); i ++) {
			int min_time_format_id = i;
			for (int j = i+1; j < current_file_info->time_field_setting->time_fields.size(); j ++) {
				if (current_file_info->time_field_setting->time_fields[j] < current_file_info->time_field_setting->time_fields[min_time_format_id])
					min_time_format_id = j;
			}
			Time_field_param *min_ind = current_file_info->time_field_setting->time_fields[min_time_format_id];
			current_file_info->time_field_setting->time_fields[min_time_format_id] = current_file_info->time_field_setting->time_fields[i];
			current_file_info->time_field_setting->time_fields[i] = min_ind;
		}
		for (int i = 0; i < current_file_info->time_field_setting->time_fields.size(); i ++)
			sprintf(current_file_info->time_field_setting->time_field_full_time_format, "%s%s", current_file_info->time_field_setting->time_fields[i]->time_format_in_data_file, current_file_info->time_field_setting->time_field_full_time_format);
	}
	else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(current_file_info->time_field_setting->specification, "file_name"), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\": the node under \"data_time_series\" node should only be named \"data_files\" or \"time_fields\". Please verify the XML file around line %d", XML_file_name, datamodel_name, time_fields_element->Row());
}


void Inout_datamodel::config_datamodel_configuration_file_info(const char *datamodel_name, bool datamodel_type_label)
{
	char datamodel_file_name[NAME_STR_SIZE];
	if (datamodel_type_label)
		sprintf(datamodel_file_name,"output_datamodel_%s.xml",datamodel_name);
	else sprintf(datamodel_file_name,"input_datamodel_%s.xml",datamodel_name);

	sprintf(this->datamodel_config_dir, "%s/CCPL_dir/datamodel/config",comp_comm_group_mgt_mgr->get_root_working_dir());
	sprintf(this->XML_file_name,"%s/%s", datamodel_config_dir, datamodel_file_name);
	sprintf(this->datamodel_input_data_dir, "%s/CCPL_dir/datamodel/data",comp_comm_group_mgt_mgr->get_root_working_dir());
	sprintf(this->datamodel_output_data_dir, "%s/%s",comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_working_dir(), datamodel_name);
}


void Inout_datamodel::find_target_root_node_of_datamodel(int host_comp_id, const char *datamodel_name, TiXmlNode *&datamodel_node, const char *keyword)
{
	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in register_output_datamodel");

	TiXmlDocument *XML_file = open_XML_file_to_read(host_comp_id, this->XML_file_name, comm, false);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, XML_file != NULL, "Can't find the Datamodel configuration file \"%s\" when registering an output datamodel handler at the model code with the annotation \"%s\". Please check.", XML_file_name, annotation);

	TiXmlNode *root_node = get_XML_first_child_of_unique_root(host_comp_id, XML_file_name, XML_file);
	TiXmlElement *root_element = XML_file->FirstChildElement();
	for (; root_node != NULL; root_node = root_node->NextSibling()) {
		if (root_node->Type() != TiXmlNode::TINYXML_ELEMENT)
			continue;
		root_element = root_node->ToElement();
		if (!words_are_the_same(root_element->Value(), keyword)) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "the first child node of the \"root\" node of a datamodel configuration file should be \"%s\", Please Verify configuration file \"%s\"", keyword, XML_file_name);
		}
		else break;
	}

	datamodel_node = root_node;
	TiXmlElement *datamodel_element = root_element;
	//check output_datamodel status
	const char *name_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, datamodel_element, "name", XML_file_name, line_number, "The \"name\" of the datamodel","output datamodel XML configuration file",true);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(name_str, datamodel_name),"The name of the output_datamodel should be \"%s\" but is \"%s\" in the XML configuration file \"%s\". Please verify the XML configuration file around line %d.", datamodel_name, name_str, XML_file_name, datamodel_element->Row());
}

Inout_datamodel::Inout_datamodel(int host_comp_id, const char *output_datamodel_name, bool datamodel_type_label, const char *annotation) //for output datamodel
{
	int line_number, num_data_files_node = 0, num_hg_grids_node = 0, num_vg_grids_node = 0, num_vd_grids_node = 0, num_out_freq_node = 0, num_nodes = 0;
	TiXmlNode *sub_node = NULL, *output_datamodel_node;
	this->host_comp_id = host_comp_id;
	this->datamodel_type = datamodel_type_label;
	this->datamodel_name = strdup(output_datamodel_name);
	this->annotation = strdup(annotation);

	config_datamodel_configuration_file_info(output_datamodel_name, datamodel_type_label);
	find_target_root_node_of_datamodel(host_comp_id, output_datamodel_name, output_datamodel_node, "output_datamodel");

	at_most_one_node_of("data_files", output_datamodel_node, sub_node, num_nodes);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_nodes == 1, "Error happens when reading the XML configuration file \"%s\" for the datamodel \"%s\": there should be one and only one XML node named \"data_files\" under XML node \"output_datamodel\" while there are %d XML nodes detected. Please verify.", XML_file_name, datamodel_name, num_nodes);
	config_data_files_for_datamodel(sub_node);
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("horizontal_grids", output_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"horizontal_grids\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL) 
		config_horizontal_grids_for_datamodel(sub_node);
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("vertical_grids", output_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"vertical_grids\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL) 
		config_vertical_grids_for_datamodel(sub_node);
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("V3D_grids", output_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"V3D_grids\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL) 
		config_v3d_grids_for_datamodel(sub_node);
	EXECUTION_REPORT(REPORT_ERROR, -1, at_most_one_node_of("fields_output_settings", output_datamodel_node, sub_node, num_nodes), "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\", %d nodes named \"fields_output_settings\" are found under node \"output_datamodel\", Please check.", XML_file_name, datamodel_name, num_nodes);
	if (sub_node != NULL) 
		config_field_output_settings_for_datamodel(sub_node);
}


bool at_most_one_node_of(const char *node_name, TiXmlNode *outer_node, TiXmlNode *&inner_node, int &num_nodes) 
{
	inner_node = NULL;
	num_nodes = 0;
	TiXmlNode *sub_node = NULL;
	for (sub_node = outer_node->FirstChild(); sub_node != NULL; sub_node = sub_node->NextSibling()) {
		if (words_are_the_same(sub_node->ToElement()->Value(), node_name)) {
			inner_node = sub_node;
			num_nodes ++;
		}
	}

	return num_nodes <= 1;
}

void Inout_datamodel::config_data_file_node_info(TiXmlNode *data_file_node, Datamodel_file_info *current_file_info)
{
	int num_stars;
	TiXmlElement *data_file_element = data_file_node->ToElement();

	const char *file_names_str = get_XML_attribute(host_comp_id, 80, data_file_element, "file_names", XML_file_name, line_number, "the \"file_names\" of an output_datamodel", "input_datamodel configuration file", true);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, file_names_str[strlen(file_names_str)-1] != '/', "Error happens in the XML configuration file \"%s\": the XML node \"file_name\" should include a right file name, Please check the XML configuration file around line %d.", XML_file_name, data_file_element->Row());
	const char *file_type_str = get_XML_attribute(host_comp_id, 80, data_file_element, "file_type", XML_file_name, line_number, "the \"file_type\" of an input_datamodel", "datamodel configuration file", true);

	generate_datamodel_files_info(file_names_str, file_type_str, num_stars, current_file_info);
	EXECUTION_REPORT(REPORT_ERROR, -1, opendir(current_file_info->file_dir) != NULL, "Error happens when reading the XML configuration file of datamodel \"%s\", the data file directory specified as \"%s\" cannot be found, Please check the XML configuration file \"%s\" around line number %d", datamodel_name, current_file_info->file_dir, XML_file_name, data_file_element->Row());

	const char *time_format_in_file_names = get_XML_attribute(host_comp_id, 80, data_file_element, "time_format_in_filename", XML_file_name, line_number, "the \"time_format\" in the file name of a datamodel", "datamodel configuration file", num_stars>=1);
	if (time_format_in_file_names != NULL) {
		current_file_info->id_time_format_in_file_name = check_time_format(time_format_in_file_names, "time_format_in_file_name");
		current_file_info->time_format_in_file_name = strdup(time_format_in_file_names);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, num_stars == 1,"Error happens when reading the XML configuration file \"%s\" around line number %d, only one * can be specified in data_file names, there are %d detected, Pleas Verify.", XML_file_name, data_file_element->Row(), num_stars);//if there are no * s in filename, do not read time_format_in_filename
	}
	else current_file_info->time_format_in_file_name = strdup("");
}


void Inout_datamodel::config_data_files_for_datamodel(TiXmlNode *data_files_node) {
	int line_number, num_stars=0;
	bool is_data_file_configured = false;

	Datamodel_file_info *current_file_info = new Datamodel_file_info();
	for (TiXmlNode *data_file_node = data_files_node->FirstChild(); data_file_node != NULL; data_file_node = data_file_node->NextSibling()) {
		TiXmlElement *data_file_element = data_file_node->ToElement();
		if (is_XML_setting_on(host_comp_id, data_file_element, XML_file_name,"the status of a data_file node of an output_datamodel", "output_datamodel XML file")) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_file_element->Value(), "data_file"), "Error happens when configuring datamodel \"%s\": any child node of a \"data_files\" node should be named \"data_file\", Please check the XML configuration file \"%s\" around line %d.", datamodel_name, XML_file_name, data_file_element->Row());
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !is_data_file_configured, "Error happens in output_datamodel configuration file \"%s\" for datamodel \"%s\": only one active XML node of \"data_file\" is enabled. Please Verify the XML file around line_number %d.", XML_file_name, datamodel_name, data_file_element->Row());
		}
		else continue;

		is_data_file_configured = true;
		config_data_file_node_info(data_file_node, current_file_info);
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, is_data_file_configured, "Error happens in output_datamodel configuration file\"%s\" for datamodel \"%s\", one \"status\" of the \"data_file\" node must be set as \"on\", now none was found, Please Verify the XML file around line_number %d.", XML_file_name, datamodel_name, data_files_node->ToElement()->Row());
	this->file_sets_info.push_back(current_file_info);
}


void Inout_datamodel::generate_datamodel_files_info(const char *file_names_str, const char *file_type_str, int &num_stars, Datamodel_file_info *current_file_info)
{
	int pos_last_star=-1, id, k;

	if (file_names_str[0] != '/') {
		current_file_info->datamodel_files_dir_name = new char [NAME_STR_SIZE];
		if (datamodel_type == OUTPUT_DATAMODEL)
			sprintf(current_file_info->datamodel_files_dir_name, "%s/%s", datamodel_output_data_dir, file_names_str);
		else sprintf(current_file_info->datamodel_files_dir_name, "%s/%s", datamodel_input_data_dir, file_names_str);
	}
	else current_file_info->datamodel_files_dir_name = strdup(file_names_str);

	for (id = strlen(current_file_info->datamodel_files_dir_name)-1; id >= 0; id --) {
		if (current_file_info->datamodel_files_dir_name[id] == '/')
			break;
	}

	current_file_info->file_dir = strdup(current_file_info->datamodel_files_dir_name);
	current_file_info->file_dir[id+1] = '\0';
	char *file_name = current_file_info->datamodel_files_dir_name + id + 1;
	create_directory(current_file_info->file_dir, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id,""), comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id,"") == 0, false, false);
	current_file_info->file_name_prefix = strdup(file_name);
	num_stars = 0;
	for (k = 0; k < strlen(file_name); k ++) {
		if (file_name[k] == '*') {
			num_stars ++;
			pos_last_star = k;
			current_file_info->file_name_prefix[k]='\0';
		}
	}
	if (pos_last_star != -1)
		current_file_info->file_name_suffix = strdup(file_name+pos_last_star+1);
	else current_file_info->file_name_suffix = strdup("");

	current_file_info->file_type = strdup(file_type_str);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_similar(file_type_str, "netcdf"), "Error happens when configuring the datamodel \"%s\" corresonding to the XML configuration file \"%s\": currently node \"file_type\" can only be specified as \"netcdf\" (case-insensitive).", datamodel_name, XML_file_name);
}


void Inout_datamodel::config_horizontal_grids_for_datamodel(TiXmlNode *hgs_node) {
	int line_number, grid_id;

	for (TiXmlNode *hg_node = hgs_node->FirstChild(); hg_node != NULL; hg_node = hg_node->NextSibling()) {
		TiXmlElement *hg_element = hg_node->ToElement();
		if (!is_XML_setting_on(host_comp_id, hg_element, XML_file_name, "the status of a \"horizontal_grid\" node of an output_datamodel", "output_datamodel xml file"))
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(hg_element->Value(), "horizontal_grid"), "Error happens when reading the XML configuration file for datamodel \"%s\": child nodes of a \"horizontal_grids\" node should be named \"horizontal_grid\". Please verify the configuration file \"%s\" around line %d.", datamodel_name, XML_file_name, hg_element->Row());

		const char *grid_name_str = get_XML_attribute(host_comp_id, 80, hg_element, "grid_name", XML_file_name, line_number, "the \"grid_name\" of an output_datamodel horizontal_grid","datamodel configuration file", true);
		char *grid_name_str2 = new char [NAME_STR_SIZE];
		rename_datamodel_grid(grid_name_str2, datamodel_name, grid_name_str);
		const char *specification_str = get_XML_attribute(host_comp_id, 80, hg_element, "specification", XML_file_name, line_number, "the \"specification\" of an output_datamodel horizontal_grid","datamodel configuration file", true);
		TiXmlNode *entry_node = hg_node->FirstChild();
		TiXmlElement *entry_element = entry_node->ToElement();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(entry_element->Value(), "entry"), "Error happens when reading the XML configuration file for datamodel \"%s\", the first child node of a \"horizontal_grid\" node should be named \"entry\", Please check file \"%s\".", datamodel_name, XML_file_name);
		if (words_are_the_same(specification_str, "CCPL_grid_file"))
			grid_id = config_horizontal_grid_via_CCPL_grid_file(entry_node, grid_name_str, grid_name_str2);
		else if (words_are_the_same(specification_str, "grid_data_file_field"))
			grid_id = config_horizontal_grid_via_grid_data_file_field(entry_node, grid_name_str, grid_name_str2);
		else if (words_are_the_same(specification_str, "uniform_lonlat_grid"))
			grid_id = config_horizontal_grid_via_uniform_lonlat_grid(entry_node, grid_name_str, grid_name_str2);
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "The \"specification\" of a horizontal_grid (currently is \"%s\") in a datamodel configuration file can only be one of \"CCPL_grid_file\", \"grid_data_file_field\", or \"uniform_lonlat_grid\". Please Verify the XML configuration file \"%s\" around line %d.", specification_str, XML_file_name, hg_element->Row());

		Field_grid_info *grid_info = new Field_grid_info();
		grid_info->grid_id = grid_id;
		grid_info->grid_name = strdup(grid_name_str2);
		grid_info->is_sigma = false;
		fields_grid_info.push_back(grid_info);
		delete [] grid_name_str2;
	}
}


int Inout_datamodel::config_horizontal_grid_via_CCPL_grid_file(TiXmlNode *grid_file_entry_node, const char *ori_grid_name, const char *grid_name) 
{
	int line_number, grid_id, ncfile_id;
	char string_annotation[NAME_STR_SIZE];

	const char *CCPL_grid_file_name_str = get_XML_attribute(host_comp_id, 80, grid_file_entry_node->ToElement(), "file_name", XML_file_name, line_number, "the \"file_name\" of an CCPL_grid_file","datamodel configuration file", true);
	char *file_name = new char [NAME_STR_SIZE];
	sprintf(file_name, "%s/%s", datamodel_input_data_dir, CCPL_grid_file_name_str);
	sprintf(string_annotation, "register H2D grid \"%s\" for datamodel \"%s\" via \"CCPL_grid_file\" at the line %d in XML configuration file \"%s\"", ori_grid_name, datamodel_name, grid_file_entry_node->ToElement()->Row(), XML_file_name);
	if (report_error_enabled) {
		int rcode = nc_open(file_name, NC_NOWRITE, &ncfile_id);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, rcode == NC_NOERR, " netcdf file \"%s\" for datamodel \"%s\" \"horizontal_grid\" does not exist or cannot be opened, Please check the XML configuration file \"%s\" around line number %d", file_name, datamodel_name, XML_file_name, grid_file_entry_node->ToElement()->Row());
		rcode = nc_close(ncfile_id);
	}
	grid_id = original_grid_mgr->register_H2D_grid_via_file(host_comp_id, grid_name, file_name, string_annotation);

	delete [] file_name;

	return grid_id;
}


int Inout_datamodel::config_horizontal_grid_via_grid_data_file_field(TiXmlNode *file_field_entry_node, const char *ori_grid_name, const char *grid_name) 
{
	int line_number;
	char full_grid_data_file_name[NAME_STR_SIZE], string_annotation[NAME_STR_SIZE];
	TiXmlElement *file_field_entry_element = file_field_entry_node->ToElement();

	const char *file_name_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "file_name", XML_file_name, line_number, "the \"file_name\" of a grid_data_file_field", "datamodel configuration file", datamodel_type);//true for output
	const char *edge_type_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "edge_type", XML_file_name, line_number, "the \"edge_type\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *coord_unit_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "coord_unit", XML_file_name, line_number, "the \"coord_unit\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *cyclic_or_acyclic_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "cyclic_or_acyclic", XML_file_name, line_number, "the \"cyclic_or_acyclic\" of a grid_data_file_field", "datamodel configuration file",true);
	const char *dim_size1_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "dim_size1", XML_file_name, line_number, "the \"dim_size1\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *dim_size2_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "dim_size2", XML_file_name, line_number, "the \"dim_size2\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *min_lon_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "min_lon", XML_file_name, line_number, "the \"min_lon\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *min_lat_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "min_lat", XML_file_name, line_number, "the \"min_lat\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *max_lon_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "max_lon", XML_file_name, line_number, "the \"max_lon\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *max_lat_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "max_lat", XML_file_name, line_number, "the \"max_lat\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *center_lon_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "center_lon", XML_file_name, line_number, "the \"center_lon\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *center_lat_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "center_lat", XML_file_name, line_number, "the \"center_lat\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *mask_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "mask", XML_file_name, line_number, "the \"mask\" of a grid_data_file_field", "datamodel configuration file", false);
	const char *area_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "area", XML_file_name, line_number, "the \"area\" of a grid_data_file_field", "datamodel configuration file", false);
	const char *vertex_lon_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "vertex_lon", XML_file_name, line_number, "the \"vertex_lon\" of a grid_data_file_field", "datamodel configuration file", false);
	const char *vertex_lat_str = get_XML_attribute(host_comp_id, 80, file_field_entry_element, "vertex_lat", XML_file_name, line_number, "the \"vertex_lat\" of a grid_data_file_field", "datamodel configuration file", false);
	if (file_name_str[0] == '/')
		strcpy(full_grid_data_file_name, file_name_str);
	else sprintf(full_grid_data_file_name, "%s/%s", datamodel_input_data_dir, file_name_str);
	sprintf(string_annotation, "register H2D grid \"%s\" for datamodel \"%s\" via grid data file \"%s\" specified at the line %d in the XML configuration file \"%s\"", ori_grid_name, datamodel_name, full_grid_data_file_name, file_field_entry_element->Row(), XML_file_name);

	return register_common_H2D_grid_for_datamodel(grid_name, full_grid_data_file_name, edge_type_str, coord_unit_str, cyclic_or_acyclic_str, dim_size1_str, dim_size2_str, min_lon_str, min_lat_str, max_lon_str, max_lat_str, center_lon_str, center_lat_str, mask_str, area_str, vertex_lon_str, vertex_lat_str, file_field_entry_element->Row(), string_annotation);
}


int Inout_datamodel::register_common_H2D_grid_for_datamodel(const char *grid_name_str, const char *file_name_str, const char *edge_type_str, const char *coord_unit_str, const char *cyclic_or_acyclic_str, const char *dim_size1_str, const char *dim_size2_str, const char *min_lon_str, const char *min_lat_str, const char *max_lon_str, const char *max_lat_str, const char *center_lon_str, const char *center_lat_str, const char *mask_str, const char *area_str, const char *vertex_lon_str, const char *vertex_lat_str, int line_number, const char *annotation) 
{
	int rcode, ncfile_id, grid_id;
	int size_center_lon=-1, size_center_lat=-1, size_mask=-1, size_area=-1, size_vertex_lon=-1, size_vertex_lat=-1;
	int *mask;
	long dim_lon_size, dim_lat_size, dim_H2D_size, dim_size1, dim_size2;
	char *center_lon, *center_lat, *vertex_lon, *vertex_lat, *area;
	//char min_lon[get_data_type_size(DATA_TYPE_DOUBLE)], min_lat[get_data_type_size(DATA_TYPE_DOUBLE)], max_lon[get_data_type_size(DATA_TYPE_DOUBLE)], max_lat[get_data_type_size(DATA_TYPE_DOUBLE)];
	char min_lon[NAME_STR_SIZE], max_lon[NAME_STR_SIZE], min_lat[NAME_STR_SIZE], max_lat[NAME_STR_SIZE];
	char data_type_for_center_lat[16], data_type_for_center_lon[16], data_type_for_vertex_lon[16], data_type_for_vertex_lat[16], data_type_for_mask[16], data_type_for_area[16];
	char data_type_temp[16];
	double min_lon_value, min_lat_value, max_lon_value, max_lat_value;

	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in register_datamodel_output_handler H2D_grid");
	bool is_root_proc = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "in register_datamodel_output_handler H2D_grid") == 0;

	char *edge_type = strdup(edge_type_str);
	char *cyclic_or_acyclic = strdup(cyclic_or_acyclic_str);
	char *coord_unit = strdup(coord_unit_str);
	IO_netcdf *netcdf_file_object = new IO_netcdf("datamodel_H2D_grid_data", file_name_str, "r", false);

	if (!is_string_a_value(dim_size1_str))
		dim_size1 = netcdf_file_object->get_dimension_size(dim_size1_str, comm, is_root_proc);
	else EXECUTION_REPORT(REPORT_ERROR, -1, is_string_a_nonnegative_integer(dim_size1_str) && sscanf(dim_size1_str, "%d", &dim_size1) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the parameter \"dim_size1\" (\"%s\") is not of integer type. Please verify the XML configuration file \"%s\".", grid_name_str, datamodel_name, dim_size1_str, XML_file_name);
	if (!is_string_a_value(dim_size2_str))
		dim_size2 = netcdf_file_object->get_dimension_size(dim_size2_str, comm, is_root_proc);
	else EXECUTION_REPORT(REPORT_ERROR, -1, is_string_a_nonnegative_integer(dim_size2_str) && sscanf(dim_size2_str, "%d", &dim_size2) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the parameter \"dim_size2\" (\"%s\") is not of integer type. Please verify the XML configuration file \"%s\".", grid_name_str, datamodel_name, dim_size2_str, XML_file_name);

	netcdf_file_object->read_file_field(center_lon_str, (void**)(&center_lon), &size_center_lon, data_type_for_center_lon, comm, is_root_proc);
	netcdf_file_object->read_file_field(center_lat_str, (void**)(&center_lat), &size_center_lat, data_type_for_center_lat, comm, is_root_proc);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, center_lon != NULL, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the longitude value for the center of each grid point (variable \"%s\" in the file) is not specified in the grid file \"%s\", Please check the XML file %s.", grid_name_str, datamodel_name, center_lon_str, file_name_str, XML_file_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, center_lat != NULL, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the latitude value for the center of each grid point (variable \"%s\" in the file) is not specified in the grid file \"%s\", Please check the XML file %s.", grid_name_str, datamodel_name, center_lat_str, file_name_str, XML_file_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_center_lat), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\", in the data file \"%s\", the data type of variables \"%s\" and \"%s\" are not the same, Please check the XML file %s.", grid_name_str, datamodel_name, file_name_str, center_lon_str, center_lat_str, XML_file_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_center_lon, DATA_TYPE_FLOAT) || words_are_the_same(data_type_for_center_lon, DATA_TYPE_DOUBLE), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\", in the data file \"%s\", the data type of variables \"%s\" is not floating-point", grid_name_str, datamodel_name, file_name_str, center_lon_str, XML_file_name);

	if (vertex_lon_str != NULL)
		netcdf_file_object->read_file_field(vertex_lon_str, (void**)(&vertex_lon), &size_vertex_lon, data_type_for_vertex_lon, comm, is_root_proc);	
	else vertex_lon = NULL;
	if (vertex_lat_str != NULL)
		netcdf_file_object->read_file_field(vertex_lat_str, (void**)(&vertex_lat), &size_vertex_lat, data_type_for_vertex_lat, comm, is_root_proc);
	else vertex_lat = NULL;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (vertex_lon_str != NULL && vertex_lon != NULL) || vertex_lon_str == NULL, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the variable \"%s\" corresponding to the XML attribute \"vertex_lon\" cannot be found in grid data file \"%s\", Please check the XML configuration file \"%s\".", grid_name_str, datamodel_name, vertex_lon_str, file_name_str, XML_file_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (vertex_lat_str != NULL && vertex_lat != NULL) || vertex_lat_str == NULL, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the variable \"%s\" corresponding to the XML attribute \"vertex_lat\" cannot be found in grid data file \"%s\", Please check the XML configuration file \"%s\".", grid_name_str, datamodel_name, vertex_lat_str, file_name_str, XML_file_name);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (vertex_lon != NULL && vertex_lat != NULL) || (vertex_lon == NULL && vertex_lat == NULL), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": in the grid data file \"%s\", the longitude and latitude values for each vertex (variables \"%s\" and \"%s\" in the grid data file) of each grid point must be specified/unspecified at the same time. Please verify the XML file \"%s\".", grid_name_str, datamodel_name, file_name_str, vertex_lon_str, vertex_lat_str, XML_file_name);
	if (vertex_lon_str != NULL && vertex_lon != NULL) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_vertex_lon), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": in the data file \"%s\", the data type of varaible \"%s\" is different from the data type of variable \"%s\", Please check the xml file %s.", grid_name_str, datamodel_name, file_name_str, center_lon_str, vertex_lon_str, XML_file_name);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_vertex_lat),"Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": in the data file \"%s\", the data type of varaible \"%s\" is different from the data type of variable \"%s\", Please check the xml file %s.", grid_name_str, datamodel_name, file_name_str, center_lon_str, vertex_lat_str, XML_file_name);
	}

	if (area_str != NULL)
		netcdf_file_object->read_file_field(area_str, (void**)(&area), &size_area, data_type_for_area, comm, is_root_proc);
	else area = NULL;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (area_str != NULL && area != NULL) || area_str == NULL, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the variable \"%s\" for parameter \"area\" cannot be found in grid file \"%s\", Please check the XML file %s.", grid_name_str, datamodel_name, area_str, file_name_str, XML_file_name);
	if (area_str != NULL && area != NULL)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_center_lon, data_type_for_area), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": in the data file \"%s\", the data type of varaible \"%s\" is different from the data type of variable \"%s\", Please check the xml file %s.", grid_name_str, datamodel_name, file_name_str, center_lon_str, vertex_lat_str, XML_file_name);

	if (mask_str != NULL)
		netcdf_file_object->read_file_field(mask_str, (void**)(&mask), &size_mask, data_type_for_mask, comm, is_root_proc);
	else mask = NULL;
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (mask_str != NULL && mask != NULL) || mask_str == NULL, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the variable \"%s\" for parameter \"mask\" cannot be found in grid file \"%s\", Please check the XML file %s.", grid_name_str, datamodel_name, mask_str, file_name_str, XML_file_name);
	if (mask_str != NULL && mask != NULL)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_mask, DATA_TYPE_INT), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": in the data file \"%s\", the data type of variable \"%s\" is not \"integer\", Please check the XML file %s.", grid_name_str, datamodel_name, file_name_str, mask_str, XML_file_name);

	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(min_lon_str, "%lf", &min_lon_value) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value of the XML attribute \"min_lon\" (\"%s\") is not of float or double type. Please verify the XML configuration file \"%s\" around the line %d.", grid_name_str, datamodel_name, min_lon_str, XML_file_name, line_number);
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(max_lon_str, "%lf", &max_lon_value) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value of the XML attribute \"max_lon\" (\"%s\") is not of float or double type. Please verify the XML configuration file \"%s\" around the line %d.", grid_name_str, datamodel_name, max_lon_str, XML_file_name, line_number);
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(min_lat_str, "%lf", &min_lat_value) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value of the XML attribute \"min_lat\" (\"%s\") is not of float or double type. Please verify the XML configuration file \"%s\" around the line %d.", grid_name_str, datamodel_name, min_lat_str, XML_file_name, line_number);
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(max_lat_str, "%lf", &max_lat_value) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value of the XML attribute \"max_lat\" (\"%s\") is not of float or double type. Please verify the XML configuration file \"%s\" around the line %d.", grid_name_str, datamodel_name, max_lat_str, XML_file_name, line_number);
	transform_datatype_of_arrays((char*)(&min_lon_value), min_lon, DATA_TYPE_DOUBLE, data_type_for_center_lon, 1);
	transform_datatype_of_arrays((char*)(&min_lat_value), min_lat, DATA_TYPE_DOUBLE, data_type_for_center_lon, 1);
	transform_datatype_of_arrays((char*)(&max_lon_value), max_lon, DATA_TYPE_DOUBLE, data_type_for_center_lon, 1);
	transform_datatype_of_arrays((char*)(&max_lat_value), max_lat, DATA_TYPE_DOUBLE, data_type_for_center_lon, 1);

	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Starting to register a H2D grid \"%s\" for datamodel \"%s\".", grid_name_str, datamodel_name);
	grid_id = original_grid_mgr->register_H2D_grid_via_global_data(host_comp_id, grid_name_str, edge_type, coord_unit_str, cyclic_or_acyclic, data_type_for_center_lon, dim_size1, dim_size2, size_center_lon, size_center_lat, size_mask, size_area, size_vertex_lon, size_vertex_lat, min_lon, max_lon, min_lat, max_lat, center_lon, center_lat, mask, area, vertex_lon, vertex_lat, annotation, API_ID_HANDLER_DATAMODEL_OUTPUT);
	delete [] center_lon;
	delete [] center_lat;
	if (vertex_lon_str != NULL && vertex_lon != NULL) {
		delete [] vertex_lon;
		delete [] vertex_lat;
	}
	if (mask_str != NULL && mask != NULL)
		delete [] mask;
	if (area_str != NULL && area != NULL)
		delete [] area;
	delete netcdf_file_object;
	delete [] edge_type;
	delete [] coord_unit; 
	delete [] cyclic_or_acyclic;

	return grid_id;
}


int Inout_datamodel::config_horizontal_grid_via_uniform_lonlat_grid(TiXmlNode *uniform_lonlat_grid_entry_node, const char *ori_grid_name, const char *grid_name) 
{
	int line_number, grid_id, num_lons, num_lats;
	char annotation[NAME_STR_SIZE], string_annotation[NAME_STR_SIZE];
	double min_lon = 0.0, max_lon = 0.0, min_lat, max_lat, new_min_lon, center_lon_array_value;
	char min_lon_buf[get_data_type_size(DATA_TYPE_DOUBLE)], min_lat_buf[get_data_type_size(DATA_TYPE_DOUBLE)], max_lon_buf[get_data_type_size(DATA_TYPE_DOUBLE)], max_lat_buf[get_data_type_size(DATA_TYPE_DOUBLE)];
	TiXmlElement *uniform_lonlat_grid_entry_element = uniform_lonlat_grid_entry_node->ToElement();
	const char *coord_unit_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "coord_unit", XML_file_name, line_number, "the \"coord_unit\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *num_lons_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "num_lons", XML_file_name, line_number, "the \"num_lons\" of a uniform_lon_lat_grid", "datamodel configuration file", true);
	const char *num_lats_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "num_lats", XML_file_name, line_number, "the \"num_lats\" of a uniform_lon_lat_grid", "datamodel configuration file", true);
	const char *cyclic_or_acyclic_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "cyclic_or_acyclic", XML_file_name, line_number, "the \"cyclic_or_acyclic\" of a grid_data_file_field", "datamodel configuration file",true);
	bool cyclic = false;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(cyclic_or_acyclic_str, "cyclic") || words_are_the_same(cyclic_or_acyclic_str, "acyclic"), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of XML attribute \"cyclic_or_acyclic\" is not \"cyclic\" or \"acyclic\". Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, cyclic_or_acyclic_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(coord_unit_str, "degrees"), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"coord_unit\" is not \"degrees\". Please check the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, coord_unit_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(num_lons_str, "%d", &num_lons) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"num_lons\" is not of integer type. Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, num_lons_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(num_lats_str, "%d", &num_lats) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"num_lats\" is not of integer type. Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, num_lats_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_lons > 0, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"num_lons\" is not larger than 0. Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, num_lons_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_lats > 0, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"num_lats\" is not larger than 0. Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, num_lats_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());

	if (words_are_the_same(cyclic_or_acyclic_str, "cyclic"))
		cyclic = true;
	const char *min_lon_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "min_lon", XML_file_name, line_number, "the \"min_lon\" of a grid_data_file_field", "datamodel configuration file", !cyclic);
	const char *min_lat_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "min_lat", XML_file_name, line_number, "the \"min_lat\" of a grid_data_file_field", "datamodel configuration file", true);
	const char *max_lon_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "max_lon", XML_file_name, line_number, "the \"max_lon\" of a grid_data_file_field", "datamodel configuration file", !cyclic);
	const char *max_lat_str = get_XML_attribute(host_comp_id, 80, uniform_lonlat_grid_entry_element, "max_lat", XML_file_name, line_number, "the \"max_lat\" of a grid_data_file_field", "datamodel configuration file", true);
	EXECUTION_REPORT(REPORT_ERROR, -1, min_lon_str == NULL || (sscanf(min_lon_str, "%lf", &min_lon) == 1), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"min_lon\" is not of float or double type. Please check the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, min_lon_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT(REPORT_ERROR, -1, max_lon_str == NULL || (sscanf(max_lon_str, "%lf", &max_lon) == 1), "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"max_lon\" is not of float or double type, which is \"%s\", Please check the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, min_lon_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !cyclic || cyclic && (min_lon == 0.0 && max_lon == 360.0) || cyclic && (min_lon_str == NULL && max_lon_str == NULL), "When the XML attribute \"cyclic_or_acyclic\" is specified as \"cyclic\", the XML attributes of \"min_lon\" and \"max_lon\" can only be specified as \"0.0\" and \"360.0\" at the same time (currently are %lf and %lf in the XML configuration file. Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, min_lon, max_lon, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(min_lat_str, "%lf", &min_lat) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"min_lat\" is not of float or double type. Please check the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, min_lon_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(max_lat_str, "%lf", &max_lat) == 1, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"max_lat\" is not of float or double type. Please check the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, min_lon_str, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, min_lat < max_lat, "Error happens when registering an H2D grid \"%s\" for output_datamodel \"%s\": the value (currently is %lf) of the XML attribute \"min_lat\" is larger than the value (currently is %lf) of the XML attribute \"max_lat\". Please verify the XML configuration file \"%s\" around line number %d.", grid_name, datamodel_name, min_lat, max_lat, XML_file_name, uniform_lonlat_grid_entry_element->Row());
	sprintf(annotation, "grid \"%s\" for datamodel \"%s\"", grid_name, datamodel_name);

	if (min_lon_str == NULL || max_lon_str == NULL) {
		min_lon = 0.0;
		max_lon = 360.0;
	}
	if (!cyclic && min_lon == 0.0 && max_lon == 360.0)
		cyclic = !cyclic;
	if (cyclic) 
		num_lons ++;

	double *center_lon_array = new double [num_lons];
	double *center_lat_array = new double [num_lats];

	for (int i = 0; i < num_lats; i ++)
		if (i == 0)
			center_lat_array[i] = min_lat;
		else center_lat_array[i] = ((max_lat - min_lat)/(num_lats-1)) * i + min_lat;
	new_min_lon = min_lon > max_lon? min_lon-360.0 : min_lon;
	for (int j = 0; j < num_lons; j++) {
		if (j == 0)
			center_lon_array_value = new_min_lon;
		else center_lon_array_value = ((max_lon-new_min_lon) / (num_lons-1)) * j + new_min_lon;
		center_lon_array[j] = center_lon_array_value < 0? center_lon_array_value + 360.0 : center_lon_array_value;
	}

	sprintf(string_annotation, "register H2D grid \"%s\" for datamodel \"%s\" via CCPL_grid_file at the line %d in the XML configuration file \"%s\"", ori_grid_name, datamodel_name, uniform_lonlat_grid_entry_element->Row(), XML_file_name);
	if (cyclic) 
		num_lons --;
	grid_id = original_grid_mgr->register_H2D_grid_via_global_data(host_comp_id, grid_name, "LON_LAT", coord_unit_str, (char*)cyclic_or_acyclic_str, DATA_TYPE_DOUBLE, num_lons, num_lats, num_lons, num_lats, -1, -1, -1, -1, (char*)&min_lon, (char*)&max_lon, (char*)&min_lat, (char*)&max_lat, (char*)center_lon_array, (char*)center_lat_array, NULL, NULL, NULL, NULL, string_annotation, API_ID_HANDLER_DATAMODEL_OUTPUT);
	delete [] center_lon_array;
	delete [] center_lat_array;
	return grid_id;
}


void Inout_datamodel::config_vertical_grids_for_datamodel(TiXmlNode *vgs_node)
{
	int line_number, dim_size2, grid_id, grid_type;
	char coord_unit_str[NAME_STR_SIZE], top_value_str[NAME_STR_SIZE];
	char data_type[16];
	//void *coord_values=NULL, *sigma_values=NULL, *Ap_values=NULL, *Bp_values=NULL;
	void *value2, *value3;
	char API_label[32], annotation[256];
	double *temp_value2, *temp_value3, top_value2=0.0;

	for (TiXmlNode *vg_node = vgs_node->FirstChild(); vg_node != NULL; vg_node = vg_node->NextSibling()) {
		TiXmlElement *vg_element = vg_node->ToElement();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(vg_element->Value(), "vertical_grid"), "Error happens when parsing the configuration of datamodel \"%s\": each child node of a \"vertical_grids\" node should be named \"vertical_grid\" but currently is \"%s\", Please verify the XML configuration file \"%s\" around line %d.", datamodel_name, vg_element->Value(), XML_file_name, vg_element->Row());
		if (!is_XML_setting_on(host_comp_id, vg_element, XML_file_name, "the status of a \"vertical_grid\" node of an output_datamodel", "output_datamodel XML configuration file"))
			continue;
		value2 = NULL, value3 = NULL;
		int API_id;
		const char *grid_name_str = get_XML_attribute(host_comp_id, 80, vg_element, "grid_name", XML_file_name, line_number, "the \"grid_name\" of an output_datamodel vertical_grid","datamodel configuration file",true);
		const char *specification_str = get_XML_attribute(host_comp_id, 80, vg_element, "specification", XML_file_name, line_number, "the \"specification\" of an output_datamodel vertical_grid","datamodel configuration file",true);
		const char *grid_type_str = get_XML_attribute(host_comp_id, 80, vg_element, "grid_type", XML_file_name, line_number, "the \"grid_type\" of an output_datamodel vertical_grid", "datamodel configuration file",true);
		MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in register_V1D_grid for output_datamodel");
		bool is_root_proc = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "in register_V1D_grid for output_datamodel") == 0;

		TiXmlNode *entry_node = vg_node->FirstChild();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, entry_node != NULL, "Error happens when reading the XML configuration file \"%s\" for datamodel \"%s\": there is no child node of the XML node \"vertical_grid\" at line %d for specifying the source of grid data. Please verify.", XML_file_name, datamodel_name, vg_node->ToElement()->Row());
		TiXmlElement *entry_element = entry_node->ToElement();
		const char *file_name_str2 = get_XML_attribute(host_comp_id, 80, entry_element, "file_name", XML_file_name, line_number, "the \"file_name\" of a output_datamodel vertical_grid", "datamodel configuration file",false);
		bool using_grid_data_file = (file_name_str2 != NULL);
		char *file_name_str = NULL;
		IO_netcdf *netcdf_file_object = NULL;
		if (using_grid_data_file) {
			file_name_str = new char [NAME_STR_SIZE];
			sprintf(file_name_str, "%s/%s", datamodel_input_data_dir, file_name_str2);
			netcdf_file_object = new IO_netcdf("V1D_grid_data", file_name_str, "r", false);
		}
		EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "start to register V1D grid \"%s\" for output_datamodel", grid_name_str);
		sprintf(annotation,"V1D grid \"%s\" for datamodel \"%s\"", grid_name_str, datamodel_name);
		char *grid_name_str2 = new char [NAME_STR_SIZE];
		rename_datamodel_grid(grid_name_str2, datamodel_name, grid_name_str);
		if (words_are_the_same(specification_str, "grid_data_file_field")) {
			if (words_are_the_same(grid_type_str, "Z")) {
				grid_type = 1;
				config_vertical_z_grid(comm, entry_node, using_grid_data_file, netcdf_file_object, coord_unit_str, &value2, dim_size2, data_type, is_root_proc);
				API_id = API_ID_GRID_MGT_REG_V1D_Z_GRID_VIA_MODEL;
			}
			else if (words_are_the_same(grid_type_str, "SIGMA")) {
				grid_type = 2;
				config_vertical_sigma_grid(comm, entry_node, using_grid_data_file, netcdf_file_object, coord_unit_str, top_value_str, &value2, dim_size2, data_type, is_root_proc);
				API_id = API_ID_GRID_MGT_REG_V1D_SIGMA_GRID_VIA_MODEL;
			}
			else if (words_are_the_same(grid_type_str, "HYBRID")) {
				grid_type = 3;
				config_vertical_hybrid_grid(comm, entry_node, using_grid_data_file, grid_name_str,  netcdf_file_object, coord_unit_str, top_value_str, &value2, &value3, dim_size2, data_type, is_root_proc);
				API_id = API_ID_GRID_MGT_REG_V1D_HYBRID_GRID_VIA_MODEL;
			}
			else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "The \"grid_type\" of a vertical_grid in a datamodel configuration file can only be one of \"Z\", \"SIGMA\", or \"HYBRID\", Please Verify the XML configuration file \"%s\".", XML_file_name);
		}
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "The \"specification\" of a vertical_grid in a datamodel configuration file can only be \"grid_data_file_field\", Please verify the XML configuration file \"%s\".", XML_file_name);
		common_checking_for_grid_registration(host_comp_id, grid_name_str2, coord_unit_str, API_id, this->annotation);

		double *temp_value2 = NULL, *temp_value3 = NULL;
		temp_value2 = new double [dim_size2];
		temp_value3 = new double [dim_size2];
		if (value3 == NULL)
			value3 = value2;
		if (grid_type != 1)
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(top_value_str, "%lf", &top_value2) == 1, "Error happens when registering an vertical grid \"%s\" for output_datamodel \"%s\": the value (\"%s\") of the XML attribute \"top_value\" is not of float or double type. Please check the XML configuration file \"%s\" around line number %d.", grid_name_str, datamodel_name, top_value_str, XML_file_name, entry_node->ToElement()->Row());

		if (words_are_the_same(data_type, DATA_TYPE_FLOAT)) {
			transform_datatype_of_arrays((float*)value2, temp_value2, dim_size2);
			transform_datatype_of_arrays((float*)value3, temp_value3, dim_size2);
		}
		else {
			transform_datatype_of_arrays((double*)value2, temp_value2, dim_size2);
			transform_datatype_of_arrays((double*)value3, temp_value3, dim_size2);
		}

		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, is_array_in_sorting_order(temp_value2, dim_size2) != 0, "Error happens when registering V1D grid \"%s\" for datamodel \"%s\": some arrays of parameters are not in a descending/ascending order. Please check the data file or the XML configuration file \"%s\" around line %d.", grid_name_str, datamodel_name, XML_file_name, vg_element->Row());
		grid_id = original_grid_mgr->register_V1D_grid_via_data(API_id, host_comp_id, grid_name_str2, grid_type, coord_unit_str, dim_size2, top_value2, temp_value2, temp_value3, true, annotation);
		//v1d_grid_ids.push_back(grid_id);
		delete [] temp_value2;
		delete [] temp_value3;
		if (value3 != value2)
			delete [] value3;
		delete [] value2;
		delete [] grid_name_str2;
		if (using_grid_data_file) {
			delete [] file_name_str;
			delete netcdf_file_object;
		}
	}
}


void Inout_datamodel::config_vertical_z_grid(MPI_Comm comm, TiXmlNode *entry_node, bool using_grid_data_file, IO_netcdf *netcdf_file_object, char *coord_unit_str, void **coord_values, int &size_coord_values, char *data_type_for_coord_values, bool is_root_proc) 
{
	int line_number;
	char report_string[NAME_STR_SIZE];
	TiXmlElement *entry_element = entry_node->ToElement();
	strcpy(coord_unit_str,get_XML_attribute(host_comp_id, 80, entry_element, "coord_unit", XML_file_name, line_number, "the \"coord_unit\" of a datamodel vertical_grid","datamodel configuration file", true));
	const char *coord_values_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, entry_element, "coord_values", XML_file_name, line_number, "the \"coord_values\" of a vertical_z_grid", "datamodel configuration file", true);
	if (using_grid_data_file)
		netcdf_file_object->read_file_field(coord_unit_str, coord_values, &size_coord_values, data_type_for_coord_values, comm, is_root_proc);
	else {
		sprintf(report_string, "Among the coord_values, one element that is not of float or double type is found. Please verify the XML configuration file \"%s\" around line number %d.", XML_file_name, entry_element->Row());
		get_values_from_string(coord_values_str, (double**)coord_values, size_coord_values, report_string, ", ", DATA_TYPE_DOUBLE);
		strcpy(data_type_for_coord_values, DATA_TYPE_DOUBLE);
	}
}


void Inout_datamodel::config_vertical_sigma_grid(MPI_Comm comm, TiXmlNode *entry_node, bool using_grid_data_file, IO_netcdf *netcdf_file_object, char *coord_unit_str, char *top_value_str, void **sigma_values, int &size_sigma_values, char *data_type_for_sigma_values, bool is_root_proc) 
{
	int line_number;
	char report_string[NAME_STR_SIZE];
	TiXmlElement *entry_element = entry_node->ToElement();
	strcpy(coord_unit_str, get_XML_attribute(host_comp_id, 80, entry_element, "coord_unit", XML_file_name, line_number, "the \"coord_unit\" of an datamodel vertical_grid", "datamodel configuration file", true));
	strcpy(top_value_str, get_XML_attribute(host_comp_id, 80, entry_element, "top_value", XML_file_name, line_number, "the \"top_value\" of an output_datmaodel vertical_grid", "datamodel configuration file", true));
	const char *sigma_values_str = get_XML_attribute(host_comp_id, 80, entry_element, "sigma_values", XML_file_name, line_number, "the \"sigma_values\" of a datamodel vertical grid", "datamodel configuration file", true);
	if (using_grid_data_file)
		netcdf_file_object->read_file_field(sigma_values_str, sigma_values, &size_sigma_values, data_type_for_sigma_values, comm, is_root_proc);
	else {
		sprintf(report_string, "Among the coord_values, one element that is not of float or double type is found. Please check the XML configuration file \"%s\" around line number %d.", XML_file_name, entry_element->Row());
		get_values_from_string(sigma_values_str, (double**)sigma_values, size_sigma_values, report_string, ", ", DATA_TYPE_DOUBLE);
		strcpy(data_type_for_sigma_values, DATA_TYPE_DOUBLE);
	}
}


void Inout_datamodel::config_vertical_hybrid_grid(MPI_Comm comm, TiXmlNode *entry_node, bool using_grid_data_file, const char *grid_name_str, IO_netcdf *netcdf_file_object, char *coord_unit_str, char *top_value_str, void **Bp_values, void **Ap_values, int &size_Ap, char *data_type_for_ap_values, bool is_root_proc) 
{
	int line_number, size_Bp;
	char data_type_for_bp_values[NAME_STR_SIZE], report_string[NAME_STR_SIZE];
	TiXmlElement *entry_element = entry_node->ToElement();

	strcpy(coord_unit_str, get_XML_attribute(host_comp_id, 80, entry_element, "coord_unit", XML_file_name, line_number, "the \"coord_unit\" of an datamodel vertical_grid", "datamodel configuration file", true));
	strcpy(top_value_str, get_XML_attribute(host_comp_id, 80, entry_element, "top_value", XML_file_name, line_number, "the \"top_value\" of an output_datmaodel vertical_grid", "datamodel configuration file", true));
	const char *Ap_str = get_XML_attribute(host_comp_id, 80, entry_element, "coef_A", XML_file_name, line_number, "the \"Ap\" values of an datamodel vertical grid", "datamodel configuration file", true);
	const char *Bp_str = get_XML_attribute(host_comp_id, 80, entry_element, "coef_B", XML_file_name, line_number, "the \"Bp\" values of an datamodel vertical grid", "datamodel configuration file", true);
	if (using_grid_data_file) {
		netcdf_file_object->read_file_field(Ap_str, Ap_values, &size_Ap, data_type_for_ap_values, comm, is_root_proc);
		netcdf_file_object->read_file_field(Bp_str, Bp_values, &size_Bp, data_type_for_bp_values, comm, is_root_proc);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(data_type_for_ap_values, data_type_for_bp_values), "Error happens when registering a V1D hybrid gird \"%s\" for datamodel \"%s\": the data type (\"%s\") for Ap (variable \"%s\" in the NetCDF data file) is different from the data type (\"%s\") for Bp (variable \"%s\" in the NetCDF data file). Please check the XML configuration file \"%s\" and the NetCDF file \"%s\".", grid_name_str, datamodel_name, data_type_for_ap_values, Ap_str, data_type_for_bp_values, Bp_str, XML_file_name, netcdf_file_object->get_file_name());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, size_Ap == size_Bp, "Error happens when registering a V1D hybrid gird \"%s\" for datamodel \"%s\": the array size (%d) for Ap (variable \"%s\" in the NetCDF data file) is different from the array size (%d) for Bp (variable \"%s\" in the NetCDF data file). Please check the XML configuration file \"%s\" and the NetCDF file \"%s\".", grid_name_str, datamodel_name, size_Ap, Ap_str, size_Bp, Bp_str, XML_file_name, netcdf_file_object->get_file_name());
	}
	else {
		sprintf(report_string, "Among the coord_values, one element that is not of float or double type found. Please check the XML configuration file \"%s\" around line number %d.", XML_file_name, entry_element->Row());
		get_values_from_string(Ap_str, (double**)Ap_values, size_Ap, report_string, ", ", DATA_TYPE_DOUBLE);
		get_values_from_string(Bp_str, (double**)Bp_values, size_Bp, report_string, ", ", DATA_TYPE_DOUBLE);
		strcpy(data_type_for_ap_values, DATA_TYPE_DOUBLE);
		strcpy(data_type_for_bp_values, DATA_TYPE_DOUBLE);
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, size_Ap == size_Bp, "Error happens when registering a V1D hybrid gird \"%s\" for datamodel \"%s\": the array size (%d) for Ap corresponding to the string \"%s\" is different form the array size (%d) for Bp corresponding to the string \"%s\". Please check the XML configuration file \"%s\".", grid_name_str, datamodel_name, size_Ap, Ap_str, size_Bp, Bp_str, XML_file_name);
	}
}


void Inout_datamodel::config_v3d_grids_for_datamodel(TiXmlNode *v3ds_node)
{
	int line_number, v3d_grid_id, v1d_subgrid_id, h2d_subgrid_id;
	const char *surface_field_type_str, *field_name_in_file_str;

	for (TiXmlNode *v3d_node = v3ds_node->FirstChild(); v3d_node != NULL; v3d_node = v3d_node->NextSibling()) {
		TiXmlElement *v3d_element = v3d_node->ToElement();
		if (!is_XML_setting_on(host_comp_id, v3d_element, XML_file_name, "the status of a \"3d_grid\" node of a datamodel", "datamodel xml file"))
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(v3d_element->Value(), "V3D_grid"), "Error happens when configuring datamodel \"%s\": the child node of a \"V3D_grids\" node should be named \"V3D_grid\". Please verify the XML configuration file \"%s\" around line %d.", datamodel_name, XML_file_name, v3d_element->Row());

		const char *grid_name_str = get_XML_attribute(host_comp_id, 80, v3d_element, "grid_name", XML_file_name, line_number, "the \"grid_name\" of a datamodel 3d_grid","datamodel configuration file",true);
		const char *mid_point_grid_name_str = get_XML_attribute(host_comp_id, 80, v3d_element, "mid_point_grid_name", XML_file_name, line_number, "the \"mid_point_grid_name\" of a datamodel vertical_grid","datamodel configuration file",false);
		const char *dimension_order_str = get_XML_attribute(host_comp_id, 80, v3d_element, "dimension_order", XML_file_name, line_number, "the \"dimension_order\" of a datamodel vertical_grid","datamodel configuration file",false);
		if (dimension_order_str != NULL)
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_similar(dimension_order_str, "H2D+V1D") || words_are_similar(dimension_order_str, "V1D+H2D"), "Error happens when configuring datamodel \"%s\": the value of the XML attribute \"dimension_order\" must be \"H2D+V1D\" or \"V1D+H2D\" (currently is \"%s\"). Please verify the XML configuration file \"%s\" around line %d.", datamodel_name, dimension_order_str, XML_file_name, v3d_element->Row());

		TiXmlNode *horizontal_sub_grid_node = NULL, *vertical_sub_grid_node = NULL;
		char horizontal_subgrid_name_str[80], vertical_subgrid_name_str[80];

		for (TiXmlNode *temp_XML_node = v3d_node->FirstChild(); temp_XML_node != NULL; temp_XML_node = temp_XML_node->NextSibling()) {
			TiXmlElement *temp_XML_element = temp_XML_node->ToElement();
			if (words_are_the_same(temp_XML_element->Value(), "horizontal_sub_grid")) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, horizontal_sub_grid_node == NULL, "Error happens when configuring datamodel \"%s\" based on the XML configuration file \"%s\": the XML node \"V3D_grid\" at line %d has more than one child of \"horizontal_sub_grid\" (only one is allowed). Please verify.", datamodel_name, XML_file_name, v3d_element->Row());
				horizontal_sub_grid_node = temp_XML_node;
			}
			else if (words_are_the_same(temp_XML_element->Value(), "vertical_sub_grid")) {
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, vertical_sub_grid_node == NULL, "Error happens when configuring datamodel \"%s\" based on the XML configuration file \"%s\": the XML node \"V3D_grid\" at line %d has more than one child of \"vertical_sub_grid\" (only one is allowed). Please verify.", datamodel_name, XML_file_name, v3d_element->Row());
				vertical_sub_grid_node = temp_XML_node;
			}
			else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "Error happens when configuring datamodel \"%s\" based on the XML configuration file \"%s\": the child node \"%s\" at line %d of the XML node \"V3D_grid\" is not \"horizontal_sub_grid\" or \"vertical_sub_grid\". Please verify.", datamodel_name, XML_file_name, temp_XML_element->Value(), temp_XML_element->Row());
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, horizontal_sub_grid_node != NULL, "Error happens when configuring datamodel \"%s\" based on the XML configuration file \"%s\": the XML node \"V3D_grid\" at line %d must have one child of \"horizontal_sub_grid\". Please verify.", datamodel_name, XML_file_name, v3d_element->Row());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, vertical_sub_grid_node != NULL, "Error happens when configuring datamodel \"%s\" based on the XML configuration file \"%s\": the XML node \"V3D_grid\" at line %d must have one child of \"vertical_sub_grid\". Please verify.", datamodel_name, XML_file_name, v3d_element->Row());
		strcpy(horizontal_subgrid_name_str, get_XML_attribute(host_comp_id, 80, horizontal_sub_grid_node->ToElement(), "grid_name",XML_file_name, line_number, "the \"grid_name\" of the horizontal_subgrid of a datamodel v3d_grid","datamodel configuration file",true));
		strcpy(vertical_subgrid_name_str, get_XML_attribute(host_comp_id, 80, vertical_sub_grid_node->ToElement(), "grid_name",XML_file_name, line_number, "the \"grid_name\" of the vertical_subgrid of a datamodel v3d_grid","datamodel configuration file", true));

		//config surface_field
		TiXmlNode *surface_field_node = vertical_sub_grid_node->FirstChild();
		if (surface_field_node != NULL) {
			TiXmlElement *surface_field_element = surface_field_node->ToElement();
			surface_field_type_str = get_XML_attribute(host_comp_id, 80, surface_field_element, "type", XML_file_name, line_number, "the surface_field \"type\" of a datamodel V3D_grid", "datamodel configuration file",false);
			field_name_in_file_str = get_XML_attribute(host_comp_id, 80, surface_field_element, "field_name_in_file", XML_file_name, line_number, "the surface_field name of a datamodel V3D_grid", "datamodel configuration file", false);
		}
		char *surface_field_name_in_file = field_name_in_file_str != NULL? strdup(field_name_in_file_str): strdup("original_name");
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, (datamodel_type == OUTPUT_DATAMODEL) || (datamodel_type == INPUT_DATAMODEL) && (surface_field_type_str != NULL), "Error happens when registering a V3D_grid \"%s\" for datamodel \"%s\": \"type\" of surface_field need not be specified when DATAMODEL is for OUTPUT, and must be specified when DATAMODEL is for INPUT. Please check the XML configuration file \"%s\" around line %d.", grid_name_str, datamodel_name, XML_file_name, v3d_element->Row());

		if (!words_are_the_same(horizontal_subgrid_name_str, "handler_output_H2D_grid") && !words_are_the_same(vertical_subgrid_name_str, "handler_output_V1D_grid")) {
			char *hg_subgrid_name = new char [NAME_STR_SIZE];
			char *vg_subgrid_name = new char [NAME_STR_SIZE];
			char *vd_grid_name = new char [NAME_STR_SIZE];
			rename_datamodel_grid(hg_subgrid_name, datamodel_name, horizontal_subgrid_name_str);
			rename_datamodel_grid(vg_subgrid_name, datamodel_name, vertical_subgrid_name_str);
			rename_datamodel_grid(vd_grid_name, datamodel_name, grid_name_str);
			h2d_subgrid_id = original_grid_mgr->get_grid_id(host_comp_id, hg_subgrid_name, "get horizontal_subgrid name for datamodel");
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, h2d_subgrid_id != -1, "Error happens when registering a V3D_grid \"%s\" for datamodel \"%s\": the horizontal subgrid with name \"%s\" has not been specified, Please check the XML configuration file \"%s\".", vd_grid_name, datamodel_name, hg_subgrid_name, XML_file_name);
			v1d_subgrid_id = original_grid_mgr->get_grid_id(host_comp_id, vg_subgrid_name, "get vertical_subgrid name for datamodel");
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, v1d_subgrid_id != -1, "Error happens when registering a V3D_grid \"%s\" for datamodel \"%s\": the verticlal subgrid with name \"%s\" has not been specified, Please check the XML configuration file \"%s\".", vd_grid_name, datamodel_name, vg_subgrid_name, XML_file_name);
			v3d_grid_id = register_datamodel_v3d_grid(dimension_order_str, vd_grid_name, h2d_subgrid_id, v1d_subgrid_id, v3d_element->Row());
			if (original_grid_mgr->search_grid_info(v3d_grid_id)->has_sigma_sub_grid())
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, surface_field_node != NULL, "Error happens when configuring \"V3D_grid\" node for datamodel \"%s\", the \"surface_field\" subnode of a V3D_grid must be set when its vertical subgrid is of \"SIGMA\" or \"HYBRID\" type, Please check XML configuration file \"%s\"", datamodel_name, XML_file_name);

			if (mid_point_grid_name_str != NULL)
				register_datamodel_mid_point_v3d_grid(mid_point_grid_name_str, v3d_grid_id);
			if (original_grid_mgr->search_grid_info(v3d_grid_id)->has_sigma_sub_grid()) {
				if (datamodel_type == OUTPUT_DATAMODEL) {
					char API_label[NAME_STR_SIZE];
					get_API_hint(-1, API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD, API_label);
					original_grid_mgr->set_3d_grid_bottom_field(host_comp_id, v3d_grid_id, -1, BOTTOM_FIELD_VARIATION_EXTERNAL, API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD, API_label, this->annotation);
				}
				else {
					V3d_grid_surface_map *new_surface_map = new V3d_grid_surface_map();
					new_surface_map->v3d_grid_id = v3d_grid_id;
					new_surface_map->surface_field_name_in_file = surface_field_name_in_file;

					EXECUTION_REPORT(REPORT_ERROR, host_comp_id, surface_field_node != NULL && field_name_in_file_str != NULL, "Error happens when registering a V3D_grid \"%s\" for datamodel \"%s\": \"field_name_in_file\" of surface_field must be specified when DATAMODEL is for INPUT. Please check the XML configuration file \"%s\" around line %d.", grid_name_str, datamodel_name, XML_file_name, v3d_element->Row());

					TiXmlNode *surface_field_data_time_node = surface_field_node->FirstChild();
					if (surface_field_data_time_node != NULL) {
						TiXmlElement *surface_field_data_time_element = surface_field_data_time_node->ToElement();
						Datamodel_file_info *current_file_info = config_input_datamodel_data_time_series(surface_field_data_time_node);
						
						Field_config_info *input_field_info = new Field_config_info();
						input_field_info->name_in_model = strdup(field_name_in_file_str);
						char surface_field_hg_grid[80];
						rename_datamodel_grid(surface_field_hg_grid, datamodel_name, horizontal_subgrid_name_str);
						input_field_info->grid_name = strdup(surface_field_hg_grid);
						input_field_info->file_set_for_use = current_file_info;
						if (fields_config_info.size() == 0) {
							std::vector<Field_config_info*> current_set_of_fields;
							current_set_of_fields.push_back(input_field_info);
							fields_config_info.push_back(current_set_of_fields);
						}
						else fields_config_info[0].push_back(input_field_info);
						new_surface_map->file_set_for_use = current_file_info;
					}
					v3d_grid_surfaces.push_back(new_surface_map);
				}
			}

			fields_grid_info[fields_grid_info.size()-1]->surface_field_name = strdup(surface_field_name_in_file);
			if (mid_point_grid_name_str != NULL)
				fields_grid_info[fields_grid_info.size()-2]->surface_field_name = strdup(surface_field_name_in_file);
			delete hg_subgrid_name, vg_subgrid_name, vd_grid_name;
		}
		else push_special_v3d_grid_info(dimension_order_str, mid_point_grid_name_str, grid_name_str, horizontal_subgrid_name_str, vertical_subgrid_name_str, surface_field_name_in_file, v3d_element->Row());

		EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "Finished registering V3D_grid for datamodel %s, grid %s", datamodel_name, grid_name_str);
	}
}

void Inout_datamodel::allocate_surface_fields_for_input_datamodel()
{
	char unit_data_type[16], field_data_type[16], field_unit[16], API_label[NAME_STR_SIZE];
	int pio_proc_num, io_proc_stride, io_proc_mark;
	MPI_Comm io_comm;
	get_API_hint(-1, API_ID_GRID_MGT_SET_3D_GRID_EXTERNAL_BOT_FLD, API_label);

#ifdef USE_PARALLEL_IO
	pio_proc_num = datamodel_mgr->get_comp_PIO_proc_setting(host_comp_id, io_proc_stride, io_proc_mark, io_comm);
#endif

	for (int i = 0; i < v3d_grid_surfaces.size(); i++) {
		V3d_grid_surface_map *v3d_grid_surface_map = v3d_grid_surfaces[i];

		const char *example_input_data_file = randomly_match_a_data_file_in_datamodel(v3d_grid_surface_map->file_set_for_use);
		IO_netcdf *random_file_object = new IO_netcdf(example_input_data_file, example_input_data_file, "r", true);

#ifdef USE_PARALLEL_IO
		Decomp_info *new_decomp = decomps_info_mgr->generate_parallel_decomp_for_parallel_IO(original_grid_mgr->search_grid_info(v3d_grid_surface_map->v3d_grid_id), io_proc_stride, pio_proc_num, io_proc_mark);
#else
		Decomp_info *new_decomp = decomps_info_mgr->generate_default_parallel_decomp_serial(original_grid_mgr->search_grid_info(v3d_grid_surface_map->v3d_grid_id));
#endif
		random_file_object->get_field_datatype(v3d_grid_surface_map->surface_field_name_in_file, field_data_type);
		bool find_unit = random_file_object->get_file_field_attribute(v3d_grid_surface_map->surface_field_name_in_file, "unit", field_unit, unit_data_type);
		if (!find_unit) sprintf(field_unit, "%s", "Pa");
		EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "surface_field data_type:%s", field_data_type);
		datamodel_mgr->add_handlers_field_mem_buf_mark();
		Field_mem_info *surface_field_instance = memory_manager->alloc_mem(v3d_grid_surface_map->surface_field_name_in_file, new_decomp->get_decomp_id(), original_grid_mgr->search_grid_info(v3d_grid_surface_map->v3d_grid_id)->get_H2D_sub_grid()->get_grid_id(), BUF_MARK_IO_FIELD_MIRROR ^ datamodel_mgr->get_handlers_field_mem_buf_mark(), field_data_type, field_unit, this->annotation, false, false);
		surface_field_instance->set_usage_tag(REG_FIELD_TAG_NONE);
		original_grid_mgr->set_3d_grid_bottom_field(host_comp_id, v3d_grid_surface_map->v3d_grid_id, surface_field_instance->get_field_instance_id(), BOTTOM_FIELD_VARIATION_DYNAMIC, API_ID_GRID_MGT_SET_3D_GRID_DYN_BOT_FLD, API_label, this->annotation);
		input_surface_mems.push_back(surface_field_instance);
		delete example_input_data_file, random_file_object;
	}
}

void Inout_datamodel::push_special_v3d_grid_info(const char *dimension_order_str, const char *mid_point_grid_name, const char *vd_grid_name, char *h2d_subgrid_name, char *v1d_subgrid_name, char *surface_field_name_in_file, int line_number) 
{
	V3d_grid_info *v3d_grid = new V3d_grid_info();
	v3d_grid->dimension_order = strdup(dimension_order_str == NULL? "H2D+V1D": dimension_order_str);
	v3d_grid->grid_name = strdup(vd_grid_name);
	v3d_grid->h2d_subgrid_name = strdup(h2d_subgrid_name);
	v3d_grid->v1d_subgrid_name = strdup(v1d_subgrid_name);
	v3d_grid->surface_field_name = strdup(surface_field_name_in_file);
	if (mid_point_grid_name != NULL)
		v3d_grid->mid_point_grid_name = strdup(mid_point_grid_name);
	else v3d_grid->mid_point_grid_name = strdup("");
	v3d_grid->line_number = line_number;
	special_v3ds.push_back(v3d_grid);
}


bool Inout_datamodel::is_special_v3d_grid(const char *v3d_grid_name, V3d_grid_info *&v3d_grid_info) {
	v3d_grid_info = NULL;
	bool return_val = false;
	for (int i = 0; i < special_v3ds.size(); i ++) {
		char *new_grid_name1 = new char [NAME_STR_SIZE];
		rename_datamodel_grid(new_grid_name1, datamodel_name, special_v3ds[i]->grid_name);
		char *new_grid_name2 = new char [NAME_STR_SIZE];
		rename_datamodel_grid(new_grid_name2, datamodel_name, special_v3ds[i]->mid_point_grid_name);
		if (words_are_the_same(new_grid_name1, v3d_grid_name) || words_are_the_same(new_grid_name2, v3d_grid_name)) {
			v3d_grid_info = special_v3ds[i];
			return_val = true;
			break;
		}
		delete [] new_grid_name1;
		delete [] new_grid_name2;
	}
	return return_val;
}


void Inout_datamodel::register_datamodel_mid_point_v3d_grid(const char *mid_point_grid_name_str, int v3d_grid_id) 
{
	int mid_3d_grid_id, mid_1d_grid_id;
	char API_label[256];
	Field_grid_info *grid_info = new Field_grid_info();

	get_API_hint(-1, API_ID_GRID_MGT_REG_MID_POINT_GRID, API_label);
	char *mid_point_grid_name = new char [NAME_STR_SIZE];
	rename_datamodel_grid(mid_point_grid_name, datamodel_name, mid_point_grid_name_str);
	original_grid_mgr->register_mid_point_grid(v3d_grid_id, &mid_3d_grid_id, &mid_1d_grid_id, -1, NULL, "register mid point grid for datamodel", API_label);
	original_grid_mgr->search_grid_info(mid_3d_grid_id)->rename_grid_name(mid_point_grid_name);

	grid_info->grid_id = mid_3d_grid_id;
	grid_info->grid_name = strdup(mid_point_grid_name);
	//grid_info->is_sigma = false;//set is_sigma to false to avoid double registering in export_interface
	grid_info->is_sigma = original_grid_mgr->search_grid_info(mid_3d_grid_id)->has_sigma_sub_grid();
	fields_grid_info.push_back(grid_info);
	delete mid_point_grid_name;
}


int Inout_datamodel::register_datamodel_v3d_grid(const char *dimension_order_str, char *vd_grid_name, int h2d_subgrid_id, int v1d_subgrid_id, int line_number) 
{
	int v3d_grid_id = -1;
	Field_grid_info *grid_info = new Field_grid_info();

	if (dimension_order_str == NULL || words_are_similar(dimension_order_str, "H2D+V1D"))
		v3d_grid_id = original_grid_mgr->register_md_grid_via_multi_grids(host_comp_id, vd_grid_name, h2d_subgrid_id, v1d_subgrid_id, -1, -1, NULL, true, "register v3d_grid for datamodel");
	else if (words_are_similar(dimension_order_str, "V1D+H2D"))
		v3d_grid_id = original_grid_mgr->register_md_grid_via_multi_grids(host_comp_id, vd_grid_name, v1d_subgrid_id, h2d_subgrid_id, -1, -1, NULL, true, "register v3d_grid for datamodel");
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Error happens when registering a V3D_grid \"%s\" for datamodel \"%s\": the \"dimension_order\" of a V3D_grid must be one of \"H2D+V1D\" or \"V1D+H2D\" (case-insensitive). Please check the XML configuration file \"%s\" around line %d.", vd_grid_name, datamodel_name, XML_file_name, line_number);

	grid_info->grid_id = v3d_grid_id;
	grid_info->grid_name = strdup(vd_grid_name);
	grid_info->is_sigma = original_grid_mgr->search_grid_info(v3d_grid_id)->has_sigma_sub_grid();
	fields_grid_info.push_back(grid_info);

	return v3d_grid_id;
}


void Inout_datamodel::config_field_output_settings_for_datamodel(TiXmlNode *field_output_settings_node) 
{
	int line_number, fields_output_setting_index = 0;

	for (TiXmlNode *fields_output_setting_node = field_output_settings_node->FirstChild(); fields_output_setting_node != NULL; fields_output_setting_node = fields_output_setting_node->NextSibling()) {
		TiXmlElement *field_output_setting_element = fields_output_setting_node->ToElement();
		if (!is_XML_setting_on(host_comp_id, field_output_setting_element, XML_file_name, "the status of a \"fields_output_setting\" node of an output_datamodel", "output_datamodel xml file"))
			continue;
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(field_output_setting_element->Value(), "fields_output_setting"), "Error happens when configuring datamodel \"%s\": the child node of a \"fields_output_settings\" node should be named \"fields_output_setting\" (currently is \"%s\"). Please check the XML configuration file \"%s\" around line %d.", datamodel_name, field_output_setting_element->Value(), XML_file_name, field_output_setting_element->Row());

		config_fields_output_common_setting(fields_output_setting_node);
		int num_output_time_setting_configed = 0, num_fields_configed = 0;
		for (TiXmlNode *sub_node = fields_output_setting_node->FirstChild(); sub_node != NULL; sub_node = sub_node->NextSibling()) {
			TiXmlElement *sub_element = sub_node->ToElement();
			if (words_are_the_same(sub_element->Value(), "time_setting")) {
				if ((num_output_time_setting_configed == 0))
					config_output_time_setting(sub_node);
				num_output_time_setting_configed ++;
			}
			if (words_are_the_same(sub_element->Value(), "fields")) {
				if ((num_fields_configed == 0))
					config_field_info(sub_node, fields_output_setting_index);
				num_fields_configed ++;
			}
		}
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_output_time_setting_configed == 1, "Error happens when configuring datamodel \"%s\": the XML node \"time_setting\" must be set and can only be set once in an XML node \"fields_output_setting\", Please verify the XML file \"%s\" around line number %d.", datamodel_name, XML_file_name, fields_output_setting_node->ToElement()->Row());
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, num_fields_configed == 1, "Error happens when configuring datamodel \"%s\": the XML node \"fields\" must be set and can only be set once in an XML node \"fields_output_setting\". Please verify the XML file \"%s\" around line number %d.", datamodel_name, XML_file_name, fields_output_setting_node->ToElement()->Row());
		fields_output_setting_index ++;
	}
}


void Inout_datamodel::config_fields_output_common_setting(TiXmlNode *fields_output_setting_node) 
{
	int line_number, grid_index;
	Fields_output_common_setting *default_setting = new Fields_output_common_setting();
	TiXmlElement *fields_output_setting_element = fields_output_setting_node->ToElement();
	const char *file_mid_name = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "file_mid_name", XML_file_name, line_number, "The \"file_mid_name\" of the output_datamodel", "output datamodel XML configuration file",false);
	const char *time_format_in_datafile_str = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "time_format_in_data_file", XML_file_name, line_number, "The \"time_format_in_data_file\" of the output_datamodel","output datamodel xml file",false);
	int id_time_format_in_data_file = -1;

	if (time_format_in_datafile_str != NULL)
		id_time_format_in_data_file = check_time_format(time_format_in_datafile_str, "time_format in output_datamodel data file");
	default_setting->id_time_format_in_data_file = id_time_format_in_data_file;

	const char *field_specification_str = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "field_specification", XML_file_name, line_number, "The \"field_specification\" of the output_datamodel","output datamodel xml file",true);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(field_specification_str, "default") || words_are_the_same(field_specification_str, "fields") || words_are_the_same(field_specification_str, "hybrid"), "Error happens when configuring the XML node \"fields_output_settings\": the value (currently is \"%s\") of the XML attribute \"field_specification\" can only be one of \"default\" \"fields\" or \"hybrid\". Please verify the XML configuration file \"%s\" around line number %d", field_specification_str, XML_file_name, fields_output_setting_element->Row());

	const char *default_operation_str = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "default_operation", XML_file_name, line_number, "The \"default_operation\" of the output_datamodel","output datamodel xml file", true);
	const char *default_output_grid_str = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "default_output_grid", XML_file_name, line_number, "The \"default_v3d_grid\" of the output_datamodel","output datamodel xml file", false);
	const char *default_float_type_str = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "default_float_type", XML_file_name, line_number, "The \"default_float_type\" of the output_datamodel","output datamodel xml file", false);
	const char *default_integer_type_str = get_XML_attribute(host_comp_id, 80, fields_output_setting_element, "default_integer_type", XML_file_name, line_number, "The \"default_integer_type\" of the output_datamodel","output datamodel xml file", false);
	if (file_mid_name != NULL)
		default_setting->file_mid_name = strdup(file_mid_name);
	else default_setting->file_mid_name = strdup("");

	default_setting->field_specification = strdup(field_specification_str);
	default_setting->default_operation = strdup(default_operation_str);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, words_are_the_same(default_operation_str, "inst") || words_are_the_same(default_operation_str, "aver"), "Error happens when configuring the XML node \"fields_output_settings\": the value (currently is \"%s\") of the XML attribute \"default_operation\" can only be one of \"inst\" or \"aver\". Please verify the XML configuration file \"%s\" around line number %d", default_operation_str, XML_file_name, fields_output_setting_element->Row());

	if (default_float_type_str != NULL) {
		if (words_are_the_same(default_float_type_str, DATA_TYPE_FLOAT) || words_are_the_same(default_float_type_str, DATA_TYPE_DOUBLE))
			default_setting->default_float_type = strdup(default_float_type_str);
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Error happens when configuring the XML node \"fields_output_settings\" for datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"default_float_type\" must be one of \"float\" or \"double\". Please verify the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, fields_output_setting_element->Row());
	} 
    else default_setting->default_float_type = strdup("");

	if (default_integer_type_str != NULL) {
		if (words_are_the_same(default_integer_type_str, DATA_TYPE_INT) || words_are_the_same(default_integer_type_str, DATA_TYPE_SHORT) || words_are_the_same(default_integer_type_str, DATA_TYPE_LONG))
			default_setting->default_integer_type = strdup(default_integer_type_str);
		else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false , "Error happens when configuring the XML node \"fields_output_settings\" for datamodel \"%s\": the value (currently is \"%s\") of the XML attribute \"default_integer_type_str\" must be one of \"int\", \"long\" or \"short\". Please verify the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, fields_output_setting_element->Row());
	} 
	else default_setting->default_integer_type = strdup("");

	default_setting->default_output_grid_id = get_datamodel_grid_id(default_output_grid_str, grid_index, false, true);
	char default_output_grid_name1[NAME_STR_SIZE];
	if (default_output_grid_str != NULL) {
		rename_datamodel_grid(default_output_grid_name1, datamodel_name, default_output_grid_str);
		default_setting->default_output_grid_name = strdup(default_output_grid_name1);
	} 
	else default_setting->default_output_grid_name = strdup("");
	default_settings.push_back(default_setting);
}


int Inout_datamodel::get_datamodel_grid_id(const char *grid_name, int &index, bool has_grid_name_prefix, bool report_error) 
{
	int default_grid_id = -1, i;
	char *default_name;
	V3d_grid_info *v3d_grid_info;


	if (grid_name == NULL || words_are_the_same(grid_name, "")) 
		return -1;

	index = -1;
	default_name = new char[NAME_STR_SIZE];
	if (!has_grid_name_prefix)
		rename_datamodel_grid(default_name, datamodel_name, grid_name);
	else strcpy(default_name, grid_name);
	for (i = 0; i < fields_grid_info.size(); i ++) {
		if (words_are_the_same(default_name, fields_grid_info[i]->grid_name)) {
			default_grid_id = fields_grid_info[i]->grid_id;
			index = i;
			break;
		}
	}
	if (report_error)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, default_grid_id != -1 || is_special_v3d_grid(default_name, v3d_grid_info), "Error happens when reading the XML configuration file \"%s\": grid name \"%s\" is invalid. Please verify the XML configruation file.", XML_file_name, grid_name);
	delete [] default_name;
	return default_grid_id;
}


void Inout_datamodel::config_output_time_setting(TiXmlNode *time_setting_node) 
{
	TiXmlElement *time_setting_element = time_setting_node->ToElement();
	int file_freq_count_value, recurse_level = 0;
	const char *file_freq_unit_str = get_XML_attribute(host_comp_id, 80, time_setting_element, "file_freq_unit", XML_file_name, line_number, "the \"file_freq_unit\" of time_setting", "output datamodel xml file", false);
	const char *file_freq_count_str = get_XML_attribute(host_comp_id, 80, time_setting_element, "file_freq_count", XML_file_name, line_number, "the \"file_freq_count\" of time_setting", "output datamodel xml file", false);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (file_freq_unit_str != NULL && file_freq_count_str != NULL) || (file_freq_unit_str == NULL && file_freq_count_str == NULL), "Error happens in configuring a \"time_setting\" node of the datamodel \"%s\": the values of \"file_freq_unit\" and \"file_freq_count\" must be set or unset at the same time. Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, time_setting_element->Row());
	if (file_freq_unit_str != NULL || file_freq_count_str != NULL) {
		EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(file_freq_count_str, "%d", &file_freq_count_value) == 1, "Error happens in configuring a \"time_setting\" node of the datamodel \"%s\": the value (\"%s\") of \"file_freq_count\" is not of integer type. Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, file_freq_count_str, XML_file_name, time_setting_element->Row());
		file_timer_ids.push_back(timer_mgr->define_timer(host_comp_id, file_freq_unit_str, file_freq_count_value, 0, 0, this->annotation));
	}
	else file_timer_ids.push_back(timer_mgr->define_timer(host_comp_id, "steps", 1, 0, 0, this->annotation));

	Coupling_timer *time_setting_timer = config_period_slots_points_node(time_setting_node, NULL, ++recurse_level);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, time_setting_timer->get_num_children() > 0, "Error happens in configuring a \"time_setting\" node of the datamodel \"%s\": the \"time_setting\" node does not have active children nodes that should specify details of the time setting. Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, time_setting_element->Row());

	io_timer_ids.push_back(time_setting_timer->get_timer_id());
}


Coupling_timer *Inout_datamodel::config_period_slots_points_node(TiXmlNode *period_slots_points_node, const char *parent_frequency_unit, int recurse_level)
{
	TiXmlElement *period_slots_points_element = period_slots_points_node->ToElement();
	Coupling_timer *time_slot_period_periodical_timer = NULL;
	const char *current_frequency_unit = NULL;
	char annotation_hint[NAME_STR_SIZE];
	if (words_are_the_same(period_slots_points_element->Value(), "time_period") || words_are_the_same(period_slots_points_element->Value(), "time_slots")) {
		int freq_count_value;
		const char *freq_unit_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, period_slots_points_element, "freq_unit", XML_file_name, line_number, "the \"freq_unit\" of time_slots", "output datamodel xml file", true);
		const char *freq_count_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, period_slots_points_element, "freq_count", XML_file_name, line_number, "the \"freq_count\" of time_slots", "output datamodel xml file", true);
		EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(freq_count_str, "%d", &freq_count_value) == 1, "Error happens when configuring \"freq_count\" of a \"time_period\" or \"time_slots\" node for the datamodel \"%s\": the value (\"%s\") of \"freq_count\" is not of integer type. Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, freq_count_str, XML_file_name, period_slots_points_element->Row());
		current_frequency_unit = strdup(freq_unit_str);
		if (words_are_the_same(period_slots_points_element->Value(), "time_period")) {
			sprintf(annotation_hint, "check the Coupling_timer with freq_unit of \"%s\", freq_count of \"%s\"", freq_unit_str, freq_count_str);
			time_slot_period_periodical_timer = new Coupling_timer(host_comp_id, -1, freq_unit_str, freq_count_value, 0, 0, annotation_hint);
		}
		else {
			const char *time_format_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, period_slots_points_element, "time_format", XML_file_name, line_number, "the \"time_format\" of time_slots", "output datamodel xml file", true);
			const char *slots_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, period_slots_points_element, "slots", XML_file_name, line_number, "the \"slots\" of time_slots", "output datamodel xml file", true);
			std::vector<std::pair<long, long> > slots_or_points;
			get_time_slots_points_from_str(slots_or_points, slots_str, period_slots_points_element, true);
			sprintf(annotation_hint, "check the time_slot with slot values of \"%s\"", slots_str);
			Time_slots_points *new_time_slot = new Time_slots_points(slots_or_points, true, time_format_str, parent_frequency_unit, annotation_hint);
			sprintf(annotation_hint, "check the periodical timer with freq_unit of \"%s\", freq_count of \"%s\"", freq_unit_str, freq_count_str);
			Periodical_timer *time_slot_freq = new Periodical_timer(host_comp_id, -1, freq_unit_str, freq_count_value, 0, 0, annotation_hint);
			sprintf(annotation_hint, "check the Coupling_timer with time_slots values of \"%s\" and periodical timer with freq_unit \"%s\" and freq_count \"%s\"", slots_str, freq_unit_str, freq_count_str);
			time_slot_period_periodical_timer = new Coupling_timer(host_comp_id, -1, new_time_slot, time_slot_freq, annotation_hint);
		}
	}
	else if (words_are_the_same(period_slots_points_element->Value(), "time_points")) {
		const char *time_points_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, period_slots_points_element, "time_points", XML_file_name, line_number, "the \"time_points\" of time_points", "output datamodel xml file", true);
		const char *time_format_str = get_XML_attribute(host_comp_id, NAME_STR_SIZE, period_slots_points_element, "time_format", XML_file_name, line_number, "the \"time_format\" of time_points", "output datamodel XML file", true);

		std::vector<std::pair<long, long> > slots_or_points;
		get_time_slots_points_from_str(slots_or_points, time_points_str, period_slots_points_element, false);
		sprintf(annotation_hint, "check the time_points with point values of \"%s\"", time_points_str);
		Time_slots_points *new_time_point = new Time_slots_points(slots_or_points, false, time_format_str, parent_frequency_unit, annotation_hint);
		sprintf(annotation_hint, "check the Coupling_timer with time_points values of \"%s\"", time_points_str);
		time_slot_period_periodical_timer = new Coupling_timer(host_comp_id, -1, new_time_point, NULL, annotation_hint);
	}
	else if (words_are_the_same(period_slots_points_element->Value(), "time_setting")) {
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, recurse_level == 1, "Error happens when configuring a time setting for the datamodel \"%s\": the \"time_setting\" node should not appear in the inner level (level %d) of time setting node. Please verify the XML configuration file \"%s\" around line %d", datamodel_name, recurse_level, XML_file_name, period_slots_points_element->Row());
		time_slot_period_periodical_timer = timer_mgr->get_timer(timer_mgr->define_timer(host_comp_id, NULL, NULL, "time_setting timer"));
	}
	else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, false, "Error happens when configuring a time setting for the datamodel \"%s\": the attribute of a XML node is incorrect (currently is \"%s\" but not \"time_period\", \"time_slots\", or \"time_points\"). Please verify the XML configuration file \"%s\" around line %d", datamodel_name, period_slots_points_element->Value(), XML_file_name, period_slots_points_element->Row());

	for (TiXmlNode *sub_node = period_slots_points_node->FirstChild(); sub_node != NULL; sub_node = sub_node->NextSibling()) {
		TiXmlElement *sub_element = sub_node->ToElement();
		if (!is_XML_setting_on(host_comp_id, sub_element, XML_file_name, "the status of a \"time_period\" or \"time_slots\" or \"time_points\" node of a datamodel", "datamodel XML file"))
			continue;
		Coupling_timer *child_timer = config_period_slots_points_node(sub_node, current_frequency_unit, ++recurse_level);
		sprintf(annotation_hint, "add child coupling timer for current timer with frequency unit \"%s\"", current_frequency_unit);
		time_slot_period_periodical_timer->add_child_coupling_timer(child_timer, annotation_hint);
	}
	if (current_frequency_unit != NULL)
		delete current_frequency_unit;
	return time_slot_period_periodical_timer;
}


void Inout_datamodel::get_time_slots_points_from_str(std::vector<std::pair<long, long> > &slots_or_points, const char *time_slots_str, TiXmlElement *period_slots_points_element, bool slots_nor_points) {
	int left_pos[NAME_STR_SIZE], right_pos[NAME_STR_SIZE], lsize = 0, rsize = 0, array_buffer_size;
	char str_buffer[NAME_STR_SIZE], report_hint[NAME_STR_SIZE];
	long *array_buffer;
	left_pos[0] = 0;
	right_pos[0] = 0;

	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, strlen(time_slots_str) < NAME_STR_SIZE, "Software error in Inout_datamodel::get_time_slots_points_from_str");

	if (slots_nor_points) {
		for (int i = 0; i < strlen(time_slots_str); i ++) {
			if (time_slots_str[i] == '[')
				left_pos[lsize++] = i;
			else if (time_slots_str[i] == ']')
				right_pos[rsize++] = i;
			else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (time_slots_str[i] >= '0' && time_slots_str[i] <= '9') || time_slots_str[i] == ',' || time_slots_str[i] == ' ', "Error happens in configuring the \"slots\" string (\"%s\") of a \"time_slots\" node of datamodel \"%s\": the character \'%c\' beyond \'[\' , \']\',  \' \' and \',\' is not allowed in the string. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, time_slots_str[i], XML_file_name, period_slots_points_element->Row());
		}
		for (int j = -1; j < lsize || j < rsize; j ++) {
			std::pair<long, long> slot_pair;
			int upper_bound = strlen(time_slots_str);
			int lower_bound = 0;
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, lsize == rsize && (j == -1 || left_pos[j] < right_pos[j]), "Error happens in configuring the \"slots\" string (\"%s\") of a \"time_slots\" node of datamodel \"%s\": \'[\' and \']\' are not correctly paired in the string. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, XML_file_name, period_slots_points_element->Row());
			if (j < lsize -1) 
				upper_bound = left_pos[j+1];
			if (j >= 0)
				lower_bound = right_pos[j] + 1;
			for (int pos = lower_bound; pos < upper_bound; pos ++)
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, time_slots_str[pos] == ' ' || time_slots_str[pos] == ',', "Error happens in configuring the \"slots\" string (\"%s\") of a \"time_slots\" node of datamodel \"%s\": the character \'%c\' is not allowed between \']\' and \'[\' or after the last \']\'. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, time_slots_str[pos], XML_file_name, period_slots_points_element->Row());
			if (j == -1)
				continue;
			strncpy(str_buffer, time_slots_str + left_pos[j] + 1, right_pos[j] - left_pos[j] - 1);
			str_buffer[right_pos[j]-left_pos[j]-1] = '\0';
			sprintf(report_hint, "Error happens in configuring the \"time_points\" string (\"%s\") of a \"time_points\" node of datamodel \"%s\": found special character between \'[\' and \']\'. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, XML_file_name, period_slots_points_element->Row());

			get_values_from_string(str_buffer, &array_buffer, array_buffer_size, report_hint, ", ", DATA_TYPE_LONG);
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, array_buffer_size == 2 || array_buffer[0] <= array_buffer[1], "Error happens in configuring the \"slots\" string (\"%s\") of a \"time_slots\" node of datamodel \"%s\": the sub-string (\"%s\") between a pair of \"[\" and \"]\" does not contain two integers in non-descending order. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, str_buffer, XML_file_name, period_slots_points_element->Row());
			slot_pair.first = array_buffer[0];
			slot_pair.second = array_buffer[1];
			slots_or_points.push_back(slot_pair);
			if (array_buffer != NULL)
				delete [] array_buffer;
		}
	}
	else {
		sprintf(report_hint, "Error happens in configuring the \"time_points\" string (\"%s\") of a \"time_points\" node of datamodel \"%s\": found special character. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, XML_file_name, period_slots_points_element->Row());
		get_values_from_string(time_slots_str, &array_buffer, array_buffer_size, report_hint, ", ", DATA_TYPE_LONG);
		for (int i = 0; i < array_buffer_size; i ++) {
			std::pair<long, long> slot_pair;
			slot_pair.first = array_buffer[i];
			slot_pair.second = NULL;
			slots_or_points.push_back(slot_pair);
		}
		if (array_buffer != NULL)
			delete [] array_buffer;
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, slots_or_points.size() > 0, "Error happens in configuring the \"slots\" string (\"%s\") of a \"time_slots\" node of datamodel \"%s\":the \"slots\" string is an empty string without values. Please verify the XML configuration file \"%s\" around line_number %d", time_slots_str, datamodel_name, str_buffer, XML_file_name, period_slots_points_element->Row());
}

void Inout_datamodel::config_field_info(TiXmlNode *fields_node, int index) 
{
	std::vector<Field_config_info*> current_set_of_fields;
	int temp_index;

	for (TiXmlNode *field_node = fields_node->FirstChild(); field_node != NULL; field_node = field_node->NextSibling()) {
		TiXmlElement *field_element = field_node->ToElement();
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, words_are_the_same(field_element->Value(), "field"), "Error happens when configuring datamodel \"%s\": the child node of a \"fields\" node should be named \"field\", Please check the XML configuration file \"%s\" around line %d.", datamodel_name, XML_file_name, field_element->Row());

		Field_config_info *output_field_info = new Field_config_info();
		const char *name_in_model_str = strdup(get_XML_attribute(host_comp_id, 80, field_element, "name_in_model", XML_file_name, line_number, "The \"name_in_model\" of the datamodel field","output datamodel xml file",true));
		const char *grid_name_str = get_XML_attribute(host_comp_id, 80, field_element, "grid_name", XML_file_name, line_number, "The \"grid_name\" of the datamodel field","datamodel xml file",false);
		const char *name_in_file_str = get_XML_attribute(host_comp_id, 80, field_element, "name_in_file", XML_file_name, line_number, "The \"name_in_file\" of the datamodel field","datamodel xml file",false);
		const char *float_datatype_str = get_XML_attribute(host_comp_id, 80, field_element, "float_type", XML_file_name, line_number, "The \"float_datatype\" of the datamodel field","datamodel xml file",false);
		const char *integer_datatype_str = get_XML_attribute(host_comp_id, 80, field_element, "integer_type", XML_file_name, line_number, "The \"integer_datatype\" of the datamodel field","datamodel xml file",false);
		//const char *operation_str = get_XML_attribute(host_comp_id, 80, field_element, "operation", XML_file_name, line_number, "The \"operation\" of the datamodel field","datamodel xml file",false);
		const char *unit_str = get_XML_attribute(host_comp_id, 80, field_element, "unit", XML_file_name, line_number, "The \"unit\" of the datamodel field","datamodel xml file",false);
		const char *add_offset_str = get_XML_attribute(host_comp_id, 80, field_element, "add_offset", XML_file_name, line_number, "The \"add_offset\" of the datamodel field","datamodel xml file",false);
		const char *scale_factor_str = get_XML_attribute(host_comp_id, 80, field_element, "scale_factor", XML_file_name, line_number, "The \"scale_factor\" of the datamodel field","datamodel xml file",false);
		const char *value_min_bound_str = get_XML_attribute(host_comp_id, 80, field_element, "value_min_bound", XML_file_name, line_number, "The \"value_min_bound\" of the datamodel field","datamodel xml file",false);
		const char *value_max_bound_str = get_XML_attribute(host_comp_id, 80, field_element, "value_max_bound", XML_file_name, line_number, "The \"value_max_bound\" of the datamodel field","datamodel xml file",false);
		const char *min_max_data_type_str = get_XML_attribute(host_comp_id, 80, field_element, "min_max_data_type", XML_file_name, line_number, "The \"min_max_data_type\" of the datamodel field","datamodel xml file",false);

		char grid_name_str2[NAME_STR_SIZE];
		if (grid_name_str != NULL && !words_are_the_same(grid_name_str, "original_grid")) {
			rename_datamodel_grid(grid_name_str2, datamodel_name, grid_name_str);
			output_field_info->grid_name = strdup(grid_name_str2);
			get_datamodel_grid_id(grid_name_str2, temp_index, true, true);
		}
		else if (grid_name_str == NULL && default_settings[index]->default_output_grid_id != -1)
			output_field_info->grid_name = strdup(default_settings[index]->default_output_grid_name);
		else if (grid_name_str != NULL)
			output_field_info->grid_name = strdup("original_grid");
		else output_field_info->grid_name = strdup("");

		output_field_info->name_in_model = strdup(name_in_model_str);
		if (name_in_file_str != NULL)
			output_field_info->name_in_file = strdup(name_in_file_str);
		else output_field_info->name_in_file = strdup(output_field_info->name_in_model);
		if (value_min_bound_str != NULL) {
			EXECUTION_REPORT(REPORT_ERROR, host_comp_id, min_max_data_type_str != NULL, "\"min_max_data_type\" must be set when \"value_min_bound\" is offered in datamodel configuration \"%s\". Please check XML configuration file \"%s\".", datamodel_name, XML_file_name);
			output_field_info->have_min_bound = true;
			memcpy((char*)&(output_field_info->min_bound_value), value_min_bound_str, get_data_type_size(min_max_data_type_str));
		}
		else output_field_info->have_min_bound = false;

		if (value_max_bound_str != NULL) {
			EXECUTION_REPORT(REPORT_ERROR, host_comp_id, min_max_data_type_str != NULL, "\"min_max_data_type\" must be set when \"value_max_bound\" is offered in datamodel configuration \"%s\". Please check XML configuration file \"%s\".", datamodel_name, XML_file_name);
			output_field_info->have_max_bound = true;
			memcpy((char*)&(output_field_info->max_bound_value), value_max_bound_str, get_data_type_size(min_max_data_type_str));
		}
		else output_field_info->have_max_bound = false;

		for (int i = 0; i < current_set_of_fields.size(); i ++) {
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !words_are_the_same(output_field_info->name_in_file, current_set_of_fields[i]->name_in_file), "Error happens when configuring fields information for datamodel \"%s\": the name of the field \"%s\" (\"s\") is the same with another field under the same XML node \"fields_output_setting\", which is not allowed. Please check the XML configuration file \"%s\" around line_number %d", datamodel_name, name_in_model_str, output_field_info->name_in_file, XML_file_name, field_element->Row());
			EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, !words_are_the_same(name_in_model_str, current_set_of_fields[i]->name_in_model), "Error happens when configuring fields information for datamodel \"%s\": the field \"%s\" has been specified more than once under the same XML node \"fields_output_setting\". Please check the XML configuration file \"%s\" around line_number %d", datamodel_name, name_in_model_str, XML_file_name, field_element->Row());
		}

		if (float_datatype_str != NULL) {
			if (words_are_the_same(float_datatype_str, DATA_TYPE_FLOAT) || words_are_the_same(float_datatype_str, DATA_TYPE_DOUBLE))
				output_field_info->float_datatype = strdup(float_datatype_str);
			else EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Error happens when configuring \"field\" node for datamodel \"%s\": parameter \"float_type\" must be one of \"float\" or \"double\". Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, field_element->Row());
		}
		else if (!words_are_the_same(default_settings[index]->default_float_type, ""))
			output_field_info->float_datatype = strdup(default_settings[index]->default_float_type);
		else output_field_info->float_datatype = strdup("");

		if (integer_datatype_str != NULL) {
			if (words_are_the_same(integer_datatype_str, DATA_TYPE_INT) || words_are_the_same(integer_datatype_str, DATA_TYPE_SHORT) || words_are_the_same(integer_datatype_str, DATA_TYPE_LONG))
				output_field_info->integer_datatype = strdup(integer_datatype_str);
			else 
				EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false , "Error happens when configuring \"field\" node for datamodel \"%s\": parameter \"integer_type\" must be one of \"int\", \"short\" or \"long\", Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, field_element->Row());
		}
		else if (!words_are_the_same(default_settings[index]->default_integer_type, ""))
			output_field_info->integer_datatype = strdup(default_settings[index]->default_integer_type);
		else output_field_info->integer_datatype = strdup("");

		if (!words_are_the_same(default_settings[index]->default_operation, ""))
			output_field_info->operation = strdup(default_settings[index]->default_operation);
		else output_field_info->operation = strdup("inst");

		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, (add_offset_str != NULL && scale_factor_str != NULL) || (add_offset_str == NULL && scale_factor_str == NULL), "Error happens when configuring \"field\" node for datamodel \"%s\": parameter \"add_offset\" and \"scale_factor\" must be specified/unspecified at the same time, Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, XML_file_name, field_element->Row());
		if (add_offset_str != NULL) {
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(add_offset_str, "%lf", &(output_field_info->add_offset)) == 1, "Error happens when configuring \"field\" node for datamodel \"%s\": parameter \"add_offset\" is not of float or double type (which is \"%s\"), Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, add_offset_str, XML_file_name, field_element->Row());
			EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(scale_factor_str, "%lf", &(output_field_info->scale_factor)) == 1, "Error happens when configuring \"field\" node for datamodel \"%s\": the parameter \"scale_factor\" is not of float or double type (which is \"%s\"), Please check the XML configuration file \"%s\" around line number %d.", datamodel_name, scale_factor_str, XML_file_name, field_element->Row());
		}
		else {
			output_field_info->add_offset = 0.0;
			output_field_info->scale_factor = 1.0;
		}
		current_set_of_fields.push_back(output_field_info);
	}
	fields_config_info.push_back(current_set_of_fields);
}


Field_grid_info *Inout_datamodel::search_datamodel_field_grid_info(int dst_grid_id) 
{
	for (int i = 0; i < fields_grid_info.size(); i ++) {
		if (fields_grid_info[i]->grid_id == dst_grid_id)
			return fields_grid_info[i];
	}
	return NULL;
}


Inout_datamodel::~Inout_datamodel()
{
	delete datamodel_name;
	delete annotation;

	for (int i = 0; i < file_sets_info.size(); i++)
		delete file_sets_info[i];

	for (int i = 0; i <default_settings.size(); i ++)
		delete default_settings[i];

	for (int i = 0; i < fields_grid_info.size(); i ++)
		delete fields_grid_info[i];

	for (int i = 0; i < special_v3ds.size(); i ++)
		delete special_v3ds[i];
}


Datamodel_mgt::~Datamodel_mgt() 
{
	for (int i = 0; i < output_handlers.size(); i ++)
		delete output_handlers[i];

	for (int i = 0; i < output_datamodels.size(); i ++)
		delete output_datamodels[i];

	for (int i = 0; i < input_datamodels.size(); i ++)
		delete input_datamodels[i];

	for (int i = 0; i < input_instances.size(); i++)
		delete input_instances[i];

	for (int i = 0; i < input_handler_controllers.size(); i++)
		delete input_handler_controllers[i];

	for (int i = 0; i < input_handler_operators.size(); i++)
		delete input_handler_operators[i];
}


void Datamodel_mgt::handle_normal_output(int handler_id, bool bypass_timer, int API_id, const char *handler_annotation) 
{
	Output_handler *output_handler;
	char API_label[NAME_STR_SIZE];
	get_API_hint(-1, API_id, API_label);
	if (!is_legal_output_handler_id(handler_id))
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Error happens when handling a normal output through the API \"%s\": the given interface ID 0x%x is illegal. Please check the model code with the annotation \"%s\"", API_label, handler_id, handler_annotation);
	output_handler = get_output_handler(handler_id);
	EXECUTION_REPORT_LOG(REPORT_LOG, output_handler->get_host_comp_id(), true, "Begin to execute handler \"%s\" (model code annotation is \"%s\")", output_handler->get_handler_name(), handler_annotation);
	output_handler->execute_handler(bypass_timer, API_id, NULL, NULL, handler_annotation);
	EXECUTION_REPORT_LOG(REPORT_LOG, output_handler->get_host_comp_id(), true, "Finish executing handler \"%s\" (model code annotation is \"%s\")", output_handler->get_handler_name(), handler_annotation);
}

void Datamodel_mgt::handle_explicit_output(int handler_id, bool bypass_timer, int API_id, const char *handler_annotation)
{
	Output_handler *output_handler = get_output_handler(handler_id);
	char API_label[NAME_STR_SIZE];
	get_API_hint(-1, API_id, API_label);
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, output_handler->get_host_comp_id(), output_handler->get_handler_implicit_or_explicit(), "Error happens when executing handler \"%s\", Should not use explicit handler API \"%s\" to execute an implicit hanler, Please check.", output_handler->get_handler_name(), API_label);
	handle_normal_output(handler_id, bypass_timer, API_id, handler_annotation);
}


void Datamodel_mgt::handle_implicit_output(const char *handler_annotation) {
	get_all_implicit_handler_ids();
	for (int i = 0; i < implicit_output_handler_ids.size(); i ++) {
		handle_normal_output(implicit_output_handler_ids[i], false, API_ID_HANDLE_NORMAL_IMPLICIT_OUTPUT, handler_annotation);
	}
}


void Datamodel_mgt::get_all_implicit_handler_ids() 
{
	for (int i = 0; i < output_handlers.size(); i ++) {
		if (!output_handlers[i]->get_handler_implicit_or_explicit())//implicit
			implicit_output_handler_ids.push_back(output_handlers[i]->get_handler_id());
	}
}

void Output_handler::execute_handler_average_procedure(bool bypass_timer)
{
	for (int i = 0; i < output_procedures.size(); i++) {
		for (int j = 0; j < output_procedures[i]->outer_level_averaging_algorithm.size(); j ++) {
			if (output_procedures[i]->outer_level_averaging_algorithm[j] != NULL)
				output_procedures[i]->outer_level_averaging_algorithm[j]->run(bypass_timer || timer_mgr->get_timer(output_procedures[i]->interface_timer_id)->is_timer_on());
		}
	}
}


void Output_handler::execute_handler_coupling_procedure(bool bypass_timer, const char *execute_annotation) 
{
	int *field_update_status = new int [num_fields];
	int local_field_index;
	handler_export_interface->execute(bypass_timer, API_ID_INTERFACE_EXECUTE_WITH_ID, field_update_status, num_fields, execute_annotation);
	handler_export_interface->dump_active_coupling_connections();
	delete field_update_status;
	for (int i = 0; i < output_procedures.size(); i ++) {
		if (output_procedures[i]->interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT) {
			std::vector<const char*> imported_fields_name;
			int field_size = output_procedures[i]->interface->get_num_dst_fields();
			int *import_field_update_status = new int [field_size];
			output_procedures[i]->interface->execute(bypass_timer, API_ID_INTERFACE_EXECUTE_WITH_ID, import_field_update_status, field_size, execute_annotation); 
			output_procedures[i]->interface->dump_active_coupling_connections();
			delete [] import_field_update_status;
		}
	}
}

void Output_handler::execute_handler(bool bypass_timer, int API_id, const char *unique_external_field_io_name, void *external_file_object, const char *execute_annotation) 
{
	//if (level == 1)		
		//comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_IO, TIMING_OUTPUT_ALL, -1, "output datamodel");
	if (level == 1 && timer_mgr->get_timer(export_timer_id)->is_timer_on())
		execute_handler_average_procedure(bypass_timer);
	
	if (level == 2)
		execute_handler_coupling_procedure(bypass_timer, execute_annotation);
	//begin to output files
	std::vector<int> default_field_ids;
	int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "output_handler");
	for (int i = 0; i < output_procedures.size(); i ++) {
		//if (output_procedures[i]->interface->get_interface_type() == COUPLING_INTERFACE_MARK_IMPORT) {
#ifdef USE_PARALLEL_IO
			IO_pnetcdf *netcdf_file_object = NULL;
#else
			IO_netcdf *netcdf_file_object = NULL;
#endif
			if (bypass_timer || timer_mgr->get_timer(output_procedures[i]->interface_timer_id)->is_timer_on()) {
				if (external_file_object != NULL) {
#ifdef USE_PARALLEL_IO
					if (io_proc_mark == 1)
						output_procedures[i]->netcdf_file_object = (IO_pnetcdf*)external_file_object;
#else
					if (local_proc_id == 0)
						output_procedures[i]->netcdf_file_object = (IO_netcdf*)external_file_object;
#endif
				}
				else if (bypass_timer || output_procedures[i]->netcdf_file_object == NULL || timer_mgr->get_timer(output_procedures[i]->file_timer_id)->is_timer_on()) {
					char time_string[NAME_STR_SIZE], full_file_name[NAME_STR_SIZE];
					char *file_operation_str = output_procedures[i]->inst_or_aver == 1? strdup("aver"): strdup("inst");
					get_time_string(time_string, inout_datamodel->file_sets_info[0]->time_format_in_file_name, false, NULL, time_mgr->get_current_year(), time_mgr->get_current_month(), time_mgr->get_current_day(), time_mgr->get_current_second());
					if (strlen(output_procedures[i]->file_mid_name) == 0)
						sprintf(full_file_name, "%s/%s.%s.%s%s.nc", inout_datamodel->file_sets_info[0]->file_dir, inout_datamodel->file_sets_info[0]->file_name_prefix, file_operation_str, time_string, inout_datamodel->file_sets_info[0]->file_name_suffix);
					else sprintf(full_file_name, "%s/%s.%s_%s.%s%s.nc", inout_datamodel->file_sets_info[0]->file_dir, inout_datamodel->file_sets_info[0]->file_name_prefix, output_procedures[i]->file_mid_name, file_operation_str, time_string, inout_datamodel->file_sets_info[0]->file_name_suffix);
					EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "full_file_name: %s\n", full_file_name);
					if (output_procedures[i]->netcdf_file_object != NULL)
						delete output_procedures[i]->netcdf_file_object;

#ifdef USE_PARALLEL_IO
					if (io_proc_mark == 1)
						output_procedures[i]->netcdf_file_object = new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, full_file_name, full_file_name, "w", !bypass_timer);
#else
					if (local_proc_id == 0)
						output_procedures[i]->netcdf_file_object = new IO_netcdf(full_file_name, full_file_name, "w", !bypass_timer);
#endif
					delete file_operation_str;
				}

				int local_field_index;
				for (int m = 0; m < output_procedures[i]->fields_name.size(); m++) {
					int field_info_inner_id = -1;
					char *output_field_name;
					Field_mem_info *out_field_instance;
					if (unique_external_field_io_name != NULL)
						output_field_name = strdup(unique_external_field_io_name);
					else if (fields_config_info.size() != 0 && find_field_info_from_datamodel_config(output_procedures[i]->fields_name[m], fields_config_info[i], field_info_inner_id))
						output_field_name = strdup(fields_config_info[i][field_info_inner_id]->name_in_file);
					else output_field_name = strdup(output_procedures[i]->fields_name[m]);
                    netcdf_file_object = output_procedures[i]->netcdf_file_object;

					out_field_instance = memory_manager->get_field_instance(output_procedures[i]->dst_field_instance_ids[m]);
					Field_mem_info *src_field_instance = memory_manager->get_field_instance(output_procedures[i]->src_field_instance_ids[m]);
					out_field_instance->check_field_sum(report_internal_log_enabled, true, "before writing data into a file");
					strcpy(out_field_instance->get_field_data()->get_grid_data_field()->field_name_in_IO_file, output_field_name);
					if (level == 2) {
						comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_performance_timing_mgr()->performance_timing_start(TIMING_TYPE_IO, TIMING_IO_OUTPUT, -1, "output datamodel");
#ifdef USE_PARALLEL_IO
						if (io_proc_mark == 1) {
							EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "about to write file :%s, filed:%s, io_comm:%d", netcdf_file_object->get_file_name(), output_field_name, netcdf_file_object->get_io_comm());
							if (netcdf_file_object->get_io_with_time_info())
								netcdf_file_object->write_grided_data(host_comp_id, out_field_instance, true, time_mgr->get_current_date(), time_mgr->get_current_second(), false);
							else netcdf_file_object->write_grided_data(host_comp_id, out_field_instance, true, -1, -1, true);
							EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, true, "after write file :%s, filed:%s", netcdf_file_object->get_file_name(), output_field_name);
						}
						else if (out_field_instance->get_decomp_id() != -1)
							EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, out_field_instance->get_size_of_field() == 0, "Software error in Output_handler::execute_handler");
#else
						if (local_proc_id == 0) {
							if (netcdf_file_object->get_io_with_time_info())
								netcdf_file_object->write_grided_data(out_field_instance->get_field_data(), true, time_mgr->get_current_date(), time_mgr->get_current_second(), false);
							else netcdf_file_object->write_grided_data(out_field_instance->get_field_data(), true, -1, -1, true);
						}
#endif
						comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_IO, TIMING_IO_OUTPUT, -1, "output datamodel");
					}
					else {
						if (output_procedures[i]->averaging_field_mems[m] != NULL) {
							fields_gather_scatter_mgr->gather_field(output_procedures[i]->averaging_field_mems[m], out_field_instance, (void*)netcdf_file_object, output_field_name);
						}
						else {
                        	fields_gather_scatter_mgr->gather_field(src_field_instance, out_field_instance, (void*)netcdf_file_object, output_field_name);
						}
					}
					delete output_field_name;
				}

			}
		//}
	}
	//if (level == 1)		
		//comp_comm_group_mgt_mgr->search_global_node(host_comp_id)->get_performance_timing_mgr()->performance_timing_stop(TIMING_TYPE_IO, TIMING_OUTPUT_ALL, -1, "output datamodel");
}


bool Output_handler::check_handler_field_instance_match(Field_mem_info *field_instance, Field_mem_info *&restart_field_instance) {
	std::vector<const char*> interface_fields_name;
	int field_local_index;
	output_procedures[0]->interface->get_fields_name(&interface_fields_name);
	Field_mem_info *interface_field = output_procedures[0]->interface->search_registered_field_instance(interface_fields_name[0], field_local_index);
	if (interface_field->get_grid_id() == field_instance->get_grid_id() && interface_field->get_decomp_id() == field_instance->get_decomp_id() && words_are_the_same(interface_field->get_data_type(), field_instance->get_data_type())) {
		restart_field_instance = interface_field;
		return true;
	}
	restart_field_instance = NULL;
	return false;
}

Output_handler *Datamodel_mgt::get_output_handler(int handler_id) 
{
	if (!is_legal_output_handler_id(handler_id))
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Wrong output handler_id:%d, Please check.", handler_id);
	return output_handlers[handler_id&TYPE_ID_SUFFIX_MASK];
}


bool Datamodel_mgt::is_legal_output_handler_id(int handler_id) 
{
	if ((handler_id & TYPE_ID_PREFIX_MASK) != TYPE_OUTPUT_HANDLER_ID_PREFIX)
		return false;
	return (handler_id&TYPE_ID_SUFFIX_MASK) < output_handlers.size();
}


Output_handler::~Output_handler() 
{
	delete datamodel_name, handler_name, annotation;
	delete inout_datamodel;

	for (int i = 0; i < output_procedures.size(); i ++)
		delete output_procedures[i];

	for (int i = 0; i < import_fields_info.size(); i ++)
		delete import_fields_info[i];

	for (int i = 0; i < fields_config_info.size(); i++) {
		for (int j = 0; j < fields_config_info[i].size(); j++)
			delete fields_config_info[i][j];
	}
}


int Datamodel_mgt::get_comp_PIO_proc_setting(int comp_id, int &io_proc_stride, int &io_proc_mark, MPI_Comm &io_comm) 
{
	Comp_comm_group_mgt_node *comp_node = comp_comm_group_mgt_mgr->search_global_node(comp_id);

	comp_node->generate_PIO_comm_info();
	io_proc_stride = comp_node->get_PIO_proc_ID_stride();
	io_proc_mark = comp_node->get_PIO_proc_mark();
	io_comm = comp_node->get_PIO_sub_comm();

	return comp_node->get_num_PIO_procs();
}


void Datamodel_mgt::read_in_field_from_datafile(int field_instance_id, const char *data_file_name, const char *file_type, const char *field_name_in_file, int necessity, int grid_id_for_file, char *value_min_bound, char *value_max_bound, bool have_min_bound, bool have_max_bound, const char *min_max_data_type, int &successful, const char *annotation)
{
	int API_id = API_ID_HANDLER_DATA_FILE_INPUT;
	char API_label[256];

	get_API_hint(-1, API_ID_HANDLER_DATA_FILE_INPUT, API_label);
	char *input_handler_name = strdup("read_in_field_handler_controller");

	Input_handler_controller *new_input_handler_controller = new Input_handler_controller(input_handler_name, NULL, get_next_input_handler_controller_id(), &field_instance_id, 1, 1, -1, file_type, &necessity, NULL, data_file_name, field_name_in_file, grid_id_for_file, value_min_bound, value_max_bound, have_min_bound, have_max_bound, min_max_data_type, annotation);
	new_input_handler_controller->execute_input_handler_controller(annotation);
	//sucessful
	successful = 1;
	delete new_input_handler_controller;
}


Input_handler_controller::Input_handler_controller(const char *handler_name, const char *instance_name, int handler_id, int *field_instance_ids, int num_fields, int field_instance_ids_size, int input_timer_id, const char *file_type, int *necessity, int *field_connected_status, const char *data_file_name, const char *field_name_in_file, int grid_id_for_file, char *value_min_bound, char *value_max_bound, bool have_min_bound, bool have_max_bound, const char *min_max_data_type, const char *annotation)
{
	char full_file_name[NAME_STR_SIZE], datamodel_input_data_dir[NAME_STR_SIZE];
	char tmp_min_bound[8], tmp_max_bound[8];
	char *io_field_name;

	host_comp_id = memory_manager->get_field_instance(field_instance_ids[0])->get_comp_id();
	this->handler_type = INPUT_HANDLER;
	this->input_timer_id = input_timer_id;
	this->time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);
	this->annotation = strdup(annotation);
	this->handler_id = handler_id;

	if (instance_name != NULL) {
		input_instance = datamodel_mgr->search_input_instance(instance_name);
		inout_datamodel = datamodel_mgr->search_input_datamodel(host_comp_id, input_instance->get_datamodel_name());
		this->fields_config_info = inout_datamodel->fields_config_info;
	}
	else {
		inout_datamodel = NULL;
		config_input_field_info(value_min_bound, value_max_bound, have_min_bound, have_max_bound, min_max_data_type);
	}

	common_checking_for_read_in_field_handler(handler_name, host_comp_id, field_instance_ids, num_fields, field_instance_ids_size, data_file_name, field_name_in_file, necessity, grid_id_for_file, annotation);

	//get model_field_data_type and model_field_unit in previlige
	if (data_file_name[0] != '/') {
		sprintf(datamodel_input_data_dir, "%s/CCPL_dir/datamodel/data",comp_comm_group_mgt_mgr->get_root_working_dir());
		sprintf(full_file_name, "%s/%s", datamodel_input_data_dir, data_file_name);
	}
	else strcpy(full_file_name, data_file_name);

#ifdef USE_PARALLEL_IO
	pio_proc_num = datamodel_mgr->get_comp_PIO_proc_setting(host_comp_id, io_proc_stride, io_proc_mark, io_comm);
#endif
	//get data_type in file
#ifdef USE_PARALLEL_IO
	IO_pnetcdf *netcdf_file_object = new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, full_file_name, full_file_name, "r", true);
#else
	IO_netcdf *netcdf_file_object = new IO_netcdf(full_file_name, full_file_name, "r", true);
#endif

	//config field mapping for each field
	for (int i = 0; i < inout_datamodel->input_surface_mems.size(); i++)
		config_input_fields_model_io_info((void*)netcdf_file_object, inout_datamodel->input_surface_mems[i]->get_field_instance_id(), full_file_name, instance_name, grid_id_for_file, field_name_in_file);
	for (int i = 0; i < num_fields; i++) {
		config_input_fields_model_io_info((void*)netcdf_file_object, field_instance_ids[i], full_file_name, instance_name, grid_id_for_file, field_name_in_file);
	}
	delete netcdf_file_object;
	initialize_input_file_time_info_structure(data_file_name);
}


void Input_handler_controller::config_input_fields_model_io_info(void *example_netcdf_file_object, int model_field_instance_id, const char *full_file_name, const char *instance_name, int grid_id_for_file, const char *field_name_in_file)
{
	char unit_data_type[16], io_field_data_type[16], io_field_unit[16];
	Decomp_info *new_decomp = NULL;
	int io_field_decomp_id = -1, io_field_grid_id = -1, grid_index = -1, field_index_in_xml = -1;

#ifdef USE_PARALLEL_IO
	IO_pnetcdf *netcdf_file_object = (IO_pnetcdf*) example_netcdf_file_object;
#else
	IO_netcdf *netcdf_file_object = (IO_netcdf*) example_netcdf_file_object;
#endif

	Import_field_info *input_field_info = new Import_field_info();
	input_field_info->model_field_instance_id = model_field_instance_id;
	Field_mem_info *model_field_instance = memory_manager->get_field_instance(model_field_instance_id);

	//get io_field_grid_id & io_field_decomp_id
	if (instance_name == NULL) {
		if (grid_id_for_file != -1)
			io_field_grid_id = grid_id_for_file;
		else io_field_grid_id = model_field_instance->get_grid_id();
	}
	else {
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, find_field_info_from_datamodel_config(model_field_instance->get_field_name(), fields_config_info[0], field_index_in_xml), "Can't find field %s in XML configuration file \"%s\" with the annotation \"%s\"", model_field_instance->get_field_name(), inout_datamodel->get_XML_file_name());
		EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "field_index_in_xml:%d", field_index_in_xml);
		io_field_grid_id = inout_datamodel->get_datamodel_grid_id(fields_config_info[0][field_index_in_xml]->grid_name, grid_index, true, true);
	}

#ifdef USE_PARALLEL_IO
	new_decomp = decomps_info_mgr->generate_parallel_decomp_for_parallel_IO(original_grid_mgr->search_grid_info(io_field_grid_id), io_proc_stride, pio_proc_num, io_proc_mark);
#else
	new_decomp = decomps_info_mgr->generate_default_parallel_decomp_serial(original_grid_mgr->search_grid_info(io_field_grid_id));
#endif

	if (new_decomp != NULL)
		io_field_decomp_id = new_decomp->get_decomp_id();

	if (io_field_decomp_id != -1) {
		int new_io_field_decomp_id = decomps_info_mgr->generate_empty_decomp(io_field_decomp_id);
		io_field_decomp_id = new_io_field_decomp_id;
	}
	if (instance_name != NULL)
		input_field_info->field_name_in_file = strdup(input_instance->get_field_name_in_file(model_field_instance->get_field_name()));
	else {
		if (field_name_in_file != NULL)
			input_field_info->field_name_in_file = strdup(field_name_in_file);
		else input_field_info->field_name_in_file = strdup(model_field_instance->get_field_name());//delete
	}

	netcdf_file_object->get_field_datatype(input_field_info->field_name_in_file, io_field_data_type);
#ifdef USE_PARALLEL_IO
	bool find_unit = netcdf_file_object->get_file_field_attribute(input_field_info->field_name_in_file, "unit", io_field_unit, unit_data_type, false);
#else
	bool find_unit = netcdf_file_object->get_file_field_attribute(input_field_info->field_name_in_file, "unit", io_field_unit, unit_data_type);
#endif
	if (!find_unit) sprintf(io_field_unit, "%s", model_field_instance->get_unit());
	EXECUTION_REPORT_LOG(REPORT_LOG, host_comp_id, !find_unit, "Did not find attribute \"unit\" for field \"%s\" in data file \"%s\". Using field unit in model instead", input_field_info->field_name_in_file, full_file_name);

	//allocate null io_field_instance
	datamodel_mgr->add_handlers_field_mem_buf_mark();
	Field_mem_info *io_field_instance = memory_manager->alloc_mem(model_field_instance->get_field_name(), io_field_decomp_id, io_field_grid_id, BUF_MARK_IO_FIELD_MIRROR ^ datamodel_mgr->get_handlers_field_mem_buf_mark(), io_field_data_type, io_field_unit, annotation, false, false);
	io_field_instance->set_usage_tag(REG_FIELD_TAG_NONE);
	input_field_info->io_field_instance_id = io_field_instance->get_field_instance_id();
	input_field_info->import_field_operation = NULL;
	if (instance_name != NULL) {
		input_field_info->file_set_for_use = fields_config_info[0][field_index_in_xml]->file_set_for_use;
		for (int j = 0; j < 2; j ++) {
			int new_buf_mark = datamodel_mgr->get_handlers_field_mem_buf_mark() | ((((get_input_handler_controller_id() + 16)&TYPE_ID_SUFFIX_MASK)) << 12);
			input_field_info->averaging_field_mems.push_back(memory_manager->alloc_mem(model_field_instance, BUF_MARK_AVERAGED_INNER, new_buf_mark, NULL, false, false));
		}
	}
	input_fields_info.push_back(input_field_info);
}


void Input_handler_controller::config_input_field_info(char *value_min_bound, char *value_max_bound, bool have_min_bound, bool have_max_bound, const char *min_max_data_type)
{
	std::vector<Field_config_info*> input_fields_config_info;
	Field_config_info *input_field_config_info = new Field_config_info();
	input_fields_config_info.push_back(input_field_config_info);
	fields_config_info.push_back(input_fields_config_info);

	input_instance = NULL;

	if (have_min_bound) {
		memcpy((char*)&(input_field_config_info->min_bound_value), value_min_bound, get_data_type_size(min_max_data_type));
		check_API_parameter_double(host_comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in Input_handler_controller::Input_handler_controller::"), "registering an Input_handler_controller", input_field_config_info->min_bound_value, "value_min_bound", annotation);
	}
	if (have_max_bound) {
		memcpy((char*)&(input_field_config_info->max_bound_value), value_max_bound, get_data_type_size(min_max_data_type));
		check_API_parameter_double(host_comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in Input_handler_controller::Input_handler_controller"), "registering an Input_handler_controller", input_field_config_info->max_bound_value, "value_max_bound", annotation);
	}
}

void Input_handler_controller::execute_input_handler_controller(const char *execute_annotation) {
	if (input_instance != NULL) {
		config_next_input_file_time_info(execute_annotation);
	}
	for (int i = 0; i < input_fields_info.size(); i++) {
		Field_mem_info *model_field_instance = memory_manager->get_field_instance(input_fields_info[i]->model_field_instance_id);
		Field_mem_info *io_field_instance = memory_manager->get_field_instance(input_fields_info[i]->io_field_instance_id);
		if (input_fields_info[i]->file_set_for_use != NULL) {
			//update values in averaging_field_mems
			if (input_fields_info[i]->file_set_for_use->move_flag) {
				memory_manager->copy_field_data_values(input_fields_info[i]->averaging_field_mems[0], input_fields_info[i]->averaging_field_mems[1]);
			}
			for (int j = 0; j < 2; j ++) {
				if (input_fields_info[i]->file_set_for_use->next_input_file_times[j]->read_flag) {
					fields_gather_scatter_mgr->scatter_field(input_fields_info[i]->averaging_field_mems[j], io_field_instance, input_fields_info[i]->file_set_for_use->next_input_file_times[j], input_fields_info[i]->field_name_in_file);
				}
			}
			//set value for model_field_instance
			if (input_fields_info[i]->file_set_for_use->perfect_match)
				memory_manager->copy_field_data_values(model_field_instance, input_fields_info[i]->averaging_field_mems[0]);
			else interpolate_read_in_value_on_target_time(model_field_instance, input_fields_info[i]);
		}
		else {//readin_single
			fields_gather_scatter_mgr->scatter_field(model_field_instance, io_field_instance, next_input_file_times[0], input_fields_info[i]->field_name_in_file);
		}
	}
	//check_min_max_bound_for_readin_field();
}


void Input_handler_controller::interpolate_read_in_value_on_target_time(Field_mem_info *model_field_instance, Import_field_info *input_field_info) {
	int delta1, delta2, delta, left_time, right_time;
	long end_point_value;
	std::vector<float> fracs;

	Datamodel_file_info *current_file_info = input_field_info->file_set_for_use;
	int left_ind = current_file_info->next_input_file_times[0]->global_ind, right_ind = current_file_info->next_input_file_times[1]->global_ind;
	
	if (words_are_the_same(current_file_info->time_field_setting->specification,"file_field")) {
		left_time = current_file_info->time_fields_info->time_fields[left_ind];
		right_time = current_file_info->time_fields_info->time_fields[right_ind];
	}
	else {
		left_time = current_file_info->file_name_times[left_ind];
		right_time = current_file_info->file_name_times[right_ind];
	}
	if (right_ind == 0)
		while (file_time_value > right_time) right_time += period_time_value_add_count;

	delta1 = file_time_value - left_time;
	delta2 = right_time - file_time_value;
	delta = right_time - left_time;

	fracs.push_back((float)delta2/delta);
	fracs.push_back((float)delta1/delta);
	if (words_are_the_same(model_field_instance->get_data_type(), DATA_TYPE_FLOAT))
		linear_time_interpolation<float>(model_field_instance, fracs, input_field_info->averaging_field_mems);
	else if (words_are_the_same(model_field_instance->get_data_type(), DATA_TYPE_DOUBLE))
		linear_time_interpolation<double>(model_field_instance, fracs, input_field_info->averaging_field_mems);
	else if (words_are_the_same(model_field_instance->get_data_type(), DATA_TYPE_INT))
		linear_time_interpolation<int>(model_field_instance, fracs, input_field_info->averaging_field_mems);
	else if (words_are_the_same(model_field_instance->get_data_type(), DATA_TYPE_SHORT))
		linear_time_interpolation<short>(model_field_instance, fracs, input_field_info->averaging_field_mems);
	else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "");
}


template<typename T> void Input_handler_controller::linear_time_interpolation(Field_mem_info *target_field_mem, std::vector<float> fracs, std::vector<Field_mem_info*> averaging_field_mems)
{
	int num_chunks = averaging_field_mems[0]->get_num_chunks();
	for (int j = 0; j < num_chunks; j++) {
		T*data_buf = (T*) target_field_mem->get_chunk_buf(j);
		for (int k = 0; k < target_field_mem->get_chunk_data_buf_size(j); k++) {
			data_buf[k] = 0;
			for (int i = 0; i < averaging_field_mems.size(); i++) {
				T *averaging_field_value = (T*) averaging_field_mems[i]->get_chunk_buf(j);
				data_buf[k] += averaging_field_value[k] * fracs[i];
			}
		}
	}
}


void Input_handler_controller::initialize_input_file_time_info_structure(const char *data_file_name)
{
	if (input_instance == NULL)
		config_read_in_single_file_time(data_file_name);
	else {
		if (input_instance->offset_info != NULL) calculate_unit = input_instance->offset_info->id_offset_unit;
		else if (input_instance->period_info != NULL) calculate_unit = input_instance->period_info->id_period_unit;
		else calculate_unit = 0;

		switch (calculate_unit) {//YYYYMMDDSSSSS
			case 1 : 
				period_time_value_add_count = 10^9;
				break;
			case 2 : 
				period_time_value_add_count = 10^7;
				break;
			case 3 : 
				period_time_value_add_count = 10^5;
				break;
			case 4 : 
				period_time_value_add_count = 1;
				break;
			default: 0;
		}

		for (int i = 0; i < input_instance->input_datamodel->file_sets_info.size(); i++) {
			Datamodel_file_info *current_file_info = input_instance->input_datamodel->file_sets_info[i];
			current_file_info->next_input_file_times.push_back(new Input_file_time_info());
			current_file_info->next_input_file_times.push_back(new Input_file_time_info());

			if (strlen(current_file_info->time_format_in_file_name) == 0) {
				EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(current_file_info->time_field_setting->specification, "file_field"), "time_fields specification must be set as \"file_field\" when only on file is specified in \"data_files\" node. Please verify XML configuration file \"%s\"", inout_datamodel->get_XML_file_name());
				sprintf(current_file_info->next_input_file_times[0]->full_file_name, "%s/%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, current_file_info->file_name_suffix);
				sprintf(current_file_info->next_input_file_times[1]->full_file_name, "%s/%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, current_file_info->file_name_suffix);
			}
		}
	}
}


void Input_handler_controller::config_next_input_file_time_info(const char *execute_annotation)
{
	//function: 确定文件名、时次（文件名/时次：0个，1个或2个）
	int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "Input_handler_controller");
	config_readin_multiple_file_time_infos();
}


void Input_handler_controller::config_read_in_single_file_time(const char *data_file_name)
{
	Input_file_time_info *file_time_info = new Input_file_time_info();
	char datamodel_input_data_dir[NAME_STR_SIZE];

	if (data_file_name[0] != '/') {
		sprintf(datamodel_input_data_dir, "%s/CCPL_dir/datamodel/data",comp_comm_group_mgt_mgr->get_root_working_dir());
		sprintf(file_time_info->full_file_name, "%s/%s", datamodel_input_data_dir, data_file_name);
	}
	else strcpy(file_time_info->full_file_name, data_file_name);
#ifdef USE_PARALLEL_IO
	file_time_info->netcdf_file_object = new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, file_time_info->full_file_name, file_time_info->full_file_name, "r", false);
#else
	file_time_info->netcdf_file_object = new IO_netcdf(file_time_info->full_file_name, file_time_info->full_file_name, "r", false);
#endif
	next_input_file_times.push_back(file_time_info);
}


void Input_handler_controller::config_readin_multiple_file_time_infos()
{
	//(model_current_time - model_start_time + period_start_time + offset)%period_count+period_start_time
	//at most two next_input_file_times, one for each readin file
	int calculate_count = 0, period_add_count = 0;
	int period_start_year = 0, period_start_month = 0, period_start_day = 0, period_start_second = 0;
	int file_current_year = 0, file_current_month = 0, file_current_day = 0, file_current_second = 0;
	int elapsed_year = 0, elapsed_month = 0, elapsed_day = 0, elapsed_second = 0;
	this->time_mgr = components_time_mgrs->get_time_mgr(host_comp_id);

	if (input_instance->offset_info == NULL && input_instance->period_info == NULL) {
		file_current_year = time_mgr->get_current_year();
		file_current_month = time_mgr->get_current_month();
		file_current_day = time_mgr->get_current_day();
		file_current_second = time_mgr->get_current_second();
	}
	else {
		if (calculate_unit == 3 || calculate_unit == 4) {
			time_mgr->get_elapsed_days_from_start_date(&elapsed_day, &elapsed_second);
			calculate_count = calculate_unit == 3? elapsed_day: elapsed_day*3600*24+elapsed_second;
		}
		else {//1, 2
			get_num_elapsed_time(host_comp_id, elapsed_year, elapsed_month, elapsed_day, elapsed_second);
			EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "elapsed_year:%d, elapsed_month:%d, elapsed_day:%d, elapsed_second:%d", elapsed_year, elapsed_month, elapsed_day, elapsed_second);
			calculate_count = generate_time_on_target_unit(elapsed_year, elapsed_month, elapsed_day, elapsed_second, calculate_unit);
		}

		if (input_instance->offset_info != NULL)
			calculate_count += input_instance->offset_info->offset_count;

		if (input_instance->period_info != NULL) {
			Period_setting *period_set = input_instance->period_info;
			get_time_from_string(host_comp_id, period_set->period_data_start_time, period_set->period_time_format, period_set->id_period_time_format, false, 0, period_start_year, period_start_month, period_start_day, period_start_second);
			period_add_count = calculate_count%period_set->period_count;
			add_time_offset(host_comp_id, period_start_year, period_start_month, period_start_day, period_start_second, period_add_count, calculate_unit);
			if (calculate_unit != 4)
				add_time_offset(host_comp_id, period_start_year, period_start_month, period_start_day, period_start_second, elapsed_second, 4);
			if (calculate_unit == 2 || calculate_unit == 1)
				add_time_offset(host_comp_id, period_start_year, period_start_month, period_start_day, period_start_second, elapsed_day, 3);
			if (calculate_unit == 1)
				add_time_offset(host_comp_id, period_start_year, period_start_month, period_start_day, period_start_second, elapsed_month, 2);
			file_current_year = period_start_year;
			file_current_month = period_start_month;
			file_current_day = period_start_day;
			file_current_second = period_start_second;
		}
	}

	EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "file_current_year:%d, file_current_month:%d, file_current_day:%d, file_current_second:%d, current_year:%d, current_month:%d, current_day:%d, current_second:%d", file_current_year, file_current_month, file_current_day, file_current_second, time_mgr->get_current_year(), time_mgr->get_current_month(), time_mgr->get_current_day(), time_mgr->get_current_second());
	find_readin_files_for_target_time(file_current_year, file_current_month, file_current_day, file_current_second);
}


void Input_handler_controller::find_readin_files_for_target_time(int file_current_year, int file_current_month, int file_current_day, int file_current_second)
{
	char file_target_time[32], time_string_left[8], time_string_right[8];
	int k, time_pos_left = -1, time_pos_right = -1, time_index_left = -1, time_index_right = -1;
	bool is_root_proc = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "in  Input_handler_controller") == 0;
	MPI_Comm comm = comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, "in Input_handler_controller");

	for (int i = 0; i < input_instance->input_datamodel->file_sets_info.size(); i ++) {
		Datamodel_file_info *current_file_info = input_instance->input_datamodel->file_sets_info[i];
		get_time_string(file_target_time, "YYYYMMDDSSSSS", true, current_file_info->time_format_in_file_name, file_current_year, file_current_month, file_current_day, file_current_second);
		time_string_and_format_match(host_comp_id, file_target_time, "YYYYMMDDSSSSS", TIME_FORMAT_YYYYMMDDSSSSS, false, 0, file_time_value);
		if (words_are_the_same(current_file_info->time_field_setting->specification, "file_field")) {
			if (current_file_info->next_input_file_times[0]->time_position == -1) {
				EXECUTION_REPORT(REPORT_ERROR, host_comp_id, file_time_value >= current_file_info->time_fields_info->time_fields[0], "Error happens when matching current time(%ld) with time_field time, current time is smaller than the smallest time_field time(%ld). Please verify.", file_time_value);
				int left = 0, right = current_file_info->time_fields_info->time_fields.size(), mid = 0;
				while (left < right) {
					mid = (left+right)/2;
					if (file_time_value < current_file_info->time_fields_info->time_fields[mid]) right = mid;
					else left = mid + 1;
				}
				if (file_time_value == current_file_info->time_fields_info->time_fields[left-1]) config_time_field_perfect_match(left-1, current_file_info);//perfect match
				else if (left == current_file_info->time_fields_info->time_fields.size())
					config_time_field_left_right_match(left-1, 0, current_file_info);
				else config_time_field_left_right_match(left-1, left, current_file_info);//need time interpolation
			}
			else {//not the first time step
				int global_ind = current_file_info->time_fields_info->time_fields[current_file_info->next_input_file_times[0]->global_ind];
				while (file_time_value > current_file_info->time_fields_info->time_fields[global_ind]) global_ind ++;
				if (file_time_value == current_file_info->time_fields_info->time_fields[global_ind])//perfect match
					config_time_field_perfect_match(global_ind, current_file_info);
				else if (global_ind == current_file_info->time_fields_info->time_fields.size())
					config_time_field_left_right_match(global_ind-1, 0, current_file_info);
				else config_time_field_left_right_match(global_ind-1, global_ind, current_file_info);
			}
		}
		else {//file_name
			if (current_file_info->next_input_file_times[0]->time_position == -1) {//first time_step
				int left = 0, right = current_file_info->file_name_times.size(), mid = 0;
				while (left < right) {
					mid = (left+right)/2;
					if (file_time_value < current_file_info->file_name_times[mid]) right = mid;
					else left = mid + 1;
				}
				if (file_time_value == current_file_info->file_name_times[left-1])
					config_file_name_time_perfect_match(left-1, current_file_info);
				else if (left == current_file_info->file_name_times.size()) 
					config_file_name_time_left_right_match(left-1, 0, current_file_info);
				else config_file_name_time_left_right_match(left-1, left, current_file_info);
			}
			else {
				int global_ind = current_file_info->next_input_file_times[0]->time_ind_in_dir;
				while (global_ind < current_file_info->file_name_times.size() && file_time_value > current_file_info->file_name_times[global_ind]) global_ind ++;
				if (file_time_value == current_file_info->file_name_times[global_ind])
					config_file_name_time_perfect_match(global_ind, current_file_info);
				else if (global_ind == current_file_info->file_name_times.size())
					config_file_name_time_left_right_match(global_ind-1, 0, current_file_info);
				else config_file_name_time_left_right_match(global_ind-1, global_ind, current_file_info);
			}
		}
	}
}


void Input_handler_controller::get_full_file_name_from_time_value_on_target_format(Datamodel_file_info *current_file_info, int time_ind_in_dir, char *full_file_name_buff)
{	char time_string[32];
	get_time_string_from_time_value_and_format(host_comp_id, current_file_info->file_name_times[time_ind_in_dir], "YYYYMMDDSSSSS", current_file_info->time_format_in_file_name, time_string);
	sprintf(full_file_name_buff, "%s/%s%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, time_string, current_file_info->file_name_suffix);
}


void Input_handler_controller::config_file_name_time_perfect_match(int global_ind, Datamodel_file_info *current_file_info)
{
	bool current_read_flag = false;
	if (current_file_info->next_input_file_times[0]->global_ind != global_ind) {
		current_read_flag = true;
		if (current_file_info->next_input_file_times[1]->global_ind == global_ind) {
			strcpy(current_file_info->next_input_file_times[0]->full_file_name, current_file_info->next_input_file_times[1]->full_file_name);
			current_file_info->move_flag = true;
		}
		else get_full_file_name_from_time_value_on_target_format(current_file_info, global_ind, current_file_info->next_input_file_times[0]->full_file_name);
		delete current_file_info->next_input_file_times[0]->netcdf_file_object;
#ifdef USE_PARALLEL_IO
		IO_pnetcdf *netcdf_file_object = new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, current_file_info->next_input_file_times[0]->full_file_name, current_file_info->next_input_file_times[0]->full_file_name, "r", true);
#else
		IO_netcdf *netcdf_file_object = new IO_netcdf(current_file_info->next_input_file_times[0]->full_file_name, current_file_info->next_input_file_times[0]->full_file_name, "r", true);
#endif
		current_file_info->next_input_file_times[0]->set_value((void*)netcdf_file_object, global_ind, -1, global_ind, current_read_flag);
	}
	current_file_info->next_input_file_times[1]->read_flag = false;
	current_file_info->perfect_match = true;
}

void Input_handler_controller::config_file_name_time_left_right_match(int left, int right, Datamodel_file_info *current_file_info)
{
	std::vector<bool> read_flags(2);
	std::vector<int> global_index(2);
	global_index[0] = left;
	global_index[1] = right;

	if (current_file_info->next_input_file_times[0]->global_ind != left || current_file_info->next_input_file_times[1]->global_ind != right) {
		current_file_info->move_flag = false;
		if (current_file_info->next_input_file_times[0]->global_ind == left || current_file_info->next_input_file_times[1]->global_ind == left) {
			read_flags[0] = false;
			read_flags[1] = true;
			if (current_file_info->next_input_file_times[1]->global_ind == left) {//mem_cpy
				current_file_info->move_flag = true;
			}
		}
		else {
			for (int i = 0; i < read_flags.size(); i++) 
				read_flags[i] = true;
		}

		for (int i = 0; i < current_file_info->next_input_file_times.size(); i++) {
			void *netcdf_file_object;
			get_full_file_name_from_time_value_on_target_format(current_file_info, global_index[i], current_file_info->next_input_file_times[i]->full_file_name);
			delete current_file_info->next_input_file_times[i]->netcdf_file_object;
#ifdef USE_PARALLEL_IO
			netcdf_file_object = (void*) new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, current_file_info->next_input_file_times[i]->full_file_name, current_file_info->next_input_file_times[i]->full_file_name, "r", false);
#else
			netcdf_file_object = (void*) new IO_netcdf(current_file_info->next_input_file_times[i]->full_file_name, current_file_info->next_input_file_times[i]->full_file_name, "r", false);
#endif
			current_file_info->next_input_file_times[i]->set_value(netcdf_file_object, global_index[i], -1, global_index[i], read_flags[i]);
		}
	}
	else {
		for (int i = 0; i < current_file_info->next_input_file_times.size(); i++)
			current_file_info->next_input_file_times[i]->read_flag = false;
	}
	current_file_info->perfect_match = false;
}

void Input_handler_controller::config_time_field_perfect_match(int global_ind, Datamodel_file_info *current_file_info) 
{
	int ind = 0;
	while (global_ind > current_file_info->time_fields_info->last_time_global_ind[ind]) ind ++;
	int time_position = global_ind - current_file_info->time_fields_info->last_time_global_ind[ind-1];
	if (current_file_info->next_input_file_times[0]->global_ind != global_ind) {
		if (current_file_info->next_input_file_times[1]->global_ind == global_ind) {//mem_cpy
			current_file_info->next_input_file_times[0]->copy_info(current_file_info->next_input_file_times[1]);
			current_file_info->move_flag = true;
			current_file_info->next_input_file_times[0]->read_flag = false;
		}
		else {//update full_file_name
			//next_input_file_times[0]->netcdf_file_object
			if (current_file_info->next_input_file_times[0]->time_ind_in_dir != ind && current_file_info->next_input_file_times[1]->time_ind_in_dir != ind) {//0 !=, 1 != : set full_file_name from time_value on target format
				//set full_file_name for next_input_file_time_info[0]
				get_full_file_name_from_time_value_on_target_format(current_file_info, ind, current_file_info->next_input_file_times[0]->full_file_name);
			}
			else if (current_file_info->next_input_file_times[1]->time_ind_in_dir == ind) {//1=: copy
				strcpy(current_file_info->next_input_file_times[0]->full_file_name, current_file_info->next_input_file_times[1]->full_file_name);
			}
		}
		//update netcdf_file_object
		delete current_file_info->next_input_file_times[0]->netcdf_file_object;
#ifdef USE_PARALLEL_IO
		IO_pnetcdf *netcdf_file_object = new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, current_file_info->next_input_file_times[0]->full_file_name, current_file_info->next_input_file_times[0]->full_file_name, "r", true);
#else
		IO_netcdf *netcdf_file_object = new IO_netcdf(current_file_info->next_input_file_times[0]->full_file_name, current_file_info->next_input_file_times[0]->full_file_name, "r", true);
#endif
		current_file_info->next_input_file_times[0]->set_value((void*)netcdf_file_object, ind, time_position, global_ind, true);
	}
	else {//unlikely to happen
		for (int i = 0; i < current_file_info->next_input_file_times.size(); i++)
			current_file_info->next_input_file_times[i]->read_flag = false;
	}
	//next_input_file_times[1]
	current_file_info->next_input_file_times[1]->read_flag = false;
	current_file_info->perfect_match = false;
}

void Input_handler_controller::config_time_field_left_right_match(int left, int right, Datamodel_file_info *current_file_info) 
{
	int ind = 0;
	std::vector<int> time_ind(2), time_position(2), global_index(2);
	std::vector<bool> read_flags;

	while (left > current_file_info->time_fields_info->last_time_global_ind[ind]) ind ++;
	time_ind[0] = ind;
	if (right > current_file_info->time_fields_info->last_time_global_ind[ind]) ind ++;
	time_ind[1] = ind;
	time_position[0] = left - current_file_info->time_fields_info->last_time_global_ind[ind-1];
	time_position[1] = right - current_file_info->time_fields_info->last_time_global_ind[ind-1];
	global_index[0] = left;
	global_index[1] = right;

	if (current_file_info->next_input_file_times[0]->global_ind != left || current_file_info->next_input_file_times[1]->global_ind != right) {
		if (current_file_info->next_input_file_times[0]->global_ind == left || current_file_info->next_input_file_times[1]->global_ind == left) {
			read_flags[0] = false;
			read_flags[1] = true;
			if (current_file_info->next_input_file_times[1]->global_ind == left) {//mem_cpy
				current_file_info->move_flag = true;
			}
		}
		else {
			for (int i = 0; i < read_flags.size(); i++)
				read_flags[i] = true;
		}

		for (int i = 0; i < current_file_info->next_input_file_times.size(); i++) {
			//move_flag also need to new netcdf_file_object
			void *netcdf_file_object;
			if (current_file_info->next_input_file_times[i]->time_ind_in_dir != time_ind[i]) {get_full_file_name_from_time_value_on_target_format(current_file_info, time_ind[i], current_file_info->next_input_file_times[i]->full_file_name);
				delete current_file_info->next_input_file_times[i]->netcdf_file_object;
#ifdef USE_PARALLEL_IO
				netcdf_file_object = (void*) new IO_pnetcdf(host_comp_id, pio_proc_num, io_proc_mark, this->io_comm, current_file_info->next_input_file_times[i]->full_file_name, current_file_info->next_input_file_times[i]->full_file_name, "r", true);
#else
				netcdf_file_object = (void*) new IO_netcdf(current_file_info->next_input_file_times[i]->full_file_name, current_file_info->next_input_file_times[i]->full_file_name, "r", true);
#endif
			}
			else netcdf_file_object = (void*) current_file_info->next_input_file_times[i]->netcdf_file_object;
			current_file_info->next_input_file_times[i]->set_value(netcdf_file_object, time_ind[i], time_position[i], global_index[i], read_flags[i]);
		}
	}
	else {
		for (int i = 0; i < current_file_info->next_input_file_times.size(); i++)
			current_file_info->next_input_file_times[i]->read_flag = false;
	}
	current_file_info->perfect_match = true;
}


void Input_handler_controller::common_checking_for_read_in_field_handler(const char *handler_name, int comp_id, int *field_instance_ids, int num_fields, int field_instance_ids_size, const char *data_file_name, const char *field_name_in_file, int *necessity, int grid_id_for_file, const char *annotation)
{
	int API_id = API_ID_HANDLER_DATA_FILE_INPUT;
	char API_label[256];
	get_API_hint(-1, API_ID_HANDLER_DATA_FILE_INPUT, API_label);

	int temp_int = (grid_id_for_file == -1)? 0 : 1;
	check_API_parameter_int(comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Input_handler_controller::Intput_handler_controller"), NULL, temp_int, "grid_id_for_file", annotation);

	/*if (necessity[0] != -1) {
		for (int i = 0; i < num_fields; i++)
		check_API_parameter_int(comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Input_handler_controller::Input_handler_controller"), "registering an Input_handler_controller", necessity[i], "necessity", annotation);
	}*/

	check_API_parameter_string(comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Input_handler_controller::Input_handler_controller"), "registering an input datamodel handler", handler_name, "handler_name", annotation);
	if (strlen(data_file_name) != 0)
		check_API_parameter_string(comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Input_handler_controller::Input_handler_controller"), "registering an input datamodel handler", data_file_name, "data_file_name", annotation);
	if (field_name_in_file != NULL)
		check_API_parameter_string(comp_id, API_ID_HANDLER_DATA_FILE_INPUT, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(comp_id, "in Input_handler_controller::Input_handler_controller"), "registering an input datamodel handler", field_name_in_file, "field_name_in_file", annotation);
	datamodel_mgr->common_checking_for_datamodel_handler_registration(num_fields, field_instance_ids, 1, grid_id_for_file, -1, -1, -1, field_instance_ids_size, handler_name, annotation, false);
}


Input_handler_controller::~Input_handler_controller()
{
	for (int i = 0; i < next_input_file_times.size(); i++)
		delete next_input_file_times[i];

	for (int i = 0; i < input_fields_info.size(); i++)
		delete input_fields_info[i];

	for (int i = 0; i < fields_config_info.size(); i++) {
		for (int j = 0; j < fields_config_info[i].size(); j++)
			delete fields_config_info[i][j];
	}
}


int Datamodel_mgt::register_input_handler_operator(int field_instance_id, int io_grid_id, const char *io_data_type, const char *io_unit)
{
	int API_id = API_ID_HANDLER_DATA_FILE_INPUT;
	char API_label[256];

	get_API_hint(-1, API_ID_HANDLER_DATA_FILE_INPUT, API_label);
	char *handler_name = strdup("read_in_field_handler_operator");
	EXECUTION_REPORT_LOG(REPORT_LOG, -1, true, "about to register an Input_handler_operator");

	int handler_id = get_next_input_handler_operator_id();
	Input_handler_operator *input_operator = new Input_handler_operator(handler_name, handler_id, field_instance_id, io_grid_id, io_data_type, io_unit, handler_name);
	input_handler_operators.push_back(input_operator);
	return handler_id;
}

Input_handler_operator::Input_handler_operator(const char *handler_name, int handler_id, int field_instance_id, int io_grid_id, const char *io_data_type, const char *io_unit, const char *annotation)
{
	int field_necessity = 1;
	char handler_import_interface_name[NAME_STR_SIZE];

	this->handler_type = INPUT_HANDLER;
	this->annotation = strdup(annotation);
	this->handler_name = strdup(handler_name);
	Field_mem_info *model_field_instance = memory_manager->get_field_instance(field_instance_id);
	this->host_comp_id = model_field_instance->get_comp_id();

#ifdef USE_PARALLEL_IO
	pio_proc_num = datamodel_mgr->get_comp_PIO_proc_setting(host_comp_id, this->io_proc_stride, this->io_proc_mark, this->io_comm);
#endif
	input_model_timer_id = timer_mgr->define_timer(host_comp_id, "steps", 1, 0, 0, annotation);

	sprintf(handler_import_interface_name, "%s_import_interface", handler_name);
	handler_import_interface = new Inout_interface(handler_import_interface_name, TYPE_INPUT_HANDLER_ID_PREFIX|datamodel_mgr->get_num_input_handler_operators(), COUPLING_INTERFACE_MARK_IMPORT, 1, &field_instance_id, 1, input_model_timer_id, 0, "field_instance_IDs", annotation, API_ID_HANDLER_DATA_FILE_INPUT, INTERFACE_SOURCE_IO_INPUT, true);
	handler_import_interface->set_fields_necessity(&field_necessity, 1, this->annotation);

	config_input_handler_operator_export_interface(field_instance_id, io_grid_id, io_data_type, io_unit);

	inout_interface_mgr->do_coupling_generation_between_interfaces_intra_one_component_model(handler_import_interface, handler_export_interface, -1, false);
}


void Input_handler_operator::config_input_handler_operator_export_interface(int field_instance_id, int io_grid_id, const char *io_data_type, const char *io_unit)
{
	int io_field_decomp_id;
	char interface_name[NAME_STR_SIZE];

	export_interface_field_info = new Import_field_info();
	export_interface_field_info->import_field_operation = NULL;
	Field_mem_info *model_field_instance = memory_manager->get_field_instance(field_instance_id);
	int grid_id_for_io = determine_field_io_grid(model_field_instance, io_grid_id);//can be skipped?already done in controller?
	//decomposition
	if (grid_id_for_io != -1) {
#ifdef USE_PARALLEL_IO
		Decomp_info *new_decomp = decomps_info_mgr->generate_parallel_decomp_for_parallel_IO(original_grid_mgr->search_grid_info(grid_id_for_io), this->io_proc_stride, this->pio_proc_num, this->io_proc_mark);
#else
		Decomp_info *new_decomp = decomps_info_mgr->generate_default_parallel_decomp_serial(original_grid_mgr->search_grid_info(grid_id_for_io));
#endif
		if (new_decomp != NULL)
			io_field_decomp_id = new_decomp->get_decomp_id();
	}
	else io_field_decomp_id = host_comp_id;

	datamodel_mgr->add_handlers_field_mem_buf_mark();
	Field_mem_info *io_field_instance = memory_manager->alloc_mem(model_field_instance->get_field_name(), io_field_decomp_id, grid_id_for_io, BUF_MARK_IO_FIELD_MIRROR ^ datamodel_mgr->get_handlers_field_mem_buf_mark(), io_data_type, io_unit, this->annotation, false, false);
	io_field_instance->set_usage_tag(REG_FIELD_TAG_NONE);

	export_interface_field_info->model_field_instance_id = field_instance_id;
	export_interface_field_info->io_field_instance_id = io_field_instance->get_field_instance_id();
	//interface
	int export_field_instance_id = io_field_instance->get_field_instance_id();
	sprintf(interface_name, "%s_export_interface", this->handler_name);
	input_io_timer_id = timer_mgr->define_timer(host_comp_id, "steps", 1, 0, 0, annotation);
	handler_export_interface = new Inout_interface(interface_name, TYPE_OUTPUT_HANDLER_ID_PREFIX|datamodel_mgr->get_num_input_handler_operators(), COUPLING_INTERFACE_MARK_EXPORT, 1, &export_field_instance_id, 1, input_io_timer_id, 0, "field_instance_IDs", this->annotation, API_ID_HANDLER_DATA_FILE_INPUT, INTERFACE_SOURCE_IO_INPUT, true);
}


int Input_handler_operator::determine_field_io_grid(Field_mem_info *model_field_instance, int io_grid_id)
{
	if (model_field_instance->get_grid_id() == -1) {
		return -1;
	}
	if (io_grid_id != -1)//can be skipped?
		return io_grid_id;
	else return model_field_instance->get_grid_id();
	//check_grid_dim?
}


Field_mem_info *Input_handler_operator::get_unique_IO_field_mem()
{
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, handler_export_interface->get_num_fields_mem_registered() == 1, "Software error in Input_handler_operator::get_unique_IO_field_mem, num_fields_mem_registered:%d", handler_export_interface->get_num_fields_mem_registered());
	return handler_export_interface->get_field_mem(0);
}


void Input_handler_operator::execute_input_handler_operator(const char *field_name, Input_file_time_info *next_input_file_time, const char *annotation)
{
	execute_read_in_field(field_name, next_input_file_time, annotation);
	//execute export/import interface
	execute_input_handler_operator_interface(true, annotation);
}


void Input_handler_operator::execute_read_in_field(const char *field_name, Input_file_time_info *next_input_file_time, const char *annotation)
{
	//no need to specify field_name
	Field_mem_info *io_field_instance = get_unique_IO_field_mem();
	int local_proc_id = comp_comm_group_mgt_mgr->get_current_proc_id_in_comp(host_comp_id, "Input_handler_operator");
	strcpy(io_field_instance->get_field_data()->get_grid_data_field()->field_name_in_IO_file, field_name);

#ifdef USE_PARALLEL_IO
	if (io_proc_mark == 1) {
		if (next_input_file_time->netcdf_file_object->get_io_with_time_info())
			next_input_file_time->netcdf_file_object->read_field_data(host_comp_id, io_field_instance, next_input_file_time->time_position, true);//temp
		else next_input_file_time->netcdf_file_object->read_field_data(host_comp_id, io_field_instance, -1, true);
	}
	else if (io_field_instance->get_decomp_id() != -1)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, io_field_instance->get_size_of_field() == 0, "Software error in Output_handler::execute_handler");
#else
	if (local_proc_id == 0) {
		if (next_input_file_time->netcdf_file_object->get_io_with_time_info())
			next_input_file_time->netcdf_file_object->read_data(io_field_instance->get_field_data()->get_grid_data_field(), next_input_file_time->time_position, true);//temp
		else next_input_file_time->netcdf_file_object->read_data(io_field_instance->get_field_data()->get_grid_data_field(), -1, true);
	}
#endif
}


void Input_handler_operator::execute_input_handler_operator_interface(bool bypass_timer, const char *execute_annotation)
{
	int export_field_update_status, import_field_update_status;

	get_unique_IO_field_mem()->define_field_values(false);
	handler_export_interface->execute(bypass_timer, API_ID_INTERFACE_EXECUTE_WITH_ID, &export_field_update_status, 1, execute_annotation);
	handler_export_interface->dump_active_coupling_connections();

	handler_import_interface->get_field_mem(0)->define_field_values(false);
	handler_import_interface->execute(bypass_timer, API_ID_INTERFACE_EXECUTE_WITH_ID, &import_field_update_status, 1, execute_annotation);
	handler_import_interface->dump_active_coupling_connections();
}


Input_handler_operator::~Input_handler_operator() 
{
	delete handler_export_interface;
	delete handler_import_interface;
	delete export_interface_field_info;
	delete annotation, handler_name;
}

int Datamodel_mgt::register_datamodel_input_handler(const char *handler_name, int *field_instance_ids, int num_fields, int field_instance_ids_size, const char *config_input_instance_name, int input_timer_id, int *necessity, int *field_connected_status, const char *annotation)
{
	Input_instance *new_input_instance;
	Input_handler_controller *new_input_handler_controller;
	int API_id = API_ID_HANDLER_DATAMODEL_INPUT, i;
	char API_label[256];
	get_API_hint(-1, API_ID_HANDLER_DATAMODEL_INPUT, API_label);
	char *datamodel_name;

	int host_comp_id = memory_manager->get_field_instance(field_instance_ids[0])->get_comp_id();

	for (i = 0; i < input_instances.size(); i ++) {
		if (words_are_the_same(input_instances[i]->get_input_instance_name(), config_input_instance_name))
			new_input_instance = input_instances[i];
	}

	if (i == input_instances.size()) {
		new_input_instance = new Input_instance(host_comp_id, config_input_instance_name, annotation);
		input_instances.push_back(new_input_instance);
		input_datamodels.push_back(new_input_instance->get_input_datamodel());
	}

	for (int m = 0; m < input_handler_controllers.size(); i++)
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, host_comp_id, !words_are_the_same(handler_name, input_handler_controllers[m]->get_handler_name()) && host_comp_id == input_handler_controllers[m]->get_host_comp_id(), "Error happens when calling \"%s\" at the model code with the annotation \"%s\", a handler with the name \"%s\" has already been registered at the model code with the annotation \"%s\", Please Verify.", API_label, annotation, handler_name, input_handler_controllers[m]->get_handler_annotation());

	const char *example_input_data_file = new_input_instance->get_input_datamodel()->randomly_match_a_data_file_in_datamodel(NULL);
	EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "example_input_data_file:%s", example_input_data_file);
	new_input_handler_controller = new Input_handler_controller(handler_name, new_input_instance->get_input_instance_name(), get_next_input_handler_controller_id(), field_instance_ids, num_fields, field_instance_ids_size, input_timer_id, "netcdf", necessity, field_connected_status, example_input_data_file, NULL, -1, NULL, NULL, false, false, NULL, annotation);
	input_handler_controllers.push_back(new_input_handler_controller);
	delete example_input_data_file;
	EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "Input_handler_controller_id:%d", new_input_handler_controller->get_input_handler_controller_id());
	return new_input_handler_controller->get_input_handler_controller_id();
}


const char *Inout_datamodel::randomly_match_a_data_file_in_datamodel(Datamodel_file_info *current_file_info)
{
	char time_string[8];
	char *full_file_name = new char[NAME_STR_SIZE];
	if (current_file_info == NULL) current_file_info = file_sets_info[0];

	if (strlen(current_file_info->time_format_in_file_name) != 0) {
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, current_file_info->file_name_times.size() != 0, "file_time size:%d", current_file_info->file_name_times.size());
		long time_in_file = current_file_info->file_name_times[0];
		get_time_string_from_time_value_and_format(host_comp_id, time_in_file, "YYYYMMDDSSSSS", current_file_info->time_format_in_file_name, time_string);
		EXECUTION_REPORT(REPORT_LOG, host_comp_id, true, "random time_string:%s", time_string);
		sprintf(full_file_name, "%s/%s%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, time_string, current_file_info->file_name_suffix);
	}
	else sprintf(full_file_name, "%s/%s%s.nc", current_file_info->file_dir, current_file_info->file_name_prefix, current_file_info->file_name_suffix);
	return full_file_name;
}



Input_instance::Input_instance(int host_comp_id, const char *config_input_instance_name, const char *annotation)
{
	this->host_comp_id = host_comp_id;
	this->annotation = strdup(annotation);
	XML_file_name = new char [strlen(comp_comm_group_mgt_mgr->get_config_root_dir())+256];
	char *datamodel_config_dir = new char [strlen(comp_comm_group_mgt_mgr->get_config_root_dir())+256];
	sprintf(datamodel_config_dir, "%s/CCPL_dir/datamodel/config",comp_comm_group_mgt_mgr->get_root_working_dir());

	Comp_comm_group_mgt_node *node = comp_comm_group_mgt_mgr->get_global_node_of_local_comp(host_comp_id, true, "in Input_instance::load_coupling_connection_config_file");
	const char *current_comp_full_name = node->get_full_name();
	sprintf(XML_file_name, "%s/%s.input_instances.xml", datamodel_config_dir, current_comp_full_name);
	config_input_instance_from_xml_file(XML_file_name, comp_comm_group_mgt_mgr->get_comm_group_of_local_comp(host_comp_id, ""), config_input_instance_name);
}


void Input_instance::config_input_instance_from_xml_file(char *XML_file_name, MPI_Comm comm, const char *config_input_instance_name)
{
	int line_number;
	bool find_target_input_instance = false;

	TiXmlDocument *XML_file = open_XML_file_to_read(host_comp_id, XML_file_name, comm, false);
    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, XML_file != NULL, "the \"config_input_instance_name\" (\"%s\") for specifying a coupling_connection configuration file (\"%s\") is wrong", config_input_instance_name, XML_file_name);

    TiXmlElement *root_XML_element, *XML_element;
    TiXmlNode *input_instance_node;

    TiXmlNode *root_XML_element_node = get_XML_first_child_of_unique_root(host_comp_id, XML_file_name, XML_file);

    for (; root_XML_element_node != NULL; root_XML_element_node = root_XML_element_node->NextSibling()) {
    	if (root_XML_element_node->Type() != TiXmlNode::TINYXML_ELEMENT)
    		continue;
    	root_XML_element = root_XML_element_node->ToElement();
    	if (words_are_the_same(root_XML_element->Value(), "input_instances"))
    		break;
    }
    EXECUTION_REPORT(REPORT_ERROR, host_comp_id, root_XML_element_node != NULL, "Couldn't find node \"input_instances\" in XML file \"%s\"", XML_file_name);
    for (input_instance_node = root_XML_element_node->FirstChild(); input_instance_node != NULL; input_instance_node = input_instance_node->NextSibling()) {
    	TiXmlElement *input_instance_element = input_instance_node->ToElement();
    	if (!is_XML_setting_on(host_comp_id, input_instance_element, XML_file_name,"the status of a data_file node of an input_instances node", "coupling_connection XML file"))
    		continue;
    	const char *instance_name_str = get_XML_attribute(host_comp_id, -1, input_instance_element, "instance_name", XML_file_name, line_number, "the instance_name of the input_instance", "the coupling connection configuration file", true);
    	if (!words_are_the_same(instance_name_str, config_input_instance_name))
    		continue;

    	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, !find_target_input_instance, "Error happens when configuring Input_instance \"%s\", another Input_instance with the same name has already been found. Please verify xml configuration file \"%s\"", instance_name_str, XML_file_name);
    	this->instance_name = strdup(instance_name_str);

    	find_target_input_instance = true;
    	this->datamodel_name = get_XML_attribute(host_comp_id, -1, input_instance_element, "datamodel_name", XML_file_name, line_number, "the instance_name of the input_instance", "the coupling connection configuration file", true);
    	input_datamodel = new Inout_datamodel(host_comp_id, this->datamodel_name, this->annotation);

    	TiXmlNode *time_mapping_node = input_instance_node->FirstChild();
    	config_time_mapping_info_from_xml(time_mapping_node);

		TiXmlNode *field_renaming_node = time_mapping_node->NextSibling();
    	config_input_instance_field_renamings(field_renaming_node);
    }
}

void Input_instance::config_time_mapping_info_from_xml(TiXmlNode *time_mapping_node)
{
	TiXmlElement *time_mapping_element = time_mapping_node->ToElement();
	bool config_period = false, config_offset = false;
	period_info = NULL, offset_info = NULL;

	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, words_are_the_same(time_mapping_element->Value(), "time_mapping_configurations"), "The first node under \"input_instance\" node should name \"time_mapping_configurations\". Please verify the XML configuration file \"%s\"", XML_file_name);

	for (TiXmlNode *setting_node = time_mapping_node->FirstChild(); setting_node != NULL; setting_node = setting_node->NextSibling()) {
		TiXmlElement *setting_element = setting_node->ToElement();
		if (!is_XML_setting_on(host_comp_id, setting_element, XML_file_name,"the status of a period/offset_setting node of an input_instances node", "coupling_connection XML file"))
			continue;
		if (words_are_the_same(setting_element->Value(), "period_setting")) {
			EXECUTION_REPORT(REPORT_ERROR, host_comp_id, !config_period, "Error happens when configuring input_instance \"%s\", \"period_setting\" node has been set twice. Please verify the XML configuration file \"%s\"", XML_file_name);
			config_period_setting(setting_element);
			config_period = true;
		}
		else if (words_are_the_same(setting_element->Value(), "offset_setting")) {
			EXECUTION_REPORT(REPORT_ERROR, host_comp_id, !config_offset, "Error happens when configuring input_instance \"%s\", \"offset_setting\" node has been set twice. Please verify the XML configuration file \"%s\"", instance_name, XML_file_name);
			config_offset_setting(setting_element);
			config_offset = true;
		}
		else EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "Error happens when configuring input_instance \"%s\", \"time_mapping_configurations\" node should only be conposed of \"period_setting\" node or \"offset_setting\" node. Please verify the XML configuration file \"%s\"", instance_name, XML_file_name);
	}
	if (config_period && config_offset) {
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, period_info->id_period_unit == offset_info->id_offset_unit, "Error happens in Input_instance \"time_mapping_configurations\": \"period_unit\" and \"period_count\" should be the same. Please verify XML configuration file \"%s\" around line %d", XML_file_name, time_mapping_element->Row());
		EXECUTION_REPORT(REPORT_ERROR, host_comp_id, period_info->period_count >= offset_info->offset_count,"Error happens in Input_instance \"time_mapping_configurations\": \"period_count\" should be largher than \"offset_count\". Please verify XML configuration file \"%s\" around line %d", XML_file_name, time_mapping_element->Row());
	}
}

void Input_instance::config_offset_setting(TiXmlElement *setting_element)
{
	offset_info = new Offset_setting();
	const char *offset_unit_str = get_XML_attribute(host_comp_id, -1, setting_element, "offset_unit", XML_file_name, line_number, "the offset_unit of the input_instance", "the coupling connection configuration file", true);
	offset_info->id_offset_unit = set_unit(offset_unit_str, "offset_unit");
	offset_info->offset_unit = strdup(offset_unit_str);

	const char *offset_count_str = get_XML_attribute(host_comp_id, -1, setting_element, "offset_count", XML_file_name, line_number, "the offset_count of the input_instance", "the coupling connection configuration file", true);
	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, sscanf(offset_count_str, "%d", &(offset_info->offset_count)) == 1, "Error happens when configuring Input_instance \"%s\", \"offset_count\" is not an integer. Please verify xml configuration file \"%s\"", instance_name, XML_file_name);
}

void Input_instance::config_period_setting(TiXmlElement *setting_element)
{
	period_info = new Period_setting();
	const char *period_data_start_time_str = get_XML_attribute(host_comp_id, -1, setting_element, "period_data_start_time", XML_file_name, line_number, "the period_data_start_time of the input_instance", "the coupling connection configuration file", true);
	EXECUTION_REPORT(REPORT_ERROR, -1, sscanf(period_data_start_time_str, "%d", &(period_info->period_data_start_time)) == 1, "Error happens when configuring Input_instance \"%s\", \"period_data_start_time\" is not an integer. Please verify xml configuration file \"%s\"",instance_name, XML_file_name);

	const char *period_time_format_str = get_XML_attribute(host_comp_id, -1, setting_element, "period_time_format", XML_file_name, line_number, "the period_time_format of the input_instance", "the coupling connection configuration file", true);
	period_info->id_period_time_format = check_time_format(period_time_format_str, "period_time_format");
	period_info->period_time_format = strdup(period_time_format_str);

	const char *period_unit_str = get_XML_attribute(host_comp_id, -1, setting_element, "period_unit", XML_file_name, line_number, "the period_unit of the input_instance", "the coupling connection configuration file", true);
	period_info->id_period_unit = set_unit(period_unit_str, "period_unit");
	period_info->period_unit = strdup(period_unit_str);

	const char *period_count_str = get_XML_attribute(host_comp_id, -1, setting_element, "period_count", XML_file_name, line_number, "the period_count of the input_instance", "the coupling connection configuration file", true);
	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, sscanf(period_count_str, "%d", &(period_info->period_count)) == 1, "Error happens when configuring Input_instance \"%s\", \"period_count\" is not an integer. Please verify xml configuration file \"%s\"", instance_name, XML_file_name);
	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, get_time_format_smallest_unit(period_info->id_period_time_format) == period_info->id_period_unit, "The smallest unit of \"period_start_time\" (which is\"%d\") should be the same with \"period_unit\" (which is \"%s\"). Please verify XML configuration file \"%s\" under node \"period_setting\" around line %d", period_info->period_data_start_time, period_info->period_unit, XML_file_name, setting_element->Row());
}


void Input_instance::config_input_instance_field_renamings(TiXmlNode *field_renaming_node)
{
	TiXmlElement *field_renameing_element = field_renaming_node->ToElement();
	for (TiXmlNode *entry_node = field_renaming_node->FirstChild(); entry_node != NULL; entry_node = entry_node->NextSibling()) {
		TiXmlElement *entry_element = entry_node->ToElement();
		Field_rename_map *rename = new Field_rename_map();

		rename->name_in_model = get_XML_attribute(host_comp_id, -1, entry_element, "name_in_model", XML_file_name, line_number, "the name_in_model of the input_instance", "the coupling connection configuration file", true);
		rename->name_in_file = get_XML_attribute(host_comp_id, -1, entry_element, "name_in_file", XML_file_name, line_number, "the name_in_file of the input_instance", "the coupling connection configuration file", true);
		field_renamings.push_back(rename);
	}
}


const char *Input_instance::get_field_name_in_file(const char *field_name_in_model)
{
	for (int i = 0; i < field_renamings.size(); i++) {
		if (words_are_the_same(field_renamings[i]->name_in_model, field_name_in_model))
			return field_renamings[i]->name_in_file;
	}
	EXECUTION_REPORT(REPORT_ERROR, host_comp_id, false, "Software error::Input_instance: field name_in_model:%s for Input_instance \"%s\" not configured in input_instance configuration file \"%s\"", field_name_in_model, instance_name, XML_file_name);
	return NULL;
}

Input_instance::~Input_instance()
{
	delete XML_file_name, annotation;
}

Input_instance *Datamodel_mgt::search_input_instance(const char *target_instance_name) {//no need to consider host_comp_id, for models have their own input_instance configuration file
	for (int i = 0; i < input_instances.size(); i++) {
		if (words_are_the_same(input_instances[i]->get_input_instance_name(), target_instance_name))
			return input_instances[i];
	}
	EXECUTION_REPORT(REPORT_ERROR, -1, false, "Software error:target input_instance name \"%s\" can't be found", target_instance_name);
}

Inout_datamodel *Datamodel_mgt::search_input_datamodel(int host_comp_id, const char *datamodel_name)//host_comp_id is needed because input_datamodels can be shared by different comps
{
	for (int i = 0; i < input_datamodels.size(); i++) {
		if (input_datamodels[i]->get_host_comp_id() == host_comp_id && words_are_the_same(input_datamodels[i]->get_datamodel_name(), datamodel_name))
			return input_datamodels[i];
	}
	EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Error happens when searching for datamodel named \"%s\", Please verify the corresponding XML configuration file.", datamodel_name);
	return NULL;
}

void Datamodel_mgt::handle_explicit_input(int handler_id, const char *annotation)
{
	Input_handler_controller *current_input_handler_controller = get_input_handler_controller(handler_id);
	current_input_handler_controller->execute_input_handler_controller(annotation);
}


bool Datamodel_mgt::is_legal_input_handler_controller_id(int handler_id)
{
	if ((handler_id & TYPE_ID_PREFIX_MASK) != TYPE_INPUT_HANDLER_ID_PREFIX)
		return false;
	return (handler_id&TYPE_ID_SUFFIX_MASK) < input_handler_controllers.size();
}


Input_handler_controller *Datamodel_mgt::get_input_handler_controller(int handler_id) 
{
	if (!is_legal_input_handler_controller_id(handler_id))
		EXECUTION_REPORT_ERROR_OPTIONALLY(REPORT_ERROR, -1, false, "Wrong Input_handler_controller_id:%d, Please check.", handler_id);
	return input_handler_controllers[handler_id&TYPE_ID_SUFFIX_MASK];
}
