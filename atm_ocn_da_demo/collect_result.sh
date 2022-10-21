#!/bin/bash
#set -x
ROOT=/public/home/zhangxz/dafcc_test/sample_da_dafcc1/
da_demo_root=$ROOT/atm_ocn_da_demo/

atm_grid_info="20km"
ocn_grid_info="20km"
atm_proc="10"
ocn_proc="200"
#ocn_proc="$1"
#atm_proc="$3"
#atm_grid_info="${2}km"
#ocn_grid_info="${2}km"
new="true"
#new="false"
#grep_new="true"
grep_new="false"

result_path=$da_demo_root/results/$ocn_grid_info/${ocn_proc}_${atm_proc}
echo $result_path
if [ ! -e $result_path ]; then
    mkdir -p $result_path
fi
if [ $new == true  ]; then
    mv $da_demo_root/sample_da_dafcc1.* $result_path
    mv $da_demo_root/CCPL_dir/run/CCPL_logs/by_executables/ocn_demo  $result_path
    grep -rin "ccpl time\|total time" $result_path/sample_da_dafcc1.output > $result_path/model_output_time.out
    grep -rin "TIME in" $result_path/ocn_demo > $result_path/ccpl_output_time.out
fi
if [ $grep_new == true  ]; then
    grep -rin "ccpl time\|total time" $result_path/sample_da_dafcc1.output > $result_path/model_output_time.out
    grep -rin "TIME in" $result_path/ocn_demo > $result_path/ccpl_output_time.out
fi

ls $result_path

date>> ${ocn_grid_info}.results
echo "${ocn_grid_info} ${ocn_proc}_${atm_proc}"                        >> ${ocn_grid_info}.results
python get_tims.py ${ocn_grid_info} ${ocn_proc}_${atm_proc}            >> ${ocn_grid_info}.results
python get_ccpl_run_tims.py ${ocn_grid_info} ${ocn_proc}_${atm_proc}   >> ${ocn_grid_info}.results
#


