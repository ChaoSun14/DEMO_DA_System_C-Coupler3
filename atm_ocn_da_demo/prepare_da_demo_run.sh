#!/bin/bash
#set -x
ROOT=/public/home/zhangxz/dafcc_test/sample_da_dafcc1
da_demo_root=$ROOT/inputdata/da_algorithm_demo
echo $ROOT

ensemble_num=10
atm_grid_info="20km"
ocn_grid_info="20km"
num_3d_vars=3
#ensemble_num=$1
#atm_grid_info="${2}km"
#ocn_grid_info="${2}km"
echo $ensemble_num 
echo $ocn_grid_info

ens_run_dir=$ROOT/atm_ocn_da_demo/da_demo
if [ ! -e $ens_run_dir ]; then
    mkdir -p $ens_run_dir
fi
sed -i '6s/NANALS    =.*/NANALS    = '${ensemble_num}'/' $da_demo_root/$ocn_grid_info/ocn_da_algorithm_demo.nml
sed -i '7s/NVARS     =.*/NVARS     = '${num_3d_vars}'/' $da_demo_root/$ocn_grid_info/ocn_da_algorithm_demo.nml
cp $da_demo_root/$ocn_grid_info/ocn_da_algorithm_demo.nml $ens_run_dir/
cp $da_demo_root/$ocn_grid_info/ocn_demo.nc $ens_run_dir

sed -i '6s/NANALS    =.*/NANALS    = '${ensemble_num}'/' $da_demo_root/$atm_grid_info/atm_da_algorithm_demo.nml
sed -i '7s/NVARS     =.*/NVARS     = '${num_3d_vars}'/' $da_demo_root/$atm_grid_info/atm_da_algorithm_demo.nml
cp $da_demo_root/$atm_grid_info/atm_da_algorithm_demo.nml $ens_run_dir/
cp $da_demo_root/$atm_grid_info/atm_demo.nc $ens_run_dir

for ((i = 1; i <= ensemble_num ; i++))
do
    ens_index=$(printf "%d" "$i")
    echo "ensemble index: " ${ens_index}
    atm_run_dir=$ROOT/atm_ocn_da_demo/run/atm/atm_demo/${ens_index}/data
    ocn_run_dir=$ROOT/atm_ocn_da_demo/run/ocn/ocn_demo/${ens_index}/data
    if [ ! -e $ocn_run_dir/record_time ]; then
        touch $ocn_run_dir/record_time
    fi
    cp $da_demo_root/$ocn_grid_info/ocn_da_algorithm_demo.nml $ocn_run_dir/
    cp $da_demo_root/$atm_grid_info/atm_da_algorithm_demo.nml $atm_run_dir/
done
#


