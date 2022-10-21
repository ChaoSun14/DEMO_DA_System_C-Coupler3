#!/bin/csh -f

# === Note by Li Ruizhe ===
# Required paramters:
#   DATA_DIR
#   DATA_DIR
#   MODEL_DATA_DIR
#
#   RUNTYPE
#   GRID
#   RUNNING_SETTINGS_RUN_REFCASE
#   RUNNING_SETTINGS_RUN_START_DATE
#   RAMP_CO2_START_YMD
#   DOUT_L_MSNAME
# =========================

set called=($_)
if ("$called" != "") then
    set me = "$called[2]"    # the script was sourced from this location
endif
if ("$0" != "csh") then
    set me = "$0"                # the script was run from this location
endif
set me = `readlink -f $me`
set MYPATH = `dirname "$me"`

cd $DATA_DIR


set start_year = `echo $RUNNING_SETTINGS_RUN_START_DATE | awk '{print substr($0,1,4)}'`
set start_month = `echo $RUNNING_SETTINGS_RUN_START_DATE | awk '{print substr($0,5,2)}'`
set start_day = `echo $RUNNING_SETTINGS_RUN_START_DATE | awk '{print substr($0,7,2)}'`
@ start_hour = $RUNNING_SETTINGS_RUN_START_SECOND / 3600
@ start_minute = ( $RUNNING_SETTINGS_RUN_START_SECOND / 60 ) % 60
@ start_second = $RUNNING_SETTINGS_RUN_START_SECOND % 60
set start_hour = "$start_hour"
set start_minute = "$start_minute"
set start_second = "$start_second"
if ( $start_hour < 10 ) set start_hour = "0$start_hour"
if ( $start_minute < 10 ) set start_minute = "0$start_minute"
if ( $start_second < 10 ) set start_second = "0$start_second"

set stop_year = `echo $RUNNING_SETTINGS_RUN_STOP_DATE | awk '{print substr($0,1,4)}'`
set stop_month = `echo $RUNNING_SETTINGS_RUN_STOP_DATE | awk '{print substr($0,5,2)}'`
set stop_day = `echo $RUNNING_SETTINGS_RUN_STOP_DATE | awk '{print substr($0,7,2)}'`
@ stop_hour = ( $RUNNING_SETTINGS_RUN_STOP_SECOND) / 3600
@ stop_minute = ( ( $RUNNING_SETTINGS_RUN_STOP_SECOND) / 60 ) % 60
@ stop_second = ( $RUNNING_SETTINGS_RUN_STOP_SECOND) % 60
set stop_hour = "$stop_hour"
set stop_minute = "$stop_minute"
set stop_second = "$stop_second"
if ( $stop_hour < 10 ) set stop_hour = "0$stop_hour"
if ( $stop_minute < 10 ) set stop_minute = "0$stop_minute"
if ( $stop_second < 10 ) set stop_second = "0$stop_second"

set restart = .false.
if ($RUN_TYPE == 'restart') then
    set restart = .true.
endif
echo
if ( $CUSTOM_SETTINGS_GRID_INFO == "5km" ) then
    cp -rf $MODEL_DATA_DIR/atm*5km.nc ./atm_demo.nc
    echo "5KM Grid used"
else if ( $CUSTOM_SETTINGS_GRID_INFO == "10km" ) then
    cp -rf $MODEL_DATA_DIR/atm*10km.nc ./atm_demo.nc
    echo "10KM Grid used"
else if ( $CUSTOM_SETTINGS_GRID_INFO == "20km" ) then
    cp -rf $MODEL_DATA_DIR/atm*20km.nc ./atm_demo.nc
    echo "20KM Grid used"
else if ( $CUSTOM_SETTINGS_GRID_INFO == "50km" ) then
    cp -rf $MODEL_DATA_DIR/atm*50km.nc ./atm_demo.nc
    echo "50KM Grid used"
else
    echo "No grid file corresponding to the grid_info("$CUSTOM_SETTINGS_GRID_INFO"), please check!"
endif

cp -rf $MODEL_DATA_DIR/atm_demo.nml atm_demo.nml
$MYPATH/set_namelist time_length  $CUSTOM_SETTINGS_TIME_LENGTH
$MYPATH/set_namelist time_step  $CUSTOM_SETTINGS_TIME_STEP
$MYPATH/set_namelist coupling_freq $CUSTOM_SETTINGS_COUPLING_FREQ


