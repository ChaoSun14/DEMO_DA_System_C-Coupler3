#!/bin/bash

ensemble_idx=${1}
CASEROOT=${2}
RUN_DIR="${CASEROOT}/run"
if [ $ensemble_idx -gt 1 ]; then
    cp "${RUN_DIR}/ensemble_1/atm/wrf/exe/wrf" "${RUN_DIR}/ensemble_${ensemble_idx}/atm/wrf/exe/"
fi
