#!/bin/bash
export OUTPUT_FILENAME=test4
export RESTART_FILE_IN=test3.rs
if [[ $string == *"My long"* ]]; then
  echo "It's there!"
fi
RESTART=`cat makefile | awk '/^CPPFLAGS/ {print substr($0, index($0,$3))}' `
if [ -z "${RESTART##*dynRestart*}" ] && [ ! -f "$RESTART_FILE_IN" ] 
then
  echo "Clef cpp dynRestart mais le fichier '$RESTART_FILE_IN' n existe pas"
  exit
fi
export MODEL_TYPE=my_project_model
export OPA_NC_AUTHOR=`whoami`
#
export PATH=.:$PATH
export OPA_NC_SOURCE=libcdfpack4_3
export OPA_SOURCE=mf210_F90.0
export OPA_NC_ID=4
export OPA_CPP_KEYS=`cat makefile | awk '/^CPPFLAGS/ {print substr($0, index($0,$3))}' `
export NETCDF_OUTPUT_FILE=${OUTPUT_FILENAME}.nc
export DIAG_OUTPUT_FILE=${OUTPUT_FILENAME}.dg
export RESTART_FILE_OUT=${OUTPUT_FILENAME}.rs
rm -f $NETCDF_OUTPUT_FILE $DIAG_OUTPUT_FILE
trap '{ export EXITCODE="$?" ; }' INT
./main.exe
export EXITCODE="$?" 
echo '['`date -u +"%d/%m/%y %H:%M:%S UTC" `']' "$EXITCODE $MODEL_TYPE $NETCDF_OUTPUT_FILE $DIAG_OUTPUT_FILE $OPA_CPP_KEYS" >> simul.log
