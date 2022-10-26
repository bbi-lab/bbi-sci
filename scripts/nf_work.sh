#!/bin/bash

#
# Script: nf_work.sh 
#
# Description:
#   Extract lines from files in the work directory.
#   Steps:
#     o  run Nextflow with the command line option '-with-trace <trace_file_name>'.
#     o  edit this file to set
#          o  trace file path
#          o  work directory created by Nextflow
#          o  process block to examine
#          o  the name of the files to search
#          o  the (grep) search string
#


TRACE_FILE="/net/bbi/vol1/data/bge/bbi/tests/nobackup/nf.RNA3-014-a.centos_7.test_10/demux_out/analyze.trace.tsv"
WORK_PARENT="/net/bbi/vol1/data/bge/bbi/tests/nobackup/nf.RNA3-014-a.centos_7.test_10/demux_out/work"


PROCESS_BLOCK="make_cds"
SEARCH_FILE=".command.sh"
FILTER_STRING="barnyard"


if [ "${PROCESS_BLOCK}" ]; then
  WORK_DIR_LIST=`grep -v "^task_id" ${TRACE_FILE} | awk 'BEGIN{FS="\t"}{if($4 ~ /'${PROCESS_BLOCK}'/){print$2}}'`
else
  WORK_DIR_LIST=`grep -v "^task_id" ${TRACE_FILE} | awk 'BEGIN{FS="\t"}{print$2}'`
fi


for DIR in ${WORK_DIR_LIST}
do
  DPATH=`compgen -o dirnames -- ${WORK_PARENT}/$DIR`
  echo "== $DPATH  file: ${SEARCH_FILE}"
  grep "${FILTER_STRING}" "${DPATH}/${SEARCH_FILE}"
  echo
done

