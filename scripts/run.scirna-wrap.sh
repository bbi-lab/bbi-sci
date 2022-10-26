#!/bin/bash

#
# Path to the Nextflow processing run configuration file.
#
PWD=`pwd`
CONFIG_FILE="$PWD/wrap.config"

#
# Nextflow executable and pipeline script locations.
#
NEXTFLOW="/net/gs/vol1/home/bge/bin/nextflow"
NF_MAIN="/net/gs/vol1/home/bge/git/bbi-wrap/main.nf"

#
# Current date and time.
#
NOW=`date '+%Y%m%d_%H%M%S'`

#
# Get the path to the analyze output directory from
# the configuration file and set the Nextflow work
# directory to be in the analyze output directory.
#
ROOT_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.root_path"){print$2}}' | sed 's/"//g'`
WORK_DIR="$ROOT_DIR/work_wrap"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
# REPORT_FIL=$ROOT_DIR/wrap.report.${NOW}.html
TRACE_FIL=$ROOT_DIR/wrap.trace.${NOW}.tsv
# TIMELINE_FIL=$ROOT_DIR/wrap.timeline.${NOW}.html

#
# Nextflow run parameters.
#
# PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FIL -resume"

pushd $ROOT_DIR

date > ./run_start.${NOW}.txt

#
# Run Nextflow sci-RNA wrap pipeline.
#
$NEXTFLOW run $NF_MAIN $PARS

date > ./run_finish.${NOW}.txt

popd
