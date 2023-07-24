#!/bin/bash

#
# Path to the Nextflow processing run configuration file.
#
PWD=`pwd`
CONFIG_FILE="$PWD/experiment.config"

#
# Nextflow executable and pipeline script locations.
# Note: set the paths in the three variables below.
#
NEXTFLOW="$HOME/bin/nextflow"
NF_HOME="$HOME/git/bbi-sci"
NF_MAIN="${NF_HOME}/main.nf"

#
# Current date and time.
#
NOW=`date '+%Y%m%d_%H%M%S'`

#
# Get the path to the analyze output directory from
# the configuration file and set the Nextflow work
# directory to be in the analyze output directory.
#
OUTPUT_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.output_dir"){print$2}}' | sed 's/"//g'`
ANALYZE_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.demux_out"){print$2}}' | sed 's/"//g'`
WORK_DIR="$ANALYZE_DIR/work_analyze"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
# REPORT_FIL=$ANALYZE_DIR/analyze.report.${NOW}.html
TRACE_FIL=$ANALYZE_DIR/analyze.trace.${NOW}.tsv
# TIMELINE_FIL=$ANALYZE_DIR/analyze.timeline.${NOW}.html

#
# Nextflow run parameters.
#
# PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FIL -resume"

mkdir -p $ANALYZE_DIR
pushd $ANALYZE_DIR

date > ./run_start.${NOW}.txt

#
# Run Nextflow sci-RNA analyze pipeline.
#
$NEXTFLOW run $NF_MAIN $PARS

date > ./run_finish.${NOW}.txt

popd

#
# Set run directory file and directory permissions.
#
${NF_HOME}/scripts/set_run_permissions.sh ${OUTPUT_DIR}

