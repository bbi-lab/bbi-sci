#!/bin/bash

#
# The master copy of this script is in bbi-dmux/scripts. Please make changes to
# that copy and copy that file to all repositories that use it:
#   bbi-dmux
#   bbi-sci
#   bbi-wrap
#   bbi-sciatac-demux
#   bbe-sciatac-analyze
#   bbi-sciatac-wrap
#

#
# Executable files.
# This is a regex alternation expression.
#
exec_files='(run.scirna-analyze.sh|run.scirna-demux.sh|run.sciatac-demux.sh|run.sciatac-analyze.sh|set_run_permissions.sh)'

run_path="$1"
echo "Change permissions in path $run_path"
echo

# Set directory permissions.
find "$run_path" -type d -print0 | xargs -0 chmod -R 770 

# Set file permissions.
find "$run_path" -type f -print0 | egrep -z -v "${exec_files}" | xargs -0 chmod -R 660
find "$run_path" -type f -print0 | egrep -z    "${exec_files}" | xargs -0 chmod -R 770


