

#
# run_meta-velvetg.sh
#

# Authors: Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This script provides a wrapper around meta-velvetg assembler pipeline
# GRANT: ANR BIOADAPT CYTOPHENO 2012

# based on the preliminary work of 2 interns (03-07/2014): Alexandrina Brodrug and Myriam Shafie

# This software is governed by the CeCILL license, Version 2.0 (the "License"), under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license, Version 2.0 (the "License"), as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license, Version 2.0 (the "License"), and that you accept its terms.

# Date: 2014-11-02

declare -r VERSION="dev"

########################
# SECTION CONFIGURATION
#######################

### SESSION VARIABLES ###

NAMESPACE="MITOLOGY"

WORKING_DIR=$(pwd)
DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%m%d%H%M%S')
EXECUTED_COMMAND="${BASH_SOURCE[0]} $*"
PROG_NAME="$(basename ${BASH_SOURCE[0]})"
SESSION_TAG=${NAMESPACE}_${USER}_${SESSION_ID}

LOG_DIR="log"
DEBUGFILE=${SESSION_TAG}.log
ERROR_TMP_MODEL="/tmp/${PROG_NAME%.*}_error_${SESSION_TAG}.XXXXXX"
ERROR_TMP=$(mktemp "$ERROR_TMP_MODEL")
[[ $? -eq 0 ]] && echo "Execute ${PROG_NAME%.*}"| tee -a $ERROR_TMP || exit 1

[[ $VERSION -eq "dev" ]] && PROG_PATH=$(realpath $(dirname ${BASH_SOURCE[0]}));PIPELINE_USER_CONFIG=$(find ${PROG_PATH} -iname "mitology-pipeline_user.config") || PIPELINE_USER_CONFIG=/usr/local/share/mitology-pipeline/etc/mitology-pipeline_user.config

PIDS_ARR=()
WAITALL_TIMEOUT=259200
WAITALL_INTERVAL=60
WAITALL_DELAY=60

### LOGGING CONFIGURATION ###

# load log4sh (disabling properties file warning) and clear the default
# configuration
LOG4SH_CONFIGURATION='none' . /usr/local/share/log4sh/build/log4sh 2>/dev/null
[[ $? != 0 ]] && $(echo "Error loading log4sh lib" >&2; exit 1)
log4sh_resetConfiguration

# set the global logging level
# add and configure a FileAppender that outputs to STDERR
logger_addAppender stderr
appender_setType stderr FileAppender
appender_file_setFile stderr STDERR
appender_setLevel stderr FATAL
appender_setLayout stderr PatternLayout
appender_setPattern stderr '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions stderr
# add and configure console appender that outputs to standard output
logger_addAppender console
appender_setType console ConsoleAppender
appender_setLevel console INFO
appender_setLayout console PatternLayout
appender_setPattern console '%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] %t - %m'
appender_activateOptions console
appender_exists console && logger_debug "Console appender is enabled." || logger_warn "Console appender was not enabled. Maybe a log4sh error occured."

### LOAD LIB ###

# bash-common lib
[[ $VERSION == "dev" ]] && LIB_PATH=$(realpath $(dirname ${BASH_SOURCE[0]}))/../../../../bash-common/share/bash-common/lib/bash-common_lib.inc || LIB_PATH=/usr/local/share/bash-common/lib/bash-common_lib.inc

logger_debug "[Library] Loading $LIB_PATH"
. $LIB_PATH
if [[ $? -ne 0 ]]; then
    logger_fatal "Error loading bash common lib: $LIB_PATH"
    exit 1
fi

# mitology-pipeline lib
[[ $VERSION == "dev" ]] && LIB_PATH=$(realpath $(dirname $0))/../share/mitology-pipeline/lib/mitology-pipeline_lib.inc || LIB_PATH=/usr/local/share/mitology-pipeline/lib/mitology-pipeline_lib.inc

logger_debug "[Library] Loading $LIB_PATH"
. $LIB_PATH
if [[ $? -ne 0 ]]; then
    logger_fatal "Error loading mitology pipeline lib: $LIB_PATH"
    exit 1
fi

### USAGE ###
Usage()
{
printf %s "\
Program: ${PROG_NAME}
Version: $VERSION

Copyright 2014 Joseph Tran <Joseph.Tran@versailles.inra.fr>

Licensed under the CeCILL License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Usage: ${PROG_NAME} -c|--configfile CONFIG_FILE -o|--out_dir OUTPUT_DIR [-N|--namespace NAMESPACE] [-P|--paired_end PAIRED_END] [-S|--singletons SINGLETONS] [--scaffold no|yes] [--skip_config] [-d|--debug] [-e|--email_address VALID_EMAIL_ADDR] [-h|--help]

Mandatory:
-c|--config_file CONFIG_FILE            The user configuration file, CONFIG_FILE, listing the pipeline parameters by section.
                                        You can get a copy there: $PIPELINE_USER_CONFIG.
                                        Only the sections (and following parameters) listed here are mandatory:
                                        - contig_assembler
                                        - scaffolder
                                        - velveth
                                        - velvetg
                                        - meta_velvetg 
-o|--out_dir OUTPUT_DIR                 The output directory.

Options:
-N|--namespace NAMESPACE                The namespace to use for the pipeline section/parameters.
-P|--paired_end PAIRED_END              The interleaved paired end sequences.
-S|--singletons SINGLETONS              The singletons sequences file.
--scaffold no|yes                       Disable or enable the scaffolding process.
--skip_config                           Will skip loading the mandatory config file.
                                        Useful when the config was pre-loaded by a caller script.
-d|--debug                              Enable debugging mode in the console.
-e|--email_address VALID_EMAIL_ADDR     An optional but valid email address to send pipeline job/error status notifications
-h|--help                               Displays this message.

"
}











#=====
# END
#=====

logger_info "[End] Run successfully the $PROG_NAME pipeline."
logger_info "[End] Will exit now."

# close all appenders
appender_exists stderr && appender_close stderr
appender_exists console && appender_close console
appender_exists debuggerF && appender_close debuggerF

exit 0

