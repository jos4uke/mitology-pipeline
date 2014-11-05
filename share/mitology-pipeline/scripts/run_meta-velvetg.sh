#!/usr/bin/env bash

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

# DEFAULTS
NAMESPACE_DEFAULT="MITOLOGY"
SCAFFOLD_DEFAULT="no"
CONFIG_SECTIONS=("contig_assembler" "scaffolder" "velveth" "velvetg" "meta_velvetg")

WORKING_DIR=$(pwd)
DATE=$(date '+%F_%Hh%Mm%Ss')
SESSION_ID=$(date '+%Y%m%d%H%M%S')
EXECUTED_COMMAND="${BASH_SOURCE[0]} $@"
PROG_NAME="$(basename ${BASH_SOURCE[0]})"
SESSION_TAG=${NAMESPACE_DEFAULT}_${USER}_${SESSION_ID}

DEBUGFILE=${SESSION_TAG}.log
ERROR_TMP_MODEL="/tmp/${PROG_NAME%.*}_error_${SESSION_TAG}.XXXXXX"
ERROR_TMP=$(mktemp "$ERROR_TMP_MODEL")
[[ $? -eq 0 ]] && echo -e "Execute ${PROG_NAME%.*} pipeline. Version: $VERSION"| tee -a $ERROR_TMP || exit 1

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
logger_addAppender ${PROG_NAME}.stderr2
appender_setType ${PROG_NAME}.stderr2 FileAppender
appender_file_setFile ${PROG_NAME}.stderr2 STDERR
appender_setLevel ${PROG_NAME}.stderr2 FATAL
appender_setLayout ${PROG_NAME}.stderr2 PatternLayout
appender_setPattern ${PROG_NAME}.stderr2 "%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] ${PROG_NAME} - %m"
appender_activateOptions ${PROG_NAME}.stderr2
# add and configure console appender that outputs to standard output
logger_addAppender ${PROG_NAME}.console2
appender_setType ${PROG_NAME}.console2 ConsoleAppender
appender_setLevel ${PROG_NAME}.console2 INFO
appender_setLayout ${PROG_NAME}.console2 PatternLayout
appender_setPattern ${PROG_NAME}.console2 "%d{HH:mm:ss,SSS} %-4rs [%F:%-5p] ${PROG_NAME} - %m"
appender_activateOptions ${PROG_NAME}.console2
appender_exists ${PROG_NAME}.console2 && logger_debug "Console appender is enabled." || logger_warn "Console appender was not enabled. Maybe a log4sh error occured."

### PRINT PROG VERSION AND COMMAND
echo "$PROG_NAME pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1 | logger_info
echo "Executed command: ${EXECUTED_COMMAND}" | tee -a $ERROR_TMP 2>&1 | logger_info

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

Usage: ${PROG_NAME} -c|--configfile CONFIG_FILE -o|--out_dir OUTPUT_DIR -P|--paired_end PAIRED_END [-S|--singletons SINGLETONS] [--scaffold no|yes] [-N|--namespace NAMESPACE] [--skip_config] [-d|--debug] [-e|--email_address VALID_EMAIL_ADDR] [-h|--help]

Mandatory:
-c|--config_file CONFIG_FILE            The user configuration file, CONFIG_FILE, listing the pipeline parameters by section.
                                        You can get a copy there: $PIPELINE_USER_CONFIG.
                                        Only the sections (and following parameters) listed here are mandatory:
                                        $(echo "${CONFIG_SECTIONS[@]}" | sed -e 's/ /, /g')
-o|--out_dir OUTPUT_DIR                 The output directory.
-P|--paired_end PAIRED_END              The interleaved paired end sequences.

Options:
-N|--namespace NAMESPACE                The namespace to use for the pipeline section/parameters.
-S|--singletons SINGLETONS              The singletons sequences file.
--scaffold no|yes                       Disable or enable the scaffolding process.
--skip_config                           Will skip loading the mandatory config file.
                                        Useful when the config was pre-loaded by a caller script.
-d|--debug                              Enable debugging mode in the console.
-e|--email_address VALID_EMAIL_ADDR     An optional but valid email address to send pipeline job/error status notifications
-h|--help                               Displays this message.

"
}

# Configure opts
### NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately;
CONFIGURE_OPTS=`getopt -o hc:N:P:S:o:e:d --long help,config_file:,namespace:,paired_end:,singletons:,out_dir:,debug,skip_config:,email_address: \
    -n 'run_meta-velvetg.sh' -- "$@"`

if [[ $? != 0 ]] ; then Usage >&2 ; exit 1 ; fi

# Note the quotes around `$CONFIGURE_OPTS'
eval set -- "$CONFIGURE_OPTS"

while true; do
    case "$1" in
        -h | --help ) Usage >&2; exit 1;;
        -c | --config_file ) CONFIGFILE="$2"; shift 2 ;;
        -o | --out_dir ) OUTPUT_DIR="$2"; shift 2 ;;
        -N | --namespace ) NAMESPACE="$2"; shift 2 ;;
        -P | --paired_end ) PAIRED_END="$2"; shift 2 ;;
        -S | --singletons ) SINGLETONS="$2"; shift 2;;
		--scaffold) $SCAFFOLD="$2"; shift 2;;
		--skip-config) $SKIP_CONFIG=true; shift;; 
        -d | --debug )
                    appender_setLevel console DEBUG;
                    appender_activateOptions console;
                    shift 1 ;;
        -e | --email_address ) EMAIL="$2"; shift 2 ;;
        -- ) shift; break ;;
        * ) break ;;
    esac
done

### VALIDATION ###

# mandatory
if [[ ! -s $CONFIGFILE ]]; then
    if [[ $SKIP_CONFIG ]]; then
		logger_info "Set skipping loading config file, $([[ -n $CONFIGFILE ]] && echo $CONFIGFILE)."
	else
		logger_fatal "Config file, $CONFIGFILE, does not exist or is empty. See Usage with --help option.";
    	exit 1;
	fi
else
	if [[ $SKIP_CONFIG ]]; then
		logger_info "Set skipping loading config file, $CONFIGFILE."
	fi
fi

if [[ -z $OUTPUT_DIR ]]; then
    logger_fatal "Output directory must be not null. See Usage with --help option.";
    exit 1;
fi

if [[ ! -s $PAIRED_END ]]; then
	logger_fatal "Interleaved paired end file, $PAIRED_END, does not exist or is empty. See Usage with --help option."
	exit 1;
fi

# not mandatory
if [[ -z $NAMESPACE ]]; then
	NAMESPACE=$NAMESPACE_DEFAULT
	logger_info "Use default namespace: $NAMESPACE"
fi

if [[ ! -s $SINGLETONS ]]; then
	logger_warn "Singletons file, SINGLETONS, does not exist or is empty. See Usage with --help option."
fi

if [[ "$SCAFFOLD" -eq "no" || "$SCAFFOLD" -eq "yes" ]]; then
	logger_info "The user scaffold value is now, ${SCAFFOLD}."
else
	logger_warn "The user scaffold value, $SCAFFOLD, is not recognized. The authorized values are no or yes."
	SCAFFOLD=$SCAFFOLD_DEFAULT
	logger_warn "Will use the default scaffold value, ${SCAFFOLD}."
fi





if [[ -z $ASSEMBLER_CLI ]]; then
    logger_warn "Assembler string must be not null. See Usage with --help option."
else
    # check if assembler is supported else use default
    if [[ ! $(elementIn "$ASSEMBLER_CLI" "${SUPPORTED_ASSEMBLERS[@]}") ]]; then
        logger_warn "Given assembler, $ASSEMBLER_CLI, is not currently supported. See Usage with --help option."
        exit 1
    else
        logger_info "Replace the default assembler, $ASSEMBLER_DEFAULT, by the user provided assembler, $ASSEMBLER_CLI ."
    fi
fi









#=====
# END
#=====

logger_info "[End] Run successfully the $PROG_NAME pipeline."
logger_info "[End] Will exit now."

# close all appenders
appender_exists ${PROG_NAME}.stderr2 && appender_close ${PROG_NAME}.stderr
appender_exists ${PROG_NAME}.console2 && appender_close ${PROG_NAME}.console2
appender_exists ${PROG_NAME}.debuggerF2 && appender_close ${PROG_NAME}.debuggerF2

exit 0

