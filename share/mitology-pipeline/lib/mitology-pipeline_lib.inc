#!/usr/bin/env bash 
#
# MITOLOGY PIPELINE FUNCTIONS LIB
#
# author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
# date: 2014-10-13

# dependencies: bash-common_lib.inc 

# VERSION="dev"

#
# Process management
#

#----------------------------------------------------------
# Function: run_cli
# Purpose: wrap a command line and manage its execution
run_cli() {
	# usage
	run_cli_usage() { 
		printf %s "\
		Usage: run_cli [-h|--help] -c|--cmd COMMAND -t|--tag TAG -e|--cli_err COMMAND_ERR_FILE -E|--eval_err EVAL_ERR_FILE
		
		Mandatory:
		-c|--cmd COMMAND				The command line in COMMAND to execute
		-t|--tag TAG					The tag in TAG to report messages for log
		-e|--cli_err COMMAND_ERR_FILE	The standard error file in COMMAND_ERR_FILE to catch errors when executing the command line
		-E|--eval_err EVAL_ERR_FILE		The standard error file in EVAL_ERR_FILE to catch errors when evaluating the command line	

		Options:
		-h|--help						Displays this message.

	"
	}

	# opts
	local my_opts=`getopt -o hc:t:e:E: --long help,cmd:,tag:,cmd_err:,eval_err: \
		-n 'run_cli' -- "$@"`

	if [[ $? != 0 ]]; then run_cli_usage >&2; exit 1; fi

	# eval opts
	eval set -- "$my_opts"

	while true; do
		case "$1" in
			-h | --help ) run_cli_usage >&2; exit 1;;
			-c | --cmd ) local my_cli=$2; shift 2;;
			-t | --tag ) local my_tag=$2; shift 2;;
			-e | --cli_err ) local my_cli_err=$2; shift 2;;
			-E | --eval_err ) local my_eval_err=$2; shift 2;;
			-- ) shift; break;;
			* ) break;;
		esac
	done

	# validate
	if [[ -z $my_cli ]]; then logger_fatal "[run_cli] Command line string, $my_cli, must be not null. See Usage with --help option."; exit 1; fi
	if [[ -z $my_tag ]]; then logger_fatal "[run_cli] Tag string, $my_tag, must be not null. See Usage with --help option."; exit 1; fi
	if [[ ! -d $(dirname $my_cli_err) ]]; then logger_fatal "[run_cli] Parent directory of the command line error file, $my_cli_err, must exist. See Usage with --help option."; exit 1; fi
	if [[ ! -d $(dirname $my_eval_err) ]]; then logger_fatal "[run_cli] Parent directory of the eval command line error file, $my_eval_err, must exist. See Usage with --help option."; exit 1; fi
	
	# cli prefix to remove
	local prefix="cd*;[[:space:]]"
	
	#
	# run
	#
	# logging calls resolution depend on parent script which is supposed to load log4sh
	logger_debug "[$my_tag] $my_cli"
	eval  "$my_cli" 2>$my_eval_err
	local pid=$!
	local rtrn=$?
	local eval_failed_msg="[$my_tag] An error occured while eval $my_cli command line."
	# GLOBALS: $OUTPUT_DIR, $LOG_DIR, $DEBUGFILE, $SESSION_TAG, $EMAIL, $PIDS_ARR
	exit_on_error "$my_eval_err" "$eval_failed_msg" "$rtrn" "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
	
	# add pid to array
	local PIDS_ARR=()
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
	# wait until command line process finish and reinit pid array
	for p in "${PIDS_ARR[@]}"; do
		local pid_list_failed_msg="[$my_tag] Failed getting process status for process $p."
		logger_trace "[$my_tag] $(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>$my_eval_err)"
		local rtrn=$?
		exit_on_error "$my_eval_err" "$pid_list_failed_msg" "$rtrn" "$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
	done
	logger_info "[$my_tag] Wait for all $(echo "${my_cli##$prefix}" | cut -d" " -f1) processes to finish before proceeding to next step."
	waitalluntiltimeout "${PIDS_ARR[@]}" 2>/dev/null
	if [[ -s $my_cli_err ]]; then
		logger_warn "[$my_tag] Some messages were thrown to standard error while executing $my_cli . See $my_cli_err for more details."
	fi
	logger_info "[$my_tag] All $(echo "$my_cli" | cut -d" " -f1) processes finished."
	PIDS_ARR=()
}


