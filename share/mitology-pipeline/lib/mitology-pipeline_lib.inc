#!/usr/bin/env bash 
#
# MITOLOGY PIPELINE FUNCTIONS LIB
#
# author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
# date: 2014-10-13

# dependencies: bash-common_lib.inc 

# VERSION="dev"

#
# Array utils
#

#----------------------------------------------------------
# Function: elementIn 
# Purpose: Tells if the given string is in the array and returns the matched value 
elementIn() {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && echo "$1"; return 0; done
    return 1
}

#----------------------------------------------------------
# Function: shortenAndSortOptions
# Purpose: Will shorten options (-option) output from buildCommandLineOptions function (--option)
#          as velveth/velvetg/meta-velvetg required it (not ANSI compliant), and sort them alphabetically
#          before returning the array of options
shortenAndSortOptions() {
	echo "${@}" | sed -e 's/^--/-/' -e 's/ --/\n-/g' | sort
}

#
# File utils
#

#----------------------------------------------------------
# Function: createDir
# Purpose: will create a directory given a name
createDir() {
	# usage 
	run_cli_usage() {
		printf %s "\
		Usage: createDir [-h|--help] -n|name DIRNAME -t|--tag TAG -e|--err_file ERROR_FILE [-d|--debug] 

		Mandatory:
		-n|--name DIRNAME               The name of directory to create.
		-t|--tag TAG                    The tag in TAG to report messages for log
		-e|--err_file ERROR_FILE        The error file in ERROR_FILE to catch errors while creating the directory

		Options:
		-d|--debug                      Enable debug logging by piping log4sh logger_debug to the stdout
		-h|--help                       Displays this message.

		"
	}

	# opts
	local my_opts=`getopt -o hn:t:e:d --long help,name:,tag:,err_file:,debug \
					-n 'createDir' -- "$@"`	
	
	if [[ $? != 0 ]]; then run_cli_usage >&2; exit 1; fi

	# eval opts
	eval set -- "$my_opts"

	while true; do
		case "$1" in
			-h | --help ) run_cli_usage >&2; exit 1;;
			-n | --name ) local my_dirname=$2; shift 2;;
			-t | --tag ) local my_tag="run_$2"; shift 2;;
			-e | --err_file ) local my_errfile=$2; shift 2;;
			-d | --debug ) local my_DEBUG=true; shift ;;
			-- ) shift; break;;
			* ) break;;
		esac
	done
	
	# validate
	if [[ -z $my_dirname ]]; then logger_fatal "[createDir] The given DIRNAME, ${my_dirname}, must be not null. See Usage with --help option."; exit 1; fi
	if [[ -z $my_tag ]]; then logger_fatal "[createDir] Tag string, $my_tag, must be not null. See Usage with --help option."; exit 1; fi
	if [[ ! -d $(dirname $my_errfile) ]]; then logger_fatal "[createDir] Parent directory of the command line error file, $my_errfile, must exist. See Usage with --help option."; exit 1; fi	

	# current directory
	local CWD=$(pwd)
	
	# errfile path
	local errDir=$(dirname $my_errfile)
    if [[ ! -d $errDir ]]; then
        my_cmd_err=$CWD/$my_errfile
    fi

	# run mkdir
	local my_cli="mkdir -p $my_dirname 2>$my_errfile"
	if [[ ! -d $my_dirname ]]; then
		eval "$my_cli"
		rtrn=$?
		mkdir_failed_msg="[$my_tag] Failed. Cannot create directory named $my_dirname"
		[[ "$rtrn" -ne 0 ]] && logger_fatal "$mkdir_failed_msg"
		exit_on_error "$my_errfile" "$mkdir_failed_msg" $rtrn "" $SESSION_TAG $EMAIL
		[[ $my_DEBUG == "true" ]] && logger_debug "[$my_tag] OK. $my_dirname directory was created successfully. Will output all $my_tag files in this directory."
	else
		logger_debug "[$my_tag] OK. $my_dirname directory already exists. Will output all $my_tag files into this directory."
	fi
}

#
# Process management
#

#----------------------------------------------------------
# Function: run_cli
# Purpose: wrap one command line and manage its execution
run_cli()
{
 # usage
    run_cli_usage() {
        printf %s "\
        Usage: run_cli [-h|--help] -c|--cmd COMMAND -t|--tag TAG -e|--cmd_err COMMAND_ERR_FILE [-d|--debug]

        Mandatory:
        -c|--cmd COMMAND                The command line in COMMAND to execute
        -t|--tag TAG                    The tag in TAG to report messages for log
        -e|--cmd_err COMMAND_ERR_FILE   The standard error file in COMMAND_ERR_FILE to catch errors while executing the command line
        -d|--debug                      enable debug logging by piping log4sh logger_debug to the COMMAND stdout

        Options:
        -h|--help                       Displays this message.

    "
    }

    # opts
    local my_opts=`getopt -o hc:t:e:d --long help,cmd:,tag:,cmd_err:,debug \
        -n 'run_cli' -- "$@"`

    if [[ $? != 0 ]]; then run_cli_usage >&2; exit 1; fi

	# eval opts
    eval set -- "$my_opts"

    while true; do
        case "$1" in
            -h | --help ) run_cli_usage >&2; exit 1;;
            -c | --cmd ) local my_cli=$2; shift 2;;
            -t | --tag ) local my_tag="run_$2"; shift 2;;
            -e | --cmd_err ) local my_cmd_err=$2; shift 2;;
            -d | --debug ) local my_DEBUG=true; shift ;;
            -- ) shift; break;;
            * ) break;;
        esac
    done

    # validate
    if [[ -z $my_cli ]]; then logger_fatal "[run_cli] Command line string, $my_cli, must be not null. See Usage with --help option."; exit 1; fi
    if [[ -z $my_tag ]]; then logger_fatal "[run_cli] Tag string, $my_tag, must be not null. See Usage with --help option."; exit 1; fi
    if [[ ! -d $(dirname $my_cmd_err) ]]; then logger_fatal "[run_cli] Parent directory of the command line error file, $my_cmd_err, must exist. See Usage with --help option."; exit 1; fi

    # current directory
    local CWD=$(pwd)

	# cmd/outDir
    local cmd=$(basename $(echo -e "${my_cli##$prefix}" | cut -d" " -f1))
    local outDir="$(dirname $my_cmd_err)"
    if [[ -d $outDir ]]; then
        local cmd_pid_file="$outDir/${cmd}.pid"
    else
        local outDir="$CWD/$outDir"
        local cmd_pid_file="$outDir/${cmd}.pid"
    fi

    local errDir=$(dirname $my_cmd_err)
    if [[ ! -d $errDir ]]; then
        my_cmd_err=$CWD/$my_cmd_err
    fi

	#
    # run
    #
    # logging calls resolution depend on parent script which is supposed to load log4sh
    logger_debug "[$my_tag] $my_cli"

    ##if [[ "$my_DEBUG" == true ]]; then
    ##   $my_cli 2>$CWD/$my_cmd_err &
    ##   pid=$!
    ##else
	#($my_cli 2>$my_cmd_err & echo $! >&3) 3>$cmd_pid_file | logger_debug 
	eval "$my_cli 2>$my_cmd_err &"
	pid=$!
	##pid=$!
    ##fi

	#pid=$(cat $cmd_pid_file)
    logger_debug "$cmd pid: $pid"
    echo $pid >$cmd_pid_file
    logger_trace "$cmd pid: $pid saved to file $cmd_pid_file"

	# GLOBALS: $OUTPUT_DIR, $LOG_DIR, $DEBUGFILE, $SESSION_TAG, $EMAIL, $PIDS_ARR
    #exit_on_error "$CWD/$my_cmd_err" "$eval_failed_msg" "$rtrn" "$CWD/$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"

    # add pid to array
    #local cmd_ps=$(ps aux | grep $USER | grep -v grep | grep "$(echo -e "${my_cli##$prefix}" | cut -d" " -f1-4)")
    #[[ -n $cmd_ps ]] && logger_trace "$cmd ps: $cmd_ps"
    #[[ -z $cmd_ps ]] && logger_trace "$cmd ps: $(ps aux | grep $USER | grep -v grep | grep $pid)"
    local PIDS_ARR=()
    PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
    ## wait until command line process finish and reinit pid array
    for p in "${PIDS_ARR[@]}"; do
        local pid_list_failed_msg="[$my_tag] Failed getting process status for process $p."
        logger_trace "[$my_tag] $(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>$my_cmd_err)"
        local rtrn=$?
        exit_on_error "$my_cmd_err" "$pid_list_failed_msg" "$rtrn" "$CWD/$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
    done
    logger_debug "[$my_tag] Wait for all $cmd processes to finish before proceeding to next step."
    local waitall_file="$outDir/${cmd}.waitall.$(date '+%Y%m%d%H%M%S').log"
    waitalluntiltimeout "${PIDS_ARR[@]}" 2>"$waitall_file"
    wait $pid
    local status=$?
    logger_debug "$cmd pid status: $status"
    cmd_status_file=${cmd_pid_file%.pid}.exit-status
    echo $status >$cmd_status_file
    logger_debug "$cmd pid status: $status saved to $cmd_status_file"
    if [[ -s $my_cmd_err ]]; then
        logger_warn "[$my_tag] Some messages were thrown to standard error while executing $my_cli . See $my_cmd_err for more details."
    fi
    logger_debug "[$my_tag] All $cmd processes finished."
    PIDS_ARR=()

	if [[ -n last_line=$(tail -1 "$waitall_file" | grep " $pid exited with non-zero exit status") ]]; then
        pid_status_error_msg="[$my_tag] $last_line"
        logger_warn "$pid_status_error_msg"
    else
        logger_debug "[$my_tag] $(tail -1 "$waitall_file" | grep " $pid exited with zero exit status")"
    fi
    return $status
}

#----------------------------------------------------------
# Function: run_one_cli
# Purpose: wrap one command line and manage its execution
run_one_cli() {
	# usage
	run_cli_usage() { 
		printf %s "\
		Usage: run_one_cli [-h|--help] -c|--cmd COMMAND -t|--tag TAG -e|--cmd_err COMMAND_ERR_FILE
		
		Mandatory:
		-c|--cmd COMMAND                The command line in COMMAND to execute
		-t|--tag TAG                    The tag in TAG to report messages for log
		-e|--cli_err COMMAND_ERR_FILE   The standard error file in COMMAND_ERR_FILE to catch errors while executing the command line

		Options:
		-h|--help						Displays this message.

	"
	}

	# opts
	local my_opts=`getopt -o hc:t:e: --long help,cmd:,tag:,cmd_err: \
		-n 'run_one_cli' -- "$@"`

	if [[ $? != 0 ]]; then run_cli_usage >&2; exit 1; fi

	# eval opts
	eval set -- "$my_opts"

	while true; do
		case "$1" in
			-h | --help ) run_cli_usage >&2; exit 1;;
			-c | --cmd ) local my_cli=$2; shift 2;;
			-t | --tag ) local my_tag="run_$2"; shift 2;;
			-e | --cmd_err ) local my_cmd_err=$2; shift 2;;
			-- ) shift; break;;
			* ) break;;
		esac
	done

	# validate
	if [[ -z $my_cli ]]; then logger_fatal "[run_cli] Command line string, $my_cli, must be not null. See Usage with --help option."; exit 1; fi
	if [[ -z $my_tag ]]; then logger_fatal "[run_cli] Tag string, $my_tag, must be not null. See Usage with --help option."; exit 1; fi
	if [[ ! -d $(dirname $my_cmd_err) ]]; then logger_fatal "[run_cli] Parent directory of the command line error file, $my_cmd_err, must exist. See Usage with --help option."; exit 1; fi
	
	# current directory
	local CWD=$(pwd)

	# cli prefix to remove
	local prefix="cd*;[[:space:]]*"

	# cmd/outDir
    local cmd=$(basename $(echo -e "${my_cli##$prefix}" | cut -d" " -f1))
    local outDir="$(dirname $my_cmd_err)"
    if [[ -d $outDir ]]; then
        local cmd_pid_file="$outDir/${cmd}.pid"
    else
        local outDir="$CWD/$outDir"
        local cmd_pid_file="$outDir/${cmd}.pid"
    fi

	local errDir=$(dirname $my_cmd_err)
	if [[ ! -d $errDir ]]; then
		my_cmd_err=$CWD/$my_cmd_err	
	fi
	
	#
	# run
	#
	# logging calls resolution depend on parent script which is supposed to load log4sh
	logger_debug "[$my_tag] $my_cli"
	
	($my_cli 2>$my_cmd_err & echo $! >&3) 3>$cmd_pid_file | logger_debug
	#pid=$!

	pid=$(cat $cmd_pid_file)
	logger_debug "$cmd pid: $pid"
	#echo $pid >$cmd_pid_file
	logger_trace "$cmd pid: $pid saved to file $cmd_pid_file"

	# GLOBALS: $OUTPUT_DIR, $LOG_DIR, $DEBUGFILE, $SESSION_TAG, $EMAIL, $PIDS_ARR
	#exit_on_error "$CWD/$my_cmd_err" "$eval_failed_msg" "$rtrn" "$CWD/$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
	
	# add pid to array
	#local cmd_ps=$(ps aux | grep $USER | grep -v grep | grep "$(echo -e "${my_cli##$prefix}" | cut -d" " -f1-4)")
	#[[ -n $cmd_ps ]] && logger_trace "$cmd ps: $cmd_ps"
	#[[ -z $cmd_ps ]] && logger_trace "$cmd ps: $(ps aux | grep $USER | grep -v grep | grep $pid)"
	local PIDS_ARR=()
	PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
	## wait until command line process finish and reinit pid array
	for p in "${PIDS_ARR[@]}"; do
		local pid_list_failed_msg="[$my_tag] Failed getting process status for process $p."
		logger_trace "[$my_tag] $(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>$my_cmd_err)"
		local rtrn=$?
		exit_on_error "$my_cmd_err" "$pid_list_failed_msg" "$rtrn" "$CWD/$OUTPUT_DIR/$LOG_DIR/$DEBUGFILE" "$SESSION_TAG" "$EMAIL"
	done
	logger_debug "[$my_tag] Wait for all $cmd processes to finish before proceeding to next step."
	#local waitall_file="$outDir/${cmd}.waitall.$(date '+%Y%m%d%H%M%S').log"
	#waitalluntiltimeout "${PIDS_ARR[@]}" 2>"$waitall_file"
	wait $pid
	local status=$?
	logger_debug "$cmd pid status: $status"
	cmd_status_file=${cmd_pid_file%.pid}.exit-status
	echo $status >$cmd_status_file
	logger_debug "$cmd pid status: $status saved to $cmd_status_file"
	if [[ -s $my_cmd_err ]]; then
		logger_warn "[$my_tag] Some messages were thrown to standard error while executing $my_cli . See $my_cmd_err for more details."
	fi
	logger_debug "[$my_tag] All $cmd processes finished."
    PIDS_ARR=()

#	last_line=$(tail -1 "$waitall_file" | grep " $pid exited with non-zero exit status")
#	if [[ -n last_line=$(tail -1 "$waitall_file" | grep " $pid exited with non-zero exit status") ]]; then
#		pid_status_error_msg="[$my_tag] $last_line"
#		logger_warn "$pid_status_error_msg"	
#	else
#		logger_debug "[$my_tag] $(tail -1 "$waitall_file" | grep " $pid exited with zero exit status")"
#	fi
	return $status
}


