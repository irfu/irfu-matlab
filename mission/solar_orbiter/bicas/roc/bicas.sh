#!/bin/bash
#====================================================================================================
# Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
# First created 2016-04-xx
#
# The primary purpose of this bash script is to call a MATLAB script.
# See the MATLAB script called for more documentation.
#
#
# 
# PARAMETERS: [--developer-bash-settings] <<--log <dir>; Other arguments>>
#
# --developer-bash-settings : Unofficial flag for using particular settings.
# --log <dir>               : Path to directory where log file will be created. 
# Other arguments           : Argument parsed in the MATLAB code.
#
#
#
# NOTE: The script reads the "--log" parameter value and pipes "all" stdout to a log file. 
# It also uses a prefix to filter out some lines that are piped to the actual stdout.
# In addition, a symlink (with a constant filename) to the current log file is created.
# Overwriting the old symlink might technically not be compliant with the RCS ICD.
# 
# NOTE: The script is designed to handle important basic errors, such as
# 1) Can not find MATLAB executable (not explicitly, but implicitly)
# 2) MATLAB can not find the MATLAB script.
# 
# NOTE: The script is designed to return the same exit codes as the MATLAB code.
# 
# 
# 
# IMPLEMENTATION NOTE: Reasons to pipe MATLAB's stdout to the log file
# from this bash script rather than from MATLAB code:
# 1) can handle stdout from the MATLAB startup process
# 2) can handle stdout from library MATLAB code that prints to stdout (code that should not be altered, e.g. in irfu-matlab, e.g. irf.log)
# 3) avoid manually having to pipe everything to log file in MATLAB code.
# 
# IMPLEMENTATION NOTE: The same arguments that are passed to this script are passed on to the MATLAB code, including
# arguments which are unnecessary there, including --log <dir> and --developer-bash-settings.
# PRO (Reasons for):
# 1) It is best to let as much parameter-parsing (incl. syntax checking) code to be in one place as possible (the MATLAB code).
# 2) Enable the MATLAB code to print all the command-line parameters in error messages and to the log file.
# 3) Enable the MATLAB code to print its syntax which is then identical to the syntax of this script.
# CON: The MATLAB code needs some knowledge of the parameter parsing in this script.
#
# IMPLEMENTATION NOTE: The name of this file may only contain alphanumeric characters. Can therefore not use 
# file extension ".sh" to indicate that it is a shell script.
# /ROC-TST-GSE-ICD-00023-LES, iss02rev01, Section 3.1.
# 
#====================================================================================================
# TODO: Redirect (copy) all stderr to stdout and hence to the log file.
# QUESTION: If log file already exists, overwrite, amend, or error?
#====================================================================================================



# Search pattern that determines which stdout (in particular
# from MATLAB) that is actually piped to stdout.
# NOTE: This pattern must match the corresponding prefix used in MATLAB.
declare -r stdout_pattern='^STDOUT: '    # NOTE: Includes "^" = beginning of line.

# Name of script such as it is called from inside MATLAB, i.e. without file extension.
declare -r matlab_script_name="bicas"







#======================================
# Find argument for the log directory.
#======================================
for ((i=1; i<$#; i++))                   # NOTE: Omit the last argument.
do
    if [[ "${!i}" == "--log" ]]
    then
        let j=i+1
        declare -r log_dir="${!j}"
        break
    fi
done



#=====================================
# Construct path and name of log file
#=====================================
if [[ ! -z "$log_dir" ]]
then
    if [[ ! -d "$log_dir" ]]
    then
        echo "Can not find specified log directory \"$log_dir\"." >> /dev/stderr
        exit 1
    fi
    base=$(basename "$0")
    timestamp=$(date +%Y-%m-%d_%H.%M.%S)   # LOCAL time, not UTC.
    log_file="${log_dir}/${base}_${timestamp}.log"
    log_file_symlink="${log_dir}/${base}_latest.log"
    
    
    
    #=====================================================
    # Test if can create (write) to the proposed log file
    #=====================================================
    touch "$log_file"
    if [[ $? != 0 ]]
    then
        echo "Can not write to log file \"$log_file\"." >> /dev/stderr
        exit 1
    fi



    #=================================================================
    # Create symlink to latest/current log.
    # NOTE: This may technically be in violation of the RSC ICD 
    # when removing/overwriting an old symlink.
    # 
    # test options:
    #        -h FILE
    #               FILE exists and is a symbolic link (same as -L)
    # 
    # ln options:
    #        -r, --relative
    #               create symbolic links relative to link location
    #=================================================================
    if [[ -h "$log_file_symlink" ]]
    then
        rm -v "$log_file_symlink" >> "$log_file"
    fi
    ln -srv "$log_file" "$log_file_symlink" >> "$log_file"

else
    log_file='/dev/null'
fi



# DEBUG
# echo log_dir  = "$log_dir"
# echo log_file = "$log_file"
# echo log_file_symlink = "$log_file_symlink"



echo "Command-line arguments: $@" >> "$log_file"



launch_matlab_code() {
    
    # Find (absolute) paths by using the path to this executable + relative paths.
    local temp=$(dirname "$0")/..
    local root_path=$(readlink -f "$temp")
    local matlab_src_path="${root_path}/src/"
    local irfumatlab_path="${root_path}/lib/irfu-matlab/"
    local matlab_command="matlab"
#     local matlab_command="matlabXXX"   # DEBUG



    if [[ "$1" == "--developer-bash-settings" ]]
    then
        #============================================================================
        # This option exists to enable working with customized settings on 
        # a development machine (at IRF-U).
        #
        # Explanation for setting irfumatlab_path:
        # (1) irfu-matlab is a general git repository that contains
        #     the RPW BIAS RCS-specific code in a subdirectory.
        # (2) When the software is run in the ROC pipeline, irfu-matlab has to be
        #     a subdirectory due to this software's prescribed tree structure.
        # (3) I want to avoid using symbolic links inside a git repository.
        #============================================================================
        #shift   # Remove the flag from the argument list!
        echo "===================================="
        echo "WARNING!   USING DEVELOPER SETTINGS!"
        echo "===================================="
        irfumatlab_path=$(readlink -f "$root_path"/../../../)
        local matlab_command="/home/data/MATLAB/R2016a/bin/matlab"
#         local matlab_command="/home/data/MATLAB/R2016a/bin/matlabXXX"    # DEBUG
    fi
    
    # DEBUG
    echo "temp            = \"$temp\""
    echo "root_path       = \"$root_path\""
    echo "matlab_src_path = \"$matlab_src_path\""
    echo "irfumatlab_path = \"$irfumatlab_path\""
    echo "matlab_command  = \"$matlab_command\""
    
    

    #=======================================================================
    # Convert arguments into one long string that MATLAB code
    # can interpret as a list of arguments.
    # NOTE: Can not handle arguments containing the character single quote.
    #=======================================================================
    matlab_arg_list=""
    for arg in "$@"
    do
        if [[ -z "$matlab_arg_list" ]]   # if string is empty
        then
            matlab_arg_list="'$arg'"
        else
            matlab_arg_list="$matlab_arg_list, '$arg'"
        fi
    done



    #=============================================================================
    # MATLAB options:
    #
    #     -r MATLAB_command    - Start MATLAB and execute the MATLAB_command.
    #     -nodesktop           - Do not start the MATLAB desktop. Use the current
    #                            terminal for commands. The Java virtual machine
    #                            will be started.
    #     -nojvm               - Shut off all Java support by not starting the
    #                            Java virtual machine. In particular the MATLAB
    #                            desktop will not be started.
    #     -nosplash            - Do not display the splash screen during startup.
    #
    # NOTE: The MATLAB code supplied here is to make sure that the code works well in case of failure.
    # 1) MATLAB will NOT EXIT unless the code itself calls "exit".
    # 2) MATLAB hangs if it can not find the MATALB script.
    # NOTE: It is very important to end lines with semicolons and backslash for
    # the command line-supplied MATLAB code to not yield syntax error AND NOT HANG.
    # 
    # NOTE: Do not change directory, so that the current working directory is the same
    # in MATLAB as here. This is necessary for MATLAB code to correctly interpret relative paths (arguments).
    #=============================================================================
    echo "Launching MATLAB, using \"${matlab_command}\"."
    "${matlab_command}" -nodisplay -r "\
        try;\
            addpath(genpath('${irfumatlab_path}'), '-begin');\
            addpath(genpath('${matlab_src_path}'), '-begin');\
            error_code = $matlab_script_name($matlab_arg_list);\
            quit(error_code);\
        catch e;\
            msg = sprintf('Could not launch the MATLAB script. The reason MIGHT be that MATLAB could not find it.\nException message: %s\n', e.message);\
            fprintf(1, msg);\
            fprintf(2, msg);\
            quit(1);\
        end"
    matlab_error_code=$?
    echo "Bash launch function: matlab_error_code = $matlab_error_code"    # NOTE: Printed to stdout ==> Piped to log file.
    
    
    
    #=======================================================================================================
    # NOTE: "exit" inside a bash function when the function is being piped
    # 1) does NOT EXIT the entire script and
    # 2) does NOT MAKE THE EXIT CODE AVAILABLE the usual way.
    # Instead, the script does not quit, and the exit code ends up in builtin variable $PIPESTATUS (array).
    #=======================================================================================================
    exit $matlab_error_code
}



if [[ -z "$log_dir" ]]
then
    # CASE: No log file specified.
    launch_matlab_code "$@"                      | grep "$stdout_pattern" | sed "s/${stdout_pattern}//"
    matlab_error_code=${PIPESTATUS[0]}    # Get exit code of the first component of piping.
else
    # CASE: Log file specified.
    
    # tee options:
    #        -a, --append
    #               append to the given FILEs, do not overwrite
    launch_matlab_code "$@" | tee -a "$log_file" | grep "$stdout_pattern" | sed "s/${stdout_pattern}//"
    matlab_error_code=${PIPESTATUS[0]}    # Get exit code of the first component of piping.
fi



exit $matlab_error_code
