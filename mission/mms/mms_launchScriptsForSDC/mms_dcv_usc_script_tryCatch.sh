#!/bin/bash

# This is used for starting the MATLAB process of irfu-matlab code on SDC for the MMS mission.
# Author: T. Nilsson, IRFU
# Date: 2014/02/28
#
# Usage: place script in the same folder as has irfu-matlab as a subfolder, then run
#  "./script.sh <mmsX_dcv_filename> <mmsX_101_filename>", with the following
#  input arguments:
#    <mmsX_***_dcv_filename.cdf> = Filename of DC V data to be processed for Usc. Including path and extention.
#    <mmsX_***_101_filename.cdf> = Filename of HK 101 data (with sunpulse) to be processed for Usc. Including path and extention.
#  output files created:
#    <mmsX_***_l2_dcv_yyyymmddHHMMSS_vX.Y.Z.cdf>         = File placed in $DROPBOX_ROOT
#    <DATE_IRFU.log>                                     = Logfile of run, placed in $LOG_PATH_ROOT/mmsX/sdp/.
#    <mmsX_***_l2_dcv_yyyymmddHHMMSS_vX.Y.Z_runTime.log> = File to identify output file, as per e-mail of 2013/11/22. Also placed in $LOG_PATH_ROOT/mmsX/sdp/.
#    return code 0, if everything is ok
#    return code 100-198 if MMS Matlab specific code error. <<-- FIXME: These are still to be defined and specified.
#
# Note: The script assumes it is located in the folder which has irfu-matlab as a subfolder.

# make sure that the correct number of arguments are provided
if [ ${#} -lt 2 ] || [ ${#} -gt 2 ] ; then
  exit 166  # SDC-defined error code for "incorrect usage"
fi

# For debug: display input arguments to terminal
#echo $1 # First argument should be DCV file

# SET ENVIRONMENT MATLABPATH by finding all dir and subdir of 'irfu-matlab' with full path.
#find `pwd` -type d \( -name '@*' -o -name '+*' -o -name '.git' \) -prune -o \( -path "*irfu-matlab*" -type d \) -printf %p:
export MATLABPATH="$(find `pwd` -type d \( -name '@*' -o -name '+*' -o -name '.git' \) -prune -o \( -path "*irfu-matlab*" -type d \) -printf %p:)$MATLABPATH"

# RUN THIS ONLY IF DCE and Sunpulse file inputs are given
if [ ${#} -eq 2 ] ;  then
   matlab -nodesktop -nosplash -nodisplay -nojvm -r "try, mms_sdc_sdp_proc('usc','$1','$2'), catch err, if(strcmp(err.identifier,'MATLAB:SDCcode')) irf.log('critical',['Bash error catch worked: ', err.identifier, '. With message: ', err.message]); exit(str2num(err.message)); else irf.log('critical',['Bash error catch: ', err.identifier, '. With message: ', err.message]); exit(199); end; end, exit(0)"
fi

# For debug: display message 'all done'.
#echo $? # Should return 0 if ok, or 116 if input arguments wrong, 199 if error occured calling built in Matlab function, 100-198 if Matlab MMS specific code had errors...
#echo 'Back in bash, all done.'
