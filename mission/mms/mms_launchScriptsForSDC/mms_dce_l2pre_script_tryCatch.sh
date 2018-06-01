#!/bin/bash

#
# This is used for starting the MATLAB process of irfu-matlab code on SDC for the MMS mission.
# Author: T. Nilsson, IRFU
# Date: 2016/02/11
# Updated: 2016/02/18, added support for processing Brst L1b dce segments directly to L2Pre without going via L2A brst. (Note L2A Fast dce2d is requried for corresponding day).
# Updated: 2017/01/09, reading XML files (containing information about manuevers) require jvm.
#
# Usage: place script in the same folder as has irfu-matlab as a subfolder, then run
#  "./script.sh <mmsX_dce_filename> <mmsX_dfg_l2pre_filename>", with the following
#  input arguments: (order is irrelevant).
#
#    <mmsX_edp_*_l2a_filename.cdf> = Filename of DC E l2a data to be processed for 'L2pre'. Including path and extension.
#    <mmsX_dfg_*_l2pre_filename.cdf> = Filename of DFG srvy L2pre covering the same interval as the "dce" file,
#                                      for Fast/Slow dce l2a this means the dfg srvy l2pre file of same day.
#
# OR if processing BRST segments directly from L1b:
#
#    <mmsX_***_dce_filename.cdf> = Filename of DC E Brst L1b data to be processed for L2Pre. Including path and extension.
#    <mmsX_***_105_filename.cdf> = Filename of HK 105 data to be processed for L2Pre. Including path and extention.
#    <mmsX_***_10e_filename.cdf> = Filename of HK 10E data (with guard settings) to be processed for L2Pre. Including path and extention.
#    <mmsX_dfg_brst_l2pre_filename.cdf> = Filename of corresponding DFG Brst L2Pre segments covering the same time interval as the DCE Brst L1b file. Including path and extention.
#    <mmsX_aspoc_l2_srvy_***_filename.cdf> = Filename of ASPOC L2 srvy data (with aspoc status) to be processed for L2Pre. Including path and extention. (optional, if not included script will go looking for it).
#    <mmsX_DEFATT_***> = Filename of DEFATT data (with phase) to be processed for L2Pre. Including path and extention. (optional, if not included script will go looking for it)
#    <mmsX_***_101_filename.cdf> = Filename of HK 101 data (with sunpulse) to be processed for L2Pre. Including path and extention. (if no DEFATT exist, this will be used otherwise not required)
#    <mmsX_***_l2a_fast_yyyymmddHHMMSS_vX.Y.Z.cdf> = Corresponding L2A DCE2d Fast mode file created previously. (optional, if not included script will go looking for it).
#
# if using multiple HK/Aspoc/Defatt/Dfg input files when processing Brst segments, separate these by a colon (:) without additional spaces,
# ie. two HK 101 data files would be the following
# <mmsX_***_101_filename.cdf>:<mmsX_***_101_filename2.cdf>
#
#  output files created:
#    <mmsX_***_l2pre_yyyymmddHHMMSS_vX.Y.Z.cdf>         = File placed in $DROPBOX_ROOT
#    <DATE_IRFU.log>                                  = Logfile of run, placed in $LOG_PATH_ROOT/mmsX/edp/.
#    <mmsX_***_l2pre_yyyymmddHHMMSS_vX.Y.Z_runTime.log> = File to identify output file, as per e-mail of 2013/11/22. Also placed in $LOG_PATH_ROOT/mmsX/edp/.
#
#  return code 0, if ok.
#  return code 166, if error caused by incorrect usage.
#  return code 196, if error empty l1b dce file. (too early processing?)
#  return code 197, if error I/O DEFATT ascii file.
#  return code 198, if error I/O cdf file (mainly zlib compressed aspoc)
#  return code 199, if error during Matlab process.
#
# Note: The script assumes it is located in the folder which has irfu-matlab as a subfolder.

# User definable constants
MATLAB_EXE=/tools/matlab/R2017a/bin/matlab # SDC location of installed Matlab. # XXX: change this to whereever matlab is located.
MATLAB_FLAGS="-nodesktop -nosplash -nodisplay"
IRFU_MATLAB=/mms/itfhome/mms-sdp/software/irfu-matlab # SDC location of irfu-matlab. # XXX: change this to whereever irfu-matlab is located.

# No need to edit after this line
# add IRFU_MATLAB and IRFU_MATLAB/mission/mms to path used by Matlab.
export MATLABPATH=$IRFU_MATLAB:$IRFU_MATLAB/mission/mms

PROCESS_NAME=
case "$0" in
	*mms_dce_ql_script_tryCatch*) PROCESS_NAME=ql ;;
	*mms_dcv_usc_script_tryCatch*) PROCESS_NAME=scpot ;;
	*mms_dce_l2pre_script_tryCatch*) PROCESS_NAME=l2pre ;;
	*mms_dce_l2a_script_tryCatch*) PROCESS_NAME=l2a ;;
	*)
	echo "ERROR: urecognized name of the caller routine"
	exit 166
	;;
esac
echo $PROCESS_NAME

# make sure that the correct number of arguments are provided
if [ ${#} -lt 2 ] || [ ${#} -gt 8 ] ; then
	echo "ERROR: Wrong number of input parameters: min: 2, max: 8"
	exit 166  # SDC-defined error code for "incorrect usage"
fi

# test that Matlab binary (startup script) is executable
if [ ! -x $MATLAB_EXE ] ; then 
	echo "ERROR: Matlab [$MATLAB_EXE] not found/not executable"
	exit 166  # SDC-defined error code for "incorrect usage"
fi

# RUN Matlab and try to run mms_sdc_sdp_proc, and if any errors are caught
# get file name of log file created, check if that file exist,
# if so attach it (8 Bit ASCII encoded) to a mail sent to mms-ops@irfu.se.
# then exit with 199 (if errors occurred) or with 0 (if no errors).
# exit with 197 if error reading DEFATT files (incorrect times of start/stop-> Epoch error).
# exit with 198 if error reading cdf file (mostly related to ASPOC files),

$MATLAB_EXE $MATLAB_FLAGS -r\
  "\
  try;\
  mms_sdc_sdp_proc('$PROCESS_NAME','$1','$2','$3','$4','$5','$6','$7','$8');\
  catch ME;\
    logFile=irf.log('log_out');\
    errStr=[];\
    for k=1:length(ME.stack),\
      errStr=[errStr,' Error in file: ',ME.stack(k).file,' in function: ',ME.stack(k).name,' at line: ',num2str(ME.stack(k).line),'. '];\
    end;\
    errStr=strrep(errStr,'(','_'); errStr=strrep(errStr,')','_');\
    mess=strrep(ME.message,'(','_'); mess=strrep(mess,')','_');\
    iden=strrep(ME.identifier,'(','_'); iden=strrep(iden,')','_');\
    if(exist(logFile,'file')),\
      unix(['echo ''''Error ',iden,' with message: ',mess,errStr,''''' | mail -a', logFile,' -s MMS_SDC_Error mms-ops@irfu.se']);\
    else,\
      unix(['echo ''''No log found, error: ',iden,' with message: ',mess,errStr,''''' | mail -s MMS_SDC_Error mms-ops@irfu.se']);\
    end;\
    if(~isempty(regexpi(iden,'spdfcdfreadc')) || ~isempty(regexpi(iden,'spdfcdfinfoc'))),exit(198);\
    elseif(regexpi(errStr,'list_ancillary')),exit(197);\
    end;\
    exit(199);\
 end,\
 exit(0)"

