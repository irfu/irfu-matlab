#!/bin/bash

# This is used for starting the MATLAB process of irfu-matlab code on SDC for the MMS mission.
# Author: T. Nilsson, IRFU
# Date: 2014/03/04
# Updated: 2015/01/14, change comments to reflect change SDP->EDP and add new input argument HK 10E file.
# Updated: 2015/03/20, automatically send mail to mms-ops@irfu.se if error occurs in Matlab.
# Updated: 2015/11/18, add optional ASPOC L2 srvy as input, remove "(" or ")" from errStr (ie anonumous function error ensure bash can still send mail).
# Updated: 2015/12/03, new return codes 198 if I/O error reading zlib compressed cdf file (mainly aspoc), and 197 if error reading DEFATT ascii file.
# Updated: 2015/12/16, new return code 196 if L1b dce file in empty ("Epoch" variable contains zero Written records).
# Updated: 2017/01/09, reading XML files (containing information about manuevers) require jvm.
#
# Usage: place script in the same folder as has irfu-matlab as a subfolder, then run
#  "./script.sh <mmsX_dce_filename> <mmsX_dcv_filename> <mmsX_101_filename> <mmsX_10e_filename> <mmsX_105_filename>", with the following
#  input arguments: (order is irrelevant as long as the OptionalDataDescriptor is "_dce_", "_dcv_","_10e_", "_105_" or "_101_").
#    <mmsX_***_dce_filename.cdf> = Filename of DC E data to be processed for 'xyz'. Including path and extension.
#    <mmsX_***_dcv_filename.cdf> = Filename of DC V data to be processed for 'xyz'. Including path and extention.
#    <mmsX_***_105_filename.cdf> = Filename of HK 105 data to be processed for 'xyz'. Including path and extention.
#    <mmsX_***_10e_filename.cdf> = Filename of HK 10E data (with guard settings) to be processed for 'xyz'. Including path and extention.
#    <mmsX_aspoc_l2_srvy_***_filename.cdf> = Filename of ASPOC L2 srvy data (with aspoc status) to be processed for 'xyz'. Including path and extention. (optional, if not included script will go looking for it).
#    <mmsX_DEFATT_***> = Filename of DEFATT data (with phase) to be processed for 'xyz'. Including path and extention. (optional, if not included script will go looking for it)
#    <mmsX_***_101_filename.cdf> = Filename of HK 101 data (with sunpulse) to be processed for 'xyz'. Including path and extention. (if no DEFATT exist, this will be used otherwise not required)
#    <mmsX_***_l2a_fast_yyyymmddHHMMSS_vX.Y.Z.cdf> = IF processing Brst mode, include the corresponding L2A DCE2d Fast mode file created previously.
# if using multiple HK input files, separate these by a colon (:) without additional spaces,
# ie. two HK 101 data files would be the following
# <mmsX_***_101_filename.cdf>:<mmsX_***_101_filename2.cdf>
#  output files created:
#    <mmsX_***_dce2d_yyyymmddHHMMSS_vX.Y.Z.cdf>         = File placed in $DROPBOX_ROOT
#    <DATE_IRFU.log>                                  = Logfile of run, placed in $LOG_PATH_ROOT/mmsX/edp/.
#    <mmsX_***_dce2d_yyyymmddHHMMSS_vX.Y.Z_runTime.log> = File to identify output file, as per e-mail of 2013/11/22. Also placed in $LOG_PATH_ROOT/mmsX/edp/.
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
if [ ${#} -lt 3 ] || [ ${#} -gt 8 ] ; then
	echo "ERROR: Wrong number of input parameters: min: 3, max: 8"
	exit 166  # SDC-defined error code for "incorrect usage"
fi

# test that Matlab binary (startup script) is executable
if [ ! -x $MATLAB_EXE ] ; then 
	echo "ERROR: Matlab [$MATLAB_EXE] not found/not executable"
	exit 166  # SDC-defined error code for "incorrect usage"
fi

# Verify awk and cdfdump exists as commands in this environment.
command -v awk >/dev/null 2>&1 && awkExist=true || { awkExist=false; }
if [ -x $CDF_BASE/bin/cdfdump ]; then cdfdumpExist=true; else cdfdumpExist=false; fi

# If awk & cdfdump exist, try to identify L1b dce file and check number of records written. If zero records in "Epoch" variable, return exit code 196.
if [ "$awkExist"==true ] && [ "$cdfdumpExist"==true ]; then
   for var in "$@"
    do
      if [ -f $var ]; then # its a single file, not one of the combined files separated by ":".
         if [[ $var =~ mms[1-4]_edp_(fast|slow|brst|comm)_l1b_dce[0-9]{0,3}_20[0-9]{6,12}_v[0-9]+.[0-9]+.[0-9]+.cdf ]]; then
            # File matches the MMS L1b dce file naming convention.
            # Dump one record of variable "Epoch" in the file along with some metadata
	    # awk for line with "Written records  12345/12345(max)" and print only the 
	    # number (before the "/" sign).
            nrec=`$CDF_BASE/bin/cdfdump -dump data -recordrange "1,1" -vars "Epoch" $var | awk -F"[ /]+" '/Written/ {print $3}'`
            if [ "$nrec" == 0 ]; then
               exit 196
            fi
         fi # if l1b dce
      fi # if -f var
   done # for var
fi # if awk & cdfdump exist

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

