#!/bin/bash

# This is a template: CHANGE xyz to whichever process you are running.
#
# This is used for starting the MATLAB process of irfu-matlab code on SDC for the MMS mission.
# Author: T. Nilsson, IRFU
# Date: 2014/03/04
# Updated: 2015/01/14, change comments to reflect change SDP->EDP and add new input argument HK 10E file.
#
# Usage: place script in the same folder as has irfu-matlab as a subfolder, then run
#  "./script.sh <mmsX_dce_filename> <mmsX_dcv_filename> <mmsX_101_filename> <mmsX_10e_filename>", with the following
#  input arguments: (order is irrelevant as long as the OptionalDataDescriptor is "_dce_", "_dcv_","_10e_" or "_101_").
#    <mmsX_***_dce_filename.cdf> = Filename of DC E data to be processed for 'xyz'. Including path and extension.
#    <mmsX_***_dcv_filename.cdf> = Filename of DC V data to be processed for 'xyz'. Including path and extention.
#    <mmsX_***_101_filename.cdf> = Filename of HK 101 data (with sunpulse) to be processed for 'xyz'. Including path and extention.
#    <mmsX_***_10e_filename.cdf> = Filename of HK 10E data (with guard settings) to be processed for 'xyz'. Including path and extention.
#  output files created:
#    <mmsX_***_xyz_yyyymmddHHMMSS_vX.Y.Z.cdf>         = File placed in $DROPBOX_ROOT
#    <DATE_IRFU.log>                                  = Logfile of run, placed in $LOG_PATH_ROOT/mmsX/edp/.
#    <mmsX_***_xyz_yyyymmddHHMMSS_vX.Y.Z_runTime.log> = File to identify output file, as per e-mail of 2013/11/22. Also placed in $LOG_PATH_ROOT/mmsX/edp/.
#
#  return code 0 if ok.
#  return code 199 if built in Matlab error.
#  return code 100-198 if MMS Matlab specific code error.
#
# Note: The script assumes it is located in the folder which has irfu-matlab as a subfolder.

# User definable constants
MATLAB_EXE=/tools/matlab/R2013b/bin/matlab # SDC location of installed Matlab. # XXX: change this to whereever matlab is located.
MATLAB_FLAGS="-nodesktop -nosplash -nodisplay -nojvm"
IRFU_MATLAB=/mms/itfhome/mms-sdp/software/irfu-matlab # SDC location of irfu-matlab. # XXX: change this to whereever irfu-matlab is located.

# No need to edit after this line
# add IRFU_MATLAB and IRFU_MATLAB/mission/mms to path used by Matlab.
export MATLABPATH=$IRFU_MATLAB:$IRFU_MATLAB/mission/mms

PROCESS_NAME=
case "$0" in
	*mms_dce_ql_script_tryCatch*) PROCESS_NAME=ql ;;
	*mms_dce_sitl_script_tryCatch*) PROCESS_NAME=sitl ;;
	*mms_dcv_usc_script_tryCatch*) PROCESS_NAME=scpot ;;
	*)
	echo "ERROR: urecognized name of the caller routine"
	exit 166
	;;
esac
echo $PROCESS_NAME

# make sure that the correct number of arguments are provided
if [ ${#} -lt 2 ] || [ ${#} -gt 4 ] ; then
	echo "ERROR: Wrong number of input parameters: min: 2, max: 4"
	exit 166  # SDC-defined error code for "incorrect usage"
fi

# test that Matlab binary (startup script) is executable
if [ ! -x $MATLAB_EXE ] ; then 
	echo "ERROR: Matlab [$MATLAB_EXE] not found/not executable"
	exit 166  # SDC-defined error code for "incorrect usage"
fi

# SET ENVIRONMENT LD_LIBRARY_PATH in order for the linking to CDF.h and cdflib.so to properly work.
if [ ! -e $CDF_BASE/lib/libcdf.so ]; then
	echo "ERROR: no libcdf.so in CDF_BASE ($CDF_BASE)"
	exit 166  # SDC-defined error code for "incorrect usage"
fi
if [ "X$LD_LIBRARY_PATH" = "X" ]; then
	export LD_LIBRARY_PATH="$CDF_BASE/lib"
else
	export LD_LIBRARY_PATH="$CDF_BASE/lib:$LD_LIBRARY_PATH"
fi

# RUN THIS IF ONLY ONE FILE EXISTS (DCE)
$MATLAB_EXE $MATLAB_FLAGS -r "try, mms_sdc_sdp_proc('$PROCESS_NAME','$1','$2','$3','$4'), catch err, if(strcmp(err.identifier,'MATLAB:SDCcode')) irf.log('critical',['Bash error catch worked: ', err.identifier, '. With message: ', err.message]); exit(str2num(err.message)); else irf.log('critical',['Bash error catch: ', err.identifier, '. With message: ', err.message]); exit(199); end; end, exit(0)"
