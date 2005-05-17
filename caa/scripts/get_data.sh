#!/bin/sh
#
# Script for running c_get_batch using job_file as input.
# Job file has two columns: start_time [ISO time] and dt [sec]
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$


if [ "X$1" = "X" ]
then
	echo Usage: get_data.sh job-id
	exit 1
fi

if ! [ -e "$1.dat" ]
then
	echo file $1.dat not found
	exit 1
fi

matlab_setup='TMP=/tmp LD_LIBRARY_PATH=$IS_MAT_LIB:$LD_LIBRARY_PATH'
matlab_cmd='/usr/local/matlab/bin/matlab -c 1712@flexlmtmw1.uu.se:1712@flexlmtmw2.uu.se:1712@flexlmtmw3.uu.se -nojvm -nodisplay'
#matlab_cmd=/bin/cat
out_dir=/data/caa/raw/$1
log_dir=/data/caa/log-raw/$1

echo Starting job $1 
if ! [ -d $out_dir ]
then
	echo creating out_dir: $out_dir
	mkdir $out_dir
fi
if ! [ -d $log_dir ]
then
	echo creating log_dir: $log_dir
	mkdir $log_dir
fi
jobs_def=`cat $1.dat`

count=0
for job in $jobs_def
do
	count=$(($count+1))
	if [ $count = 1 ]
	then
		start_time=$job
		dt=''
	else
		dt=$job
		echo Processing $start_time dt=$dt sec ...
		#echo $matlab_setup
		#echo $matlab_cmd
		export $matlab_setup
		echo "irf_log('log_out','$log_dir/$start_time.log'); caa_get_batch('$start_time',$dt,'$out_dir'); exit" | $matlab_cmd > $log_dir/$start_time-get_data.log 2>&1	
	fi

	if [ $count = 2 ]
	then
		count=0
	fi

done

echo done with job $1
