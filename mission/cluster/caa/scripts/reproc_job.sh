#!/bin/sh
#
# Script for reprocessing E in the data directories
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$


if [ "X$1" = "X" ]
then
	echo Usage: reproc_job_arch.sh job-id
	exit 1
fi

if ! [ -d $1 ]
then
	echo directory $1 not found
	exit 1
fi

matlab_setup='TMP=/tmp LD_LIBRARY_PATH=$IS_MAT_LIB:$LD_LIBRARY_PATH'
matlab_cmd='/usr/local/matlab/bin/matlab -c 1712@flexlmtmw1.uu.se:1712@flexlmtmw2.uu.se:1712@flexlmtmw3.uu.se -nojvm -nodisplay'

log_dir=/data/caa/log-raw/$1
if ! [ -d $log_dir ]
then
	echo creating log_dir: $log_dir
	mkdir $log_dir
fi


echo Starting job $1 

events=`(cd $1;find . -depth 1 -type d -name 200\*_\*)`

if [ "X$events" = "X" ]
then
    echo no events found
    exit 1
fi

for event in $events
do
	echo Processing $event
	export $matlab_setup
	cluster="1\
	2\
	3\
	4"
	for cli in $cluster 
	do
		if [ -d "$1/$event/C$cli" ]
		then
			ints=`(cd $1/$event/C$cli;find . -depth 1 -type d -name 200\*_\*)`
			for int in $ints
			do
				donef=".done_reproc"
				rm -f "$1/$event/C$cli/$int/$donef"
				
				echo -n Processing "$1/$event/C$cli/$int"

				(cd "$1/$event/C$cli/$int";\
				echo "disp(sprintf('\nLOG %s : %s \n',datestr(now),pwd));\
				irf_log('log_out','$log_dir/reproc.log');\
	   		c_get_batch(0,0,'sc_list',$cli,'varsproc',{'dief'});\
				[s,w] = unix('touch $donef');\
				exit" | $matlab_cmd >> $log_dir/reproc_comm.log 2>&1)

				if ! [ -f "$1/$event/C$cli/$int/$donef" ]; then
					echo " Error!"
					#printf '\n-----------ERROR------------\n\n'
					#tail -22 ${log_dir}/reproc_comm.log
					#printf '\n------------END-------------\n\n... '
				else
					echo " Done."
				fi
			done
		fi
	done
done
echo done with job $1
