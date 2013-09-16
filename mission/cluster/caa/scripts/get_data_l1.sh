#!/bin/sh
#
# Script for running c_get_batch_l1 using job_file as input.
#
# Usage:_l1.sh  [options] YYYY MM DD [NDAYS]"
#		-d | --data                 Process data [default]
#		-nd | --no-data             Do not process data
#		-sth | --start-hour HOUR    Process only one interval starting at HOUR
#		-j | --job-id ID            Give Job ID, otherwise computed from date
#		-sp | --splot               Make summary plots
#		-fs | --full-scale          Plot full scale, not only 0..12.5 Hz
#		-cpdf | --com-pdf | --common-pdf 
#                               Create a common PDF for the whole job
#		-de | --disp-err | --display-errors
#                               Display error messages on the screen				
#		-dd | --disp-date | --display-date
#                               Display date/time for timing the production	
#		-h | --help | -help | '-?'  Display usage
#		-v | --version              Display version
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$

error()
{
	echo "$@" 1>&2
	usage_and_exit 1
}

usage()
{
	echo "Usage: $PROGRAM [-h] [-v] [-d|--data] [-nd|--no-data] [-sp|--splot] [-sth HOUR] YYYY MM DD [NDAYS]"
}

usage_and_exit()
{
	usage
	exit $1
}

version()
{
	echo "$PROGRAM version $VERSION \$Revision$"
}

get_one_int()
{
	cdir=${YYYY}${MM}${DAY}_${HOURS}00
	donef=$out_dir/$cdir/.done_get_data_l1_data
	rm -f $donef

	echo "irf_log('log_out','$log_dir/$start_time.log');\
 	caa_get_batch_l1('$start_time',$dt,'$out_dir');\
	[s,w] = unix('touch $donef');\
 	exit" | $MATLAB " -nodisplay" >> $log_dir/$start_time-get_data.log 2>&1	

	#Empty dir means NO DATA
	if [ -d $out_dir/$cdir ]
	then
		if ! [ -f $donef ]; then 
			if [ "X$disperr" = "Xyes" ]
			then
				printf '\n-----------ERROR------------\n\n'
				tail -22 ${log_dir}/${start_time}-get_data.log
				printf '\n------------END-------------\n\n... '
			else
				echo -n " Error"; 
			fi
		fi
	else
		echo -n " No data"
	fi
}

do_one_splot()
{
	cdir=${YYYY}${MM}${DAY}_${HOURS}00
	if [ -d $out_dir/$cdir ]
	then
		donef=$out_dir/$cdir/.done_get_data_l1_splot
		rm -f $donef
		
		xtraops="'saveps'"
		[ "X$fullscale" = "Xyes" ] && xtraops="'save','fullscale'"

		echo "irf_log('log_out','$log_dir/$start_time-splot.log');\
 		caa_pl_summary_l1('$start_time',$dt,'$out_dir/$cdir',${xtraops});\
		[s,w] = unix('touch $donef');\
 		exit" | $MATLAB ' -nosplash' >> $log_dir/$start_time-get_data.log 2>&1	

		if ! [ -f $donef ]; then 
			if [ "X$disperr" = "Xyes" ]
			then
				printf '\n-----------ERROR------------\n\n'
				tail -22 ${log_dir}/${start_time}-get_data.log
				printf '\n------------END-------------\n\n... '
			else
				echo -n " Error"; 
			fi
		fi
	else
		echo -n " No data"
	fi
}

PROGRAM=`basename $0`
VERSION=1.0
MATLABSETUP='TMP=/tmp LD_LIBRARY_PATH=$IS_MAT_LIB:$LD_LIBRARY_PATH'
MATLAB='/usr/local/matlab/bin/matlab -c 1712@flexlmtmw1.uu.se:1712@flexlmtmw2.uu.se:1712@flexlmtmw3.uu.se -nojvm'
#MATLAB=/bin/cat
INT_HOURS=3
MAXHOURS=24
HOURS=
data=yes
splot=no
fullscale=no
cpdf=no
disperr=no;
dispdate=no;
while [ $# -gt 0 ]
do
	case $1 in
		-d | --data )
		data=yes
		;;
		-j | --job-id )
		JOBNAME="$2"
		shift
		;;
		-nd | --no-data )
		data=no
		;;
		-h | --help | -help | '-?' )
		usage_and_exit 0
		;;
		-sp | --splot )
		splot=yes
		;;
		-fs | --full-scale )
		fullscale=yes
		;;
		-sth | --start-hour )
		HOURS="$2"
		shift
		;;
		-cpdf | --com-pdf | --common-pdf )
		cpdf=yes
		;;
		-de | --disp-err | --display-errors )
		disperr=yes;
		;;
		-dd | --disp-date | --display-date )
		dispdate=yes;
		;;
		-v | --version )
		version
		exit 0
		;;
		-*)
		error "Unrecognized option $1"
		;;
		*)
		break
		;;
	esac
	shift
done

if [ $# -lt 3 ]; then usage_and_exit 1; fi

if [ "X$HOURS" = "X" ]
then
	HOURS=0
else
	MAXHOURS=$(($HOURS+$INT_HOURS))
	if [ $MAXHOURS -gt 24 ]; then MAXHOURS=24; fi
fi

if [ -z $4 ]
then
	NDAYS=1
else
	NDAYS=$4
fi

YYYY=$1
MM=$2
DD=$3

JOBNAME=${JOBNAME:-"L1-${YYYY}${MM}${DD}"}

out_dir=/data/caa/raw/$JOBNAME
log_dir=/data/caa/log-raw/$JOBNAME

if [ "X`hostname -s`" = "Xamanda" ]
then
	out_dir="/export${out_dir}"
	log_dir="/export${log_dir}"
fi

[ "X$dispdate" = "Xyes" ] && echo "Starting at `date`"
echo "Job ID      $JOBNAME" 
echo "Get data    $data"
echo "Summ plot   $splot" 
[ "X$splot" = "Xyes" ] && echo "Full scale  $fullscale" 

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

export $MATLABSETUP
DAYS=0
while [ $DAYS -lt $NDAYS ]
do
	DAY=$(($DD+$DAYS))
	if [ $DAY -lt 10 ] && [ ${#DAY} -lt 2 ]; then DAY="0$DAY"; fi

	echo -n Processing ${YYYY}-${MM}-${DAY}... 

	while [ $HOURS -lt $MAXHOURS ]
	do
		if [ $HOURS -lt 10 ] && [ ${#HOURS} -lt 2 ]; then HOURS="0${HOURS}";	fi
		echo -n " $HOURS"

		start_time=${YYYY}-${MM}-${DAY}T${HOURS}:00:00.000Z
		dt="$INT_HOURS*60*60"

		if [ "X$data" = "Xyes" ]; then get_one_int; fi
		if [ "X$splot" = "Xyes" ]; then (cd $out_dir; do_one_splot); fi

		HOURS=$(($HOURS+$INT_HOURS))
	done
	echo " Done."

	if [ "X$splot" = "Xyes" ]
	then
		echo -n Joining PDFs...
		xtraops=
		[ "X$fullscale" = "Xyes" ] && xtraops="FULL"
		fmask=EFW_SPLOT_L1${xtraops}__${YYYY}${MM}${DAY}_*00.pdf
		if ! [ -z "`find $out_dir -name $fmask`" ]
		then
			(cd $out_dir; pdfjoin --outfile EFW_SP_COMM_L1${xtraops}__${YYYY}${MM}${DAY}.pdf $fmask |grep Finished)
		else
			echo No files.
		fi
	fi
  
	DAYS=$(($DAYS+1))
	HOURS=0
done

if [ "X$cpdf" = "Xyes" ]
then
	echo -n Joining PDFs...
	xtraops=
	[ "X$fullscale" = "Xyes" ] && xtraops="FULL"
	fmask=EFW_SPLOT_L1${xtraops}__${YYYY}${MM}*_*00.pdf
	(cd $out_dir; pdfjoin --outfile EFW_SP_COMM_L1${xtraops}__${YYYY}${MM}${DD}_${YYYY}${MM}${DAY}.pdf $fmask |grep Finished) 
fi

echo Done with job $JOBNAME
[ "X$dispdate" = "Xyes" ] && echo "Finished at `date`"
