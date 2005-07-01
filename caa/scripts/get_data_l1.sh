#!/bin/sh
#
# Script for running c_get_batch_l1 using job_file as input.
# Job file has two columns: start_time [ISO time] and dt [sec]
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
	echo "Usage: $PROGRAM [-h] [-v] [-d|--data] [-nd|--no-data] [-sp|--splot] YYYY MM DD [NDAYS]"
}

usage_and_exit()
{
	usage
	exit $1
}

version()
{
	echo "$PROGRAM version $VERSION"
}

get_one_int()
{
	echo Processing $start_time dt=$dt sec ...
	echo "irf_log('log_out','$log_dir/$start_time.log'); caa_get_batch_l1('$start_time',$dt,'$out_dir'); exit" | $MATLAB " -nodisplay" >> $log_dir/$start_time-get_data.log 2>&1	
}
do_one_splot()
{
	echo Summary plot $start_time dt=$dt sec ...
	echo "irf_log('log_out','$log_dir/$start_time-splot.log'); caa_pl_summary_l1('$start_time',$dt,'$out_dir/$@','save'); exit" | $MATLAB ' -nosplash' >> $log_dir/$start_time-get_data.log 2>&1	
}

PROGRAM=`basename $0`
VERSION=1.0
MATLABSETUP='TMP=/tmp LD_LIBRARY_PATH=$IS_MAT_LIB:$LD_LIBRARY_PATH'
MATLAB='/usr/local/matlab/bin/matlab -c 1712@flexlmtmw1.uu.se:1712@flexlmtmw2.uu.se:1712@flexlmtmw3.uu.se -nojvm'
#MATLAB=/bin/cat
INT_HOURS=3
data=yes
splot=no
while test $# -gt 0
do
	case $1 in
		-d | --data )
		data=yes
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

if test $# -lt 3
then
	usage_and_exit 1
fi

if test -z $4
then
	NDAYS=1
else
	NDAYS=$4
fi

YYYY=$1
MM=$2
DD=$3

JOBNAME="L1-$YYYY$MM$DD"

out_dir=/data/caa/raw/$JOBNAME
log_dir=/data/caa/log-raw/$JOBNAME

echo Starting job $JOBNAME 
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
HOURS=0
while test $DAYS -lt $NDAYS
do
	DAY=$(($DD+$DAYS))
	if test $DAY -lt 10
	then
		DAY="0$DAY"
	fi
	while test $HOURS -lt 24
	do
		if test $HOURS -lt 10
		then
			HOURS="0$HOURS"
		fi
		start_time=$YYYY-$MM-$DAY'T'$HOURS':00:00.000Z'
		dt="$INT_HOURS*60*60"
		
		if test "$data" = "yes"
		then
			get_one_int
		fi
		if test "$splot" = "yes"
		then
			cdir=$YYYY$MM$DAY'_'$HOURS'00'
			if test -d $out_dir/$cdir
			then
				(cd $out_dir; do_one_splot $cdir)
			else
				echo Skipping $start_time
			fi
		fi

		HOURS=$(($HOURS+$INT_HOURS))
	done
	DAYS=$(($DAYS+1))
done

if test "$splot" = "yes"
then
	echo Joining PDFs...
	(cd $out_dir; pdfjoin --outfile EFW_SP_COMM_L1__$YYYY$MM$DD.pdf EFW_SPLOT_L1__*.pdf )
fi
	
echo done with job $JOBNAME
