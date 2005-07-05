#!/bin/sh
#
# Script for running c_get_batch_l1 using job_file as input.
#
# Usage:_l1.sh  [options] YYYY MM DD [NDAYS]"
#		-d | --data                 Process data [default]
#		-nd | --no-data             Do not process data
#		-sp | --splot               Make summary plots
#		-cpdf | --com-pdf | --common-pdf 
#                               Create a common PDF for the whole job
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
	echo "Usage: $PROGRAM [-h] [-v] [-d|--data] [-nd|--no-data] [-sp|--splot] YYYY MM DD [NDAYS]"
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
	donef=$out_dir/${YYYY}${MM}${DAY}_${HOURS}00/.done_get_data_l1_data
	rm -f $donef

	echo "irf_log('log_out','$log_dir/$start_time.log');\
 	caa_get_batch_l1('$start_time',$dt,'$out_dir');\
	[s,w] = unix('touch $donef');\
 	exit" | $MATLAB " -nodisplay" >> $log_dir/$start_time-get_data.log 2>&1	

	if ! [ -f $donef ]; then echo -n " Error"; fi
}

do_one_splot()
{
	donef=$out_dir/${YYYY}${MM}${DAY}_${HOURS}00/.done_get_data_l1_splot
	rm -f $donef

	echo "irf_log('log_out','$log_dir/$start_time-splot.log');\
 	caa_pl_summary_l1('$start_time',$dt,'$out_dir/$@','save');\
	[s,w] = unix('touch $donef');\
 	exit" | $MATLAB ' -nosplash' >> $log_dir/$start_time-get_data.log 2>&1	

	if ! [ -f $donef ]; then echo -n " Error"; fi
}

PROGRAM=`basename $0`
VERSION=1.0
MATLABSETUP='TMP=/tmp LD_LIBRARY_PATH=$IS_MAT_LIB:$LD_LIBRARY_PATH'
#MATLAB='/usr/local/matlab/bin/matlab -c 1712@flexlmtmw1.uu.se:1712@flexlmtmw2.uu.se:1712@flexlmtmw3.uu.se -nojvm'
MATLAB=/bin/cat
INT_HOURS=3
data=yes
splot=no
cpdf=no
while [ $# -gt 0 ]
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
		-cpdf | --com-pdf | --common-pdf )
		cpdf=yes
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

if [ -z $4 ]
then
	NDAYS=1
else
	NDAYS=$4
fi

YYYY=$1
MM=$2
DD=$3

JOBNAME="L1-${YYYY}${MM}${DD}"

out_dir=/data/caa/raw/$JOBNAME
log_dir=/data/caa/log-raw/$JOBNAME

echo Starting job $JOBNAME 
echo "Data  = $data"
echo "SPlot = $splot" 

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
while [ $DAYS -lt $NDAYS ]
do
	DAY=$(($DD+$DAYS))
	if [ $DAY -lt 10 ]; then DAY="0$DAY"; fi

	echo -n Processing ${YYYY}-${MM}-${DAY}... 

	while [ $HOURS -lt 24 ]
	do
		if [ $HOURS -lt 10 ]; then HOURS="0${HOURS}";	fi
		echo -n " $HOURS"

		start_time=$YYYY-$MM-$DAY'T'$HOURS':00:00.000Z'
		dt="$INT_HOURS*60*60"

		if [ "X$data" = "Xyes" ]; then get_one_int; fi

		if [ "X$splot" = "Xyes" ]
		then
			cdir=$YYYY$MM$DAY'_'$HOURS'00'
			if [ -d $out_dir/$cdir ]
			then
				(cd $out_dir; do_one_splot $cdir)
			else
				echo -n " No data"
			fi
		fi

		HOURS=$(($HOURS+$INT_HOURS))
	done
	echo " Done."

	if [ "X$splot" = "Xyes" ]
	then
		echo -n Joining PDFs...
		fmask=EFW_SPLOT_L1__${YYYY}${MM}${DAY}_*00.pdf
		if ! [ -z "`find $out_dir -name $fmask`" ]
		then
			(cd $out_dir; pdfjoin --outfile EFW_SP_COMM_L1__${YYYY}${MM}${DAY}.pdf $fmask |grep Finished)
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
	fmask=EFW_SPLOT_L1__${YYYY}${MM}*_*00.pdf
	(cd $out_dir; pdfjoin --outfile EFW_SP_COMM_L1__${YYYY}${MM}${DD}_${YYYY}${MM}${DAY}.pdf $fmask |grep Finished) 
fi

echo Done with job $JOBNAME
