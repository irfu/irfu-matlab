#!/bin/sh
#
# Script for running c_get_batch_l1 using job_file as input.
#
# Usage:mv_problem_l1.sh  [options] -j JOB-ID YYYY MM DD HH [HH ...]"
#		-d | --data                 Process data [default]
#		-j | --job-id ID            Give Job ID, otherwise computed from date
#		-c | --c-list               Give sc list
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
	echo "Usage: $PROGRAM [-h] [-v] [-j|--job-id] [-c|--c-list] YYYY MM DD HH [HH ...]"
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

mv_one_c()
{
	cdir=${YYYY}${MM}${DD}_${HOURS}00
	if [ ! -d "${problem_dir}/$JOBNAME/$cdir" ]; then mkdir -p ${problem_dir}/$JOBNAME/$cdir; fi
	echo -n "Moving $JOBNAME/$cdir/C$cli to problems ... "
	mv ${data_dir}/$JOBNAME/$cdir/C$cli ${problem_dir}/$JOBNAME/$cdir	
	echo Done.
}
mv_all_c()
{
	cdir=${YYYY}${MM}${DD}_${HOURS}00
	if [ ! -d "${problem_dir}/$JOBNAME" ]; then mkdir -p ${problem_dir}/$JOBNAME; fi
	echo -n "Moving $JOBNAME/$cdir to problems ... "
	mv ${data_dir}/$JOBNAME/$cdir ${problem_dir}/$JOBNAME	
	echo Done.
}

PROGRAM=`basename $0`
VERSION=1.0
data_dir=/data/caa/raw
problem_dir=${data_dir}/problem
clist=1234

while [ $# -gt 0 ]
do
	case $1 in
		-d | --data )
		data=yes
		;;
		-j | --job-id )
		JOBNAME="$2";
		shift
		;;
		-nd | --no-data )
		data=no
		;;
		-h | --help | -help | '-?' )
		usage_and_exit 0
		;;
		-c | --c-list )
		clist="$2";
		shift
		;;
		-v | --version )
		version
		exit 0
		;;
		-*)
		error "Unrecognized option $1";
		;;
		*)
		break
		;;
	esac
	shift
done

if [ $# -lt 4 ]; then usage_and_exit 1; fi

if [ "X$JOBNAME" = "X" ]; then echo "-j JOBNAME is missing"; usage_and_exit 1; fi
if [ ! -d "${data_dir}/$JOBNAME" ]; then error "Bad JOB-ID. ${data_dir}/$JOBNAME does not exist"; fi

YYYY=$1
MM=$2
DD=$3
shift;shift;shift

while [ $# -gt 0 ]
do
	HOURS=`echo $1|tr -d ,`
	shift
	if [ "X$clist" = "X1234" ]
	then
		mv_all_c
	else
		tmp_list=$clist
		while [ "X$tmp_list" != "X" ]
		do
			cli=`echo $tmp_list|cut -c 1`
			mv_one_c
			tmp_list=`echo $tmp_list|tr -d $cli`
		done
	fi
done
