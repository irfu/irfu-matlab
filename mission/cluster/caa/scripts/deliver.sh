#!/bin/sh
#
# Script for data delivery to the CAA
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
echo "Usage: $PROGRAM [-h] [-v] [-i|--in-dir IN_DIR] [-o|--out-dir OUT_DIR] [-l|--log-dir LOG_DIR] [-e|--log-ext LOG_EXT]"
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

PROGRAM=`basename $0`
VERSION=1.0
DELAY=3
in_dir=/usr/caa/q
out_dir=/data/caa/delivered
log_dir=/data/caa/log
log_ext=


while [ $# -gt 0 ]
do
	case $1 in
		-i | --in-dir )
		in_dir="$2"
		shift
		;;
		-o | --out-dir )
		out_dir="$2"
		shift
		;;
		-l | --log-dir )
		log_dir="$2"
		shift
		;;
		-e | --log-ext )
		log_ext="$2"
		shift
		;;
		-h | --help | -help | '-?' )
		usage_and_exit 0
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

if [ ! -d $in_dir ]; then error "IN_DIR : $in_dir does not exist!"; fi
if [ ! -d $out_dir ]
then
	  echo creating out_dir: $out_dir
		mkdir -p $out_dir
		if [ ! -d $out_dir ]; then error "OUT_DIR : $out_dir does not exist!"; fi
fi
if [ ! -w $out_dir ]; then error "Cannot write to OUT_DIR: ${out_dir}!" ; fi
if [ ! -d $log_dir ]
then
	echo creating log_dir: $log_dir
	mkdir -p $log_dir
	if [ ! -d $log_dir ]; then error "LOG_DIR : $log_dir does not exist!"; fi
fi

cd $in_dir
files=`find . -depth 1 -name \*.cef.gz|sed -e 's=^[.]/=='|sort`
if [ "X$files" = "X" ]
then
	echo no CEF.GZ files in $in_dir
	exit 1
fi

if [ "X$log_ext" = "X" ]
then 
	log_file=${log_dir}/deliver.log
else
	log_file=${log_dir}/deliver-${log_ext}.log
fi
if [ ! -e ${log_file} ]; then touch ${log_file}; fi
if [ ! -w ${log_file} ]; then error "Cannot write to LOG_FILE: ${log_file}!" ; fi

count=0
for fname in $files
do
	#Every 10 steps we pause for DELAY sec in case 
	#we need to abort the operation
	count=$(($count+1))
	if [ $count = 10 ]
	then
		echo -n "Waiting $DELAY sec for ctrl-c... "
		sleep $DELAY
		echo "Done."
		count=0
	fi

	fsize=`du -hs $fname`
	echo -n "$fsize: delivering..."
	scp $fname efw@caa1.estec.esa.int:/c/data-20/EFW/
	if ! [ $? = 0 ]; then echo "Error!"; exit 1; fi
	md=`md5sum $fname| awk '{print $1}'`
	echo "$fname  $md `LANG= date`">>${log_file}
	echo -n " moving..."
	mv $fname $out_dir
	if ! [ $? = 0 ]; then echo "Error!"; exit 1; fi
	echo " Done."
done
