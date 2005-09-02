#!/bin/sh
#
# Script for data delivery
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$

in_dir=/usr/ttt/q
if ! [ -d $log_dir ]
then
	echo IN_DIR : $in_dir does not exist!
	exit 1
fi

log_dir=/data/caa/log
if ! [ -d $log_dir ]
then
	echo creating log_dir: $log_dir
	mkdir $log_dir
fi

out_dir=/data/caa/delivered
if ! [ -d $log_dir ]
then
	  echo creating out_dir: $out_dir
		mkdir $log_dir
fi

cd $in_dir
#files=`ls *.cef.gz`
files=`find . -depth 1 -name \*.cef.gz|sed -e 's=^[.]/=='`
if [ "X$files" = "X" ]
then
	echo no CEF.GZ files in $in_dir
	exit 1
fi

for fname in $files
do
	fsize=`du -hs $fname`
	echo -n "$fsize: delivering..."
	scp $fname efw@caa1.estec.esa.int:/c/data-20/EFW/
	echo "$fname  `LANG= date`">>${log_dir}/deliver.log
	echo -n " moving..."
	mv $fname $out_dir
	echo " Done."
done
