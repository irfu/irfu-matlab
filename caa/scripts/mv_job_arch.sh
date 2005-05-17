#!/bin/sh
#
# Script for cleaning the data directories and saving to archive.
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$


if [ "X$1" = "X" ]
then
	echo Usage: mv_job_arch.sh job-id
	exit 1
fi

if ! [ -d $1 ]
then
	echo directory $1 not found
	exit 1
fi

out_dir=/data/caa/raw-arch/$1

echo Starting job $1 
if ! [ -d $out_dir ]
then
	echo creating out_dir: $out_dir
	mkdir $out_dir
fi
events=`(cd $1;find . -depth 1 -type d -name 200\*_\*)`

if [ "X$events" = "X" ]
then
    echo no events found
    exit 1
fi

for event in $events
do
	echo Processing $event
	if ! [ -d $out_dir/$event ]
	then
		echo creating out_dir: $out_dir/$event
	    mkdir $out_dir/$event
	else
		mv $out_dir/$event $out_dir/$event.bak
		mkdir $out_dir/$event
	fi
	(cd $1/$event;rm -f *.ps *.pdf *.png mBr.mat mEDSIf.mat mEdB.mat tB_*.0*)
	mv $1/$event/*.mat $1/$event/.version $out_dir/$event
	rm -rf $1/$event
	if [ -d $out_dir/$event.bak ]
	then
		rm -rf $out_dir/$event.bak
	fi
done
rmdir $1
echo done with job $1
