#!/bin/sh
#
# Download data from the CAA
#
# $Id$
#
# ----------------------------------------------------------------------------
# SPDX-License-Identifier: Beerware
# "THE BEER-WARE LICENSE" (Revision 42):
# <yuri@irfu.se> wrote this file.  As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
# ----------------------------------------------------------------------------


usage()
{
	echo "Usage: $PROGRAM DATASETS ISO_START ISO_STOP"
}

usage_and_exit()
{
	usage
	exit $1
}

error()
{
	echo "$@" 1>&2
	usage_and_exit 1
}

cleanup()
{
	echo -n "Removing temporary data... "
	rm -rf $TMPFILE
	if [ $? -ne 0 ]; then
		echo "$0: error"
	fi
	echo Done.
	exit $1
}

PROGRAM=`basename $0`

[ "x$1" = "x" ] && usage_and_exit 1 
[ "x$2" = "x" ] && usage_and_exit 1 
[ "x$3" = "x" ] && usage_and_exit 1 

TMPFILE=`mktemp -q -d /tmp/caa_get.XXXXXX`
if [ $? -ne 0 ]; then
	error "Can't create temp file, exiting..."
fi

echo ""
echo "Downloading DATA ..."
echo ""

(cd $TMPFILE && wget "http://caa.estec.esa.int/caa_query/?uname=Khotyaintsev&pwd=xxx&dataset_id=$1&time_range=$2/$3&format=cdf")

echo ""
if [ $? -ne 0 ]; then
	echo "$0: error"
	exit 1
fi

ERRF=`find $TMPFILE -name CAA_Error.log\*`
if test $ERRF; then
	echo "ERROR :"
	echo ""
	cat $ERRF
	echo ""
	cleanup 1
fi

INFF=`find $TMPFILE -name CAA_Info.log\*`
if test $INFF; then
	echo "Download scheduled:"
	echo ""
	cat $INFF
	echo ""
	cleanup 0
fi


FILES=`find $TMPFILE -name \*.zip`
if test $FILES; then
	echo "Unzipping downloaded data : "
	echo ""
	unzip $FILES
	if [ $? -ne 0 ]; then
		echo "$0: error"
		exit 1
	fi
	echo ""
	cleanup 0
fi

echo "Nothing was downloaded :("
echo ""
cleanup 0
