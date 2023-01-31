#!/bin/sh
#
# Script to check checksums for delivered files
#
# usage : chksum.sh LOG_FILE
#
# ----------------------------------------------------------------------------
# SPDX-License-Identifier: Beerware
# "THE BEER-WARE LICENSE" (Revision 42):
# <yuri@irfu.se> wrote this file.  As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
# ----------------------------------------------------------------------------


[ "x$1" = "x" ] && echo "chksum.sh LOG_FILE" && exit 1


[ ! -f $1 ] && echo "Cannot find $1" && exit 1

logf=`cat $1|awk '{print $1 " " $2}'`
pro=

for s in $logf
do
	if [ "x$fn" = "x" ]
	then
		fn=$s
	else
		if test -f $fn
		then
			md=`md5sum $fn|awk '{print $1}'`
			if [ "x$s" != "x$md" ]
			then
				echo "Bad sum : $fn"
				pro='yes'
			fi
		else
			echo "No file : $fn"
			pro='yes'
		fi
		fn=
	fi
done

[ "x$pro" = "x" ] && echo "All is fine"
