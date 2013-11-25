#!/bin/sh
#
# Usage : ls *.cef | maarble_ingest.sh

CEFMERGE=/usr/local/bin/cefmerge
if [ ! -x $CEFMERGE ]; then ("$CEFMERGE does not exist/not executable" && exit 1); fi
QTRAN=/usr/local/bin/Qtran
if [ ! -x $QTRAN ]; then ("$QTRAN does not exist/not executable" && exit 1); fi

BASEDIR=/data/caa/MAARBLE
DELIVERY_DIR="$BASEDIR/Delivery"
if [ ! -d $DELIVERY_DIR ]; then 
	mkdir $DELIVERY_DIR || (echo "Cannot create DELIVERY dir" && exit 1)
	mkdir $DELIVERY_DIR/CEF || (echo "Cannot create CEF DELIVERY dir" && exit 1)
	mkdir $DELIVERY_DIR/CDF || (echo "Cannot create CDF DELIVERY dir" && exit 1)
	echo Created DELIVERY_DIR : $DELIVERY_DIR
fi
LOGDIR="$BASEDIR/Log"
if [ ! -d $LOGDIR ]; then 
	mkdir $LOGDIR || (echo "Cannot create LOG dir" && exit 1)
	echo Created LOGDIR : $LOGDIR
fi
FAILED="$BASEDIR/Failed"
if [ ! -d $FAILED ]; then 
	mkdir $FAILED || (echo "Cannot create FAILED dir" && exit 1)
	echo Created FAILED_DIR : $FAILED
fi

TMPDIR=`mktemp -d -t MAARBLE.XXXXXX` || (echo "Cannot mktemp" && exit 1)
INCLUDES=$TMPDIR/include
mkdir $INCLUDES
for inst in irf noa uofa iap; do
	cp $BASEDIR/Upload/$inst/HEADERS/*.ceh $INCLUDES > /dev/null 2>&1 
done

while read fname; do
	NAME=`echo "$fname" | cut -d'.' -f1`
	echo -n Processing $NAME ...
	LOG=$LOGDIR/$NAME.log
	EXTENSION=`echo "$fname" | cut -d'.' -f2`
	if [ ! -e $fname ]; then
		echo Cannot find $fname
 		continue
	fi   
	$CEFMERGE -I $INCLUDES -O $TMPDIR $fname >> $LOG 2>&1 || (cp $fname $FAILED && continue) 
	newfile=`find $TMPDIR -name \*.cef` 
	$QTRAN $newfile >> $LOG 2>&1
	echo moving $newfile $DELIVERY_DIR/CEF >> $LOG
	mv $newfile $DELIVERY_DIR/CEF 
	newfile=`find $TMPDIR -name \*.cdf`	
	if [ -n "$newfile" ]; then
		echo moving $newfile to $DELIVERY_DIR/CDF >> $LOG
		mv $newfile $DELIVERY_DIR/CDF
	fi
	echo Done.
done
rm -rf $INCLUDES 
rmdir $TMPDIR
