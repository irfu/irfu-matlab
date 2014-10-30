#!/bin/sh
#
# Usage : ls *.cef | maarble_ingest.sh

CEFMERGE=/usr/local/bin/cefmerge
if [ ! -x $CEFMERGE ]; then 
	echo "$CEFMERGE does not exist/not executable" && exit 1; 
fi
QTRAN=/usr/local/bin/Qtran
if [ ! -x $QTRAN ]; then 
	echo "$QTRAN does not exist/not executable" && exit 1;
fi

umask 002

BASEDIR=/data/caa/MAARBLE
DBDIR="$BASEDIR/WaveDatabase"
DELIVERY_DIR="$BASEDIR/Delivery"
if [ ! -d $DELIVERY_DIR/CEF ]; then 
	mkdir -p $DELIVERY_DIR/CEF || exit 1
	echo Created DELIVERY_DIR : $DELIVERY_DIR
fi
LOGDIR="$BASEDIR/Log"
if [ ! -d $LOGDIR ]; then 
	mkdir $LOGDIR || exit 1 
	echo Created LOGDIR : $LOGDIR
fi
FAILED="$BASEDIR/Failed"
if [ ! -d $FAILED ]; then 
	mkdir $FAILED || exit 1
	echo Created FAILED_DIR : $FAILED
fi

TMPDIR=`mktemp -d -t MAARBLE.XXXXXX` || exit 1
INCLUDES=$TMPDIR/include
mkdir $INCLUDES
for inst in irf noa uofa iap; do
	cp $BASEDIR/Upload/$inst/HEADERS/*.ceh $INCLUDES > /dev/null 2>&1 
done

STATUS=''
while read fname; do
	case $STATUS in
		OK|Failed) echo $STATUS;;
		*) ;;
	esac
	STATUS=Failed
	STATUS_CODE=0
	NAME=`echo "$fname" | cut -d'.' -f1`
	echo -n "$NAME  "
	LOG=$LOGDIR/$NAME.log
	if [ ! -e $fname ]; then
		echo "Cannot find $fname" && continue
	fi   
	rm -f $LOG
	$CEFMERGE -I $INCLUDES -O $TMPDIR $fname >> $LOG 2>&1 || (cp $fname $FAILED && continue) 
	newfile=`find $TMPDIR -name \*.cef` 
	echo -n "$QTRAN $newfile ... " >> $LOG
	QTRAN_OUT=`$QTRAN $newfile 2>&1`
	STATUS_CODE=`echo $QTRAN_OUT | grep Done`
	if [ -z "$STATUS_CODE" ]; then
		echo Failed >> $LOG && echo $QTRAN_OUT >> $LOG && cp $fname $FAILED && continue
	else echo $STATUS_CODE >> $LOG
	fi
	
	# Get destination directory for the data
	DATASET_NAME=`echo $NAME|awk -F'__' '{print $1}'`
	case "$DATASET_NAME" in
  		C[1-4]_CP_AUX_MAARBLE_*)
  			SHORT_NAME=`echo $DATASET_NAME|awk -F'CP_AUX_MAARBLE' '{print $2}'`
			PROJ=Cluster
			MEMBER=C`echo $DATASET_NAME|cut -c2`
		;;
		CC_CP_AUX_MAARBLE_*)
  			SHORT_NAME=`echo $DATASET_NAME|awk -F'CC_CP_AUX_MAARBLE_' '{print $2}'`
  			case "$SHORT_NAME" in
				CHAMP*) PROJ=CHAMP;;
  				TH[A-E]_*) PROJ=THEMIS;;
  				G1[1-2]_*) PROJ=GOES;;
  				DOB_*|HOR_*|KEV_*|KIR_*|NUR_*|OUJ_*|RVK_*|SOD_*|UPS_*|TRO_*)
  					PROJ=IMAGE;;
  				FCHU_*|GILL_*|ISLL_*|MCMU_*|PINA_*|RANK_*)
  					PROJ=CARISMA;;
  				*) echo "Unknown non-Cluster project" && exit 1;;
  			esac
  			MEMBER=`echo $SHORT_NAME|cut -d'_' -f1`
  		;;
  		*) echo "Unknown non-Cluster parameter" && exit 1;;
	esac
	case "$SHORT_NAME" in 
		*_ULF_PC1) DB=ULF; DSET=PC1;;
		*_ULF_PC12) DB=ULF; DSET=PC12;;
		*_ULF_PC35) DB=ULF; DSET=PC35;;
		*_ULF_FACMATR) DB=ULF; DSET=FACMATR;;
		*_VLF*) DB=VLF; DSET=VLF;;
		*) echo "Unknown Cluster parameter" && exit 1;;
	esac
	DEST=$DBDIR/$DB/$PROJ/$MEMBER/$DSET
	if [ ! -d $DEST/CEF ]; then
		mkdir -p $DEST/CEF || exit 1
	fi
	if [ ! -d $DEST/CDF ]; then
		mkdir -p $DEST/CDF || exit 1
	fi
	if [ $DSET != "FACMATR" ] && [ ! -d $DEST/PNG ]; then
		mkdir -p $DEST/PNG || exit 1
	fi

	isGzipped=`file $fname | grep gzip`
	if [ -z "$isGzipped" ]; then
		echo compressing $fname >> $LOG
		gzip $fname || exit 1	
		echo moving $fname.gz $DELIVERY_DIR >> $LOG
		mv $fname.gz $DELIVERY_DIR/CEF/ || exit 1
	else
		echo moving $fname $DELIVERY_DIR >> $LOG
		mv $fname $DELIVERY_DIR/CEF/ || exit 1
	fi

	echo compressing $newfile >> $LOG
	gzip $newfile || exit 1	
	echo moving $newfile.gz $DEST/CEF >> $LOG
	mv $newfile.gz $DEST/CEF/ || exit 1 
	 
	newfile=`find $TMPDIR -name \*.cdf`	
	if [ -n "$newfile" ]; then
		echo moving $newfile to $DEST/CDF >> $LOG
		mv $newfile $DEST/CDF/ || exit 1
	else
		continue
	fi
	STATUS=OK
done
case $STATUS in
	OK|Failed) echo $STATUS;;
	*) ;;
esac
rm -rf $INCLUDES 
rmdir $TMPDIR
