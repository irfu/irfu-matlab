#!/bin/sh
#
# Script to produce gap files for C1..C4

NOW=`date -u +20%y-%m-%dT%H:%M:%SZ`

CLI=1
while test $CLI -le 4
do
	case $CLI in
		1)
			string='2001-12-29T00:00:00Z_XXX_L1_P1_L1_P12 2009-10-15T00:00:00Z_XXX_L1_P4_L1_P34 XXX_2003-03-26T00:00:00Z_L1_P32'
			;;
		2)
			string='2007-05-14T00:00:00Z_XXX_1_L1_P12 XXX_2007-11-23T00:00:00Z_L1_P32'
			;;
		3)
			string='2014-11-04T00:00:00Z_XXX_L2_E_L3_E_L3_SFIT_L3_DER 2014-11-04T00:00:00Z_XXX_L1_P2 2002-07-30T00:00:00Z_XXX_L1_P1_L1_P12 2011-06-01T00:00:00Z_XXX_L1_P3_L1_P32_L1_P34 XXX_2003-03-26T00:00:00Z_L1_P32'
			;;
		4)
			string='2013-07-02T00:00:00Z_XXX_L1_P4_L1_P34'
			;;
		*)
			echo "should not be here"
			;;
	esac
	START=""
	STOP=""
	for word in $string; do
		(IFS='_'; for subword in $word; do 
			if   [ "$START" == "" ]; then START=$subword; continue;
			elif [ "$STOP" == "" ]; then STOP=$subword; continue;
			else
				PROBE=$subword
				if [ "$PROBE" == "L1" ] || [ "$PROBE" == "L2" ] || [ "$PROBE" == "L3" ]; then PREF=$PROBE; continue; fi
				PROBE=${PREF}_${PROBE}
				
				if [ "$START" == "XXX" ]; then START="2000-09-01T00:00:00Z"; fi
				if [ "$STOP" == "XXX" ]; then STOP="2020-01-01T00:00:00Z"; fi
				STARTSHORT=`echo $START | cut -c 1-4,6,7,9,10`
				FILE_TIME_SPAN=$START/$STOP
				fout=C${CLI}_CP_EFW_${PROBE}__${STARTSHORT}_V00
				
                echo "$fout.cef : $FILE_TIME_SPAN"
                cat gap_skeleton.cef| \
                sed -e "s=__GENERATION_DATE__=${NOW}="| \
                sed -e "s=ZZZ=${STARTSHORT}="| \
                sed -e "s=XXX=${CLI}="| \
                sed -e "s=YYY=${PROBE}="| \
                sed -e "s=YYY=${PROBE}="| \
                sed -e "s=__LOGICAL_FILE_ID__=${fout}="| \
                sed -e "s=__FILE_TIME_SPAN__=${FILE_TIME_SPAN}=" > $fout.cef
			fi
		done)
	done
	CLI=$(($CLI+1))
done
