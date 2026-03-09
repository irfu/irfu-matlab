#!/bin/sh
#
# Script to produce headers for C1..C4

NOW=`date -u +20%y-%m-%dT%H:%M:%SZ\ UTC`

f=CX_CH_EFW_INST.ceh
CLI=1
while test $CLI -le 4
do
   fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Writing $fout"
   cat $f|sed -e "s=DDD=${NOW}="|sed -e "s=XXX=${CLI}=" > $fout
   cat C${CLI}_CH_EFW_INST_CAV.ceh >> $fout
   CLI=$(($CLI+1))
done

flist="CX_CH_EFW_L1_P.ceh \
CX_CH_EFW_L1_E.ceh \
CX_CH_EFW_L2_P.ceh \
CX_CH_EFW_L3_P.ceh \
CX_CH_EFW_L2_E.ceh \
CX_CH_EFW_L3_E.ceh \
CX_CH_EFW_L3_SFIT.ceh \
CX_CH_EFW_L3_DER.ceh \
CX_CH_EFW_L2_HK.ceh \
CX_CQ_EFW_INST.ceh \
CX_CH_EFW_L1_IB.ceh \
CX_CH_EFW_L2_PB.ceh \
CX_CH_EFW_L2_EB.ceh \
CX_CH_EFW_L2_BB.ceh"


for f in $flist
do
   CLI=1
   while test $CLI -le 4
   do
      fout=`echo $f| sed -e "s=X=${CLI}="`
      echo "Writing $fout"
      cat $f|sed -e "s=DDD=${NOW}="|sed -e "s=XXX=${CLI}=" > $fout
      CLI=$(($CLI+1))
   done
done

f="CX_CH_EFW_L1_PY.ceh"
CLI=1
while test $CLI -le 4
do
	PROBE=1
	while test $PROBE -le 4
	do
		fout=`echo $f| sed -e "s=X=${CLI}=" | sed -e "s=Y=${PROBE}="`
		echo "Writing $fout"
		cat $f|sed -e "s=DDD=${NOW}="|sed -e "s=XXX=${CLI}="|sed -e "s=YYY=${PROBE}="|sed -e "s=YYY=${PROBE}=" > $fout
		PROBE=$(($PROBE+1))
	done
	CLI=$(($CLI+1))
done

f="CX_CH_EFW_L1_PYY.ceh"
CLI=1
while test $CLI -le 4
do
	PROBES="12 \
	32 \
	34"
	for PROBE in $PROBES
	do
		# No p32 on C4 yet
		[ "$CLI" = "4" ] && [ "$PROBE" = "32" ] && continue
		fout=`echo $f| sed -e "s=X=${CLI}=" | sed -e "s=YY=${PROBE}="`
		echo "Writing $fout"
		PROBE1=`echo $PROBE|cut -c 1`
		PROBE2=`echo $PROBE|cut -c 2`
		cat $f|sed -e "s=DDD=${NOW}="|sed -e "s=XXX=${CLI}="|sed -e "s=TTT=${PROBE}="|sed -e "s=YYY=${PROBE1}="|sed -e "s=ZZZ=${PROBE2}=" > $fout
		PROBE=$(($PROBE+1))
	done
	CLI=$(($CLI+1))
done

#fix dataset end date
flist="CX_CH_EFW_L1_E.ceh \
CX_CH_EFW_L2_E.ceh \
CX_CH_EFW_L3_E.ceh \
CX_CH_EFW_L3_SFIT.ceh \
CX_CH_EFW_L3_DER.ceh \
CX_CH_EFW_L2_EB.ceh"
CLI=1
ENDDATE="2018-12-10"
for f in $flist
do
	fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Processing $fout"
   sed -i.bak "s/MMM/${ENDDATE}/g" $fout
done
CLI=2
ENDDATE="2022-08-23"
for f in $flist
do
	fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Processing $fout"
   sed -i.bak "s/MMM/${ENDDATE}/g" $fout
done

CLI=3
ENDDATE="2011-03-05"
flist="CX_CH_EFW_L1_E.ceh \
CX_CH_EFW_L2_E.ceh"
for f in $flist
do
	fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Processing $fout"
   sed -i.bak "s/MMM/${ENDDATE}/g" $fout
done
ENDDATE="2014-11-03"
flist="CX_CH_EFW_L3_E.ceh \
CX_CH_EFW_L3_SFIT.ceh \
CX_CH_EFW_L3_DER.ceh \
CX_CH_EFW_L2_EB.ceh"
for f in $flist
do
	fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Processing $fout"
   sed -i.bak "s/MMM/${ENDDATE}/g" $fout
done
ENDDATE="2024-04-28"
flist="CX_CH_EFW_L2_PB.ceh \
CX_CH_EFW_L1_P.ceh \
CX_CH_EFW_L2_P.ceh \
CX_CH_EFW_L3_P.ceh" 
for f in $flist
do
	fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Processing $fout"
   sed -i.bak "s/MMM/${ENDDATE}/g" $fout
done

#default end of mission
ENDDATE="2026-12-31" 
flist="CX_CH_EFW_L1_P.ceh \
CX_CH_EFW_L2_P.ceh \
CX_CH_EFW_L3_P.ceh \
CX_CH_EFW_L1_E.ceh \
CX_CH_EFW_L2_E.ceh \
CX_CH_EFW_L3_E.ceh \
CX_CH_EFW_L3_SFIT.ceh \
CX_CH_EFW_L3_DER.ceh \
CX_CH_EFW_L2_PB.ceh \
CX_CH_EFW_L2_EB.ceh"
CLI=1
while test $CLI -le 4
do
	for f in $flist
	do
		fout=`echo $f| sed -e "s=X=${CLI}="`
   	echo "Processing $fout"
   	sed -i.bak "s/MMM/${ENDDATE}/g" $fout
	done
	CLI=$(($CLI+1))
done
