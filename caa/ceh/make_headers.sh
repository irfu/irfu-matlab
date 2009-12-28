#!/bin/sh
#
# $Id$

f=CX_CH_EFW_INST.ceh
CLI=1
while test $CLI -le 4
do
   fout=`echo $f| sed -e "s=X=${CLI}="`
   echo "Writing $fout"
   cat $f|sed -e "s=XXX=${CLI}=" > $fout
   cat C${CLI}_CH_EFW_INST_CAV.ceh >> $fout
   CLI=$(($CLI+1))
done

flist="CX_CH_EFW_L2_P.ceh \
CX_CH_EFW_L3_P.ceh \
CX_CH_EFW_L2_E.ceh \
CX_CH_EFW_L3_E.ceh \
CX_CH_EFW_L3_DER.ceh"


for f in $flist
do
   CLI=1
   while test $CLI -le 4
   do
      fout=`echo $f| sed -e "s=X=${CLI}="`
      echo "Writing $fout"
      cat $f|sed -e "s=XXX=${CLI}=" > $fout
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
		cat $f|sed -e "s=XXX=${CLI}="|sed -e "s=YYY=${PROBE}=" > $fout
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
		fout=`echo $f| sed -e "s=X=${CLI}=" | sed -e "s=YY=${PROBE}="`
		echo "Writing $fout"
		cat $f|sed -e "s=XXX=${CLI}="|sed -e "s=YYY=${PROBE}=" > $fout
		PROBE=$(($PROBE+1))
	done
	CLI=$(($CLI+1))
done
