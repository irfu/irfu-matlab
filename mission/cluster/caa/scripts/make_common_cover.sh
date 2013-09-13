#!/bin/sh
#
# Make a common file from COVERAGE_##_01.TXT
#
# (c) 2005, Yuri Khotyaintsev
#
# $Id$

sc_list="1 \
2 \
3 \
4"
tmp_f=common_tmp.txt

rm -rf $tmp_f; touch $tmp_f

for cli in $sc_list
do
	cat COVERAGE_$cli'_01.TXT'|grep Normal |awk "{print \$1 \"T\" \$2 \".000Z \" \$4\"T\" \$5 \".000Z $cli 0\"}" >> $tmp_f
	cat COVERAGE_$cli'_01.TXT'|grep Burst  |awk "{print \$1 \"T\" \$2 \".000Z \" \$4\"T\" \$5 \".000Z $cli 1\"}" >> $tmp_f
done

cat $tmp_f| sort > COVERAGE_COMM.dat
rm $tmp_f
