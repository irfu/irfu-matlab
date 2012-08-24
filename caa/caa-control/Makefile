#
# $Id$
#
XSLTPROC=xsltproc

all: check ns_ops_epo.xml txt cef
	${XSLTPROC} html.xsl ns_ops_epo.xml > ns_ops.html
	rm -f ns_ops_epo.xml

txt: ns_ops_epo.xml
	${XSLTPROC} text.xsl ns_ops_epo.xml

cef: ns_ops_epo.xml
	for cln in 1 2 3 4 ; do \
		${XSLTPROC} --param clusterno $$cln caveatscef.xsl ns_ops_epo.xml > C$${cln}_CQ_EFW_INST_V00.cef ; \
	done

ns_ops_epo.xml: ns_ops.xml
	${XSLTPROC} addepoch.xsl ns_ops.xml > ns_ops_epo.xml

clean:
	rm -f ns_ops_epo.xml ns_ops_c[1-4].dat

test: test_epoch.xml test_epoch.xsl
	echo "Tesing iso2epoch..."
	${XSLTPROC} test_epoch.xsl test_epoch.xml

update: txt
	scp ns_ops_c[1-4].dat ns_ops.xml yuri@db:/export/data/cluster/caa-control

check:
	@echo "Checking for incorrect ISO time..."
	@echo `grep "start=" ns_ops.xml| grep -v 'Z"'`
