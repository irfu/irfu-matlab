#
# Script to process ns_ops.xml
#
XSLTPROC=xsltproc

all: html txt cef

html: check ns_ops_epo.xml
	${XSLTPROC} --stringparam current-date `date +%Y-%m-%d` html.xsl ns_ops_epo.xml > ns_ops.html

txt: ns_ops_epo.xml
	${XSLTPROC} text.xsl ns_ops_epo.xml

cef: ns_ops_epo.xml
	for cln in 1 2 3 4 ; do \
		${XSLTPROC} --param clusterno $$cln caveatscef.xsl ns_ops_epo.xml > C$${cln}_CQ_EFW_INST__20010128_V02.cef ; \
	done

ns_ops_epo.xml: ns_ops.xml
	${XSLTPROC} addepoch.xsl ns_ops.xml > ns_ops_epo.xml

clean:
	rm -f ns_ops.html ns_ops_epo.xml ns_ops_c[1-4].dat C[1-4]_CQ_EFW_INST__20010128_V[0-9][0-9].cef

test: test_epoch.xml test_epoch.xsl
	echo "Tesing iso2epoch..."
	${XSLTPROC} test_epoch.xsl test_epoch.xml

check:
	@echo "Checking for incorrect ISO time..."
	@echo `grep "start=" ns_ops.xml| grep -v 'Z"'`
