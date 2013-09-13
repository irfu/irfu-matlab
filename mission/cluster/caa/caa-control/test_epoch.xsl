<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>
		 
<xsl:import href="iso2epoch.xsl"/> 

<!-- 
This stylesheet adds ISDAT epoch for start time as an extra attribute.

(c) 2004-2005, Yuri Khotyaintsev

$Id$
 -->

<xsl:output method="text" version="1.0" encoding="UTF-8"/>

<xsl:template match="operation">
	<xsl:variable name="epo">
		<xsl:call-template name="iso2epoch">
			<xsl:with-param name="isotime" select="@start"/>
		</xsl:call-template>
	</xsl:variable>
	<xsl:value-of select="@start"/>
	<xsl:text> : </xsl:text>
	<xsl:choose>
		<xsl:when test="$epo = @cst"><xsl:text>Passed</xsl:text></xsl:when>
		<xsl:otherwise><xsl:text>Failed</xsl:text></xsl:otherwise>
	</xsl:choose>
</xsl:template>
</xsl:stylesheet>
