<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>
		 
<xsl:import href="iso2epoch.xsl"/> 

<!-- 
This stylesheet adds ISDAT epoch for start time as an extra attribute.

(c) 2004-2005, Yuri Khotyaintsev

$Id$
 -->

<xsl:output method="xml" version="1.0" encoding="UTF-8"/>

<xsl:template match="node() | @*">
	<xsl:copy>
		<xsl:apply-templates select="@* | node()"/>
	</xsl:copy>
</xsl:template>

<xsl:template match="operation">
	<xsl:element name="{local-name(.)}" namespace="{namespace-uri(..)}">
		<xsl:attribute name="startepoch">
			<xsl:call-template name="iso2epoch">
				<xsl:with-param name="isotime" select="@start"/>
			</xsl:call-template>
		</xsl:attribute>
		<xsl:copy-of select="@*"/>
		<xsl:apply-templates/>
	</xsl:element>
</xsl:template>
</xsl:stylesheet>
