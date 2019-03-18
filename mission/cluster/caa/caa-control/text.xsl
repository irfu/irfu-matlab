<xsl:stylesheet version="1.0"
     xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	 xmlns:saxon="http://icl.com/saxon"
	 extension-element-prefixes="saxon">
<xsl:output method="text"/>

<!-- 
This stylesheet converts list XML file into four text files 
ns_ops_c[1-4].dat which can be loaded into Matlab

Substitutions for SDESC:
no_tm			0
bad_tm			1
bad_data		2
no_p[1-4]		1[1-4]
no_p[12/32/34]		1[12/32/34]
no_p23                  132 (alias for no_p32)
hxonly			15
bad_bias		16
bad_hx                  17
bad_lx                  18
high_bias		19
no_10Hz_filt	255
no_spin_fits	260
spec_bias		261
TDB				999
not listed		9999

(c) 2004-2005, Yuri Khotyaintsev

$Id$
 -->

<xsl:template match="/">
	<saxon:output href="ns_ops_c1.dat" method="text">
		<xsl:apply-templates select="//operation[@c1='yes']">
			<xsl:sort order="ascending" data-type="number" select="@startepoch"/>
		</xsl:apply-templates>
	</saxon:output>
	<saxon:output href="ns_ops_c2.dat" method="text">
		<xsl:apply-templates select="//operation[@c2='yes']">
			<xsl:sort order="ascending" data-type="number" select="@startepoch"/>
		</xsl:apply-templates>
	</saxon:output>
	<saxon:output href="ns_ops_c3.dat" method="text">
		<xsl:apply-templates select="//operation[@c3='yes']">
			<xsl:sort order="ascending" data-type="number" select="@startepoch"/>
		</xsl:apply-templates>
	</saxon:output>
	<saxon:output href="ns_ops_c4.dat" method="text">
		<xsl:apply-templates select="//operation[@c4='yes']">
			<xsl:sort order="ascending" data-type="number" select="@startepoch"/>
		</xsl:apply-templates>
	</saxon:output>
</xsl:template>

<xsl:template match="operation">
	<xsl:value-of select="@startepoch"/><xsl:text> </xsl:text>
	<xsl:value-of select="@dt"/><xsl:text> </xsl:text>
	<xsl:call-template name="yn2n">
		<xsl:with-param name="s" select="@plan"></xsl:with-param>
	</xsl:call-template><xsl:text> </xsl:text>
	<xsl:choose>
		<xsl:when test="@sdesc='no_tm'">0</xsl:when>
		<xsl:when test="@sdesc='bad_tm'">1</xsl:when>
		<xsl:when test="@sdesc='bad_data'">2</xsl:when>
		<xsl:when test="@sdesc='no_p1'">11</xsl:when>
		<xsl:when test="@sdesc='no_p2'">12</xsl:when>
		<xsl:when test="@sdesc='no_p3'">13</xsl:when>
		<xsl:when test="@sdesc='no_p4'">14</xsl:when>
		<xsl:when test="@sdesc='hxonly'">15</xsl:when>
		<xsl:when test="@sdesc='bad_bias'">16</xsl:when>
		<xsl:when test="@sdesc='bad_hx'">17</xsl:when>
		<xsl:when test="@sdesc='bad_lx'">18</xsl:when>
		<xsl:when test="@sdesc='high_bias'">19</xsl:when>
		<xsl:when test="@sdesc='no_p12'">112</xsl:when>
		<xsl:when test="@sdesc='no_p23'">132</xsl:when>
		<xsl:when test="@sdesc='no_p32'">132</xsl:when>
		<xsl:when test="@sdesc='no_p34'">134</xsl:when>
		<xsl:when test="@sdesc='no_10Hz_filt'">255</xsl:when>
		<xsl:when test="@sdesc='no_spin_fits'">260</xsl:when>
		<xsl:when test="@sdesc='spec_bias'">261</xsl:when>
		<xsl:when test="@sdesc='TBD'">999</xsl:when>
		<xsl:otherwise>9999</xsl:otherwise>
	</xsl:choose>
	
	<xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template name="yn2n">
	<xsl:param name="s"/>
	<xsl:choose>
		<xsl:when test="starts-with($s,'no')">
			<xsl:value-of select="0"/>
		</xsl:when>
		<xsl:otherwise>
			<xsl:value-of select="1"/>
		</xsl:otherwise>
	</xsl:choose>
</xsl:template>

</xsl:stylesheet>
