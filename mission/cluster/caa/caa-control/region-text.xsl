<xsl:stylesheet version="1.0"
     xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	 xmlns:saxon="http://icl.com/saxon"
	 extension-element-prefixes="saxon">
<xsl:output method="text"/>

<!-- 
This stylesheet converts regions XML file into a text file 
which can be loaded into Matlab

Substitutions for DESC:
mp			1
sh			2
bs			3
sw			4
psh			5
psp			6
lo			7
cusp		8
az      9
nodata		0
undef		999

(c) 2004-2005, Yuri Khotyaintsev

$Id$
 -->

<xsl:template match="/">
	<xsl:apply-templates select="//region"/>
</xsl:template>
		
<xsl:template match="region">
	<xsl:value-of select="@start"/><xsl:text> </xsl:text>
	<xsl:choose>
		<xsl:when test="@desc='nodata'">0</xsl:when>
		<xsl:when test="@desc='mp'">1</xsl:when>
		<xsl:when test="@desc='sh'">2</xsl:when>
		<xsl:when test="@desc='bs'">3</xsl:when>
		<xsl:when test="@desc='sw'">4</xsl:when>
		<xsl:when test="@desc='psh'">5</xsl:when>
		<xsl:when test="@desc='psp'">6</xsl:when>
		<xsl:when test="@desc='lo'">7</xsl:when>
		<xsl:when test="@desc='cusp'">8</xsl:when>
		<xsl:when test="@desc='az'">9</xsl:when>
		<xsl:when test="@desc='undef'">999</xsl:when>
		<xsl:otherwise>9999</xsl:otherwise>
	</xsl:choose>
	
	<xsl:text>&#xa;</xsl:text>
</xsl:template>

</xsl:stylesheet>
