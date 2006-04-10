<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>

<!-- 
(c) 2004-2005, Yuri Khotyaintsev

$Id$
 -->
 
<xsl:template name="iso2epoch">
	<!-- 
	This template converts ISO time string yyyy-mm-ddTHH:MM:ss.wwwwwwZ
	to ISDAT epoch 
	-->
	<xsl:param name="isotime"/>
	<xsl:variable name="y" select="number(substring($isotime,1,4))"/>
	<xsl:variable name="m" select="number(substring($isotime,6,2))"/>
	<xsl:variable name="d" select="number(substring($isotime,9,2))"/>
	<xsl:variable name="h" select="number(substring($isotime,12,2))"/>
	<xsl:variable name="min" select="number(substring($isotime,15,2))"/>
	<xsl:variable name="s" select="number(substring($isotime,18,string-length($isotime)-18))"/>
	
	<xsl:variable name="daym">
		<xsl:choose>
			<xsl:when test="($y mod 4) = 0">000031060091121152182213244274305335</xsl:when>
			<xsl:otherwise>000031059090120151181212243273304334</xsl:otherwise>
		</xsl:choose>
	</xsl:variable>
	<xsl:variable name="dy" select="$y - 1969"/>
	<!-- number of leap years -->
	<xsl:variable name="nly">
		<xsl:choose>
			<xsl:when test="($dy mod 4) = 0">
				<xsl:value-of select="($dy div 4)"/>
			</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="floor($dy div 4)"/>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:variable>
	
	<xsl:value-of select="($dy - $nly - 1)*31536000 + $nly*31622400 + ((number(substring($daym,($m - 1)*3 + 1,3) + $d - 1)*24 + $h)*60 + $min)*60 + $s"/>

</xsl:template>
</xsl:stylesheet>
