<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>
<xsl:output method="html"/>

<!-- 
This stylesheet converts list XML file into 
simple HTML page.

(c) 2004-2005, Yuri Khotyaintsev

$Id$
 -->

<xsl:template match="/">
	<HTML>
	<HEAD>
		<META http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>
		<META http-equiv="Content-Language" content="en-us"/>
		<TITLE><xsl:value-of select="//title"/></TITLE>
		<LINK rel="StyleSheet" href="https://www.cluster.irfu.se/efw.css" type="text/css" media="screen"/>
	</HEAD>
	<BODY>
	<H1 ALIGN="center">
		EFW Non-Standard Operations Created (branch devel): <xsl:value-of select="$current-date"/>
	</H1>
	<H2>
		<xsl:value-of select="//description"/>
	</H2>
	<TABLE width="99%" cellpadding="2" cellspacing="2" border="1">
		<TR>
			<TH width="16%">Start</TH>
			<TH>dt (sec)</TH>
			<TH>Cluster ID</TH>
			<TH>Description</TH>
			<TH  width="20%">Resources</TH>
			<TH>Plan</TH>
			<TH>short Desc</TH>
		</TR>
		<xsl:apply-templates select="//operation">
			<xsl:sort order="ascending" data-type="number" select="@startepoch"/>
		</xsl:apply-templates>
	</TABLE>
	</BODY>
	</HTML>
</xsl:template>

<xsl:template match="res">
	<DIV>
	<xsl:for-each select="child::link">
		<xsl:number count="*" format="1"/>
		<xsl:text>. </xsl:text>
		<A HREF="https://www.cluster.irfu.se/efw/ops/ops_files/{@href}">
			<xsl:apply-templates/>
		</A><BR/>
	</xsl:for-each>
	</DIV>
</xsl:template>

<xsl:template match="@c1">
	<xsl:if test="contains('yes',.)">
		<xsl:text>1 </xsl:text>
	</xsl:if>
</xsl:template>
<xsl:template match="@c2">
	<xsl:if test="contains('yes',.)">
		<xsl:text>2 </xsl:text>
	</xsl:if>
</xsl:template>
<xsl:template match="@c3">
	<xsl:if test="contains('yes',.)">
		<xsl:text>3 </xsl:text>
	</xsl:if>
</xsl:template>
<xsl:template match="@c4">
	<xsl:if test="contains('yes',.)">
		<xsl:text>4 </xsl:text>
	</xsl:if>
</xsl:template>

<xsl:template match="operation">
	<TR>
		<TD>
			<xsl:call-template name="nicedate">
				<xsl:with-param name="isotime" select="@start"/>
			</xsl:call-template>
			<xsl:text> </xsl:text>
			<xsl:call-template name="nicetime">
				<xsl:with-param name="isotime" select="@start"/>
			</xsl:call-template>
		</TD>
		<TD>
			<xsl:value-of select="@dt"/>
		</TD>
		<TD>
			<xsl:apply-templates select="@c1"/>
			<xsl:apply-templates select="@c2"/>
			<xsl:apply-templates select="@c3"/>
			<xsl:apply-templates select="@c4"/>
		</TD>
		<TD>
			<xsl:value-of select="child::desc"/>
		</TD>
		<TD>
			<xsl:apply-templates select="child::res"/>
		</TD>
		<TD ALIGN="center">
			<xsl:value-of select="@plan"/>
		</TD>
		<TD ALIGN="center">
			<xsl:value-of select="@sdesc"/>
		</TD>
	</TR>
</xsl:template>

<xsl:template name="nicedate">
	<xsl:param name="isotime"/>
	<xsl:value-of select="substring($isotime,1,10)"/>
</xsl:template>
<xsl:template name="nicetime">
	<xsl:param name="isotime"/>
	<xsl:value-of select="substring($isotime,12,8)"/>
</xsl:template>

</xsl:stylesheet>
