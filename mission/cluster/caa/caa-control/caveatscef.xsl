<?xml version='1.0' encoding='UTF-8'?>
<xsl:stylesheet version='1.0' xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
  xmlns:date='http://exslt.org/dates-and-times' xmlns:dyn='http://exslt.org/dynamic'>
<xsl:output method="text"/>
<xsl:param name="clusterno"></xsl:param>

<!--
This stylesheet converts xml_op_epo.xml file into separate
caveats .cef files.

Valid parameter clusterno: values 1-4

Uses exslt date & time and dynamic functions. Works with xsltproc.

Copyright (c) 2011-23 Jan Karlsson

Updated to exclude infinity (-1) stop time items

$Id$
 -->

<xsl:template match="/">!-------------------- CEF ASCII FILE --------------------|
! created on <xsl:value-of select="substring(date:date(),1,10)"/>
  <xsl:text> </xsl:text>
  <xsl:value-of select="substring(date:time(),1,8)"/>
!--------------------------------------------------------|
FILE_NAME = "C<xsl:value-of select="$clusterno"/>_CQ_EFW_INST__20010128_V02.cef"
FILE_FORMAT_VERSION = "CEF-2.0"
END_OF_RECORD_MARKER = "$"
include = "CL_CH_MISSION.ceh"
include = "C<xsl:value-of select="$clusterno"/>_CH_OBS.ceh"
include = "CL_CH_EFW_EXP.ceh"
include = "C<xsl:value-of select="$clusterno"/>_CH_EFW_INST.ceh"
include = "C<xsl:value-of select="$clusterno"/>_CQ_EFW_INST.ceh"
START_META     =   FILE_TYPE
   ENTRY       =   "cef"
END_META       =   FILE_TYPE
START_META     =   DATASET_VERSION
   ENTRY       =   "3"
END_META       =   DATASET_VERSION
START_META     =   LOGICAL_FILE_ID
   ENTRY       =   "C<xsl:value-of select="$clusterno"/>_CQ_EFW_INST__20010128_V02"
END_META       =   LOGICAL_FILE_ID
START_META     =   VERSION_NUMBER
   ENTRY       =   02
END_META       =   VERSION_NUMBER
START_META     =   FILE_TIME_SPAN
   VALUE_TYPE  =   ISO_TIME_RANGE
   ENTRY       =   2001-01-28T00:00:00Z/2022-01-01T00:00:00Z
END_META       =   FILE_TIME_SPAN
START_META     =   GENERATION_DATE
   VALUE_TYPE  =   ISO_TIME
   ENTRY       =   <xsl:value-of select="substring(date:date-time(),1,19)"/>Z
END_META       =   GENERATION_DATE
START_META     =   FILE_CAVEATS
   ENTRY       =   "Cluster <xsl:value-of select="$clusterno"/>. TBD stop time: 2022-01-01T00:00:00Z"
END_META       =   FILE_CAVEATS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       Data                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DATA_UNTIL = "END_OF_DATA"
<xsl:apply-templates select="//operation">
	<xsl:sort order="ascending" data-type="number" select="@startepoch"/>
</xsl:apply-templates>END_OF_DATA
</xsl:template>

<xsl:template match="@plan">
	<xsl:choose>
		<xsl:when test="starts-with('yes',.)">Planned</xsl:when>
		<xsl:when test="starts-with('no',.)">Not planned</xsl:when>
		<xsl:otherwise>Unknown</xsl:otherwise>
	</xsl:choose>
</xsl:template>

<xsl:template match="operation">
	<xsl:variable name="clnobool" select="dyn:evaluate(concat('@c',$clusterno))"/>
	<xsl:if test="starts-with($clnobool,'yes')">
		<xsl:variable name="dt" select="@dt"/>
		<xsl:if test="$dt != -1">
			<xsl:value-of select="@start"/>/<xsl:choose>
			<!--	<xsl:when test="$dt = -1">2025-01-01T00:00:00Z</xsl:when>        Infinity removed -->
				<xsl:when test="$dt = -157">2022-01-01T00:00:00Z</xsl:when> <!-- TBD -->
				<xsl:otherwise><xsl:value-of select="date:add(@start,concat('PT',$dt,'S'))"/></xsl:otherwise>
			</xsl:choose>,"<xsl:value-of select="@sdesc"/>","<xsl:apply-templates select="@plan"/>","<xsl:value-of select="translate(translate(child::desc,'&#x9;&#xa;',''),'&#x22;',&quot;&#x27;&quot;)"/>" $
</xsl:if> <!-- Not indented to not f**k up result file -->
	</xsl:if>
</xsl:template>

</xsl:stylesheet>
