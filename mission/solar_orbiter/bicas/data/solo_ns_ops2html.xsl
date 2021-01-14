<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">



<!--#######################
XSLT RULE FOR THE MAIN BODY
########################-->
<xsl:template match="main">

<HTML>
<HEAD>
    <TITLE>Non-Standard Operations (NSO) List</TITLE>
</HEAD>

<BODY>
    <H1 ALIGN="center">
    Non-Standard Operations (NSO) List
    </H1>

    <xsl:apply-templates select="description"/>

    <xsl:apply-templates select="notesTable"/>

    <H2>Events</H2>
    <TABLE style="width:100%" cellpadding="5" cellspacing="1" border="1">

        <TR>
            <TH>
                Start time (UTC)
            </TH>
            <TH>
                Stop time (UTC)
            </TH>
            <TH>
                NSO ID

            </TH>
            <TH>
                Description (free-form)
            </TH>
        </TR>

        <xsl:apply-templates select="eventsTable"/>
    </TABLE>

</BODY>
</HTML>

</xsl:template>



<!--####################################################
TRICK TO MAKE HTML TAG <B> in XML NOT BE REMOVED BY XLST
#####################################################-->
<xsl:template match="B">
    <B>
        <xsl:apply-templates/>
    </B>
</xsl:template>



<!--##################################
XSLT RULE FOR NOTES AS A WHOLE "TABLE"
###################################-->
<xsl:template match="notesTable">
    <H2>Notes</H2>
    <xsl:apply-templates/>
</xsl:template>



<!--####################
XSLT RULE FOR EACH NOTE
#####################-->
<xsl:template match="note">
    <P>
    <xsl:value-of select="."/>
    </P>
</xsl:template>



<!--####################
XSLT RULE FOR EACH EVENT
#####################-->
<xsl:template match="event">
    <tr>
        <td>
            <xsl:value-of select="startTimeUtc"/>
        </td>

        <td>
            <xsl:value-of select="stopTimeUtc"/>
        </td>

        <td>
            <xsl:value-of select="rcsNsoId"/>
        </td>

        <td>
            <xsl:value-of select="description"/>
        </td>
    </tr>
</xsl:template>



</xsl:stylesheet>
