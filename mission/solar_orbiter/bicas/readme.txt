
#############
 About BICAS
#############
This software, BICAS, is created for the calibration of the BIAS subsystem in
the RPW instrument on the Solar Orbiter spacecraft. The principle author of
this software is Erik P G Johansson, Swedish Institute of Space Physics (IRF),
Uppsala, Sweden. Software development began 2016-03-xx (March 2016).

The software was originally developed for processing
(1) L1R-->L2, for being run at ROC for official processing (deliveries of L2 to SOAR).
After that, the software has been extended analogously to also process more datasets at IRFU
(2) L2-->L3, for official deliveries IRFU-->ROC (and further to SOAR), and
(3) L2-->L2 (SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND; IRFU-internal unofficial dataset)
.

IMPORTANT NOTE: BICAS is designed to comply with the RCS ICD. Much documentation
on how to use this software can thus be found there. For more documentation,
see RCS ICD and RUM documents (see below).



#################
 Main executable
#################
<BICAS root dir>/roc/bicas



#############################
 CLI syntax / CLI parameters
#############################
NOTE: The official CLI parameter syntax is defined in RCS ICD, Iss02 Rev02, Section 3.2.

SYNTAX 1: ( --version | --identification | --swdescriptor | --help ) <General parameters>
SYNTAX 2: <S/W mode> <General parameters, Output parameter, Specific inputs parameters>

NOTE: In syntax 2, the position of the first arguments is important. The order
of all other (groups of) arguments is arbitrary.

--version          Print the software version.
--identification   Print the S/W descriptor release segment.
--swdescriptor     Print the S/W descriptor (not RCS ICD requirement).
--help             Print "help-ish" text



=========================
 Common input parameters
=========================
--log    <absolute path to log file>   (optional) Specifies log file.
--config <absolute path to file>       (optional) Specifies the configuration
                                       file to use.


<S/W mode>   Selects the S/W mode to use.
Available S/W modes can be found in the S/W descriptor. They are listed under
"modes". "name" specifies the string that identifies a given mode and can be
used on as CLI argument.



===========================
 Specific input parameters
===========================
Set of parameters which specify input CDF files. The exact set depends on the
exact S/W mode and can in principle be read from the S/W descriptor.
Required input parameters for a specific S/W mode can be found in the S/W
descriptor under "modes" --> (specific mode) --> "inputs"
--> Name of subsection, e.g. "input_hk", "input_sci".
Example: an input subsection "input_hk" means that there is a required parameter
"--input_hk <path_to_file>".



===============
 Example calls
===============
bicas --version --config ~/bicas.conf
bicas --identification
bicas --identification --log ~/bicas.log
bicas LFR-SURV-CWF-E
    --in_sci   L1R/2020/04/14/solo_L1R_rpw-lfr-surv-swf-e-cdag_20200414_V01.cdf
    --in_cur   BIA/2020/04/solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf
    --in_hk    HK/2020/04/14/solo_HK_rpw-bia_20200414_V02.cdf
    --out_sci  solo_L2_rpw-lfr-surv-swf-e_20200414T000000-20200415T000000_V02.cdf'
    --config   /home/erjo/bicas_batch.conf



###########################################
 Installation, set-up, system requirements
###########################################
See "install.txt".



####################################
 Known current limitations, caveats
####################################
For limitations and caveats, see the official user manual, the RUM
document.



###################################################################
 Source code conventions and dictionary of terms and abbreviations
###################################################################
See "misc_conventions.txt".



#####################################
 Non-Standard Operations (NSO) table
#####################################
BICAS contains a table/list of manually hard-coded time intervals called
"the NSO table" (NSO=Non-Standard Operations) located in
irfu-matlab/mission/solar_orbiter/bicas/data/solo_ns_ops.xml.
Every entry in the table consists of:
    (1) a start timestamp,
    (2) a stop timestamp,
    (3) a string constant that can be interpreted by BICAS, and
    (4) a human-readable comment string which is ignored.
This table is used to specify time intervals which should be treated in specific
ways in particular w.r.t. to data quality, e.g. modifying quality bits,
zVariable QUALITY_FLAG, blanking data etc. due to e.g. thruster firings and
saturation.



#########################################################
 Actions depending on quality information, special cases
#########################################################
BICAS takes many different actions depending on quality and different special
cases and it can be very hard to keep track of these. The purpose of this
section is meant to summarize this information.

NOTE: This section is very incomplete, but is intended to be amended over time.
NOTE: bicas.const contains hard-coded constants describing some of this
      behaviour. The NSO table is used for determining when certain conditions
      apply.


===============
 BIAS off (BW)
===============
L1 LFR datasets read by BICAS contain a zVariable "BW" which informs BICAS on
whether BIAS was on or off. There is no known counterpart for TDS HK. Therefore,
for TDS data, BIAS is always assumed to be on.

Condition          | Action taken when condition applies (L2 CWF,SWF,RSWF)
------------------------------------------------------------------------------
BIAS is turned off | All samples are replaced by fill values.


=====================
 BIAS current sweeps
=====================
The BIAS subunit can be used for performing bias current sweeps. This
temporarily changes the usual constant bias current which is assumed to be
constant over much longer time periods while not interrupting the sampling of
data. BICAS can determine whether BIAS is sweeping by using function
bicas.proc.L1L2.swpdet.SBDA_SCDA_with_margins().

NOTE: L1 datasets are planned to eventually contain information on sweeps being
      executed in zVariable QUALITY_BITMASK but at the time of writing
      (2025-01-16), this has still not been implemented.

Condition            | Action taken when condition applies (L2 CWF,SWF,RSWF)
------------------------------------------------------------------------------
"BIAS current sweep" | All samples are replaced by fill values.


============
 Saturation
============
BICAS can determine that L2 data is either "partially saturated" or
"fully saturated" in the following ways:
(1) time interval being labelled as "partially saturated" in the NSO table,
(2) time interval being labelled as "fully saturated" in the NSO table,
(3) algorithm examines L2 data and determines that it is "fully saturated".
NOTE: There is no algorithm for determining whether L2 data is "partially
      saturated".

When BICAS determines that one of these conditions apply (for L2 data), it takes
the following actions:

Condition             | Action taken when condition applies (L2 CWF,SWF,RSWF)
------------------------------------------------------------------------------
"partially saturated" | Set L2_QUALITY_BITMASK: "partially saturated"
                      | Cap QUALITY_FLAG<=1
------------------------------------------------------------------------------
"fully saturated"     | Set L2_QUALITY_BITMASK: "full saturation"
                      |                         AND "partial saturation"
                      | Cap QUALITY_FLAG<=0


=================
 Thruster firing
=================
BICAS can determine that data is affected by a SolO thruster firing by reading
the NSO table.
NOTE: L1 datasets are planned to eventually contain information on thruster
      firings being executed but at the time of writing (2025-01-16), this has
      still not been implemented.

Condition         | Action taken when condition applies (L2 CWF,SWF,RSWF)
--------------------------------------------------------------------------
"thruster firing" | Cap QUALITY_FLAG<=1
                  | NOTE: No quality bit is set.


===============
 "Bad density"
===============
BICAS can determine that the derived L3 density is "bad density" (less
reliable). This is determined by solo.psp2ne() returning a quality bit
("bad density") when deriving density. This applies to SOLO_L3_RPW-BIA-DENSITY
and SOLO_L3_RPW-BIA-DENSITY-10-SECONDS.

Condition     | Action taken when condition applies (DENSITY datasets)
-----------------------------------------------------------------------
"bad density" | Set L3_QUALITY_BITMASK: "bad density"
              | Cap QUALITY_FLAG<=1


=====================================================
 L2 data used for generating L3 DENSITY,EFIELD,SCPOT
=====================================================
Calculations for deriving L3 EFIELD+SCPOT from SOLO_L2_RPW-LFR-SURV-CWF-E
only uses SOLO_L2_RPW-LFR-SURV-CWF-E data for which
    QUALITY_FLAG >= PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN (=2)
or
    QUALITY_FLAG == fill value (!)
.

