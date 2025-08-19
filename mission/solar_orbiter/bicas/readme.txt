
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
    (3) a string constant (NSOID/QRCID) that can be interpreted by BICAS, and
    (4) a human-readable comment string which is ignored.
This table is used to specify time intervals which should be treated in specific
ways in particular w.r.t. to data quality, e.g. modifying quality bits,
zVariable QUALITY_FLAG & L2/L3_QUALITY_BITMASK, blanking data etc. due to e.g.
thruster firings and saturation.



#########################################################
 Special cases, actions depending on quality information
#########################################################
BICAS takes many different actions depending on quality, different special
cases, and settings (configuration), and it can be very hard to keep track of
this. The purpose of this section is meant to summarize this information.

NOTE: This section might not be complete, but is intended to become so over
      time.
NOTE: bicas.const.qrc contains hard-coded constants describing some of this
      behaviour. The NSO table is also used for determining when some of these
      conditions apply.


=======================================
 Global capping (max) for QUALITY_FLAG
=======================================
The value of zVariable QUALITY_FLAG in any CDF produced by BICAS can not exceed
the value of setting PROCESSING.ZV_QUALITY_FLAG_MAX (=3 as of 2025-02-19).

Condition                      | Action taken when condition applies
--------------------------------------------------------------------
zVariable QUALITY_FLAG         | All output CDFs (L2+L3)!:
value exceeds setting          |     QUALITY_FLAG value is capped at the value
PROCESSING.ZV_QUALITY_FLAG_MAX |     of setting
when writing CDF file.         |     PROCESSING.ZV_QUALITY_FLAG_MAX.


===============
 BIAS off (BW)
===============
L1 LFR datasets read by BICAS contain a zVariable "BW" which informs BICAS on
whether BIAS was on or off. There is no known counterpart for TDS HK. Therefore,
for TDS data, BIAS is always assumed to be on.

Condition        | Action taken when condition applies
---------------------------------------------------------------
BIAS is known to | L2 CWF, SWF/RSWF:
be turned off    |     All samples are replaced by fill values.


=====================
 BIAS current sweeps
=====================
The BIAS subunit can be used for performing bias current sweeps. This
temporarily changes the usual constant bias current which is otherwise supposed
to be constant over much longer time periods while still not interrupting
the sampling of science data. BICAS can determine whether BIAS is sweeping by
using function bicas.proc.L1L2.swpdet.SBDA_SCDA_with_margins().

Condition            | Action taken when condition applies
------------------------------------------------------------------
"BIAS current sweep" | L2 CWF, SWF/RSWF:
                     |    All samples are replaced by fill values.


------------------------------------------------------------------
Algorithm for determining whether a BIAS current sweep is underway
------------------------------------------------------------------
L1/L1R datasets are planned to eventually contain information on sweeps
being executed in zVariable QUALITY_BITMASK but at the time of writing
(2025-02-13), this has still not been implemented (is not in use) by ROC.
Until this is made use of, BICAS uses approximative algorithms for
determining whether a sweep takes place.

Before a certain timestamp specified by setting
PROCESSING.L2.SWEEP_DETECTION.SBDA_SCDA_BOUNDARY_UTC, BIAS demultiplexer mode 4
is taken as equivalent to a sweep being underway. After this timestamp, BIAS HK
currents varying outside some min-max threshold within some moving time window
is taken to be equivalent to a sweep being underway. If BIAS HK is missing, then
it is assumed that no sweep is underway.


===============
 L2 saturation
===============
BICAS can label L2 data as saturated according to one of two schemes (setting
PROCESSING.SATURATION.QUALITY_SCHEME): one old scheme (to be phased out), and
one new scheme (under evaluation; to be phased in).


-----------------------------
Old scheme: GLOBAL_SATURATION
-----------------------------
BICAS can classify saturated L2 data is either "partially saturated" or
"fully saturated" in the following ways:
(1) time interval is labelled as "partially saturated" in the NSO table,
(2) time interval is labelled as "fully saturated" in the NSO table,
(3) algorithm examines L2 data and determines that it is "fully saturated".
NOTE: There is no algorithm for determining whether L2 data is "partially
      saturated".

When BICAS determines that one of these conditions apply (for L2 data), it takes
the following actions:

Condition             | Action taken when condition applies
-------------------------------------------------------------------------
"partially saturated" | L2 CWF, SWF/RSWF:
                      |     Set L2_QUALITY_BITMASK: "partially saturated"
                      |     Cap QUALITY_FLAG<=1.
-------------------------------------------------------------------------
"fully saturated"     | L2 CWF, SWF/RSWF:
                      |     Set L2_QUALITY_BITMASK: "full saturation"
                      |                             AND "partial saturation"
                      |     Set QUALITY_FLAG=0.


------------------------------
New scheme: CHANNEL_SATURATION
------------------------------
There is one quality bit per channel, where AC and DC diffs share the same bits.
Reconstructed channels (i.e. channel samples which BICAS has derived from other
channels) and are affected by saturation on the source channels are also
labelled as saturated.

Algorithm for detecting saturation:
Step 1: For every sample, saturation is preliminarily detected using
thresholds for every sample separately.
Step 2: For a moving window on each channel (both source and reconstructed
channels), if the fraction of saturated samples (inverse sampling
frequency-weighted) from step (1) exceeds a threshold, then the entire window
counts as saturated.

Note 1: Due to step (2), individual samples may be saturated (as in exceeding
thresholds) without being labelled as saturated.
Note 2: Due to step (2), a reconstructed channel can be labelled as saturated
despite that neither of the two source channels is labelled as saturated in
unusual cases.

Condition            | Action taken when condition applies
--------------------------------------------------------------------------------
Channel is saturated | L2 CWF, SWF/RSWF:
                     |      Set QUALITY_FLAG=0.
                     |      Set L2_QUALITY_BITMASK quality bit for
                            the corresponding channel (V1/V2/V3/V12/V13/V23).


=================
 Thruster firing
=================
BICAS can determine that data is affected by a SolO thruster firing by reading
the NSO table.

NOTE: L1/L1R QUALITY_BITMASK contains a quality bit for thruster firings but at
      the time of writing (2025-08-14), this is not being used, due to covering
      too broad time intervals (multiple hours).

Condition         | Action taken when condition applies
-------------------------------------------------------
"thruster firing" | L2 CWF, SWF/RSWF:
                  |     Cap QUALITY_FLAG<=1.
                  |     NOTE: No quality bit is set.


=============
 BIAS sweeps
=============
BICAS can determine whether L2 data has been colleted during a BIAS sweep (bias
currents are changed according to a preprogrammed scheme within a short period
of time).

NOTE: Which data should be excluded due to being affected by sweeps is not
      well defined since some data before and after the actual sweep is affected
      due to commanding. Sweeps are detected(!) using algorithms which are not
      perfect: BICAS uses one of two different algorithms (SBDA or SCDA)
      depending on the time and it adds a customizable time margin in addition
      to that and therefore may remove too much or too little data. BICAS does
      not (yet) use the QUALITY_BITMASK bits to detect sweeps (as of
      2025-08-14).

Condition         | Action taken when condition applies
--------------------------------------------------------------------
Sweep is detected | L2 CWF, SWF/RSWF:
                  |     Voltages and currents are set to fill value.


==============
 ANT3 failing
==============
BICAS can determine whether ANT3 is failing by reading the NSO table.

Condition       | Action taken when condition applies
----------------------------------------------------------------------
ANT3 is failing | L2 CWF, SWF/RSWF:
                |     V3, V13_DC/AC, V23_DC/AC are set to fill values.


===============
 "Bad density"
===============
BICAS can determine that the derived L3 density is "bad density" (less
reliable). This is determined by solo.psp2ne() returning a quality bit
("bad density") when deriving density. This applies to SOLO_L3_RPW-BIA-DENSITY
and SOLO_L3_RPW-BIA-DENSITY-10-SECONDS.

Condition     | Action taken when condition applies
---------------------------------------------------------
"bad density" | L3 DENSITY:
              |     Set L3_QUALITY_BITMASK: "bad density".
              |     Cap QUALITY_FLAG<=1.


=================================================================
 L2 data used for generating L2 LFR CWF downsampled (unofficial)
=================================================================
L2 LFR CWF downsampled (SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND; unofficial dataset)
is derived from SOLO_L2_RPW-LFR-SURV-CWF-E only, but only when its quality
deemed good enough using setting PROCESSING.L2-CWF-DSR.ZV_QUALITY_FLAG_MIN (=0
as of 2025-05-23).

Condition                                    | Action taken when condition applies
----------------------------------------------------------------------------------
L2 LFR CWF QUALITY_FLAG is either            | L2 LFR CWF downsampled:
>= PROCESSING.L2-CWF-DSR.ZV_QUALITY_FLAG_MIN |     VDC, EDC, VDCSTD, EDCSTD
or fill value (!)                            |     values are set to fill
                                             |     values before downsampling.


=====================================================
 L2 data used for generating L3 DENSITY,EFIELD,SCPOT
=====================================================
L3 DENSITY+EFIELD+SCPOT is derived from SOLO_L2_RPW-LFR-SURV-CWF-E alone, but
only when its quality is deemed good enough using setting
PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN (=2 as of 2025-08-19). Saturated L2 data
as described by L2_QUALITY_BITMASK (when channel saturation is
enabled) is set to fill values before being passed on to solo.vdccal() and
psp2ne().

NOTE: PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN is *NOT* compared to the literal
      input SOLO_L2_RPW-LFR-SURV-CWF-E QUALITY_FLAG but a reconstructed
      "L2 QUALITY_FLAG" value which ignores saturation (too complicated to
      describe here; see source code).

Condition                                 | Action taken when condition applies
--------------------------------------------------------------------------------
A reconstructed variant of L2 LFR         | L3 DENSITY+EFIELD+SCPOT:
CWF QUALITY_FLAG is either                |     DENSITY+EFIELD+SCPOT values are
< PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN |     set to fill values.
or fill value (!)                         |

NOTE: Historically, this rule did not apply to DENSITY, but it does since
BICAS 8.5.0.


Condition                   | Action taken when condition applies
--------------------------------------------------------------------------
An L2 channel is saturated. | L3 DENSITY+EFIELD+SCPOT:
                            |     The affected channel(s) is not used for
                            |     deriving L3 data.


=============================================================
 Number of required actual samples per bin when downsampling
=============================================================
When downsampling data, samples are grouped into bins (short time intervals).
Each bin of samples is converted to one value (per channel), and for science
data also one additional modified standard deviation value (per channel).

Condition                               | Action taken when condition applies
------------------------------------------------------------------------------
Fewer than                              | Bin is represented by a fill value
bicas.const.N_MIN_OSR_SAMPLES_PER_BIN   | in the output data (both downsampled
non-fill value samples in a bin         | value and modified standard
when downsampling *SCIENCE* data.       | deviation).
------------------------------------------------------------------------------
Zero non-fill value samples in a bin    | Bin is represented by a fill value
when downsampling quality bitmasks      | in the output data.
(zVariables QUALITY_BITMASK,            |
L2_QUALITY_BITMASK, L3_QUALITY_BITMASK) |
and zVariable QUALITY_FLAG.             |


=========================================================
 ROC modifying BICAS output CDFs due for quality reasons
=========================================================
ROC sometimes modifies CDFs, *after* they have been produced by BICAS and before
they are officially delivered to SOAR. This may happen e.g. for time intervals
when e.g. ANT3 has failed and QUALITY_FLAG is capped to 1 by ROC. See global
attribute CAVEATS.
