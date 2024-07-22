
#############
 About BICAS
#############
BICAS = BIAS Calibration Software

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



####################
 Naming conventions
####################
The naming conventions are partly inconsistent for historical reasons.
The code should however converge on the following:
- Use the defined abbreviations in identifiers.
- Variables names: Should use camelCase.
    - Variables containing structs or instances of classes: Should have an
      uppercase initial.
    - All other variables: Should use a lowercase initial.
    - Exception: Variables which are direct analogues to zVariables are named
      as the corresponding zVariables, i.e. SCREAMING_SNAKE_CASE most of
      the time.
- Classes (not instances of classes):
  - Classes which are actually instantiated (which are not merely a collection
    of static methods): Should be CamelCase. Initial in uppercase.
  - Classes which are not instantiated: Should have short names and use
    lowercase without underscore.
- Functions should be named using snake_case but with abbreviations in
  uppercase.
- Special types of variables:
    b    = Logical (boolean) values, often logical arrays used for logical
           indexing.
    i, j = Indices into arrays (linear indexing).
    L    = Singleton instance of class bicas.Logger.
--
Unit tests (classes) have suffix "___UTEST".



###########################
 Abbreviations, dictionary
###########################
NOTE: This list primarily applies to comments and identifiers in the source
code. Most of the abbreviations are not "official" and have been created
specifically for BICAS.
--
AA, aampere
    Antenna Ampere. Calibrated ampere at the antenna.
AAPT
    Antenna Ampere/TM
ACHG
    (BIAS diff) AC High Gain.
ACLG
    (BIAS diff) AC Low Gain.
ASID
    Class bicas.proc.L1L2.AntennaSignalId.
ASR, Antenna Signal Representation.
    (1) A reference to a particular "physical antenna signal"
        (a "data channel") which BIAS-LFR/TDS is trying to measure, or
        NOTE: Class ASID represents such a reference.
    (2) a measurement of such.
    In reality, the terminology is:
    ASR         : Pointer to a specific physical antenna signal, e.g. V12_LF
                  (DC diff, antenna 1-2)
    ASR samples : Samples representing a specific ASR (as opposed to BLTS).
    NOTE: There are 9 ASRs (DC: 3 singles, 3 diffs; AC: 3 diffs), i.e. they
    can refer also to signals not represented by any single BLTS, given a
    chosen mux mode (BDM) and demultiplexer latching relay (DLR) setting.
AV, avolt, Antenna Volt
    Calibrated volt at the antennas, i.e. the final calibrated (measured)
    value, including for reconstructed signals (e.g. diffs calculated from
    singles). May also refer to offsets and values without offsets.
AVPIV
    Antenna Volt / (Per) Interface Volt
BDM, "mux mode"
    BIAS Demultiplexer Mode. Corresponds to BIAS HK ZV "HK_BIA_MODE_MUX_SET".
    Abbreviation partly just to have one consistent name (replacing "mux",
    "mux mode", "demux mode", "demultiplexer mode", "bias mux mode", etc.)
    and to distinguish it from "muxing" as a verb. Mostly used internally,
    not in the interface.
BFM
    BICAS Functionality Mode. One of the basic modes of BICAS operations:
    Print version, print identifiction, print version, print help, or process
    datasets. Exactly one of these must apply every time BICAS is run.
BIAS specification
    Document RPW-SYS-MEB-BIA-SPC-00001-IRF, "RPW Instrument -- BIAS
    Specification".
BIAS_1, ..., BIAS_5 (BIAS_i, i=1..5)
    Defined in BIAS specifications document. Equal to the physical signal at
    the physical boundary between BIAS and LFR/TDS. Unit: Interface volt.
    Should be, and mostly is, replaced by BLTS+specified unit in
    the implementation.
BLTS = BIAS-LFR/TDS SIGNAL
    Signals somewhere between the LFR/TDS ADCs and the non-antenna side of the
    BIAS demuxer including the BIAS transfer functions. Like BIAS_i, i=1..5,
    but includes various stages of calibration/non-calibration, including in
    particular
      - TM units (inside LFR/TDS),
      - Interface volt (at the physical boundary BIAS-LFR/TDS (BIAS_i)), and
      - Calibrated values inside BIAS but without demuxer addition and
        subtraction inside BIAS (i.e. including using BIAS offsets, BIAS
        transfer functions; volt).
    NOTE: This definition is partly created to avoid using term "BIAS_i" since
    it is easily confused with other things (the subsystem BIAS, bias currents),
    partly to include various stages of calibration.
BSO
    Singleton instance of class bicas.Settings which stores the settings for
    BICAS as a whole.
CA
    Cell Array.
CLI
    Command-line interface
CM3
    cm^-3
COHB
    CLI Option Header Body. The section of an CLI argument option/flag which
    excludes the prefix, e.g. "-" or "--". See bicas.utils.parse_CLI_options().
CTI
    CALIBRATION_TABLE_INDEX (L1R zVariable).
CTI2
    Element within a CDF L1R zVariable for a given CDF record,
    CALIBRATION_TABLE_INDEX(iCdfRecord, 2).
CWF
    Continuous WaveForm. Data on the form of samples over longer periods of
    time than a snapshot. Cf. SWF.
Dataset (data set)
    A CDF file on any one of a number standardized formats specified by the
    various RPW teams. All CDF files in the context of BICAS are datasets,
    except for RCTs, where it is ambiguous.
Deg
    Degrees (angle). 1 revolution=360 degrees=2*pi radians.
DLR
    Demultiplexer Latching Relay. Relay (true/false) that is part of the state
    of the demultiplexer.  See BIAS specification, section "3.4.4.14 MODE",
    "Data D3 = Diff probe 1&2(0), Diff probe 1&3(1)".
    Corresponds to BIAS HK zVariable "HK_BIA_MODE_DIFF_PROBE". The convention
    used for representing the value in BICAS is the same as in the BIAS HK ZV:
    0/false = V12
    1/true  = V13
    2023-08-29, EJ: Variable is constant for the entire mission except when it
    flips from zero to one at about 2023-08-21T12:04.
DSI, "dataset ID"
    DATASET_ID or "dataset ID". String constant which uniquely identifies a
    type of dataset.
    NOTE: DSI is defined by ROC, but the abbreviation "DSI" is defined by
    BICAS. DSI does not formally exist in SolO outside RPW, though the same
    string components (source+level+descriptor) can be used to define DSI
    analogously for SolO datasets outside RPW.
    NOTE: Must be uppercase by definition (at least when used as string
    constants by BICAS), though DSIs are represented in lowercase in the official
    dataset filenaming convention (except the "level"). Uppercase DSI has been
    used for some old ground-test datasets.
    NOTE: DSIs exclude the "-cdag" suffix. "-cdag" is defined by ROC
    for RPW-internal datasets only.
    NOTE: Variable names which refer to the GA "Dataset_ID" (and historical
    variants thereof) are not abbreviated "DSI", but to the GA name in
    analogy with other variables named after GAs.
DSR
    Downsampled/Decreased Sampling Rate. Cf. OSR.
DT
    MATLAB datetime object.
EMIDP
    (MATLAB) Error Message Identifier Part. One of the colon-separated
    parts of the MException .identifier string field (error message ID).
    NOTE: "Component" has a special meaning in the context of error
    message IDs. Therefore uses the term "part" instead.
FH
    Function Handle (Pointer).
FP
    Fill Position, in the context of FPAs.
FPA
    Class bicas.utils.FPArray.
FTF
    Forward Transfer Function = TF that describes the conversion of physical
    INPUT to OUTPUT (not the reverse). Cf. ITF.
FV
    (1) Field Value. (2) (CDF) Fill Value.
GA
    (CDF) Global Attribute(s).
GMDB
    "GA MODS DataBase", i.e. class bicas.ga.mods.Database.
GMDE
    "GA MODS DSI Entry", i.e. class bicas.ga.mods.DsiEntry.
GMVE
    "GA MODS Version Entry", i.e. class bicas.ga.mods.VersionEntry.
ICD
    Interface Control Document
ITF
    Inverse Transfer Function = TF that describes the conversion of physical
    OUTPUT to INPUT (not the reverse). Cf. FTF.
IV, ivolt, Interface Volt
    Calibrated volt at the interface between BIAS and LFR/TDS.
IVPAV
    Interface volt/antenna volt
IVPT
    Interface volt Per TM unit. IV/TM.
L2QBM
    zVariable L2_QUALITY_BITMASK.
L3QBM
    zVariable L3_QUALITY_BITMASK.
LL
    Log Level
LRX
    LFR Rx. Internal ZV-like variable (one boolean per CDF record) that
    determines whether there is data on BLTS 2-3 (LRX=1) or BLTS 4-5 (LRX=0).
    "Rx" is named after the similar LFR ZVs R0, R1 and R2 (one per LFR
    frequency, except F3). The value is constant for TDS.
LSF
    LFR Sampling Frequency (F0...F3).
    NOTE: When used as a variable (array index), 1=F0, ..., 4=F3.
LV
    Latest (Dataset) Version.
MC
    MATLAB Class. MATLAB data type represented as a string constant.
    Return value from "class()".
MSTD
    Modified STandard Deviation. Like standard deviation but using an
    arbitrary reference (e.g. median) instead of the conventional mean.
MVPM
    Milli-Volt Per Meter (mV/m).
NFP
    Non-Fill Position, in the context of FPAs.
NS
    Nanoseconds
NSO (Event)
    Non-Standard Operation. A time interval that has been manually labelled
    as special. BICAS can e.g. use this to set quality bits and remove data.
NSO table
    Separate XML file that contains a list of time intervals which have been
    manually labeled as NSOs. BICAS uses this list for potentially modifying
    processed data. Every NSO in the table is categorized using exactly one
    QRCID.
Offset
    Value (constant) that is ADDED to (not subtracted from) a measured
    value during the calibration process.
OSR
    Original Sampling Rate. Used in the context of downsampling. Cf DSR.
PFIID
    Production Function Input ID. String constant which is used to distinguish
    between different input values (datasets) to production functions. Different
    DSIs can be used for the same PFIID.
PFOID
    Production Function Output ID. Analogous with PFIID, but for output (return
    values).
Production Function
    Abstract method bicas.proc.SwmProcessing.production_function() which
    implementations in subclasses do the core of data processing.
QF
    zVariable QUALITY_FLAG.
QRC
    Quality-Related Condition. Condition that may influence quality ZVs.
    Ex: Thruster firing, saturation.
    QRCs are uniquely identified using QRCIDs.
    Note that QRCs are distinct from quality bits themselves. A QRC
    corresponds to setting zero, one, or multiple quality bits and/or cap
    QUALITY_FLAG.
QRCID
    QRC ID. String constant that represents a specific type of QRC.
    The NSO table uses QRCIDs to label NSO events.
RCS
    RPW Calibration Software. BICAS is an example of an RCS.
RCS ICD
    Originally document ROC-TST-GSE-ICD-00023-LES,
    "RPW Calibration Software ICD Documentation",
    which was later first complemented by and then later superseded by
    ROC-PRO-PIP-ICD-00037-LES,
    "RPW Calibration Software Interface Document".
    NOTE: For a period in the past, both documents were simultaneously valid
    and then referred to the ROC-SGSE and RODP pipeline separately(?).
RCT
    RPW Calibration Table. A CDF file with calibration data. See RCS ICD.
    ROC-defined term.
RCTD
    Class bicas.proc.L1L2.cal.rct.RctData.
RCTTID
    RCT Type ID. String constant that represents a *type* of RCT, not a
    particular RT file. Permitted values: "BIAS", "LFR", "TDS-CWF", "TDS-RSWF".
RCTS
    RCT CALIBRATION_TABLE (glob.attr) + CALIBRATION_TABLE_INDEX (zVariable).
    S = plural.
ROC
    RPW Operations Center.
ROC DFMD
    Document ROC-TST-GSE-NTT-00017-LES, "Data format and metadata definition for
    the ROC-SGSE data"
ROC Engineering Guidelines
    Document ROC-GEN-SYS-NTT-00019-LES, "ROC
    Engineering Guidelines for External Users"
RPS
    Radians Per Second
RSWF
    Regular Snapshot WafeForm. The term is used for TDS data, but is (at least
    for the purpose of BICAS) exactly the same as SWF.
RUM
    Document ROC-PRO-SFT-SUM-00080, "RCS User Manual"
RV
    Return Value.
sampere
    "Set current ampere". Simplified calibration value (in ampere) that is
    exactly proportional to bias current in TM.
SBDA
    Sweep BDM Detection Algorithm. BDM=4 is interpreted as equivalent to
    sweeping.
SCDA
    Sweep Current Detection Algorithm.
    A sliding window autodetection algorithm is used for detecting sweeps from
    measured HK currents. If the max-min current difference within that window
    exceeds a specified threshold, then the entire window is labeled as being
    part of a sweep.
SDID
    Class bicas.proc.L1L2.SignalDestinationId.
Sec
    Seconds
SIP
    Specific Input Parameter. BICAS CLI argument which specifies the path to
    an input or output dataset.
    The phrase (but not the abbreviation) is defined by the RCS ICD.
SKV
    bicas.SettingsKeyValue class.
SPR
    Samples (per channel) Per (CDF-like) Record. Only refers to actual data
    (currents, voltages), not metadata.
SRF
    Name of a particular s/c coordinate system. "Spacecraft Reference Frame"?
SSID
    Class bicas.proc.L1L2.SignalSourceId.
SWD, S/W descriptor
    Text on JSON format which describes among other things the S/W modes,
    including the required CLI parameters that every mode requires albeit not
    very clearly. (Defined by the RCS ICD.)
SWF
    Snapshot WaveForm. Snapshot data. Cf. CWF.
SWM, S/W mode
    (1) A "S/W mode" defines a set of required input CDF files and a set of
    output CDF files which can be derived from said input files. BICAS can
    execute only one such mode on each run. Executing such modes is the
    primary purpose of an RCS. (Defined by the RCS ICD.)
    (2) Class bicas.swm.SoftwareMode.
SWML
    Class bicas.swm.SoftwareModeList.
SWMP
    Class bicas.proc.SwmProcessing and subclasses thereof.
TBW
    To Bash Wrapper.
TF
    Transfer Function. Transfer function Z=Z(omega) in the frequency domain.
    Conventionally here on the form of a complex number (Z) as a function of
    radians/s (omega).
TM
    Telemetry units (in LFR/TDS ADC), or telecommand (TC) units. Using this term
    instead of the term "count".
TMK
    Class bicas.utils.Timekeeper.
TPIV
    TM/interface volt. TM per interface volt.
TS
    irfu-matlab's TSeries object/class.
TSF
    Threshold Saturation Flag. Bit for a sample that is set iff the sample
    value exceeds specified thresholds, i.e. it does not take any context
    into account, except for determining which thresholds to use. This is in
    contrast with (1) saturation flags for snapshots, which may summarize
    the saturation of the entire snapshot, or (2) saturation flags for CWF
    data which may take saturation before or after into account.
TTW
    TT2000 WOLS. Time format analogous to TT2000 but without leap seconds.
UFV
    Use Fill Values. Refers to CDF records for which science data should be
    overwritten with fill values. This refers to both bias currents and
    measured voltages.
VBMC
    Validate BICAS Master CDFs.
WOLS
    WithOut Leap Seconds
ZV
    CDF zVariable, or MATLAB variable that is analogous to one. First dimension
    corresponds to CDF record.
ZVS
    zVariable Struct. Struct where each field emulates the data content of a
    zVariable (either as an array or an FPA).



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
