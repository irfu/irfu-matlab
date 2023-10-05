
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
The code should however converge on the following.
- Use the defined abbreviations in identifiers.
- Variable names should use camelCase. Initial in lowercase.
- Classes:
  - Class names for classes which are actually instantiated (which are not
    merely a collection of static methods) should be CamelCase. Initial in
    uppercase.
  - Classes which are not instantiated should be named using lowercase without
    underscore.
- Functions should be named using snake_case. Abbreviations in uppercase.



###########################
 Abbreviations, dictionary
###########################
NOTE: This list also applies to comments and identifiers in the source code.
--
AA, aampere
    Antenna ampere. Calibrated ampere at the antenna.
AAPT
    Antenna ampere/TM
ASID
    bicas.proc.L1L2.AntennaSignalId
ASR, Antenna Signal Representation.
    (1) A reference to a particular "physical antenna signal"
        (a "data channel") which BIAS-LFR/TDS is trying to measure, or
    (2) a measurement of such.
    In reality, the terminology is:
    ASR         : Pointer to a specific physical antenna signal, e.g. V12_LF
                  (DC diff, antenna 1-2)
    ASR samples : Samples representing a specific ASR (as opposed to BLTS).
    NOTE: There are 9 ASRs (DC: 3 singles, 3 diffs; AC: 3 diffs), i.e. they
    can refer also to signals not represented by any single BLTS, given a
    chosen mux mode and latching relay setting.
AV, avolt, Antenna Volt
    Calibrated volt at the antennas, i.e. the final calibrated (measured)
    value, including for reconstructed signals (e.g. diffs calculated from
    singles). May also refer to offsets and values without offsets.
AVPIV
    Antenna Volt / (Per) Interface Volt
BDM, "mux mode"
    BIAS Demultiplexer Mode. Corresponds to BIAS HK ZV "HK_BIA_MODE_MUX_SET".
BIAS specification
    Document RPW-SYS-MEB-BIA-SPC-00001-IRF, "RPW Instrument -- BIAS
    Specification".
BIAS_1, ..., BIAS_5 (BIAS_i, i=1..5)
    Defined in BIAS specifications document. Equal to the physical signal at
    the physical boundary between BIAS and LFR/TDS. Unit: Interface volt.
    Mostly replaced by BLTS+specified unit in the implementation.
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
CA
    Cell Array.
CLI
    Command-line interface
CTI
    CALIBRATION_TABLE_INDEX (zVariable).
CTI2
    Second value in a CDF record of zVariable CALIBRATION_TABLE_INDEX.
Dataset (data set)
    A CDF file on any one of a number standardized formats specified by the
    various RPW teams. All CDF files in the context of BICAS are datasets.
Deg
    Degrees (angle). 1 revolution=360 degrees=2*pi radians.
DLR
    Demultiplexer Latching Relay. Relay (true/false) that is part of the state
    of the demultiplexer.  See BIAS specification, section "3.4.4.14 MODE",
    "Data D3 = Diff probe 1&2(0), Diff probe 1&3(1)".
    Corresponds to BIAS HK zVariable "HK_BIA_MODE_DIFF_PROBE".
    The convention used for representing the value is the same as in the ZV.
    0/false = V12
    1/true  = V13
    2023-08-29, EJ: Variable is constant for the entire mission except when it
    flips from zero to one at about 2023-08-21T12:04.
DSI
    DATASET_ID. String constant that uniquely identifies each type of dataset.
    NOTE: Variable names which refer to the GA "Dataset_ID" (and historical
    variants thereof) are not abbreviated "DSI", but to the GA name, in
    analogy with other variables named after GAs.
DSR
    Downsampled/Decreased Sampling Rate. Cf. OSR.
EMIDP
    (MATLAB) Error Message Identifier Part. One of the colon-separated
    parts of the MException .identifier string field (error message ID).
    NOTE: "Component" has a special meaning in the context of error
    message IDs. Therefore uses the term "part" instead.
FPA
    Class bicas.utils.FillPositionsArray
FTF
    Forward Transfer Function = TF that describes the conversion of physical
    INPUT to OUTPUT (not the reverse). Cf. ITF.
FV
    (1) Field Value. (2) (CDF) Fill Value.
GA
    (CDF) Global Attribute(s).
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
    zVar L2_QUALITY_BITMASK.
LSF
    LFR Sampling Frequency (F0...F3).
    NOTE: When used as a variable (array index), 1=F0, ..., 4=F3.
MSTD
    Modified STandard Deviation. Like standard deviation but using an
    arbitrary reference (e.g. median) instead of the conventional mean.
MVPM
    Milli-Volt Per Meter (mV/m).
NS
    Nanoseconds
NSO
    Non-Standard Operations. Functionality for making BICAS modify processed
    data based on manually compiled list of "events". Can e.g. set quality
    bits and remove data.
NSOID
    NSO ID. String constant that represents a specific type of NSO. Used in
    NSO table.
Offset
    Value (constant) that is ADDED to (not subtracted from) a measured
    value during the calibration process.
OSR
    Original sampling (rate). Used in the context of downsampling. Cf DSR.
RCS
    RPW Calibration Software. BICAS is an example of an RCS.
RCS ICD
    Originally document ROC-TST-GSE-ICD-00023-LES,
    "RPW Calibration Software ICD Documentation",
    which was later superseded by ROC-PRO-PIP-ICD-00037-LES,
    "RPW Calibration Software Interface Document".
    NOTE: "RCS ICD" does not at this time (2019-07-24) distinguish between
    these two which gives room for confusion since a later rev/iss for the old
    RCS ICD may thus be superseded by a lower rev/iss for the newer RCS ICD.
RCT
    RPW Calibration Table. CDF with calibration data. See RCS ICD. ROC-defined.
RCTID
    RCT ID. String constant that represents a type of RCT: "BIAS", "LFR",
    "TDS-CWF", "TDS-RSWF".
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
RUM
    Document ROC-PRO-SFT-SUM-00080, "RCS User Manual"
RV
    Return Value.
sampere
    "Set current ampere". Simplified calibration value (in ampere) that is
    exactly proportional to bias current in TM.
SDID
    bicas.proc.L1L2.SignalDestinationId
Sec
    Seconds
SPR
    Samples (per channel) Per (CDF-like) Record. Only refers to actual data
    (currents, voltages), not metadata.
SRF
    Name of a particular s/c coordinate system. "Spacecraft Reference Frame"?
SSID
    bicas.proc.L1L2.SignalSourceId
SWD, S/W descriptor
    Text on JSON format which describes among other things the S/W modes,
    including the required CLI parameters that every mode requires albeit not
    very clearly. (Defined by the RCS ICD.)
SWM, S/W mode
    (1) A "S/W mode" defines a set of required input CDF files and a set
    of output CDF files derived from the input files. BICAS can execute only
    one such mode on each run. Executing such modes is the primary purpose of
    an RCS. (Defined by the RCS ICD.)
    (2) Class bicas.swm.SoftwareMode.
SWML
    Class bicas.swm.SoftwareModeList.
TBW
    To Bash Wrapper.
TF
    Transfer Function. Transfer function Z=Z(omega) in the frequency domain.
    Conventionally here on the form of a complex number (Z) as a function of
    radians/s (omega).
TM
    Telemetry units (in LFR/TDS ADC), or telecommand (TC) units. Using this term
    instead of the term "count".
TPIV
    TM/interface volt. TM per interface volt.
TS
    irfu-matlab's TSeries object/class.
TTW
    TT2000 WOLS. Time format analogous to TT2000 but without leap seconds.
UFV
    Use Fill Values. Refers to CDF records which data should overwritten with
    fill values).
WOLS
    WithOut Leap Seconds
ZV
    CDF zVariable, or MATLAB variable that is analogous to one. First dimension
    corresponds to CDF record.



####################
 Naming conventions
####################
Variables are camelCase.
Structs and objects use an uppercase initial, and other variables use a
lowercase initial.
Exception: Variables which are direct analogues to zVariables are named as the
corresponding zVariables, i.e. SCREAMING_SNAKE_CASE most of the time.

"b" refers to logical values (boolean) and logical indexing
"i", "j" refers to indices into arrays.

Unit tests (classes) have suffix "___UTEST".



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
