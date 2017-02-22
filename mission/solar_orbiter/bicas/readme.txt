
#############
 About BICAS
#############
BICAS = BIAS Calibration Software

This software, BICAS, is created for the calibration of the hardware BIAS, a component in the RPW instrument on the Solar Orbiter spacecraft.
The principle author of this software is Erik P G Johansson, IRF-U, Uppsala, Sweden. Software development began 2016-03-xx (March 2016).

IMPORTANT NOTE: BICAS is designed to comply with the RCS ICD. Much documentation on how to use this software can thus be found there.
The documentation included in BICAS only covers the most important parts of the RCS ICD.



#######################
 Acronyms / Dictionary
#######################
NOTE: This list also applies to comments in the source code.
RCS                        = RPW Calibration Software. BICAS is an example of an RCS.
ICD                        = Interface Control Document
ROC DFMD                   = Document ROC-TST-GSE-NTT-00017-LES, "Data format and metadata definition for the ROC-SGSE data"
ROC Engineering Guidelines = Document ROC-OPS-PIP-NTT-00008-LES, "RPW Ground Segment - ROC Engineering Guidelines"
RCS ICD                    = Documents ROC-TST-GSE-ICD-00023-LES, "RPW Calibration Software ICD Documentation",
                             and       ROC-PRO-PIP-ICD-00037-LES, "RPW Calibration Software Interface Document"
                             which at this stage (2017-02-20) appear to be identical in their specification.
BIAS specification         = Document RPW-SYS-MEB-BIA-SPC-00001-IRF, "RPW Instrument -- BIAS Specification"
S/W mode                   = A "S/W mode" defines a set of required input CDF files and a set of output CDF files derived from the input files.
                             BICAS can execute only one such mode on each run. Executing such modes is the primary purpose of an RCS.
                             (Defined by the RCS ICD.)
S/W descriptor             = Text on JSON format which describes among other things the S/W modes, including the required CLI parameters that
                             every mode requires albeit not very clearly. (Defined by the RCS ICD.)
CLI                        = Command-line interface
dataset (data set)         = A CDF file on one of a number standardized formats specified by the various RPW teams.
                             All CDF files in the context of BICAS are datasets.



#################
 Main executable
#################
<BICAS root dir>/roc/bicas
(from Linux command-line)



#############################
 CLI syntax / CLI parameters
#############################
NOTE: The official CLI parameter syntax is defined in RCS ICD, Iss02 Rev02, Section 3.2.

SYNTAX 1: ( --identification | --version | --help ) <General parameters>
SYNTAX 2: <S/W mode> <General parameters, Output parameter, Specific inputs parameters> 

NOTE: Only the position of the first parameter is important. The order of all other parameters is arbitrary.

--version           Prints the software version.
--identification    Prints the S/W descriptor.
--help


====================
 General parameters
====================
--log    <absolute path to directory>   (optional) Specifies directory in which to put log files.
--config <absolute path to file>        (optional) Specifies the configuration file to use.

<S/W mode>   Selects the S/W mode to use.
Available S/W modes can be found in the S/W descriptor. They are listed under "modes". "name" specifies the string that identifies a given mode and
can be used on as CLI argument.



===========================
 Specific input parameters
===========================
Set of parameters which specify input CDF files. The exact set depends on the exact S/W mode and can in principle be read from the S/W descriptor.
Required input parameters for a specific S/W mode can be found in the S/W descriptor under
"modes" --> (specific mode) --> "inputs" --> Name of subsection, e.g. "input_hk", "input_sci".
Example: an input subsection "input_hk" means that there is a required parameter "--input_hk <path_to_file>".



==================
 Output parameter
==================
   --log    <absolute path to directory>   Specifies the directory in which to put output CDF files.



===============
 Example calls
===============
bicas --version --config ~/bicas.conf
bicas --identification
bicas --identification --log ~/logs
bicas LFR-SURV-CWF-E   --input_hk  ROC-SGSE_HK_RPW-BIA_7729147_CNE_V01.cdf
                       --input_sci ROC-SGSE_L2R_RPW-LFR-SURV-CWF_7729147_CNE_V01.cdf
                       --output ~/temp   --log ~/logs/   --config ~/bicas.conf



###########################################
 Installation, set-up, system requirements
###########################################
See "install.txt".
   
   
   
########################################
 S/W descriptor, roc_sw_descriptor.json
########################################
The S/W descriptor and the file "roc_sw_descriptor.json" are required and specified by the RCS ICD.
The file "roc_sw_descriptor.json" should contain the S/W descriptor. The file is NOT meant to be edited manually.
If BICAS has been updated in such a way that the software descriptor changes, then a new "roc_sw_descriptor.json"
can be generated by piping the stdout from the call "bicas --identification" to the corresponding file.



#######################
 Source code directory
#######################
<BICAS root dir>/src
<BICAS root dir>/src/+bicas/
<BICAS root dir>/src/+bicas/+tools
    Code that is not called by BICAS itself, but which is useful for development of BICAS, and for constructing master CDF files.    
<BICAS root dir>/src/+bicas/+utils
    Generic utilities used by BICAS. Could possibly be moved out of BICAS to a more generic separate library of code at some point.



####################################
 Known current limitations, caveats
####################################
This section lists known current limitations and caveats. Most of them are likely less important.
Most or all of these limitations represent features that should be implemented eventually.
Within each category, entries are sorted from most important to least important.

List last updated 2016-11-17.

* Reading CDF files:
    - BICAS might not check explicitly if an input CDF file (depends on setting), or master CDF file appears to be of the right type of dataset
      (dataset ID) before trying to process the information. NOTE: Many input datasets (i.e. NOT generated by BICAS) have not set the CDF
      Global Attributes needed to do this check.

* Processing:
    - Voltage calibration only uses proportionality factors, i.e.:
        - Transfer functions are not used.
        - No voltage offsets due to bias currents, MUX mode, diff gain etc. are added.
    - zVariables IBIAS1-IBIAS3 are set to fill values.
    - Not all S/W modes described in the S/W descriptor necessarily work in practice, due to ongoing development.
    - LFR: Algorithms can use ACQUISITION_TIME for time rather than Epoch (for comparing HK data time stamps to SCI timestamps). Depends on setting.
    - LFR: zVariables QUALITY_BITMASK, QUALITY_FLAG are set to fill values with the right size when the source dataset has empty zVariables
        (instead of copying).
    - Mux/demux modes 5-7 (calibration modes 1-2) have not been implemented.
    - The conversions between tt2000 (Epoch), ACQUISITION_TIME and UTC have not been verified.
    - Many combinations of input dataset CDF zVariables (settings) are untested even if they should work in theory.

* Writing CDF files:
    - Output CDF files only use fill values correctly for floating-point numeric zVariables.
    - Global attributes:
        - Generation_date : Date set in local time (should probably be set in e.g. UTC).
        - Some attributes which should be set (not inherited from the master CDF files) are not set,
            e.g. SPECTRAL_RANGE_MIN, SPECTRAL_RANGE_MAX, TIME_MIN, TIME_MAX.
            NOTE: Failure to set global attributes can also be due to missing global attribute values in the parent datasets.
        - Some global attributes may in principle be wrong if using a S/W mode with multiple output datasets
          (sets some global attributes to the same value for all datasets).
    - CDF Epoch's zVariable attributes SCALEMIN, SCALEMAX are not set correctly.
    (- BICAS can not produce CDF files with zVariables with zero records.)

NOTE: Input CDF files (datasets) may have incorrect, or partially incorrect Epoch values, which translate into faulty Epoch values in BICAS'
output CDF files. This is due to known limitations in LESIA's (?) early production of test datasets (and BICAS can not compensate for this).
LESIA prescribes that datasets contain an extra time zVariable ACQUISITION_TIME that is used in parallel to Epoch as a temporary solution.
New LESIA CDF files should contain the correct Epoch since ca July-August 2016.

