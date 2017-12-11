
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
* Not all S/W modes described in the S/W descriptor necessarily work in practice in all versions (commits) of the irfu-matlab git repository. This should not apply to offically delivered versions.


NOTE: Input CDF files (datasets) may have incorrect, or partially incorrect Epoch values, which translate into faulty Epoch values in BICAS' output CDF files. This is due to known limitations in LESIA's (?) early production of test datasets (and BICAS can not compensate for this). LESIA prescribes that datasets contain an extra time zVariable ACQUISITION_TIME that is used in parallel to Epoch as a temporary solution. New LESIA CDF files should contain the correct Epoch since ca July-August 2016.

For other limitations and caveats, see the official user manual (being written as of 2017-03-03),
ROC-PRO-SFT-SUM-xxxxx-LES, "BIAS Calibration Software (BICAS) User Manual"

