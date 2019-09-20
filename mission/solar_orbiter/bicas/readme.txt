
#############
 About BICAS
#############
BICAS = BIAS Calibration Software

This software, BICAS, is created for the calibration of the BIAS subsystem in the RPW instrument on the Solar Orbiter spacecraft.
The principle author of this software is Erik P G Johansson, IRF-U, Uppsala, Sweden. Software development began 2016-03-xx (March 2016).

IMPORTANT NOTE: BICAS is designed to comply with the RCS ICD. Much documentation on how to use this software can thus be found there.
Fore more documentation, see RCS ICD and RUM documents (see below).



#######################
 Acronyms / Dictionary
#######################
NOTE: This list also applies to comments in the source code.
--
RCS                        = RPW Calibration Software. BICAS is an example of an RCS.
ICD                        = Interface Control Document
ROC DFMD                   = Document ROC-TST-GSE-NTT-00017-LES, "Data format and metadata definition for the ROC-SGSE data"
ROC Engineering Guidelines = Document ROC-OPS-PIP-NTT-00008-LES, "RPW Ground Segment - ROC Engineering Guidelines"
RCS ICD                    = Originally document           ROC-TST-GSE-ICD-00023-LES, "RPW Calibration Software ICD Documentation",
                             which was later superseded by ROC-PRO-PIP-ICD-00037-LES, "RPW Calibration Software Interface Document".
                             NOTE: "RCS ICD" does not at this time (2019-07-24) distinguish between these two which gives room for confusion since
                                   a later rev/iss for the old RCS ICD may thus be superseded by a lower rev/iss for the newer RCS ICD.                             
RUM                        = Document ROC-PRO-SFT-SUM-00080, "RCS User Manual"
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



=========================
 Common input parameters
=========================
--log    <absolute path to log file>    (optional) Specifies log file.
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



===============
 Example calls
===============
bicas --version --config ~/bicas.conf
bicas --identification
bicas --identification --log ~/bicas.log
bicas LFR-SURV-CWF-E   --in_hk   ROC-SGSE_HK_RPW-BIA_7729147_CNE_V01.cdf
                       --in_sci  ROC-SGSE_L1R_RPW-LFR-SURV-CWF-E_7729147_CNE_V01.cdf
                       --out_sci ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_7729147_CNE_V01.cdf   --log ~/bicas.log   --config ~/bicas.conf



###########################################
 Installation, set-up, system requirements
###########################################
See "install.txt".



####################################
 Known current limitations, caveats
####################################
* BICAS DOES NOT YET CALIBRATE DATA. This should not apply to offically delivered versions.
  It does however read and write datasets (CDF files) in the s/w modes that it supports. It is just the output dataset content that is not entirely complete.

For other limitations and caveats, see the official user manual, the RUM document.

