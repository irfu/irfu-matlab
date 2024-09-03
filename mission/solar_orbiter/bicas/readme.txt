
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
