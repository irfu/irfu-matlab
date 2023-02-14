Directory for non-MATLAB files related to software "BIA_QL3" for generating official BIAS quicklooks/summary plots.
The software is delivered to ROC/LESIA for execution there.


USAGE
=====
Main bash script: bia_ql3
Syntax 1: --help
Syntax 2: matlab_executable YYYYMMDD input_bia_path input_lfr_wf_e_path output_dir log_dir


NOTES
=====
According to ROC (Milan Maksimovic, LESIA):
    "Summary Plot" (SP) = Plot which summarizes the content of a single dataset
    "Quicklook"         = Plot which summarizes data for entire instrument.
Hence this code generates what he calls "summary plots" and the MATLAB package is called "solo.sp".
However, ROC has named the software "BIA_QL3" (QL=Quicklooks? L3?!) so the usage is inconsistent.
--
Since the software is run by ROC/LESIA, the code has to conform to the interface which they specify.
Not a formal description, but see e.g.
https://gitlab.obspm.fr/ROC/RCS/BIA_QL3/-/issues/1


IRFU-MATLAB DIRECTORIES ASSOCIATED WITH BIA_QL3
===============================================
irfu-matlab/mission/solar_orbiter/bia_ql3/
    Non-MATLAB files
irfu-matlab/mission/solar_orbiter/+solo/+sp/
    MATLAB code
irfu-matlab/mission/solar_orbiter/+solo/+sp/bia_ql3.m
    Main MATLAB function.


OFFICIAL DELIVERIES TO ROC/LESIA
================================
irfu-matlab is the basis for separate git repo
    https://gitlab.obspm.fr/ROC/RCS/BIA_QL3
which is used for official deliveries to ROC/LESIA.


VERSION HISTORY
===============
v1.2.1: Bugfix
v1.2.2: No longer generate .matlab_datastore_* files.
