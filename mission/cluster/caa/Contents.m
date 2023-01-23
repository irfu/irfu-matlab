% Caa specific routines
%
% Data manipulation
%   caa_append_data       - concatenate two datasets
%   caa_fill_gaps         - fill gaps in the of a dataset
%   caa_filter_e          - simple high pass filter EFW E
%   caa_join_phase        - refetch missing phase
%   caa_powerfft          - compute power spectrum
%   caa_rm_blankt         - remove time intervals from data (sets to NaN)
%   remove_problems       - script REMOVE_PROBLEMS
%
% Planning
%   caa_find_mp           - find model magnetopause crossings using ACE data
%   caa_control_mp        - plot predicted MP location
%   caa_sh_plan           - identify magnetosheath/sw intervals
%
% Offsets
%   caa_corof_adc         - correct the ADC offset in data by removing average
%   caa_corof_delta       - correct delta offsets
%   caa_deltaoff_batch    - invert previously applied delta offsets
%   caa_corof_dsi         - correct offsets in DSI(ISR2)
%   caa_sh_xoff           - sunward offset and ampl correction in the sh/sw
%   caa_sh_xoff_batch     - sunward offset and ampl correction in the sh/sw
%   caa_update_deltaoff   - update delta offsets database
%
% Plotting tools
%   caa_ib_sp_batch       - produce summary plots for internal burst
%   caa_ms_pl_xoff        - visualize offset in the magnetosphere
%   caa_ms_pl_xoff_stats  - plot statistics for DdsiX1-4
%   caa_sh_pl_xoff        - visualize offset study results
%   caa_pl_efw_pea_hia    - compare E from EFW, PEACE and CIS_HIA
%   caa_pl_coverage       - plot on/off and coverage information
%   caa_pl_summary        - make summary plots of the CAA data
%   caa_pl_summary_l1     - CAA summary plot for L1 & L2 P data & EF
%   caa_spectrogram       - plot power spectrum in logarithimic scale
%   caa_comp_efw_edi_corr - compare EFW, EDI with CORROTATION
%
% Reading data from different sources
%   caa_get               - read data from caa Matlab files
%   caa_is_get            - get Cluster data from ISDAT database
%   caa_load              - load data downloaded from the CAA in CDF
%   caa_get_edi           - get EDI data
%   caa_sh_get_cis        - get CIS data
%
% CEF export tools
%   caa_get_batch_l0      - CAA get L0 data by CaaRunner
%   caa_export            - export data to CAA CEF files
%   caa_export_batch      - run CAA_EXPORT in a batch script
%   caa_export_batch_l1   - run CAA_EXPORT in a batch script
%   caa_reexport          - ReExport CEF file
%
% Batch processing tools
%   caa_ms_reproc         - Reprocess magnetospheric data
%   caa_reproc            - Reprocess CAA data
%   caa_rerun_l3e         - rerun spinfits
%   caa_sh_reproc         - Reprocess magnetosheath/sw data
%   caa_sp_batch          - produce summary plots
%
% AUX data
%   caa_efw_mode_tab      - read EFW FDM
%   caa_get_ns_ops_int    - get interval for a specific problem in NS_OPS
%   caa_ns_ops_int        - split/truncate interval according to NS_OPS
%   caa_quality           - set quality factor for the CAA products
%   caa_set_bitmask_and_quality  - sets the bitmask and quality factor
%   caa_identify_problems - identify problem areas in data,sets bitmask
%   caa_cefname2specs     - extract details from a CEF file name
%   caa_errid2str         - convert ns_ops error ID to string
%   caa_str2errid         - convert ns_ops error ID string to number
%
% Directory tools
%   caa_is_sh_interval    - check if directory contains an SH interval
%   caa_ms2sh_int         - change MS interval to SH
%   caa_is_valid_dirname  - check if directory was created by the CAA s/w
%   caa_read_interval     - read caa interval information from .interval
%   caa_read_version      - read caa file version
%   caa_sfit_probe        - set/get probe pair for spin resolution data
%   caa_sh_3h_int         - create common 3h intervals
%
% Old programs
%   caa_deliver           - deliver CEF files to the CAA
%   caa_make_jobs         - make job files
%   caa_read_coverage     - read on/off and coverage files
%   caa_reg_skel          - make skeleton XML file from coverage data
%   caa_testsfit          - test script for C_EFW_C_EFW_SPINFIT_MX
%
