% IRFU Matlab routines
%
% Time operations
%   irf_time        - Convert between time formats in epoch/vector/iso/others
%   irf_fname       - construct a filename string from time
%
% Time series operations
%   irf_abs         - estimate abs() 
%   irf_add         - add two vectors
%   irf_ang         - compute angle between vectors
%   irf_cross       - cross product
%   irf_dot         - Calculate dot product between vectors 
%   irf_filt        - filter time series
%   irf_integrate   - time integral
%   irf_minvar ...  - minimum variance analysis routines
%   irf_newxyz      - rotate vector into new coordinate system
%   irf_norm        - Normalize vector
%   irf_powerfft    - Calculate spectrograms
%   irf_resamp      - Resample time series to the time line of another times series
%   irf_tappl       - apply expression to data
%   irf_tlim        - find a subinterval defined by time limits
%   irf_wavelet     - fast wavelet spectrogram calculation
%
% Space/plasma physics specific routines
%   irf_cdf_read    - read cdf files (see also caa_load for CAA cdf files)
%   irf_disp_surf   - Interactively plot cold plasma dispersion surfaces
%   irf_e_vxb       - Compute VxB and ExB/B^2
%   irf_get_data    - Get different kind of data
%   irf.geocentric... - geocentric coordinate transformation (Hapgood 1997)
%   irf_gse2gsm     - GSE <> GSM
%   irf_jz          - Estimate the current given velocity and magnetic field
%   irf_mean        - put time series into mean field coordinates
%   irf_plasma_calc - Calculate basic plasma quantities
%   irf_plasma_wave_visualization.m - Visualize plasma waves
%   irf_shue_mp     - estimate distance to model magnetopause (Shue)
%   irf_units       - natural constants and units 
%   irf_vht         - estimate velocity of the deHoffman-Teller frame
%   lp.             - package with Langmuir probe specific routines 
%
% Plotting
%   irf_plot        - Flexible plotting routine for time series & CAA variables
%   dataobj/plot    - plot data objects. help dataobj/plot.m
%   irf_pl_mark     - Mark that time interval(s) with color backround
%   irf_pl_number_subplots - number subplots by putting in each subplot 'a', 'b', 'c', etc.
%   irf_plot_axis_align  - align x axis of different subplots
%   irf_plot_ylabels_align - left align ylabels of different subplots 
%   irf_tm          - time manager for plots created by irf_plot
%   irf_zoom        - zoom in
%   irf_figmenu     - Add to the current figure a menu with some useful commands
%   irf_timeaxis    - Add timeaxis in different formats
%
% Cluster
%   caa_download    - download data from CAA
%   caa_load        - load data downloaded from CAA 
%   caa_meta        - find metadata information for datasets
%   c_caa_var_get   - get variable from CAA files
%   c_load          - load data from MAT files (caa_load for CAA cdf files)
%   c_get_batch     - get Cluster data from ISDAT/CDF/DDS
%   c_pl_summary    - Make EFW summary plots
%
%   c_4...          - different 4 spacecraft methods
%   c_4_grad        - calculate linear gradient operators using 4 spacecraft technique
%   c_4_j           - Calculate current from using 4 spacecraft technique
%   c_4_v_gui       - timing analysis from 4 or less spacecraft 
%   c_coord_trans   - Convert vectors among GSE/GSM/DSC/DSI/ISR2 reference system
%   c_fgm_staff_combine - combine FGM and STAFF data into one timeseries
%   c_efw_scp2ne    - calculates the plasma density (Ne) for given EFW sc potential
%   c_eval          - evaluate expression for list of spacecraft
%   c_pl_flux_tube_distance - estimate the distance between flux tubes
%   c_pl_sc_conf_xyz - Plot the configuration of Cluster in GSE/GSM    
%   c_pl_sc_orient  - Plots the orientation of the Cluster, including EFW probes
%   c_pl_sc_pos_mf  - Plot spacecraft position in mean field coordinates
%   c_staff_getsa   - read STAFF data
%   c_wbd_read      - read Wideband data
%   
% Miscellaneous
%   irf_ssub         - substitute ? ! and $ in string to number or variable
%   whamp.           - different routines related to WHAMP
%   irfnotes         - notes on using irfu-matlab
%
% Bug reports, feature requests - https://github.com/irfu/irfu-matlab
