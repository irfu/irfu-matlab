% IRFU Matlab routines
%
% Time operations
%   toepoch         - convert a [YYYY MM DD hh mm ss] time specification to seconds since 1970
%   fromepoch       - Convert seconds since 1970 to [YYYY MM DD hh mm dd] time format
%   date2epoch      - converts MATLAB datenum to ISDAT epoch
%   epoch2date      - converts ISDAT epoch to MATLAB datenum
%   iso2epoch       - convert ISO time string to ISDAT epoch
%   epoch2iso       - Convert ISDAT epoch to ISO time string
%   epoch2yyyymmdd  - Convert ISDAT epoch to YYYYMMDD
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
%   irf_resamp      - Resample time series to the time line of another times series
%   irf_tappl       - apply expression to data
%   irf_tlim        - find a subinterval defined by time limits
%
% Space physics specific routines
%   irf_cdf_read    - read cdf files (see also caa_load for CAA cdf files)
%   irf_disp_surf   - Interactively plot cold plasma dispersion surfaces
%   irf_e_vxb       - Compute VxB and ExB/B^2
%   irf_gse2gsm     - GSE <> GSM
%   irf_jz          - Estimate the current given velocity and magnetic field
%   irf_mean        - put time series into mean field coordinates
%   irf_plasma_calc - Calculate basic plasma quantities
%   irf_shue_mp     - estimate distance to model magnetopause (Shue)
%   irf_units       - natural constants and units 
%   irf_vht         - estimate velocity of the deHoffman-Teller frame
%   lp_...          - Langmuir probe specific routines 
%
% Plotting
%   irf_plot        - Flexible plotting routine for time series
%   dataobj/plot    - plot data objects. help dataobj/plot.m
%   irf_pl_mark     - Mark that time interval(s) with color backround
%   irf_pl_number_subplots - number subplots by putting in each subplot 'a', 'b', 'c', etc.
%   irf_tm          - time manager for plots created by irf_plot
%   irf_zoom        - zoom in
%   irf_figmenu     - Add to the current figure a menu with some useful commands
%
% Cluster
%   c_4_v_gui       - timing analysis 
%   c_get_batch     - get Cluster data from ISDAT/CDF/DDS
%   c_load          - load data from MAT files (caa_load for CAA cdf files)
%   c_pl_summary    - Make EFW summary plots
%
%   c_4_grad        - calculate linear gradient operators using 4 spacecraft technique
%   c_4_j           - Calculate current from using 4 spacecraft technique
%   c_efw_scp2ne    - calculates the plasma density (Ne) for given EFW sc potential
%   c_eval          - evaluate expression for list of spacecraft
%   c_gse2dsi       - Convert vector from GSE to DSI reference system
%   c_pl_flux_tube_distance - estimate the distance between flux tubes
%   c_pl_sc_conf_lmn - Plot the configuration of Cluster in LMN (or arbitrary) coordinates  
%   c_pl_sc_conf_xyz - Plot the configuration of Cluster in XYZ (GSE or GSM) coordinates   
%   c_pl_sc_orient  - Plots the orientation of the EFW probes
%   c_pl_sc_pos_mf  - Plot spacecraft position in mean field coordinates
%   c_peace_plot .. - plot PEACE data (also reading and spectra routines)
%   c_staff_getsa   - read STAFF data
%   c_wbd_read      - read Wideband data
%   caa_load        - load data downloaded from CAA in CDF format
%   
% Miscellaneous
%   irf_ssub         - substitute ? ! and $ in string to number or variable
%   irf_whamp_plot_f - plot the distribution function, parameters as defined in WHAMP

%
% $Id$