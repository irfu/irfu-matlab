%
% PSP_LOAD loads PSP data into TSeries
%
% PSP_LOAD(dataDir,datatype,dateStart,dateStop)
%
% PSP_LOAD([],datatype,dateStart,dateStop) will use PSP data directory saved
% by datastore. If not defined will ask the first time.
%
% PSP_LOAD(datatype) uses saved PSP directory and global variables
% dateStart, dateStop. Data are downloaded if missing from the directory.
%
% PSP_LOAD(cdfFile,datatype) return variables into MATLAB base
%
% OUT = PSP_LOAD(..) output is TSeries if single variables is returned and
% cell array if several variables have to be returned.
%
% [OUT,DOBJ] = PSP_LOAD(..) return also the last DATABOJ
%
% dataDir:  the directory with the cdf files
% datatype: 'mag_rtn','fgm_rtn'      - FGM, RTN coordinates
%           'mag_sc','fgm_sc'        - FGM, SC coordinates
%           'wf_dvdc'                - DBF Digital Fields Board Differential Voltage Waveform
%           'wf_vdc'                 - DBF Digital Fields Board Differential Voltage Waveform
%           'wf_scm_sc'              - DBF Digital Fields Board Search Coil Magnetometer Waveform SC coordinates (usually has no valid data, only NaN's)
%           'wf_scm_sensor'          - DBF Digital Fields Board Search Coil Magnetometer Waveform sensor coordinates
%           'ac_bpf_dv34'            - DFB AC-coupled Differential Voltage, V3-V4 Antennae, Bandpass filtered, High Gain
%           'ac_bpf_SCMulfhg'        - DFB AC-coupled SCM, Bandpass filtered, Low Frequency, High Gain, u-component
%           'ac_bpf_SCMumfhg'        - DFB AC-coupled SCM, Bandpass filtered, Medium Frequency, High Gain, u-component
%           'ac_spec_dv12'           - DFB AC-coupled Differential Voltage, V1-V2 Antennae, Dipole Mode, Spectra, High Gain
%           'ac_spec_dv34'           - DFB AC-coupled Differential Voltage, V3-V4 Antennae, Dipole Mode, Spectra, High Gain
%           'ac_spec_v5'             - DFB AC-coupled V5 Antenna Voltage, Monopole Mode, Auto Spectra, High Gain
%           'ac_xspec'               - DFB AC-coupled Diff Voltage, V1-V2, V3-V4 Antennae, Dipole Mode, X-Spectra
%           'ac_spec_SCMv'           - DFB AC-coupled SCM, Low Frequency, High Gain, v-component, Sensor coordinates
%           'ac_spec_SCMu'           - DFB AC-coupled SCM, Low Frequency, High Gain, u-component, Sensor coordinates
%           'ac_spec_SCMd'           - DFB AC-coupled SCM, Low Frequency, High Gain, d-component, Sensor coordinates only Enc 1 (2018/11)
%           'ac_spec_SCMe'           - DFB AC-coupled SCM, Low Frequency, High Gain, e-component, Sensor coordinates only Enc 1 (2018/11)
%           'ac_spec_SCMf'           - DFB AC-coupled SCM, Low Frequency, High Gain, f-component, Sensor coordinates only Enc 1 (2018/11)
%           'dc_xspec_SCMvw'         - DFB DC-coupled Searchcoil Magnetometer, SCM, Cross Spectra, v-component, w-component, High Gain, Sensor coordinates
%           'dc_xspec_SCMde'         - DFB DC-coupled Searchcoil Magnetometer, SCM, Cross Spectra, d-component, e-component, High Gain, Sensor coordinates only Enc 1 (2018/11)
%           'dc_xspec_SCMdf'         - DFB DC-coupled Searchcoil Magnetometer, SCM, Cross Spectra, d-component, f-component, High Gain, Sensor coordinates only Enc 1 (2018/11)
%           'dc_xspec_SCMef'         - DFB DC-coupled Searchcoil Magnetometer, SCM, Cross Spectra, e-component, f-component, High Gain, Sensor coordinates only Enc 1 (2018/11)
%           'dc_bpf_dv34'            - DFB DC-coupled Diff Volt V34 antenas Bandpass Filtered (BPF)
%           'dc_bpf_SCMulfhg'        - DFB DC-coupled SCM, Bandpass filtered, Low Frequency, High Gain, u-component
%           'dc_bpf_SCMvlfhg'        - DFB DC-coupled SCM, Bandpass filtered, Low Frequency, High Gain, v-component
%           'dc_spec_dv12'           - DFB DC-coupled Differential Voltage V1-V2 antenas
%           'dc_spec_SCMv'           - DFB DC-coupled SCM, Low Frequency, High Gain, v-component, Sensor coordinates
%           'dc_spec_SCMu'           - DFB DC-coupled SCM, Low Frequency, High Gain, u-component, Sensor coordinates
%           'dc_spec_SCMw'           - DFB DC-coupled SCM, Low Frequency, High Gain, w-component, Sensor coordinates
%           'dc_spec_SCMd'           - DFB DC-coupled SCM, Low Frequency, High Gain, d-component, Sensor coordinates only Enc 1 (2018/11)
%           'dc_spec_SCMe'           - DFB DC-coupled SCM, Low Frequency, High Gain, e-component, Sensor coordinates only Enc 1 (2018/11)
%           'dc_spec_SCMf'           - DFB DC-coupled SCM, Low Frequency, High Gain, f-component, Sensor coordinates only Enc 1 (2018/11)
%           'l2_rfs_lfr'             - Radio Frequency Spectrometer, RFS, Low Frequency Reciever - level 2 LFR
%           'l2_rfs_hfr'             - Radio Frequency Spectrometer, RFS, High Frequency Reciever - level 2 HFR
%           'l3_rfs_lfr'             - Radio Frequency Spectrometer, RFS, Low Frequency Reciever - level 3 LFR
%           'l3_rfs_hfr'             - Radio Frequency Spectrometer, RFS, High Frequency Reciever - level 3 HFR
%           'tds_wf'                 - TDS Waveform Bursts-high-gain difference of V2(Ch2)/V5(Ch5) NOT corrected for frequency response (engineering mV).
%           'tds_wf'                 - TDS Waveform Bursts-high-gain difference of V1-V2(Ch3)/V3-V4(Ch1) NOT corrected for frequency response (engineering mV).
%           'tds_wf_scm4'            - TDS Waveform Bursts-high-gain SCM MF sensor(Ch4) NOT corrected for frequency response (engineering nT)
%           'tds_wf_sol_dist'        - TDS Waveform Bursts-provided distance to the sun at the time of this burst.
%           'spc'                    - SWEAP Solar Probe Cup data for proton moments
%           'spa_sf0'                - SWEAP SPAN-Ae SPA 3D spectra measured by SPAN-Ae sensor on the ram side of PSP
%           'spb_sf0'                - SWEAP SPAN-Be SPB 3D spectra measured by SPAN-Be sensor on the anti-ram side of PSP
%           'spa_sf1'                - SWEAP SPAN-Ae SPA Diff energy flux measured by SPAN-Ae sensor on the ram side of PSP
%           'spb_sf1'                - SWEAP SPAN-Be SPB Diff energy flux measured by SPAN-Be sensor on the anti-ram side of PSP
%           'spe'                    - SWEAP SPAN-E SPE Electron Pitch Angle Distribution
%           'spi'                    - SWEAP SPAN-I Partial Proton Moments, Differential Energy Eflux in SC and RTN coordinates
%           'spi_inst'               - SWEAP SPAN-I Partial Proton Moments, Differential Energy Eflux in instrument coordinates
%           'pl_spi'                 - Loads only density data to be used in the plasma line code
%           'pl_rfs_lfr'             - Loads V1V2 data in both auto and high-resolution for the plasma line code
%           'dbm_dvac'               - times series of electric field snapshots
%           'dbm_vac'                - times series of electric field snapshots
%           'dbm_dvdc'               - times series of electric field snapshots
%           'dbm_vdc'                - times series of electric field snapshots
%           'dbm_scm'                - times series of SCM snapshots
%           'qtn'                    - Electron density and temperature from quasi-thermal noise
%           'XXX'                    - where XXX is long or short name of variables in the list "PSP_VAR *"
%
% dateStart & dateStop: date vectors or strings for start and stop day
%            e.g. [yyyy mm dd] or 'yyyy/mm/dd' or  'yyyy mm dd' (read with datenum)
%
% Example:
%   psp_load('./','mag',[2020 01 26],[2020 02 01]); % assumes cdf files to be in the current directory
%   psp_load('psp_fld_l2_mag_RTN_20...cdf','mag');
%   rtnB = psp_load('psp_fld_l2_mag_RTN_20...cdf','mag');

function [output,pspobj] = psp_load(arg1,datatype,date_start,date_stop)
global dateStart dateEnd
persistent webOptionsSSL doGetFiles

% define webOptionsSSL during the first run of the funciton
if isnumeric(webOptionsSSL)
  userName = datastore('psp','usernameSSL');
  if isempty(userName)
    webOptionsSSL = {};
    disp('If you have access to the latest PSP data');
    disp('then from the matlab command line enter the PSP credentials.')
    disp('>datastore(''psp'',''usernameSSL'',''xxx'');')
    disp('>datastore(''psp'',''passwordSSL'',''xxx'');')
  else
    passWord = datastore('psp','passwordSSL');
    webOptionsSSL = weboptions('HeaderFields',{'Authorization',...
      ['Basic ' matlab.net.base64encode([userName ':' passWord])]});
  end
end

outputToBase  = (nargout == 0);
if nargin == 1
  datatype = arg1;
  startDatenum = datenum(dateStart); % globals
  endDatenum = datenum(dateEnd);
else
  dataDir = arg1;
  if nargin >=3
    startDatenum = datenum(date_start);
    endDatenum = datenum(date_stop);
  end
end
pspobj = []; % default output
listCdfFiles = {};
useStoredPspDirectory = false;

% If data directory not given, read in the stored value
if ~exist('dataDir','var') || isempty(dataDir)
  dataDir=datastore('psp','data_directory');
  if isempty(dataDir)
    disp('Your PSP directory is not defined!')
    disp('Please enter the location of the PSP data directory.');
    disp('In that folder should be at least subfolders "fields" and "sweap"');
    disp('that has the same structure as the ones on the FIELDS and SWEAP servers.');
    dataDir = input('full directory path without ending slash:','s');
    datastore('psp','data_directory',dataDir);
  end
  useStoredPspDirectory = true;
end

switch datatype

  case {'wf_vdc'}
    filename= 'psp_fld_l2_dfb_wf_scm';
    varnames = {'psp_fld_l2_dfb_wf_scm_hg_sensor'};
    varnamesout = {'wf_scm_sensor'};

    hourtag={'00';'06';'12';'18'};

  case {'wf_scm_sensor'}
    filename= 'psp_fld_l2_dfb_wf_scm';
    varnames = {'psp_fld_l2_dfb_wf_scm_hg_sensor'};
    varnamesout = {'wf_scm_sensor'};

    hourtag={'00';'06';'12';'18'};

  case {'wf_scm_sc'}
    filename= 'psp_fld_l2_dfb_wf_scm';
    varnames = {'psp_fld_l2_dfb_wf_scm_hg_sc'};
    varnamesout = {'wf_scm_sc'};

    hourtag={'00';'06';'12';'18'};

  case {'ac_bpf_dv34'}
    filename = 'psp_fld_l2_dfb_ac_bpf_dV34hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_bpf_dV34hg_frequency_bins';...
      'psp_fld_l2_dfb_ac_bpf_dV34hg_avg';...
      'psp_fld_l2_dfb_ac_bpf_dV34hg_peak'};
    varnamesout = {'dfb_ac_bpf_dv34_freq_bins';...
      'dfb_ac_bpf_dv34_avg';'dfb_ac_bpf_dv34_peak'};

    hourtag={''};

  case {'ac_bpf_SCMulfhg'}
    filename = 'psp_fld_l2_dfb_ac_bpf_SCMulfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_bpf_SCMulfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_bpf_SCMulfhg_avg';...
      'psp_fld_l2_dfb_ac_bpf_SCMulfhg_peak'};
    varnamesout = {'dfb_ac_bpf_SCMulfhg_freq_bins';...
      'dfb_ac_bpf_SCMulfhg_avg';'dfb_ac_bpf_SCMulfhg_peak'};

    hourtag={''};

  case {'ac_bpf_SCMumfhg'}
    filename = 'psp_fld_l2_dfb_ac_bpf_SCMumfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_bpf_SCMumfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_bpf_SCMumfhg_avg';...
      'psp_fld_l2_dfb_ac_bpf_SCMumfhg_peak'};
    varnamesout = {'dfb_ac_bpf_SCMumfhg_freq_bins';...
      'dfb_ac_bpf_SCMumfhg_avg';'dfb_ac_bpf_SCMumfhg_peak'};

    hourtag={''};

  case {'ac_spec_dv12'}
    filename = 'psp_fld_l2_dfb_ac_spec_dV12hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_dV12hg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_dV12hg'};
    varnamesout = {'ac_spec_dv12_freq_bins';'ac_spec_dv12_pw'};

    hourtag={''};

  case {'ac_spec_dv34'}
    filename = 'psp_fld_l2_dfb_ac_spec_dV34hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_dV34hg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_dV34hg'};
    varnamesout = {'ac_spec_dv34_freq_bins';'ac_spec_dv34_pw'};

    hourtag={''};

  case {'ac_spec_v5'}
    filename = 'psp_fld_l2_dfb_ac_spec_V5hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_V5hg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_V5hg'};
    varnamesout = {'ac_spec_v5_freq_bins';'ac_spec_v5_pw'};

    hourtag={''};

  case {'ac_spec_scmv'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMvlfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMvlfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_SCMvlfhg'};
    varnamesout = {'ac_spec_scmv_freq_bins';'ac_spec_scmv_pw'};

    hourtag={''};

  case {'ac_spec_scmu'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMulfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMulfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_SCMulfhg'};
    varnamesout = {'ac_spec_scmu_freq_bins';'ac_spec_scmu_pw'};

    hourtag={''};

  case {'ac_spec_scmd'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMdlfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMdlfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_SCMdlfhg'};
    varnamesout = {'ac_spec_scmd_freq_bins';'ac_spec_scmd_pw'};

    hourtag={''};

  case {'ac_spec_scme'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMelfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMelfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_SCMelfhg'};
    varnamesout = {'ac_spec_scme_freq_bins';'ac_spec_scme_pw'};

    hourtag={''};

  case {'ac_spec_scmf'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMflfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMflfhg_frequency_bins';...
      'psp_fld_l2_dfb_ac_spec_SCMflfhg'};
    varnamesout = {'ac_spec_scmf_freq_bins';'ac_spec_scmf_pw'};

    hourtag={''};


  case {'ac_xspec'}
    filename = 'psp_fld_l2_dfb_ac_xspec_dV12hg_dV34hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_xspec_dV12hg_dV34hg_frequency_bins';...
      'psp_fld_l2_dfb_ac_xspec_power_ch1_dV12hg';...
      'psp_fld_l2_dfb_ac_xspec_power_ch2_dV34hg';...
      'psp_fld_l2_dfb_ac_xspec_coh_dV12hg_dV34hg';...
      'psp_fld_l2_dfb_ac_xspec_phase_dV12hg_dV34hg';...
      'psp_fld_l2_dfb_ac_xspec_crossterm_Real_dV12hg_dV34hg';...
      'psp_fld_l2_dfb_ac_xspec_crossterm_Imag_dV12hg_dV34hg'};
    varnamesout = {'ac_xspec_dv12_dv34_freq_bins';...
      'ac_xspec_pw_dv12';'ac_xspec_pw_dv34';...
      'ac_xspec_coh_dv12_dv34';'ac_xspec_ph_dv12_dv34';...
      'ac_xspec_xterm_Re_dv12_dv34';...
      'ac_xspec_xterm_Im_dv12_dv34'};

    hourtag={''};

  case {'dc_bpf_dv34'}
    filename = 'psp_fld_l2_dfb_dc_bpf_dV34hg';
    varnames = {...
      'psp_fld_l2_dfb_dc_bpf_dV34hg_frequency_bins';...
      'psp_fld_l2_dfb_dc_bpf_dV34hg_avg';...
      'psp_fld_l2_dfb_dc_bpf_dV34hg_peak'};
    varnamesout = {'dfb_dc_bpf_dv34_freq_bins';...
      'dfb_dc_bpf_dv34_avg';'dfb_dc_bpf_dv34_peak'};

    hourtag={''};

  case {'dc_bpf_SCMulfhg'}
    filename = 'psp_fld_l2_dfb_dc_bpf_SCMulfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_bpf_SCMulfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_bpf_SCMulfhg_avg';...
      'psp_fld_l2_dfb_dc_bpf_SCMulfhg_peak'};
    varnamesout = {'dfb_dc_bpf_SCMulfhg_freq_bins';...
      'dfb_dc_bpf_SCMulfhg_avg';'dfb_dc_bpf_SCMulfhg_peak'};

    hourtag={''};

  case {'dc_bpf_SCMvlfhg'}
    filename = 'psp_fld_l2_dfb_dc_bpf_SCMvlfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_bpf_SCMvlfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_bpf_SCMvlfhg_avg';...
      'psp_fld_l2_dfb_dc_bpf_SCMvlfhg_peak'};
    varnamesout = {'dfb_dc_bpf_SCMvlfhg_freq_bins';...
      'dfb_dc_bpf_SCMvlfhg_avg';'dfb_dc_bpf_SCMvlfhg_peak'};

    hourtag={''};

  case {'dc_spec_dv12'}
    filename = 'psp_fld_l2_dfb_dc_spec_dV12hg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_dV12hg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_dV12hg'};
    varnamesout = {'dc_spec_dv12_freq_bins';'dc_spec_dv12_pw'};

    hourtag={''};

  case {'dc_spec_scmu'}
    filename = 'psp_fld_l2_dfb_dc_spec_SCMulfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_SCMulfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_SCMulfhg'};
    varnamesout = {'dc_spec_scmu_freq_bins';'dc_spec_scmu_pw'};

    hourtag={''};

  case {'dc_spec_scmv'}
    filename = 'psp_fld_l2_dfb_dc_spec_SCMvlfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_SCMvlfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_SCMvlfhg'};
    varnamesout = {'dc_spec_scmv_freq_bins';'dc_spec_scmv_pw'};

    hourtag={''};

  case {'dc_spec_scmw'}
    filename = 'psp_fld_l2_dfb_dc_spec_SCMwlfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_SCMwlfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_SCMwlfhg'};
    varnamesout = {'dc_spec_scmw_freq_bins';'dc_spec_scmw_pw'};

    hourtag={''};

  case {'dc_spec_scmd'}
    filename = 'psp_fld_l2_dfb_dc_spec_SCMdlfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_SCMdlfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_SCMdlfhg'};
    varnamesout = {'dc_spec_scmd_freq_bins';'dc_spec_scmd_pw'};

    hourtag={''};

  case {'dc_spec_scme'}
    filename = 'psp_fld_l2_dfb_dc_spec_SCMelfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_SCMelfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_SCMelfhg'};
    varnamesout = {'dc_spec_scme_freq_bins';'dc_spec_scme_pw'};

    hourtag={''};

  case {'dc_spec_scmf'}
    filename = 'psp_fld_l2_dfb_dc_spec_SCMflfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_spec_SCMflfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_spec_SCMflfhg'};
    varnamesout = {'dc_spec_scmf_freq_bins';'dc_spec_scmf_pw'};

    hourtag={''};

  case {'dc_xspec_SCMvw'}
    filename = 'psp_fld_l2_dfb_dc_xspec_SCMvlfhg_SCMwlfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_xspec_SCMvlfhg_SCMwlfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_xspec_power_ch1_SCMvlfhg';...
      'psp_fld_l2_dfb_dc_xspec_power_ch2_SCMwlfhg';...
      'psp_fld_l2_dfb_dc_xspec_coh_SCMvlfhg_SCMwlfhg';...
      'psp_fld_l2_dfb_dc_xspec_phase_SCMvlfhg_SCMwlfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Real_SCMvlfhg_SCMwlfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Imag_SCMvlfhg_SCMwlfhg'};
    varnamesout = {'dc_xspec_scmv_scmw_freq_bins';...
      'dc_xspec_pw_scmv';'dc_xspec_pw_scmw';...
      'dc_xspec_coh_scmv_scmw';'dc_xspec_ph_scmv_scmw';...
      'dc_xspec_xterm_Re_scmv_scmw';...
      'dc_xspec_xterm_Im_scmv_scmw'};

    hourtag={''};

  case {'dc_xspec_SCMde'}
    filename = 'psp_fld_l2_dfb_dc_xspec_SCMdlfhg_SCMelfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_xspec_SCMdlfhg_SCMelfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_xspec_power_ch1_SCMdlfhg';...
      'psp_fld_l2_dfb_dc_xspec_power_ch2_SCMelfhg';...
      'psp_fld_l2_dfb_dc_xspec_coh_SCMdlfhg_SCMelfhg';...
      'psp_fld_l2_dfb_dc_xspec_phase_SCMdlfhg_SCMelfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Real_SCMdlfhg_SCMelfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Imag_SCMdlfhg_SCMelfhg'};
    varnamesout = {'dc_xspec_scmd_scme_freq_bins';...
      'dc_xspec_pw_scmd';'dc_xspec_pw_scme';...
      'dc_xspec_coh_scmd_scme';'dc_xspec_ph_scmd_scme';...
      'dc_xspec_xterm_Re_scmd_scme';...
      'dc_xspec_xterm_Im_scmd_scme'};

    hourtag={''};

  case {'dc_xspec_SCMdf'}
    filename = 'psp_fld_l2_dfb_dc_xspec_SCMdlfhg_SCMflfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_xspec_SCMdlfhg_SCMflfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_xspec_power_ch1_SCMdlfhg';...
      'psp_fld_l2_dfb_dc_xspec_power_ch2_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_coh_SCMdlfhg_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_phase_SCMdlfhg_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Real_SCMdlfhg_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Imag_SCMdlfhg_SCMflfhg'};
    varnamesout = {'dc_xspec_scmd_scmf_freq_bins';...
      'dc_xspec_pw_scmd';'dc_xspec_pw_scmf';...
      'dc_xspec_coh_scmd_scmf';'dc_xspec_ph_scmd_scmf';...
      'dc_xspec_xterm_Re_scmd_scmf';...
      'dc_xspec_xterm_Im_scmd_scmf'};

    hourtag={''};

  case {'dc_xspec_SCMef'}
    filename = 'psp_fld_l2_dfb_dc_xspec_SCMelfhg_SCMflfhg';
    varnames = {...
      'psp_fld_l2_dfb_dc_xspec_SCMelfhg_SCMflfhg_frequency_bins';...
      'psp_fld_l2_dfb_dc_xspec_power_ch1_SCMelfhg';...
      'psp_fld_l2_dfb_dc_xspec_power_ch2_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_coh_SCMelfhg_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_phase_SCMelfhg_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Real_SCMelfhg_SCMflfhg';...
      'psp_fld_l2_dfb_dc_xspec_crossterm_Imag_SCMelfhg_SCMflfhg'};
    varnamesout = {'dc_xspec_scme_scmf_freq_bins';...
      'dc_xspec_pw_scme';'dc_xspec_pw_scmf';...
      'dc_xspec_coh_scme_scmf';'dc_xspec_ph_scme_scmf';...
      'dc_xspec_xterm_Re_scme_scmf';...
      'dc_xspec_xterm_Im_scme_scmf'};

    hourtag={''};

  case {'l2_rfs_lfr'}
    filename = 'psp_fld_l2_rfs_lfr';
    varnames = {...
      'psp_fld_l2_rfs_lfr_auto_averages_ch0_V1V2';...
      'frequency_lfr_auto_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_lfr_auto_averages_ch1_V3V4';...
      'frequency_lfr_auto_averages_ch1_V3V4';...
      'psp_fld_l2_rfs_lfr_auto_peaks_ch0_V1V2';...
      'psp_fld_l2_rfs_lfr_auto_peaks_ch1_V3V4';...
      'psp_fld_l2_rfs_lfr_hires_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_lfr_hires_averages_ch1_V3V4';...
      'psp_fld_l2_rfs_lfr_hires_peaks_ch0_V1V2';...
      'psp_fld_l2_rfs_lfr_hires_peaks_ch1_V3V4';...
      'psp_fld_l2_rfs_lfr_cross_im_V1V2_V3V4';...
      'psp_fld_l2_rfs_lfr_cross_re_V1V2_V3V4';...
      'psp_fld_l2_rfs_lfr_coher_V1V2_V3V4';...
      'psp_fld_l2_rfs_lfr_phase_V1V2_V3V4'};
    varnamesout = {'l2_rfs_lfr_v1v2';...
      'l2_rfs_lfr_v1v2_freq';...
      'l2_rfs_lfr_v3v4';...
      'l2_rfs_lfr_v3v4_freq';...
      'l2_rfs_lfr_auto_pks_V1V2';...
      'l2_rfs_lfr_auto_pks_V3V4';...
      'l2_rfs_lfr_hires_avs_V1V2';...
      'l2_rfs_lfr_hires_avs_V3V4';...
      'l2_rfs_lfr_hires_pks_V1V2';...
      'l2_rfs_lfr_hires_pks_V3V4';...
      'l2_rfs_lfr_cross_im_V1V2_V3V4';...
      'l2_rfs_lfr_cross_re_V1V2_V3V4';...
      'l2_rfs_lfr_coher_V1V2_V3V4';...
      'l2_rfs_lfr_phase_V1V2_V3V4'};

    hourtag={''};

  case {'l3_rfs_lfr'}
    filename = 'psp_fld_l3_rfs_lfr';
    varnames = {...
      'psp_fld_l3_rfs_lfr_auto_averages_ch0_V1V2';...
      'psp_fld_l3_rfs_lfr_auto_averages_ch0_V1V2_flux';...
      'frequency_lfr_auto_averages_ch0_V1V2';...
      'psp_fld_l3_rfs_lfr_auto_averages_ch1_V3V4';...
      'psp_fld_l3_rfs_lfr_auto_averages_ch1_V3V4_flux';...
      'frequency_lfr_auto_averages_ch1_V3V4';...
      'psp_fld_l3_rfs_lfr_hires_averages_ch0_V1V2';...
      'psp_fld_l3_rfs_lfr_hires_averages_ch0_V1V2_flux'
      'psp_fld_l3_rfs_lfr_hires_averages_ch1_V3V4';...
      'psp_fld_l3_rfs_lfr_hires_averages_ch1_V3V4_flux'
      'psp_fld_l3_rfs_lfr_cross_im_V1V2_V3V4';...
      'psp_fld_l3_rfs_lfr_cross_re_V1V2_V3V4';...
      'psp_fld_l3_rfs_lfr_coher_V1V2_V3V4';...
      'psp_fld_l3_rfs_lfr_phase_V1V2_V3V4';...
      'psp_fld_l3_rfs_lfr_PSD_FLUX';...
      'psp_fld_l3_rfs_lfr_PSD_SFU';...
      'psp_fld_l3_rfs_lfr_STOKES_V'};
    varnamesout = {...
      'l3_rfs_lfr_v1v2';...
      'l3_rfs_lfr_v1v2_flux';...
      'l3_rfs_lfr_v1v2_freq';...
      'l3_rfs_lfr_v3v4';...
      'l3_rfs_lfr_v3v4_flux';...
      'l3_rfs_lfr_v3v4_freq';...
      'l3_rfs_lfr_hires_avs_V1V2';...
      'l3_rfs_lfr_hires_avs_V1V2_flux';...
      'l3_rfs_lfr_hires_avs_V3V4';...
      'l3_rfs_lfr_hires_avs_V3V4_flux';...
      'l3_rfs_lfr_cross_im_V1V2_V3V4';...
      'l3_rfs_lfr_cross_re_V1V2_V3V4';...
      'l3_rfs_lfr_coher_V1V2_V3V4';...
      'l3_rfs_lfr_phase_V1V2_V3V4';...
      'l3_rfs_lfr_psd_flux';...
      'l3_rfs_lfr_psd_sfu';...
      'l3_rfs_lfr_stokes_v'};

    hourtag={''};

  case {'l2_rfs_hfr'}
    filename = 'psp_fld_l2_rfs_hfr';
    varnames = {...
      'psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2';...
      'frequency_hfr_auto_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4';...
      'frequency_hfr_auto_averages_ch1_V3V4'};
    varnamesout = {...
      'l2_rfs_hfr_v1v2';...
      'l2_rfs_hfr_v1v2_freq';...
      'l2_rfs_hfr_v3v4';...
      'l2_rfs_hfr_v3v4_freq'};

    hourtag={''};

  case {'l3_rfs_hfr'}
    filename = 'psp_fld_l3_rfs_hfr';
    varnames = {...
      'psp_fld_l3_rfs_hfr_auto_averages_ch0_V1V2';...
      'psp_fld_l3_rfs_hfr_auto_averages_ch0_V1V2_flux';...
      'frequency_hfr_auto_averages_ch0_V1V2';...
      'psp_fld_l3_rfs_hfr_auto_averages_ch1_V3V4';...
      'psp_fld_l3_rfs_hfr_auto_averages_ch1_V3V4_flux';...
      'frequency_hfr_auto_averages_ch1_V3V4';...
      'psp_fld_l3_rfs_hfr_cross_im_V1V2_V3V4';...
      'psp_fld_l3_rfs_hfr_cross_re_V1V2_V3V4';...
      'psp_fld_l3_rfs_hfr_coher_V1V2_V3V4';...
      'psp_fld_l3_rfs_hfr_phase_V1V2_V3V4';...
      'psp_fld_l3_rfs_hfr_PSD_FLUX';...
      'psp_fld_l3_rfs_hfr_PSD_SFU';...
      'psp_fld_l3_rfs_hfr_STOKES_V'};
    varnamesout = {...
      'l3_rfs_hfr_v1v2';...
      'l3_rfs_hfr_v1v2_flux';...
      'l3_rfs_hfr_v1v2_freq';...
      'l3_rfs_hfr_v3v4';...
      'l3_rfs_hfr_v3v4_flux';...
      'l3_rfs_hfr_v3v4_freq';...
      'l3_rfs_hfr_viv2_v3v4_cross_im';...
      'l3_rfs_hfr_viv2_v3v4_cross_re';...
      'l3_rfs_hfr_coher_v1v2_v3v4';...
      'l3_rfs_hfr_phase_v1v2_v3v4';...
      'l3_rfs_hfr_psd_flux';...
      'l3_rfs_hfr_psd_sfu';...
      'l3_rfs_hfr_stokes_v';};

    hourtag={''};

  case {'spc'}
    filename = 'psp_swp_spc_l3i';

    varnames = {...
      'DQF';...          % data quality flag
      'general_flag';...
      'np_fit';'wp_fit';'vp_fit_SC';'vp_fit_RTN';...
      'np1_fit';'wp1_fit';'vp1_fit_SC';'vp1_fit_RTN';...
      'np_moment';'wp_moment';'vp_moment_SC';'vp_moment_RTN';...
      'sc_pos_HCI';'sc_vel_HCI'};
    varnamesout = {'spc_dqf';'spc_general_flag';...
      'spc_np_fit';'spc_wp_fit';'spc_vp_fit_sc';'spc_vp_fit_rtn';...
      'spc_np1_fit';'spc_wp1_fit';'spc_vp1_fit_sc';'spc_vp1_fit_rtn';...
      'spc_np_moment';'spc_wp_moment';'spc_vp_moment_sc';...
      'spc_vp_moment_rtn';'spc_sc_pos_hci';'spc_sc_vel_hci'};

    hourtag={''};

  case {'spa_sf0'}
    filename = 'psp_swp_spa_sf0_L2_16Ax8Dx32E';
    varnames = {...
      'EFLUX';'ENERGY';...
      'THETA';'PHI';...
      'QUALITY_FLAG'};
    varnamesout = {'spa_sf0_eflux';...
      'spa_sf0_energy';'spa_sf0_theta';...
      'spa_sf0_phi';'spa_sf0_qf'};

    hourtag={''};

  case {'spb_sf0'}
    filename = 'psp_swp_spb_sf0_L2_16Ax8Dx32E';
    varnames = {...
      'EFLUX';'ENERGY';...
      'THETA';'PHI';...
      'QUALITY_FLAG'};
    varnamesout = {'spb_sf0_eflux';...
      'spb_sf0_energy';'spb_sf0_theta';...
      'spb_sf0_phi';'spb_sf0_qf'};

    hourtag={''};

  case {'spa_sf1'}
    filename = 'psp_swp_spa_sf1_L2_32E';
    varnames = {...
      'EFLUX';'ENERGY';...
      'QUALITY_FLAG'};
    varnamesout = {'spa_sf1_eflux';...
      'spa_sf1_energy';'spa_sf1_qf'};

    hourtag={''};

  case {'spb_sf1'}
    filename = 'psp_swp_spb_sf1_L2_32E';
    varnames = {...
      'EFLUX';'ENERGY';...
      'QUALITY_FLAG'};
    varnamesout = {'spb_sf1_eflux';...
      'spb_sf1_energy';'spb_sf1_qf'};

    hourtag={''};

  case {'spe'}
    filename = 'psp_swp_spe_sf0_L3_pad';
    varnames = {...
      'QUALITY_FLAG';'PITCHANGLE';...
      'EFLUX_VS_ENERGY';'ENERGY_VALS';...
      'EFLUX_VS_PA_E'};
    varnamesout = {'spe_sf0_pad_qf';'spe_sf0_pad_pa';...
      'spe_sf0_pad_eflux_vs_en';'spe_sf0_pad_en';...
      'spe_sf0_pad_eflux_vs_pa_en'};

    hourtag={''};

  case {'spi'}
    filename = 'psp_swp_spi_sf00_L3_mom';
    varnames = {...
      'QUALITY_FLAG';'DENS';...
      'VEL_INST';'VEL_SC';'VEL_RTN_SUN';...
      'T_TENSOR_INST';'TEMP';...
      'EFLUX_VS_ENERGY';'EFLUX_VS_THETA';...
      'EFLUX_VS_PHI';'ENERGY_VALS';...
      'THETA_VALS';'PHI_VALS';'SUN_DIST';...
      'VENUS_DIST';'SC_VEL_RTN_SUN';...
      'MAGF_SC';'MAGF_INST'};
    varnamesout = {'spi_sf00_qf';'spi_sf00_dens';...
      'spi_sf00_vel_inst';'spi_sf00_vel_sc';'spi_sf00_vel_rtn'
      'spi_sf00_t_tensor';'spi_sf00_temp';...
      'spi_sf00_eflux_vs_en';'spi_sf00_eflux_theta';...
      'spi_sf00_eflux_phi';'spi_sf00_energy';...
      'spi_sf00_theta';'spi_sf00_phi';'spi_sf00_sun_dist';...
      'spi_sf00_venus_dist';'spi_sf00_PSP_vel_rtn';...
      'spi_sf00_mag_sc';'spi_sf00_mag_inst'};

    hourtag={''};

  case {'spi_inst'}
    filename = 'psp_swp_spi_sf00_L3_mom_INST';
    varnames = {...
      'QUALITY_FLAG';...
      'DENS';'VEL';'T_TENSOR';'TEMP';...
      'MAGF_INST';'EFLUX_VS_ENERGY';...
      'EFLUX_VS_THETA';'EFLUX_VS_PHI';...
      'ENERGY_VALS';'THETA_VALS';'PHI_VALS'};
    varnamesout = {'spi_sf00_INST_qf';'spi_sf00_INST_dens';...
      'spi_sf00_INST_vel';'spi_sf00_INST_t_tensor';'spi_sf00_INST_temp';...
      'spi_sf00_INST_mag_inst';'spi_sf00_INST_eflux_vs_en';...
      'spi_sf00_INST_eflux_theta';'spi_sf00_INST_eflux_phi';...
      'spi_sf00_INST_energy';'spi_sf00_INST_theta';'spi_sf00_INST_phi'};

    hourtag={''};

  case 'ephem'

    filename = 'spp_fld_l1_ephem_spp_rtn_';
    varnames = {'position';'velocity'};
    varnamesout = {'R_ephem';'V_ephem'};

    hourtag={''};

  case {'pl_spi'}
    filename = 'psp_swp_spi_sf00_L3_mom_INST';
    varnames = {'DENS'};
    varnamesout = {'spi_sf00_dens'};

    hourtag={''};

  case {'pl_rfs_lfr'}
    filename = 'psp_fld_l2_rfs_lfr';
    varnames = {...
      'psp_fld_l2_rfs_lfr_auto_averages_ch1_V3V4';...
      'frequency_lfr_auto_averages_ch1_V3V4'
      'psp_fld_l2_rfs_lfr_hires_averages_ch1_V3V4';...
      'frequency_lfr_hires_averages_ch1_V3V4'};
    varnamesout = {'rfs_lfr_v3v4';'rfs_lfr_v3v4_freq';...
      'rfs_lfr_v3v4_hr';'rfs_lfr_v3v4_freq_hr'};

    hourtag={''};

  case {'qtn'}
    filename = 'psp_fld_l3_sqtn_rfs_V1V2';
    varnames = {...
      'electron_density';...
      'electron_core_temperature';...
      'electron_density_delta'};
    varnamesout = {'qtn_ne';...
      'qtn_core_Te';'qtn_ne_delta'};

    hourtag={''};

    %   case {'tds_wf_sol_dist'}
    %     filename = 'psp_fld_l2_tds_wf';
    %     varnames = {'PSP_FLD_L2_TDS_WF_SC_Solar_Distance'};
    %     varnamesout = {'tds_wf_sol_dist'};
    %
    %     hourtag={''};
    %
    %   case 'tds_wf'
    %     listCDFFiles = get_file_list('psp_fld_l2_tds_wf');
    %     output       = get_data_tds_wf(listCDFFiles);
    %     return;

  case 'tds_wf_scm4'
    listCDFFiles = get_file_list('psp_fld_l2_tds_wf');
    output       = get_data_tds_wf(listCDFFiles);
    return;

  case 'dbm_dvac'
    listCDFFiles = get_file_list('psp_fld_l2_dfb_dbm_dvac');
    output       = get_data_dbm_dvac(listCDFFiles);
    return;

  case 'dbm_vac'
    listCDFFiles = get_file_list('psp_fld_l2_dfb_dbm_vac');
    output       = get_data_dbm_vac(listCDFFiles);
    return;

  case 'dbm_dvdc'
    listCDFFiles = get_file_list('psp_fld_l2_dfb_dbm_dvdc');
    output       = get_data_dbm_dvdc(listCDFFiles);
    return;
  case 'dbm_vdc'
    listCDFFiles = get_file_list('psp_fld_l2_dfb_dbm_vdc');
    output       = get_data_dbm_vdc(listCDFFiles);
    return;
  case 'dbm_scm'
    listCDFFiles = get_file_list('psp_fld_l2_dfb_dbm_scm');
    output       = get_data_dbm_scm(listCDFFiles);
    return;

  otherwise
    if nargin == 1 || nargin == 4
      % read in variables defined in psp_var
      out = psp_var(datatype);
      if isempty(out)
        error('Data type not recognized!');
      end
      fileBaseName = out.fileName;
      varnames  = {out.varName};
      hourtag  = out.hourtag;
      varnamesout = {out.varNameShort};
      listCdfFiles = get_file_list(fileBaseName);
      nFiles = numel(listCdfFiles);

    elseif nargin==2
      % read in single variable from a given file
      nFiles = 1;
      varName = datatype;
      filesToLoadTable = dataDir;
      [filename,varName,hourtag,shortVar] = psp_var(varName);
      if isempty(filename)
        error(['Filename not known for varname: ' varName '. Consider updating psp_load().']);
      elseif (filename ~= filesToLoadTable(1:length(filename)))
        error(['varname: ' varName ' not consistent with filename: ' filesToLoadTable]);
      end
      varnames = {varName};
      varnamesout = {shortVar};
    else
      error('Data type not recognized!')
    end
end

if ~exist('filesToLoadTable','var') && isempty(listCdfFiles)

  nDays  = endDatenum-startDatenum+1; % include first day

  datenumTable   = (startDatenum:endDatenum)';
  timestampTable = datestr(datenumTable,'yyyymmdd');

  % sanity check
  if size(timestampTable,1)~=nDays
    error('Number of days is odd.');
  end

  nHourtag = length(hourtag);
  nFiles = nDays*nHourtag;

  filesToLoadTable = char(zeros(nFiles,...
    size(timestampTable,2)+length(hourtag{1})));
  iFile = 1;
  for i = 1:nDays
    for i_hh=1:nHourtag
      filesToLoadTable(iFile,:) = [timestampTable(i,:) hourtag{i_hh}];
      iFile = iFile + 1;
    end
  end
  if useStoredPspDirectory
    dataDirList = char(get_data_path(filename,datenumTable));
    filesToLoadTable = strcat(dataDirList,filesep,filename,'_',filesToLoadTable,'_v00','.cdf');
  else
    filesToLoadTable= strcat(dataDir,filesep,filename,'_',filesToLoadTable,'_v00','.cdf');
  end
end

nVar      = length(varnames);
varData   = cell(nVar,1);
varMeta   = cell(nVar,1);
epochData = cell(nVar,1);

for iFile = 1:nFiles

  if numel(listCdfFiles)>0
    fileToLoad = listCdfFiles{iFile};
  else
    fileToLoad=strtrim(filesToLoadTable(iFile,:));
    d=dir([fileToLoad(1:end-6) '*']); % list all files with different versions
    if isempty(d)
      irf.log('warning','No data, trying to download.');
      download_data(filename,datenumTable)
      d=dir([fileToLoad(1:end-6) '*']);
    end
    if numel(d) == 1
      fileToLoad = [d.folder filesep d.name];
      irf.log('warning',['Reading: ' d.name]);
    elseif numel(d) > 1
      irf.log('warning','Several version files exist!')
      fileNamesFound = sort({d(:).name}) %#ok<NOPRT>
      irf.log('warning',['Using the latest version: ' fileNamesFound{end}])
      fileToLoad = [d(1).folder filesep fileNamesFound{end}];
    else
      irf.log('warning','No file found')
      fileToLoad = [];
    end
  end

  if any(fileToLoad)

    pspobj=dataobj(fileToLoad);

    for iVar = 1:nVar

      varname    = varnames{iVar};
      varnameout = varnamesout{iVar};

      var_data=get_ts(pspobj,varname);
      disp(['Loaded: ' varnameout ' corresponding to: ' varname ' in file: ' fileToLoad  ])

      if isempty(var_data.data)
        irf.log('critical','Empty data!');

        % first data, allocate cellarray and get metadata
      elseif isempty(varData{iVar})
        varData{iVar}=cell(nFiles,1);
        epochData{iVar}=cell(nFiles,1);
        varMeta{iVar}=var_data;
      end

      varData{iVar}{iFile} = var_data.data;
      epochData{iVar}{iFile} =var_data.time.epoch;

    end

  end

  clear fileToLoad
  clear timestamp
  clear var_data

end

irf.log('notice','Done loading. Merging data, preparing output ...');

if ~outputToBase
  output = cell(1,nVar);
end

% Fix the output
for iOutputVar = 1:nVar

  if isempty(varnamesout{iOutputVar})
    varname_out = varnames{iOutputVar};
  else
    varname_out = varnamesout{iOutputVar};
  end
  indEmptyFiles = false(length(varData{iOutputVar}),1);
  for i=1:length(varData{iOutputVar})
    if isempty(varData{iOutputVar}{i})
      indEmptyFiles(i) = true;
    end
  end
  var_data = cell2mat(varData{iOutputVar}(~indEmptyFiles));
  epoch_data = cell2mat(epochData{iOutputVar}(~indEmptyFiles));
  varOut =  TSeries(EpochTT(epoch_data),var_data,...
    'metadata_from',varMeta{iOutputVar});

  if outputToBase
    assignin('base',varname_out,varOut);
  else
    output{iOutputVar} = varOut;
    if nVar==1
      output = output{1};
    end
  end

end

  function iFilesToGet = ind_dates_in_datenum_interval(fileDates)
    iFilesToGet = false(numel(fileDates),1);
    for iFileDates = 1: numel(fileDates)
      dd = datenum(fileDates{iFileDates}{:},'yyyymmdd');
      if dd >= startDatenum && dd <= endDatenum
        iFilesToGet(iFileDates) = true;
      end
    end
  end

  function out = get_data_dir(fileBaseName,datenumTable)
    % datenumTable is array of datenums
    % out is cellarray of directories
    out = psp_var(['file=' fileBaseName]);
    if (numel(out) > 1), out = out{1}; end
    dirBase = out.directory;
    out = cell(numel(datenumTable),1);
    for i  = 1:numel(datenumTable)
      MM = datestr(datenumTable(i),'mm');
      YY = datestr(datenumTable(i),'yy');
      dirFull = strrep(dirBase,'MM',MM);
      dirFull = strrep(dirFull,'YY',YY);
      %      dirFull = strrep(dirFull,'/',filesep);
      out{i} = dirFull;
    end
    out = unique(out);
  end

  function out = get_data_path(fileBaseName,datenumTable)
    % adds local dataDir path
    out = get_data_dir(fileBaseName,datenumTable);
    for iOut = 1: numel(out)
      out{iOut} = [dataDir filesep out{iOut}];
    end
  end

  function out = get_file_list(fileBaseName)
    % uses startDatenum and endDatenum
    out = psp_var(['file=' fileBaseName]);
    if (numel(out) > 1), out = out{1}; end
    dirBase = out.directory;
    out = {};

    for date = floor(startDatenum):floor(endDatenum)
      dirs = get_data_dir(fileBaseName,date);
      dirFull = [dirs{1} filesep fileBaseName '_' datestr(date,'YYYYmmDD')];
      listDirFilter = [dataDir filesep dirFull '*'];
      listDir = dir(listDirFilter);
      irf.log('debug',['Listing files in: ' listDirFilter]);
      if numel(listDir) == 0
        irf.log('critical',['No cdf files found for fileBaseName=' fileBaseName ]);
        download_data(fileBaseName,date);
        listDir = dir(listDirFilter);
      end
      out = [out fullfile({listDir.folder},{listDir.name})];
    end
  end

%   function out = get_data_tds_wf(listCDFFiles)
%     % output is structure witf fields
%     % .ts: time series with two columns tds_wf_burst_V1V2, tds_wf_burst_V3V4, tds_wf_burst_V2 and tds_wf_burst_V5
%     % .startStopMatriTT: start and stop times of snapshots in TT
%     % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots
%     %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
%     tds_wf = double([]);
%     tFinal = [];
%     for iCdfFile = 1:numel(listCDFFiles)
%       fileCDF = listCDFFiles{iCdfFile};
%       disp(['Reading: ' fileCDF]);
%       res = spdfcdfread(fileCDF,'VARIABLES', {...
%         'PSP_FLD_L2_TDS_WF_Burst_Time_Series_Times_TT2000',...
%         'PSP_FLD_L2_TDS_WF_Burst_Time_Series_V2_Engineering_mV',...
%         'PSP_FLD_L2_TDS_WF_Burst_Time_Series_V5_Engineering_mV',...
%         'PSP_FLD_L2_TDS_WF_Burst_Time_Series_V1V2_Engineering_mV',...
%         'PSP_FLD_L2_TDS_WF_Burst_Time_Series_V3V4_Engineering_mV'},...
%         'KeepEpochAsIs',true,'dataonly',true);
%       tt=res{1}; tempV2 = res{2}; tempV5 = res{3}; tempV1V2 = res{4}; tempV3V4 = res{5};
%       t=[tt(1,:);tt;tt(end,:)];t=t(:);
%       tFinal = [tFinal; t];
%       vecNaN = nan(size(tt,2),1,'single');
%       tds_wf_tempV2 = [vecNaN tempV2 vecNaN]';
%       tds_wf_tempV5 = [vecNaN tempV5 vecNaN]';
%       tds_wf_tempV1V2 = [vecNaN tempV1V2 vecNaN]';
%       tds_wf_tempV3V4 = [vecNaN tempV3V4 vecNaN]';
%       tds_wf = [tds_wf; [tds_wf_tempV1V2(:) tds_wf_tempV3V4(:) tds_wf_tempV2(:) tds_wf_tempV5(:)]];
%     end
%     tStartStopTT = reshape(tFinal(diff(tFinal)==0),2,[])';
%     tds_wf(tds_wf < -1e30) = NaN;
%     tds_wf    = TSeries(EpochTT(tFinal),tds_wf);
%     tds_wf.name = ['V1V2,' ' ' 'V3V4,' ' ' 'V2,' ' ' 'V5'];
%     tds_wf.units = 'V';
%     out = struct('ts_wf',tds_wf,'startStopTT',tStartStopTT);
%   end

  function out = get_data_dbm_dvac(listCDFFiles)
    % output is structure witf fields
    % .ts: time series with two columns dvac12 and dvac34
    % .ts_sc: time series in sc coordinate system
    % .startStopMatriTT: start and stop times of snapshots in TT
    % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots
    %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
    dbm_dvac = double([]);
    tFinal = [];
    for iCdfFile = 1:numel(listCDFFiles)
      fileCDF = listCDFFiles{iCdfFile};
      disp(['Reading: ' fileCDF]);
      res = spdfcdfread(fileCDF,'VARIABLES', {...
        'psp_fld_l2_dfb_dbm_dvac_time_series_TT2000',...
        'psp_fld_l2_dfb_dbm_dvac12',...
        'psp_fld_l2_dfb_dbm_dvac34'},...
        'KeepEpochAsIs',true,'dataonly',true);
      tt=res{1}; temp12 = res{2}; temp34 = res{3};
      t=[tt(1,:);tt;tt(end,:)];t=t(:);
      tFinal = [tFinal; t];
      vecNaN = nan(size(tt,2),1,'single');
      dbm_dvac_temp12 = [vecNaN temp12 vecNaN]';
      dbm_dvac_temp34 = [vecNaN temp34 vecNaN]';
      dbm_dvac = [dbm_dvac; [dbm_dvac_temp12(:) dbm_dvac_temp34(:)]];
    end
    tStartStopTT = reshape(tFinal(diff(tFinal)==0),2,[])';
    dbm_dvac(dbm_dvac < -1e30) = NaN;
    dbm_dvac    = TSeries(EpochTT(tFinal),dbm_dvac);
    dbm_dvac_sc = psp_coordinate_transform(dbm_dvac,'e>sc');
    out = struct('ts',dbm_dvac,'ts_sc',dbm_dvac_sc,...
      'startStopTT',tStartStopTT);
  end

  function out = get_data_dbm_vac(listCDFFiles)
    % output is structure witf fields
    % .ts: time series with five columns vac1, vac2, vac3, vac4 and vac5
    % .ts_sc: time series in sc coordinate system
    % .startStopMatriTT: start and stop times of snapshots in TT
    % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots
    %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
    dbm_vac = double([]);
    tFinal = [];
    for iCdfFile = 1:numel(listCDFFiles)
      fileCDF = listCDFFiles{iCdfFile};
      disp(['Reading: ' fileCDF]);
      res = spdfcdfread(fileCDF,'VARIABLES', {...
        'psp_fld_l2_dfb_dbm_vac_time_series_TT2000',...
        'psp_fld_l2_dfb_dbm_vac1',...
        'psp_fld_l2_dfb_dbm_vac2',...
        'psp_fld_l2_dfb_dbm_vac3',...
        'psp_fld_l2_dfb_dbm_vac4',...
        'psp_fld_l2_dfb_dbm_vac5'},...
        'KeepEpochAsIs',true,'dataonly',true);
      tt=res{1}; temp1 = res{2};
      temp2 = res{3}; temp3 = res{4};
      temp4 = res{5}; temp5 = res{6};
      t=[tt(1,:);tt;tt(end,:)];t=t(:);
      tFinal = [tFinal; t];
      vecNaN = nan(size(tt,2),1,'single');
      dbm_vac_temp1 = [vecNaN temp1 vecNaN]';
      dbm_vac_temp2 = [vecNaN temp2 vecNaN]';
      dbm_vac_temp3 = [vecNaN temp3 vecNaN]';
      dbm_vac_temp4 = [vecNaN temp4 vecNaN]';
      dbm_vac_temp5 = [vecNaN temp5 vecNaN]';
      dbm_vac = [dbm_vac; [dbm_vac_temp1(:),...
        dbm_vac_temp2(:), dbm_vac_temp3(:),...
        dbm_vac_temp4(:), dbm_vac_temp5(:)]];
    end
    tStartStopTT = reshape(tFinal(diff(tFinal)==0),2,[])';
    dbm_vac(dbm_vac < -1e30) = NaN;
    dbm_vac    = TSeries(EpochTT(tFinal),dbm_vac);
    dbm_vac_sc = psp_coordinate_transform(dbm_vac,'e>sc');
    out = struct('ts',dbm_vac,'ts_sc',dbm_vac_sc,...
      'startStopTT',tStartStopTT);
  end

  function out = get_data_dbm_dvdc(listCDFFiles)
    % output is structure witf fields
    % .ts: time series with two columns dvdc12 and dvdc34
    % .ts_sc: time series in sc coordinate system
    % .startStopMatriTT: start and stop times of snapshots in TT
    % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots
    %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
    dbm_dvdc = double([]);
    tFinal = [];
    for iCdfFile = 1:numel(listCDFFiles)
      fileCDF = listCDFFiles{iCdfFile};
      disp(['Reading: ' fileCDF]);
      res = spdfcdfread(fileCDF,'VARIABLES', {...
        'psp_fld_l2_dfb_dbm_dvdc_time_series_TT2000',...
        'psp_fld_l2_dfb_dbm_dvdc12',...
        'psp_fld_l2_dfb_dbm_dvdc34'},...
        'KeepEpochAsIs',true,'dataonly',true);
      tt=res{1}; temp12 = res{2}; temp34 = res{3};
      t=[tt(1,:);tt;tt(end,:)];t=t(:);
      tFinal = [tFinal; t];
      vecNaN = nan(size(tt,2),1,'single');
      dbm_dvdc_temp12 = [vecNaN temp12 vecNaN]';
      dbm_dvdc_temp34 = [vecNaN temp34 vecNaN]';
      dbm_dvdc = [dbm_dvdc; [dbm_dvdc_temp12(:) dbm_dvdc_temp34(:)]];
    end
    tStartStopTT = reshape(tFinal(diff(tFinal)==0),2,[])';
    dbm_dvdc(dbm_dvdc < -1e30) = NaN;
    dbm_dvdc    = TSeries(EpochTT(tFinal),dbm_dvdc);
    dbm_dvdc_sc = psp_coordinate_transform(dbm_dvdc,'e>sc');
    out = struct('ts',dbm_dvdc,'ts_sc',dbm_dvdc_sc,...
      'startStopTT',tStartStopTT);
  end

  function out = get_data_dbm_vdc(listCDFFiles)
    % output is structure witf fields
    % .ts: time series with five columns vdc1, vdc2, vdc3, vdc4 and vdc5
    % .ts_sc: time series in sc coordinate system
    % .startStopMatriTT: start and stop times of snapshots in TT
    % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots
    %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
    dbm_vdc = double([]);
    tFinal = [];
    for iCdfFile = 1:numel(listCDFFiles)
      fileCDF = listCDFFiles{iCdfFile};
      disp(['Reading: ' fileCDF]);
      res = spdfcdfread(fileCDF,'VARIABLES', {...
        'psp_fld_l2_dfb_dbm_vdc_time_series_TT2000',...
        'psp_fld_l2_dfb_dbm_vdc1',...
        'psp_fld_l2_dfb_dbm_vdc2',...
        'psp_fld_l2_dfb_dbm_vdc3',...
        'psp_fld_l2_dfb_dbm_vdc4',...
        'psp_fld_l2_dfb_dbm_vdc5'},...
        'KeepEpochAsIs',true,'dataonly',true);
      tt=res{1}; temp1 = res{2};
      temp2 = res{3}; temp3 = res{4};
      temp4 = res{5}; temp5 = res{6};
      t=[tt(1,:);tt;tt(end,:)];t=t(:);
      tFinal = [tFinal; t];
      vecNaN = nan(size(tt,2),1,'single');
      dbm_vdc_temp1 = [vecNaN temp1 vecNaN]';
      dbm_vdc_temp2 = [vecNaN temp2 vecNaN]';
      dbm_vdc_temp3 = [vecNaN temp3 vecNaN]';
      dbm_vdc_temp4 = [vecNaN temp4 vecNaN]';
      dbm_vdc_temp5 = [vecNaN temp5 vecNaN]';
      dbm_vdc = [dbm_vdc; [dbm_vdc_temp1(:),...
        dbm_vdc_temp2(:), dbm_vdc_temp3(:),...
        dbm_vdc_temp4(:), dbm_vdc_temp5(:)]];
    end
    tStartStopTT = reshape(tFinal(diff(tFinal)==0),2,[])';
    dbm_vdc(dbm_vdc < -1e30) = NaN;
    dbm_vdc    = TSeries(EpochTT(tFinal),dbm_vdc);
    dbm_vdc_sc = psp_coordinate_transform(dbm_vdc,'e>sc');
    out = struct('ts',dbm_vdc,'ts_sc',dbm_vdc_sc,...
      'startStopTT',tStartStopTT);
  end

  function out = get_data_dbm_scm(listCDFFiles)
    % output is structure witf fields
    % .ts: time series with 3 coumns - SCM u,v.w
    % .startStopMatriTT: start and stop times of snapshots in TT
    % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots
    %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
    dbm_scm_hg = double([]);
    dbm_scm_lg = double([]);
    tHGFinal = [];
    tLGFinal = [];
    for iCdfFile = 1:numel(listCDFFiles)
      fileCDF = listCDFFiles{iCdfFile};
      disp(['Reading: ' fileCDF]);
      res = spdfcdfread(fileCDF,'VARIABLES', {...
        'psp_fld_l2_dfb_dbm_scmlg_time_series_TT2000',...
        'psp_fld_l2_dfb_dbm_scmlgu',...
        'psp_fld_l2_dfb_dbm_scmlgv',...
        'psp_fld_l2_dfb_dbm_scmlgw',...
        'psp_fld_l2_dfb_dbm_scmhg_time_series_TT2000',...
        'psp_fld_l2_dfb_dbm_scmhgu',...
        'psp_fld_l2_dfb_dbm_scmhgv',...
        'psp_fld_l2_dfb_dbm_scmhgw'},...
        'KeepEpochAsIs',true,'dataonly',true);
      ttLG=res{1}; lgu = res{2}; lgv = res{3}; lgw = res{4};
      ttHG=res{5}; hgu = res{6}; hgv = res{7}; hgw = res{8};
      % create continuous timeline marking start/stop times of snapshots
      % with additional time point with the same time value
      if any(ttLG)
        tLG=[ttLG(1,:);ttLG;ttLG(end,:)];tLG=tLG(:);
        tLGFinal = [tLGFinal; tLG];
        vecNaN = nan(size(ttLG,2),1,'single');
        dbm_scm_u = [vecNaN lgu vecNaN]';
        dbm_scm_v = [vecNaN lgv vecNaN]';
        dbm_scm_w = [vecNaN lgw vecNaN]';
        dbm_scm_lg = [dbm_scm_lg; [dbm_scm_u(:) dbm_scm_v(:) dbm_scm_w(:)]];
      end
      if any(ttHG)
        tHG=[ttHG(1,:);ttHG;ttHG(end,:)];tHG=tHG(:);
        tHGFinal = [tHGFinal; tHG];
        vecNaN = nan(size(ttHG,2),1,'single');
        dbm_scm_u = [vecNaN hgu vecNaN]';
        dbm_scm_v = [vecNaN hgv vecNaN]';
        dbm_scm_w = [vecNaN hgw vecNaN]';
        dbm_scm_hg = [dbm_scm_hg; [dbm_scm_u(:) dbm_scm_v(:) dbm_scm_w(:)]];
      end
    end
    tStartStopTTLG = reshape(tLGFinal(diff(tLGFinal)==0),2,[])';
    tStartStopTTHG = reshape(tHGFinal(diff(tHGFinal)==0),2,[])';
    dbm_scm_lg(dbm_scm_lg < -1e30) = NaN;
    dbm_scm_hg(dbm_scm_hg < -1e30) = NaN;
    dbm_scmLG   = TSeries(EpochTT(tLGFinal),dbm_scm_lg);
    dbm_scmHG   = TSeries(EpochTT(tHGFinal),dbm_scm_hg);
    out = struct('ts_LG',dbm_scmLG,'ts_HG',dbm_scmHG,...
      'startStopTTLG',tStartStopTTLG,'startStopTTHG',tStartStopTTHG);
  end

  function download_data(fileBaseName,dateNum)
    if isempty(doGetFiles) ...
        || ( ~strcmp(doGetFiles,'N') && ~strcmp(doGetFiles,'Y') )
      doGetFiles = irf_ask('Shall I download Y-yes to all, y-yes? [Y/y/n/N] [%]>','doGetFiles','Y');
    end

    if ~strcmpi(doGetFiles,'y'),return;end

    dataSubDir = get_data_dir( fileBaseName,dateNum);
    dataPath   = get_data_path(fileBaseName,dateNum);

    % define data servers and weboptions
    if strcmp(dataSubDir{1}(1:5),'sweap')
      webserver = 'http://sweap.cfa.harvard.edu/pub/data/sci/';
      webOptions = {};
    elseif strcmp(dataSubDir{1}(1:6),'fields')
      webserver = 'https://research.ssl.berkeley.edu/data/psp/data/sci/';
      webOptions = webOptionsSSL;
    end

    for iDir = 1:numel(dataSubDir)
      dirPath = dataPath{iDir};
      if ~exist(dirPath,'dir')
        mkdir(dirPath);
      end
      wwwDir = [webserver dataSubDir{iDir}];
      irf.log('debug',['Getting file list from www: ' wwwDir]);
      if isempty(webOptions)
        tt = webread(wwwDir);
      else
        tt = webread(wwwDir,webOptions);
      end

      % find latest version files (sort descending by version and keep
      % unique file basenames
      ff = regexp(tt,'([\w]*)_v([\d]*).cdf','tokens');
      ffFull = arrayfun(@(x) x{1}, regexp(tt,'([\w]*_v[\d]*.cdf)','tokens'));

      filesBase = arrayfun(@(x) (x{1}{1}),ff,'UniformOutput',false);
      filesVersion = arrayfun(@(x) str2num(x{1}{2}),ff);
      [filesVersionSorted,indSort] = sort(filesVersion,'descend');
      filesBaseSorted = filesBase(indSort);
      ffFullSorted = ffFull(indSort);
      [filesBaseUnique,indUnique] = unique(filesBaseSorted);
      filesUnique = ffFullSorted(indUnique);
      fileDates = regexp(filesUnique,'.*_(20\d\d\d\d\d\d)\d*_v.*.cdf','tokens');

      % find files with dates in the interval requested
      iFilesToGet = find(ind_dates_in_datenum_interval(fileDates));
      if any(iFilesToGet)
        for ii = 1:numel(iFilesToGet)
          irf.log('warning',['Downloading: ' filesUnique{iFilesToGet(ii)}]);
          wwwLink = [wwwDir '/' filesUnique{iFilesToGet(ii)}];
          filePath = [dirPath '/' filesUnique{iFilesToGet(ii)}];
          irf.log('warning',['Downloading: ' wwwLink]);
          if isempty(webOptions)
            outFileName = websave(filePath,wwwLink);
          else
            outFileName = websave(filePath,wwwLink,webOptions);
          end
          irf.log('warning',['Downloaded to: ' outFileName]);
        end
      else
        irf.log('warning','No files to download from the server matching the date');
      end
    end

  end

end
