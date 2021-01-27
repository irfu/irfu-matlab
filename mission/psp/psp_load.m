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
%           'wf_scm'                 - DBF Digital Fields Board Search Coil Magnetometer Waveform
%           'ac_spec_dv12'           - DFB AC-coupled Differential Voltage V12 antenas
%           'ac_spec_dv34'           - DFB AC-coupled Differential Voltage V34 antenas
%           'ac_spec_v5'             - DFB AC-coupled V5 Antenna Voltage, Monopole Mode
%           'ac_spec_scmv'           - DFB AC-coupled SCM, Low Frequency, High Gain, v-component, Sensor coordinates
%           'ac_spec_scmu'           - DFB AC-coupled SCM, Low Frequency, High Gain, u-component, Sensor coordinates
%           'dc_spec_dv12'           - DFB DC-coupled Differential Voltage V1-V2 antenas
%           'dc_spec_SCMv'           - DFB DC-coupled SCM, Low Frequency, High Gain, v-component, Sensor coordinates
%           'dc_spec_SCMu'           - DFB DC-coupled SCM, Low Frequency, High Gain, u-component, Sensor coordinates
%           'dc_spec_SCMw'           - DFB DC-coupled SCM, Low Frequency, High Gain, w-component, Sensor coordinates
%           'rfs_lfr'                - Radio Frequency Spectrometer, RFS, Low Frequency Reciever, LFR
%           'rfs_hfr'                - Radio Frequency Spectrometer, RFS, High Frequency Reciever, HFR
%           'sweap','spc'            - SWEAP proton moments
%           'spe'                    - SWEAP SPE Electron Pitch Angle Distribution
%           'ephem'                  - ephemeris files
%           'dbm_dvac'               - times series of electric field snapshots
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
persistent webOptionsSSL
if isnumeric(webOptionsSSL) 
  % first run of function
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


if ~exist('dataDir') || isempty(dataDir)
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
  case {'mag_rtn', 'fgm_rtn'}
    filename= 'psp_fld_l2_mag_RTN';
    varnames = {'psp_fld_l2_mag_RTN'};
    varnamesout = {'rtnB'};
    
    hourtag={'00';'06';'12';'18'};
    
  case {'mag_sc', 'fgm_sc'}
    filename= 'psp_fld_l2_mag_SC';
    varnames = {'psp_fld_l2_mag_SC'};
    varnamesout = {'scB'};
    
    hourtag={'00';'06';'12';'18'};
    
  case {'wf_dvdc'}
    filename= 'psp_fld_l2_dfb_wf_dvdc';
    varnames = {...
      'psp_fld_l2_dfb_wf_dVdc_sensor';...
      'psp_fld_l2_dfb_wf_dVdc_sc'};
    varnamesout = {'wf_dvdc_sensor';'wf_dvdc_sc'};
    
    hourtag={'00';'06';'12';'18'};  
    
  case {'wf_scm'}
    filename= 'psp_fld_l2_dfb_wf_scm';
    varnames = {...
      'psp_fld_l2_dfb_wf_scm_hg_sensor';...
      'psp_fld_l2_dfb_wf_scm_hg_sc'};
    varnamesout = {...
      'wf_scm_sensor';'wf_scm_sc'};
    
    hourtag={'00';'06';'12';'18'};    
    
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
    
  case {'rfs_lfr'}
    filename = 'psp_fld_l2_rfs_lfr';
    varnames = {...
      'psp_fld_l2_rfs_lfr_auto_averages_ch0_V1V2';...
      'frequency_lfr_auto_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_lfr_auto_averages_ch1_V3V4';...
      'frequency_lfr_auto_averages_ch1_V3V4'};
    varnamesout = {'rfs_lfr_v1v2';'rfs_lfr_v1v2_freq';...
      'rfs_lfr_v3v4';'rfs_lfr_v3v4_freq'};
    
    hourtag={''};
    
  case {'rfs_hfr'}
    filename = 'psp_fld_l2_rfs_hfr';
    varnames = {...
      'psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2';...
      'frequency_hfr_auto_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4';...
      'frequency_hfr_auto_averages_ch1_V3V4'};
    varnamesout = {'rfs_hfr_v1v2';'rfs_hfr_v1v2_freq';...
      'rfs_hfr_v3v4';'rfs_hfr_v3v4_freq'};
    
    hourtag={''};  
    
  case {'sweap', 'spc'}    
    filename = 'psp_swp_spc_l3i';
    
    varnames = {...
      'DQF';...          % data quality flag
      'general_flag';...
      'np_fit';'wp_fit';'vp_fit_SC';'vp_fit_RTN';...
      'np1_fit';'wp1_fit';'vp1_fit_SC';'vp1_fit_RTN';...
      'np_moment';'wp_moment';'vp_moment_SC';'vp_moment_RTN';...
      'sc_pos_HCI';'sc_vel_HCI'};
    varnamesout = {'DQF';'general_flag';...
      'np_fit';'wp_fit';'vp_fit_SC';'vp_fit_RTN';...
      'np1_fit';'wp1_fit';'vp1_fit_SC';'vp1_fit_RTN';...
      'np_moment';'wp_moment';'vp_moment_SC';'vp_moment_RTN';...
      'R_HCI';'V_HCI'};
    
    hourtag={''};
    
  case {'spe'}
    filename = 'psp_swp_spe_sf0_L3_pad';
    varnames = {...
      'QUALITY_FLAG';'PITCHANGLE';...
      'EFLUX_VS_ENERGY';'ENERGY_VALS';...
      'MAGF_SC'};
    varnamesout = {'spe_qf';'spe_pitchangle';...
      'spe_dif_eflux_vs_en';'spe_energy';...
      'spe_mag_SC'};
    
    hourtag={''};
    
  case 'ephem'
    
    filename = 'spp_fld_l1_ephem_spp_rtn_';
    varnames = {'position';'velocity'};
    varnamesout = {'R_ephem';'V_ephem'};
    
    hourtag={''};
    
  case 'dbm_dvac'
    listCDFFiles = get_file_list('psp_fld_l2_dfb_dbm_dvac');
    output       = get_data_dbm_dvac(listCDFFiles);
    return;
    
  otherwise
    if nargin == 1 % psp_load(varName) assumes global dateStart dateEnd
      % read in variable, lookup up files to read from
      out = psp_var(datatype);
      fileBaseName = out.fileName;
      if isempty(fileBaseName)
        error(['Filename not known for datatype: ' datatype '. Consider updating psp_load().']);
      end
      
      varName = out.varName;
      hourtag = out.hourtag;
      shortVar = out.varNameShort;

      listCdfFiles = get_file_list(fileBaseName);
      nFiles = length(listCdfFiles);
      if nFiles == 0
        irf.log('critical',['No cdf files found for fileBaseName=' fileBaseName ]);
        doGetFiles = irf_ask('Shall I download? [y/n] [%]>','doGetFiles','y');
        if strcmp(doGetFiles,'y')
          dataSubDir = get_data_dir(fileBaseName,startDatenum:endDatenum);
          dataPath = get_data_path(fileBaseName,startDatenum:endDatenum);
          if strcmp(dataSubDir{1}(1:5),'sweap')
            webserver = 'http://sweap.cfa.harvard.edu/pub/data/sci/';
            webOptions = {};
          elseif strcmp(dataSubDir{1}(1:6),'fields')
            webserver = 'http://research.ssl.berkeley.edu/data/psp/data/sci/';
            webOptions = webOptionsSSL;
          end
          for iDir = 1:numel(dataSubDir)
            dirPath = dataPath{iDir};
            mkdir(dirPath);
            wwwDir = [webserver dataSubDir{iDir}];
            if isempty(webOptions)
              tt = webread(wwwDir);
            else
              tt = webread(wwwDir,webOptions);
            end
            ff = regexp(tt,'([\w]*.cdf)','tokens');
            files = arrayfun(@(x) (x{1}),ff);
            files = unique(files);
            fileDates = regexp(files,'.*_(20\d\d\d\d\d\d)\d*_v.*.cdf','tokens');
            iFilesToGet = find(ind_dates_in_datenum_interval(fileDates));
            if any(iFilesToGet)
              for ii = 1:numel(iFilesToGet)
                irf.log('warning',['Downloading: ' files{iFilesToGet(ii)}]);
                wwwLink = [wwwDir '/' files{iFilesToGet(ii)}];
                filePath = [dirPath '/' files{iFilesToGet(ii)}];
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
          listCdfFiles = get_file_list(fileBaseName);
          nFiles = length(listCdfFiles);
          if nFiles == 0
            irf.log('warning','No data to load');
            return;
          end
        else
          return;
        end
      end
      varnames = {varName};
      varnamesout = {shortVar};
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
    elseif nargin == 4
      % read in variable, lookup up files to read from
      nFiles = 1;
      out = psp_var(datatype);
      fileBaseName = out.fileName;
      varName = out.varName;
      hourtag = out.hourtag;
      shortVar = out.varNameShort;
      if isempty(fileBaseName)
        error(['Filename not known for datatype: ' datatype '. Consider updating psp_load().']);
      end
      listCdfFiles = get_file_list(fileBaseName);
      nFiles = length(listCdfFiles);
      if nFiles == 0
        irf.log('critical',['No cdf files found for fileBaseName = ' fileBaseName]); return;
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
  
  if ~isempty(listCdfFiles)
    fileToLoad = listCdfFiles{iFile};
  else
    fileToLoad=strtrim(filesToLoadTable(iFile,:));
    d=dir([fileToLoad(1:end-6) '*']); % list all files with different versions
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
  
  varname_out = varnamesout{iOutputVar};
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
    % uses startDatenum adn endDatenum, maybe should be made as arguments
    out = psp_var(['file=' fileBaseName]);
    if (numel(out) > 1), out = out{1}; end
    dirBase = out.directory;
    for date = floor(startDatenum):floor(endDatenum)
      MM = datestr(date,'mm');
      YY = datestr(date,'yy');
      dirFull = strrep(dirBase,'MM',MM);
      dirFull = strrep(dirFull,'YY',YY);
      dirFull = [dirFull filesep fileBaseName '_' datestr(date,'YYYYmmDD')];
      listDir = dir([dataDir filesep dirFull '*']);
      out = fullfile({listDir.folder},{listDir.name});
    end
  end

  function out = get_data_dbm_dvac(listCDFFiles)
    % output is structure witf fields 
    % .ts: time series wirt dvac12 and dvac13
    % .startStopMatriTT: start and stop times of snapshots in TT
    % .startStopLineEpoch: vector with start stop times and NaNs inbetween, to plot intervals of snapshots 
    %               irf_plot([tSnapline tSnapline*0],'-.','markersize',5);
    dbm_dvac = double([]);
    tStartEndSnapshTT = [];
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
      tStartEndSnapshTT = [tStartEndSnapshTT; reshape(t(diff(t)==0),2,[])'];
      dbm_dvac_temp = nan(numel(t),2,'single');
      vecNaN = nan(size(tt,2),1,'single');
      dbm_dvac_temp12 = [vecNaN temp12 vecNaN]';
      dbm_dvac_temp34 = [vecNaN temp34 vecNaN]';
      dbm_dvac = [dbm_dvac; [dbm_dvac_temp12(:) dbm_dvac_temp34(:)]];
    end
    dbm_dvac(dbm_dvac < -1e30) = NaN;
    tSnapLineEpoch = irf_time(tStartEndSnapshTT(:),'ttns>epoch');
    tSnapLineEpoch = reshape(tSnapLineEpoch,[],2)';
    tSnapLineEpoch(end+1,:)=NaN;
    tSnapLineEpoch = tSnapLineEpoch(:);
    dbm_dvac    = TSeries(EpochTT(tFinal),dbm_dvac);
    dbm_dvac_sc = psp_coordinate_transform(dbm_dvac,'e>sc');
    out = struct('ts_v12v34',dbm_dvac,'ts_sc',dbm_dvac_sc,...
      'startStopMatriTT',tStartEndSnapshTT,...
      'startStopLineEpoch',tSnapLineEpoch);
  end
    
end
