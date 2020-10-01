%
% PSP_LOAD loads PSP data into TSeries
%
% PSP_LOAD(dataDir,datatype,dateStart,dateStop)
%
% dataDir:  the directory with the cdf files
% datatype: 'mag','fgm'      - FGM
%           'sweap','spc'    - sweap proton moments
%           'ephem'          - ephemeris files
% dateStart & dateStop: date vectors or strings for start and stop day
%            e.g. [yyyy mm dd] or 'yyyy/mm/dd' or  'yyyy mm dd' (read with datenum)
%

function psp_load(dataDir,datatype,date_start,date_stop)


switch datatype
  case {'mag', 'fgm'}
    filename= 'psp_fld_l2_mag_RTN_';
    varnames = {'psp_fld_l2_mag_RTN'};
    varnamesout = {'rtnB'};
    
    versiontag = '_v01';
    
    hourtag={'00';'06';'12';'18'};
    
  case {'sweap', 'spc'}
    
    
    filename = 'psp_swp_spc_l3i_';
    versiontag = '_v01';
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
    
  case 'ephem'
    
    filename = 'spp_fld_l1_ephem_spp_rtn_';
    varnames = {'position';'velocity'};
    varnamesout = {'R_ephem';'V_ephem'};
    
    versiontag = '_v00';
    
    hourtag={''};
    
  otherwise
    error('Data type not recognized!')
end



nVar      = length(varnames);
varData   = cell(nVar,1);
varMeta   = cell(nVar,1);
epochData = cell(nVar,1);


tStart = datenum(date_start);
tStop  = datenum(date_stop);
nDays  = tStop-tStart+1; % include first day

datenumTable   = (tStart:tStop)';
timestampTable = datestr(datenumTable,'yyyymmdd');

% sanity check
if size(timestampTable,1)~=nDays
  error('Number of days is odd.');
end

nHourtag = length(hourtag);


for i = 1:nDays
  for i_hh=1:nHourtag
    
    timestamp = [timestampTable(i,:) hourtag{i_hh}];
    
    filename2load=[dataDir filesep filename timestamp versiontag '.cdf'];
    
    irf.log('notice',['Loading: ' filename2load]);
    
    if exist(filename2load,'file')
      
      pspobj=dataobj(filename2load);
      
      for i_var = 1:nVar
        
        varname    = varnames{i_var};
        varnameout = varnamesout{i_var};
        
        disp(['Loaded: ' filename2load ', Extracting: ' varname ', Output: ' varnameout '_' timestamp ])
        
        var_data=get_ts(pspobj,varname);
        
        if isempty(var_data.data)
          error('Empty data!');
          
          % first data, allocate cellarray and get metadata
        elseif isempty(varData{i_var})
          varData{i_var}=cell(nDays*nHourtag,1);
          epochData{i_var}=cell(nDays*nHourtag,1);
          varMeta{i_var}=var_data;
        end
        
        varData{i_var}{(i-1)*nHourtag + i_hh} =var_data.data;
        epochData{i_var}{(i-1)*nHourtag + i_hh} =var_data.time.epoch;
        
      end
      
    end
    
    clear filename2load
    clear timestamp
    clear var_data
    clear pspobj
  end
  
end

irf.log('notice','Done loading. Merging the data...');


for i_var = 1:nVar
  
  varname_out = varnamesout{i_var};
  indEmptyFiles = false(length(varData{i_var}),1);
  for i=1:length(varData{i_var})
    if isempty(varData{i_var}{i})
      indEmptyFiles(i) = true;
    end
  end
  var_data = cell2mat(varData{i_var}(~indEmptyFiles));
  epoch_data = cell2mat(epochData{i_var}(~indEmptyFiles));
  varOut =  TSeries(EpochTT(epoch_data),var_data,...
    'metadata_from',varMeta{i_var});
  assignin('base',varname_out,varOut);
  
end
