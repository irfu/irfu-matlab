function [ outFileName ] = mms_sdc_sdp_cdf_writing_2( HeaderInfo )
% MMS_SDC_SDP_CDF_WRITING writes the data to the corresponding CDF file.
%
%	filename_output = MMS_SDC_SDP_CDF_WRITING_2( HeaderInfo)
%   will write an MMS CDF file containing the data stored to a temporary 
%   output folder defined by ENVIR.DROPBOX_ROOT. HeaderInfo contains start 
%   time as well as information about source files ("Parents").
%
%   Example:
%   filename_output = mms_sdc_sdp_cdf_writing_2(HeaderInfo);
%
%	Note 1: It is assumed that other SDC processing scripts will move the
%   created output file to its final destination (from /ENIVR.DROPBOX_ROOT/
%   to /path/as/defined/by/	SDC Developer Guide, v1.7).
%
% 	See also MMS_SDC_SDP_CDF_IN_PROCESS, MMS_SDC_SDP_INIT.

% Verify that we have all information requried.
narginchk(1,1);

global ENVIR;
global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
global DATAC; % Simply recall all data from memory.

instrumentId = 'sdp';  scId = DATAC.scId;
procId = DATAC.procId; procName = MMS_CONST.SDCProcs{procId};

% NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS
% CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING
% DIRECTORY WHEN EXECUTING.
oldDir = pwd; cd(ENVIR.DROPBOX_ROOT);
outFileName = get_file_name();
irf.log('notice',['Writing to DROPBOX_ROOT/',outFileName,'.cdf']);


GATTRIB = getGlobalAttributes;
VATTRIB = getVariableAttributes;
if(HeaderInfo.numberOfSources==1)
  GATTRIB.Parents = {['CDF>',HeaderInfo.parents_1]};
elseif(HeaderInfo.numberOfSources==2)
  GATTRIB.Parents = {['CDF>',HeaderInfo.parents_1]; ['CDF>',HeaderInfo.parents_2]};
end
GATTRIB.Logical_file_id = {outFileName};

switch procId
  case {MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql}
    % Create an almost empty cdf file from skeleton, (have properly
    % formatted LABL_1 and static VATTIB and GATTRIB).
    skel = [ENVIR.CDF_BASE, filesep, 'bin', filesep, 'skeletoncdf -cdf ',...
      pwd, filesep, outFileName,' ', which('mms_sdp_sitl_dce2d.skt')];
    [status, mesg] = system(skel);
    if(status)
      % Error in creating CDF file, could be caused by existing file or
      % skeletoncdf command not found or read/write permission issues or
      % something else.
      errStr=['Error in creating CDF from skeleton. ',mesg];
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:SKELETON', errStr);
    else
      irf.log('notice', mesg);
    end
    % FIXME DUMMY DATA FOR NOW.
    % For now store data temporarly
    epochTT = num2cell(DATAC.dce.time(1:10));
    data1(:,1) = DATAC.dce.e12.data(1:10);
    data1(:,2) = DATAC.dce.e34.data(1:10);
    data1(:,3) = DATAC.dce.e56.data(1:10);
    pgse = num2cell(data1,2);
    dsl = num2cell(data1,2);
    bitmask = num2cell(uint16(DATAC.dce.e12.bitmask(1:10)));
    
    name.epoch   = sprintf('mms%i_%s_dce_epoch',scId,instrumentId);
    name.pgse    = sprintf('mms%i_%s_dce_xyz_pgse',scId,instrumentId);
    name.dsl     = sprintf('mms%i_%s_dce_xyz_dsl',scId,instrumentId);
    name.bitmask = sprintf('mms%i_%s_dce_bitmask',scId,instrumentId);
    name.label   = 'LABL_1';
    %label = {'DCE_X';'DCE_Y';'DCE_Z'};
    
    % Update VariableAttributes
    VATTRIB.CATDESC = {name.epoch, 'Time tags, UTC in TT2000'; ...
      name.pgse,    'DC E field in PGSE frame of reference';...
      name.dsl,     'DC E field in DSL frame of reference'; ...
      name.bitmask, 'Bitmask of signal.'};%...
      %name.label,   'Label'};
    VATTRIB.DEPEND_0 = {name.pgse, name.epoch; ...
      name.dsl,     name.epoch; ...
      name.bitmask, name.epoch};
    VATTRIB.DISPLAY_TYPE = {name.pgse, 'time_series'; ...
      name.dsl,     'time_series'};
    VATTRIB.FIELDNAM = {name.epoch, 'Time tags'; ...
      name.pgse,    'DC E field (pgse)'; ...
      name.dsl,     'DC E field (dsl)'; ...
      name.bitmask, 'Bitmask'};% ...
      %name.label,   'Label'};
    VATTRIB.FILLVAL = {name.epoch, int64(-9223372036854775808); ...
      name.pgse,    single(-1.0E31); ...
      name.dsl,     single(-1.0E31); ...
      name.bitmask, uint16(hex2dec('FFFF'))};
    VATTRIB.FORMAT = {name.pgse, 'F8.3'; ...
      name.dsl,     'F8.3'; ...
      name.bitmask, 'I7'};% ...
      %name.label,   'A23'};
    VATTRIB.LABL_PTR_1 = {name.pgse, name.label; ...
      name.dsl, name.label};
    VATTRIB.SI_CONVERSION = {name.pgse, '1.0e3>V/m'; ...
      name.dsl,     '1.0e3>V/m'};
    VATTRIB.UNITS = {name.pgse, 'mV/m'; ...
      name.dsl,     'mV/m'; ...
      name.bitmask, 'Bitmask'};
    VATTRIB.VALIDMIN = {name.epoch, int64(-43135816000000); ...
      name.pgse,    single(-60.14); ...
      name.dsl,     single(-60.14); ...
      name.bitmask, uint16(0)};
    VATTRIB.VALIDMAX = {name.epoch, int64(946728067183999999); ...
      name.pgse,    single(60.14); ...
      name.dsl,     single(60.14); ...
      name.bitmask, uint16(hex2dec('FFFE'))};
    VATTRIB.VAR_TYPE = {name.epoch, 'support_data'; ...
      name.pgse,    'data'; ...
      name.dsl,     'data'; ...
      name.bitmask, 'support_data'}; %...
      %name.label,   'metadata'};
    VATTRIB.MONOTON = {name.epoch, 'INCREASE'};
    
    if procId==MMS_CONST.SDCProc.sitl
      % Write to file, No QUALITY for SITL 
       spdfcdfwrite(outFileName, ...
        {name.epoch, epochTT, ...
         name.pgse, pgse, ...
         name.dsl, dsl, ...
         name.bitmask, bitmask}, 'EpochType', {name.epoch}, ...
         'GlobalAttributes', GATTRIB, 'VariableAttributes', VATTRIB, ...
         'WriteMode', 'append');
    else
      % Add Quality.
      quality = num2cell(uint16(mms_sdc_sdp_bitmask2quality('e',DATAC.dce.e12.bitmask)));
      name.quality = sprintf('mms%i_sdp_dce_quality',scId);
      VATTRIB.CATDESC  = [VATTRIB.CATDESC;  {name.quality, 'Bitmask of quality'}];
      VATTRIB.DEPEND_0 = [VATTRIB.DEPEND_0; {name.quality, name.epoch}];
      VATTRIB.FIELDNAM = [VATTRIB.FIELDNAM; {name.quality, 'Bitmask'}];
      VATTRIB.FILLVAL  = [VATTRIB.FILLVAL;  {name.quality, uint16(hex2dec('FFFF'))}];
      VATTRIB.FORMAT   = [VATTRIB.FORMAT;   {name.quality, 'I7'}];
      VATTRIB.UNITS    = [VATTRIB.UNITS;    {name.quality, 'Bitmask'}];
      VATTRIB.VALIDMIN = [VATTRIB.VALIDMIN; {name.quality, uint16(0)}];
      VATTRIB.VALIDMAX = [VATTRIB.VALIDMAX; {name.quality, uint16(hex2dec('FFFE'))}];
      VATTRIB.VAR_TYPE = [VATTRIB.VAR_TYPE; {name.quality, 'support_data'}];
      % Write to file.
      spdfcdfwrite(outFileName, ...
       {name.epoch, epochTT, ...
        name.pgse, pgse, ...
        name.dsl, dsl, ...
        name.bitmask, bitmask, ...
        name.quality, quality}, 'EpochType', {name.epoch}, ...
        'GlobalAttributes', GATTRIB, 'VariableAttributes', VATTRIB, ...
        'WriteMode', 'append');
    end
    
  case MMS_CONST.SDCProc.usc
    %% FIXME: DUMMY DATA FOR NOW.
    skel = [ENVIR.CDF_BASE, filesep, 'bin', filesep, 'skeletoncdf -cdf ',...
      pwd, filesep, outFileName,' ', which('mms_sdp_l2_usc.skt')];
    [status, mesg] = system(skel);
    if(status)
      % Error in creating CDF file, could be caused by existing file or
      % skeletoncdf command not found or read/write permission issues or
      % something else.
      errStr=['Error in creating CDF from skeleton. ',mesg];
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:SKELETON', errStr);
    else
      irf.log('notice', mesg);
    end
    % For now store data temporarly
    epochTT = num2cell(DATAC.dcv.time(1:10));
    psp_p(:,1) = DATAC.dcv.v1.data(1:10);
    psp_p(:,2) = DATAC.dcv.v2.data(1:10);
    psp_p(:,3) = DATAC.dcv.v3.data(1:10);
    psp_p(:,4) = DATAC.dcv.v4.data(1:10);
    psp_p(:,5) = DATAC.dcv.v5.data(1:10);
    psp_p(:,6) = DATAC.dcv.v6.data(1:10);
    psp_p = num2cell(psp_p,2);
    bitmask = num2cell(uint16(DATAC.dcv.v1.bitmask(1:10)));
    ESCP = num2cell(DATAC.dcv.v1.data(1:10));
    PSP = num2cell(DATAC.dcv.v2.data(1:10));
    Delta = num2cell(DATAC.dcv.v3.data(1:10));

    name.epoch   = sprintf('mms%i_%s_epoch_dcv',scId,instrumentId); % Timestamp in TT2000
    name.escp    = sprintf('mms%i_%s_escp_dcv',scId,instrumentId); % Estimated Spacecraft potential
    name.psp     = sprintf('mms%i_%s_psp_dcv',scId,instrumentId); % Probe to spacecraft potential
    name.delta   = sprintf('mms%i_%s_delta_dcv',scId,instrumentId); % Delta
    name.psp_p   = sprintf('mms%i_%s_psp_probes_dcv',scId,instrumentId); % Probe to spacecraft potential, indiv probes
    name.bitmask = sprintf('mms%i_%s_bitmask_dcv',scId,instrumentId); % Bitmask
    name.label   = 'LABL_1';
    %label = 'PSP_P1,PSP_P2,PSP_P3,PSP_P4,PSP_P5,PSP_P6'; NOT correct
    %label = '["PSP_P1","PSP_P2","PSP_P3","PSP_P4","PSP_P5","PSP_P6"]'; % NOT entirely correct either.
    %label = {'PSP_P1','PSP_P2','PSP_P3','PSP_P4','PSP_P5','PSP_P6'}; % NOT correct, 6 records..
    %label = {'"PSP_P1","PSP_P2","PSP_P3","PSP_P4","PSP_P5","PSP_P6"'}; % Cannot determine proper cdf datatype from matlab value
    
    % Update VariableAttributes
    VATTRIB.CATDESC = {name.epoch, 'Time tags, UTC in TT2000'; ...
      name.escp,    'Estimated spacecraft potential';...
      name.psp,     'Probe to spacecraft potential'; ...
      name.delta,   'Delta (between escp and psp)'; ...
      name.psp_p,   'Probe to spacecraft potential, individual probes'; ...
      name.bitmask, 'Bitmask of signal'};% ...
      %name.label,   'Label'};
    VATTRIB.DEPEND_0 = {name.escp, name.epoch; ...
      name.psp,     name.epoch; ...
      name.delta,   name.epoch; ...
      name.psp_p,   name.epoch; ...
      name.bitmask, name.epoch};
    VATTRIB.DISPLAY_TYPE = {name.escp, 'time_series'; ...
      name.psp,     'time_series'; ...
      name.delta,   'time_series'; ...
      name.psp_p,   'time_series'};
    VATTRIB.FIELDNAM = {name.epoch, 'Time tags'; ...
      name.escp,    'Spacecraft potential'; ...
      name.psp,     'Probe to spacecraft potential'; ...
      name.delta,   'Delta'; ...
      name.psp_p,   'Probe to spacecraft potential individual probe'; ...
      name.bitmask, 'Bitmask'};% ...
      %name.label,   'Label'};
    VATTRIB.FILLVAL = {name.epoch, int64(-9223372036854775808); ...
      name.escp,    single(-1.0E31); ...
      name.psp,     single(-1.0E31); ...
      name.delta,   single(-1.0E31); ...
      name.psp_p,   single(-1.0E31); ...
      name.bitmask, uint16(hex2dec('FFFF'))};
    VATTRIB.FORMAT = {name.escp, 'F8.3'; ...
      name.psp,     'F8.3'; ...
      name.delta,   'F8.3'; ...
      name.psp_p,   'F8.3'; ...
      name.bitmask, 'I7'};% ...
      %name.label,   'A23'};
    VATTRIB.LABLAXIS = {name.escp, 'Spacecraft potential';
      name.psp,     'Probe to spacecraft potential';
      name.delta,   'Delta'};
    VATTRIB.LABL_PTR_1 = {name.psp_p, name.label};
    VATTRIB.SI_CONVERSION = {name.escp, '1.0>V'; ...
      name.psp,     '1.0>V'; ...
      name.delta,   '1.0>V'; ...
      name.psp_p,   '1.0>V'};
    VATTRIB.UNITS = {name.epoch, 'ns'; ...
      name.escp,    'V'; ...
      name.psp,     'V'; ...
      name.delta,   'V'; ...
      name.psp_p,   'V'; ...
      name.bitmask, 'Bitmask'};
    VATTRIB.VALIDMIN = {name.epoch, int64(-43135816000000); ...
      name.escp,    single(-60.14); ...
      name.psp,     single(-60.14); ...
      name.delta,   single(-60.14); ...
      name.psp_p,   single(-60.14); ...
      name.bitmask, uint16(0)};
    VATTRIB.VALIDMAX = {name.epoch, int64(946728067183999999); ...
      name.escp,    single(60.14); ...
      name.psp,     single(60.14); ...
      name.delta,   single(60.14); ...
      name.psp_p,   single(60.14); ...
      name.bitmask, uint16(hex2dec('FFFE'))};
    VATTRIB.VAR_TYPE = {name.epoch, 'support_data'; ...
      name.escp,    'data'; ...
      name.psp,     'data'; ...
      name.delta,   'data'; ...
      name.psp_p,   'data'; ...
      name.bitmask, 'support_data'};% ...
      %name.label,   'metadata'};
    VATTRIB.MONOTON = {name.epoch, 'INCREASE'};
    
    % Write to file, using new spdfcdfwrite.
    spdfcdfwrite(outFileName, {name.epoch, epochTT, ...
      name.escp, ESCP, ...
      name.psp, PSP, ...
      name.delta, Delta, ...
      name.psp_p, psp_p, ...
      name.bitmask, bitmask}, 'EpochType', {name.epoch}, ...
      'GlobalAttributes', GATTRIB, 'VariableAttributes', VATTRIB, ...
      'WriteMode', 'append');
  
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

% Return to previous working directory.
cd(oldDir);

  function fileName = get_file_name
    % Generate output file name incrementing the file version if necessary
    switch procId
      case {MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql}
        subDir = procName; suf = 'dce2d';
      case MMS_CONST.SDCProc.usc
        subDir = 'l2'; suf = 'uscdcv';
      otherwise
        errStr = 'unrecognized procId';
        irf.log('critical', errStr); error(errStr)
    end
    scIdStr = sprintf('mms%d',scId);
    tmMode = DATAC.tmMode; tmModeStr = MMS_CONST.TmModes{tmMode};
    startTime =  HeaderInfo.startTime;
    verStr = sprintf('%d.%d.',MMS_CONST.Version.X,MMS_CONST.Version.Y);
    fileName = [scIdStr '_' instrumentId, '_' tmModeStr '_' subDir '_' ...
      suf '_' startTime '_v' ];
    
    % Check for preexisting files and increment file version
    dataPathPref = [ENVIR.DATA_PATH_ROOT, filesep,'science',filesep, ...
      scIdStr, filesep, instrumentId, filesep, tmModeStr, filesep, ...
      subDir, filesep, startTime(1:4), filesep, startTime(5:6), filesep];
    dataPathPref = [dataPathPref,startTime(7:8), filesep];
    
    preExistingFiles = dir([dataPathPref fileName verStr '*.cdf']);
    if numel(preExistingFiles)
      maxRev = 0;
      for iFile = 1:numel(preExistingFiles)
        rev = get_rev(preExistingFiles(iFile).name);
        if rev>maxRev, maxRev = rev; end
      end
      newVer = maxRev + 1;
    else newVer = 0;
    end
    verStr = [verStr num2str(newVer)];
    fileName = [fileName verStr];
    
    function r = get_rev(s)
      % Find revision (Z) from version string in a file name xxx_vX.Y.Z.cdf
      idxDot = find(s=='.');
      if numel(idxDot)~=3
        irf.log('warning',['Bad file name: ' s])
        r = 0; return
      end
      r = str2double(s(idxDot(2)+1:idxDot(3)-1));
    end % GET_REV
  end % get_file_name

  function GATTRIB = getGlobalAttributes
    %% Create GlobalAttribute struct
    % Source: MMS_CDF_Format_Guide, v1.5, issued 2013/11/19.
    GATTRIB=[];
    % Global Attributes REQUIRED:
%    GATTRIB.Data_type = cell(0,1); % mode_dataLevel_optionalDescriptor
%    GATTRIB.Data_version = cell(0,1); % Same as version number in filename.
%    GATTRIB.Descriptor = {'SDP>Spin-plane Double Probe'};
%    GATTRIB.Discipline = {'Space Physics>Magnetospheric Science'};
    GATTRIB.Generation_date = {datestr(now,'yyyymmdd')};
%    GATTRIB.Instrument_type = {'Electric Fields (space)'};
%    GATTRIB.Logical_file_id = cell(0,1); % Same as filename without ".cdf".
%    GATTRIB.Logical_source = cell(0,1); % Ex: mms3_sdp_fast_ql_swd (mmsSC_instrument_mode_dataLevel_optionalDescriptor)
%    GATTRIB.Logical_source_description = cell(0,1); % Description in full words..
%    GATTRIB.Mission_group = {'MMS'};
%    GATTRIB.PI_affiliation = {'SWRI, LASP, KTH'};
%    GATTRIB.PI_name = {'Burch, J, Ergun, R., Lindqvist, P.'};
%    GATTRIB.Project = {'STP>Solar-Terrestrial Physics'};
    GATTRIB.Source_name = {sprintf('MMS%i>MMS Satellite Number %i',scId,scId)}; % Or possibly 'MMS>MMS Constellation'.
%    GATTRIB.TEXT = cell(0,1);  % FIXME This attribute is an SPDF standard global attribute, which is a text description of the
      %experiment whose data is included in the CDF. A reference to a journal article(s) or to a
      %World Wide Web page describing the experiment is essential, and constitutes the
      %minimum requirement. A written description of the data set is also desirable. This
      %attribute can have as many entries as necessary to contain the desired information.
      %Typically, this attribute is about a paragraph in length and is not shown on CDAWeb.
%    GATTRIB.HTTP_LINK = {'http://mms.gsfc.nasa.gov/'}; % FIXME should point to data
%    GATTRIB.LINK_TEXT = {'Magnetospheric Multiscale (MMS) Mission - NASA'}; % FIXME as well
%    GATTRIB.LINK_TITLE = {'Magnetospheric Multiscale (MMS) Mission - NASA'}; % FIXME as well
%    GATTRIB.MODS = cell(0,1);
    % Global Attributes RECOMMENDED:
%    GATTRIB.Acknowledgement = cell(0,1);
    GATTRIB.Generated_by = {['IRFU Matlab', irf('version')]};
    % Global Attributes OPTIONAL:
%    GATTRIB.Parents = cell(0,1);
%    GATTRIB.Skeleton_version = cell(0,1);
%    GATTRIB.Rules_of_use = cell(0,1);
%    GATTRIB.Time_resolution = cell(0,1);
  end

  function VATTRIB = getVariableAttributes
    %% Create VariableAttribute struct
    VATTRIB=[];
    % Variable attributes, REQUIRED as specified.
    VATTRIB.CATDESC = cell(0,1);        % Req for data, support_data and meta_data.
    VATTRIB.DEPEND_0 = cell(0,1);       % Req for data, support_data and meta_data (if timevarying).
    VATTRIB.DEPEND_1 = cell(0,1);         % Req if variable is not using DEPEND_0.
    VATTRIB.DISPLAY_TYPE = cell(0,1);   % Req for data.
    VATTRIB.FIELDNAM = cell(0,1);       % Req for data, support_data, meta_data.
    VATTRIB.FILLVAL = cell(0,1);        % Req for data, support_data, meta_data (if timevarying).
    VATTRIB.FORMAT = cell(0,1);         % Req for data, support_data, meta_data.
    VATTRIB.FORM_PTR = cell(0,1);         % Req if variable is not using FORMAT.
    VATTRIB.LABLAXIS = cell(0,1);       % Req for data.
    VATTRIB.LABL_PTR_1 = cell(0,1);       % Req if variable is not using LABLAXIS.
    VATTRIB.SI_CONVERSION = cell(0,1);  % Req for data.
    VATTRIB.UNITS = cell(0,1);          % Req for data, support_data.
    VATTRIB.UNIT_PTR = cell(0,1);         % Req if variable is not using UNITS.
    VATTRIB.VALIDMIN = cell(0,1);       % Req for data, support_data (if timevarying).
    VATTRIB.VALIDMAX = cell(0,1);       % Req for data, support_data (if timevarying).
    VATTRIB.VAR_TYPE = cell(0,1);       % valid values: 'data', 'support_data', 'metadata' and 'ignore_data'.
    % Others, used in source CDF files from UNH
    VATTRIB.MONOTON = cell(0,1); % For time variables valid values: 'INCREASE' and 'DECREASE', can speed up reading of cdf file.
  end

end
