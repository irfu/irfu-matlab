function [ outFileName ] = mms_sdc_sdp_cdf_write( HeaderInfo )
% MMS_SDC_SDP_CDF_WRITING writes the data to the corresponding CDF file.
%
%	filename_output = MMS_SDC_SDP_CDF_WRITE( HeaderInfo)
%   will write an MMS CDF file containing the data stored to a temporary 
%   output folder defined by ENVIR.DROPBOX_ROOT. HeaderInfo contains start 
%   time as well as information about source files ("Parents").
%
%   Example:
%   filename_output = mms_sdc_sdp_cdf_write(HeaderInfo);
%
%	Note 1: It is assumed that other SDC processing scripts will move the
%   created output file to its final destination (from /ENIVR.DROPBOX_ROOT/
%   to /path/as/defined/by/	SDC Developer Guide, v1.7).
%
% 	See also MMS_SDC_SDP_LOAD, MMS_SDC_SDP_INIT.

% Verify that we have all information requried.
narginchk(1,1);

global ENVIR;
global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

INST_NAME = 'edp'; % Electric double probe
DCE_FILE = 'dce2d';
EFIELD_MAX = single(700); % Max value of E-field in mV/m with shortening factor applied
VOLTAGE_MIN = single(-120); % Min voltage
VOLTAGE_MAX = single(50); % Max voltage
QUALITY_MAX = int16(4);  % Max value of quality.
COMPRESS_LEVEL = 'gzip.6'; % Default compression level to be used for variables.

procId = mms_sdc_sdp_datamanager('procId');
if procId==MMS_CONST.Error
    errStr = 'mms_sdc_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
end
procName = MMS_CONST.SDCProcs{procId};
scId = mms_sdc_sdp_datamanager('scId');
tmMode = mms_sdc_sdp_datamanager('tmMode'); 
tmModeStr = MMS_CONST.TmModes{tmMode};
datasetPrefix = sprintf('mms%i_%s',scId,INST_NAME);

% NOTE MOVE TO DROPBOX FOLDER BEFORE TRYING TO WRITE ANYTHING AS
% CDF MAY TRY TO WRITE TEMPORARY FILES IN THE CURRENT WORKING
% DIRECTORY WHEN EXECUTING.
oldDir = pwd; cd(ENVIR.DROPBOX_ROOT);
[outFileName, verFileName] = get_file_name();
irf.log('notice',['Writing to DROPBOX_ROOT/',outFileName,'.cdf']);

GATTRIB = getGlobalAttributes;
VATTRIB = getVariableAttributes;
% Update some dynamic GlobalAttributes common to all data products.
if(HeaderInfo.numberOfSources==1)
  GATTRIB.Parents = {['CDF>',HeaderInfo.parents_1]}; % Add the only parent
elseif(HeaderInfo.numberOfSources>=2)
  GATTRIB.Parents = {['CDF>',HeaderInfo.parents_1]}; % Add first parent
  for i=2:HeaderInfo.numberOfSources % Add each extra parent source
    GATTRIB.Parents = [GATTRIB.Parents; {['CDF>',...
      HeaderInfo.(sprintf('parents_%i',i))]}];
  end
end
GATTRIB.Logical_file_id = {outFileName};    % Filename, except '.cdf'
GATTRIB.Data_version = {['v' verFileName]}; % 'vX.Y.Z'

switch procId
  case {MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql, ...
      MMS_CONST.SDCProc.l2a}
    %% SITL/QL/L2A DCE2D - get data
    % QL/SITL/L2A outputs: DCE_dsl [E_x, E_y, E_z], and for QL & L2A also
    % Bitmask and Quality.
    % (and corresponding epoch and labels)
    dataType = [tmModeStr '_' procName '_' DCE_FILE];
    dataDesc = sprintf(...
      'MMS %i dual probe %s (%s), two dimensional electric field.',...
      scId,procName,tmModeStr);

    if(procId==MMS_CONST.SDCProc.l2a)
      % Replace GATTRIB.Parents with source cdf (l2pre) parents (ie only
      % official data products listed as parents not internal (FIELDS) file.
      l2pre = mms_sdc_sdp_datamanager('l2pre');
      if( ~isstruct(l2pre) && l2pre == MMS_CONST.Error)
        errStr = 'Cannot output ''l2pre''';
        irf.log('critical', errStr);
        error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
      end
      GATTRIB.Parents = l2pre.dataObj.GlobalAttributes.Parents;
    end
    
    dce_xyz_dsl = mms_sdc_sdp_datamanager('dce_xyz_dsl');
    if ~isstruct(dce_xyz_dsl) && dce_xyz_dsl == MMS_CONST.Error
      errStr = 'Cannot output ''dce_xyz_dsl''';
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end

    epochTT = dce_xyz_dsl.time;
    dsl = dce_xyz_dsl.data;
    % Bitmask, defined as CDF_UINT1 (uint8 in Matlab)
    bitmask = uint8(dce_xyz_dsl.bitmask);
    
    name.epoch   = [datasetPrefix '_dce_epoch'];
    name.dsl     = [datasetPrefix '_dce_xyz_dsl'];
    name.bitmask = [datasetPrefix '_dce_bitmask'];
    name.label   = 'LABL_1';
    label = ['DSL_X'; 'DSL_Y'; 'DSL_Z'];
    outVars = {name.epoch, epochTT, name.label, label, name.dsl, dsl};
    % RecordBound, ie individual cdf records for each row.
    recBound = {name.epoch, name.dsl};
    %Variable Datatypes
    varDatatype = {name.epoch, 'cdf_time_tt2000', name.label, 'cdf_char', ...
      name.dsl,'cdf_real4'};
    % Compression level for each variable. Note: Not for TT2000 epoch.
    compressVars = {name.label, COMPRESS_LEVEL, name.dsl, COMPRESS_LEVEL};

    if procId~=MMS_CONST.SDCProc.sitl
      % Add Quality and Bitmask for QL or L2a
      % Quality, defined as CDF_INT2 (int16 in Matlab)
      quality = int16(mms_sdc_sdp_bitmask2quality('e',dce_xyz_dsl.bitmask));
      name.quality = [datasetPrefix '_dce_quality'];
      outVars = [outVars {name.bitmask bitmask name.quality quality}];
      recBound = [recBound {name.bitmask name.quality}];
      varDatatype = [varDatatype {name.bitmask, 'cdf_uint1', name.quality, 'cdf_int2'}];
      compressVars = [compressVars {name.bitmask, COMPRESS_LEVEL, name.quality, COMPRESS_LEVEL}];
    end

    %% Update VariableAttributes
    VATTRIB.CATDESC = {name.epoch, 'Time tags, in TT2000'; ...
      name.label,   'Label'; ...
      name.dsl,     'DC E field in DSL frame of reference'};
    VATTRIB.DEPEND_0 = {name.dsl,     name.epoch};
    VATTRIB.DISPLAY_TYPE = {name.dsl,     'time_series'};
    VATTRIB.FIELDNAM = {name.epoch, 'Time tags'; ...
      name.label,   'Label'; ...
      name.dsl,     'DC E field (dsl)'};
    VATTRIB.FILLVAL = {name.epoch, int64(-9223372036854775808); ...
      name.dsl,     single(-1.0E31)};
    VATTRIB.FORMAT = {name.label,   'A23'; ...
      name.dsl,     'F8.3'};
    VATTRIB.LABL_PTR_1 = {name.dsl, name.label};
    VATTRIB.SI_CONVERSION = {name.dsl,     '1.0e3>V/m'};
    VATTRIB.UNITS = {name.dsl,     'mV/m'};
    VATTRIB.VALIDMIN = {name.epoch, spdfcomputett2000([1990,01,01,0,0,0,0,0,0]); ...
      name.dsl,     -EFIELD_MAX};
    VATTRIB.VALIDMAX = {name.epoch, spdfcomputett2000([2100,01,01,0,0,0,0,0,0]); ...
      name.dsl,     EFIELD_MAX};
    VATTRIB.VAR_TYPE = {name.epoch, 'support_data'; ...
      name.label,   'metadata'; ...
      name.dsl,     'data'};
    VATTRIB.MONOTON = {name.epoch, 'INCREASE'};
    if procId~=MMS_CONST.SDCProc.sitl
      % Add Quality and Bitmask for QL or L2PRE
      VATTRIB.CATDESC  = [VATTRIB.CATDESC; {...
        name.bitmask, 'Status bitmask';...
        name.quality, 'Quality indicator'}];
      VATTRIB.DEPEND_0 = [VATTRIB.DEPEND_0; {...
        name.bitmask, name.epoch;...
        name.quality, name.epoch}];
      VATTRIB.FIELDNAM = [VATTRIB.FIELDNAM; {...
        name.bitmask, 'Status bitmask';...
        name.quality, 'Quality indicator'}];
      VATTRIB.FILLVAL  = [VATTRIB.FILLVAL;  {...
        name.bitmask, uint8(255);...
        name.quality, int16(-32768)}];
      VATTRIB.FORMAT   = [VATTRIB.FORMAT;   {...
        name.bitmask, 'I7';...
        name.quality, 'I7'}];
      VATTRIB.UNITS    = [VATTRIB.UNITS;    {...
        name.bitmask, 'unitless';...
        name.quality, 'unitless'}];
      VATTRIB.VALIDMIN = [VATTRIB.VALIDMIN; {...
         name.bitmask, uint8(0);...
         name.quality, int16(-32767)}];
      VATTRIB.VALIDMAX = [VATTRIB.VALIDMAX; {...
        name.bitmask, uint8(254);...
        name.quality, QUALITY_MAX}];
      VATTRIB.VAR_TYPE = [VATTRIB.VAR_TYPE; {...
        name.bitmask, 'support_data';...
        name.quality, 'support_data'}];
    end

  case MMS_CONST.SDCProc.l2pre
    dataType = [tmModeStr '_' procName '_' DCE_FILE];
    dataDesc = sprintf(...
      'MMS %i dual probe %s (%s), two dimensional electric field. L2Pre file',...
      scId,procName,tmModeStr);

    % Verify output data exist. L2Pre outputs: DCE [e12, e34, e56],
    % Spinfits [e12, e34], bitmask of DCE [e12, e34, e56] and a quality.
    % (and corresponding epoch and labels).
    dce = mms_sdc_sdp_datamanager('dce');
    if ~isstruct(dce) && dce == MMS_CONST.Error
      errStr = 'Cannot output ''dce''';
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end
    spinfits = mms_sdc_sdp_datamanager('spinfits');
    if ~isstruct(spinfits) && spinfits == MMS_CONST.Error
      errStr = 'Cannot output ''spinfits''';
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end

    epochTT = dce.time;
    dcedata = [dce.e12.data, dce.e34.data, dce.e56.data];
    % Bitmask, defined as CDF_UINT1 (uint8 in Matlab)
    bitmask = [uint8(dce.e12.bitmask), uint8(dce.e34.bitmask), uint8(dce.e56.bitmask)];
    quality = int16(mms_sdc_sdp_bitmask2quality('e',bitmask(:,1)));

    name.epoch   = [datasetPrefix '_dce_epoch'];
    name.dce     = [datasetPrefix '_dce_data'];
    name.bitmask = [datasetPrefix '_dce_bitmask'];
    name.quality = [datasetPrefix '_dce_quality'];
    name.label   = 'LABL_1';
    name.sfitsEpoch = [datasetPrefix '_dce_spinfits_epoch'];
    name.label2 = 'LABL_2';
    label = ['E_12'; 'E_34'; 'E_56'];

    outVars = {name.epoch, epochTT, name.label, label, name.dce, dcedata, ...
      name.bitmask, bitmask, name.quality, quality, ...
      name.sfitsEpoch, spinfits.time};
    % RecordBound, ie individual cdf records for each row.
    recBound = {name.epoch, name.dce, name.bitmask, name.quality, name.sfitsEpoch};
    %Variable Datatypes
    varDatatype = {name.epoch, 'cdf_time_tt2000', name.label, 'cdf_char', ...
      name.dce,'cdf_real4', name.bitmask,'cdf_uint1', name.quality,'cdf_int2',...
      name.sfitsEpoch, 'cdf_time_tt2000'};
    % Compression level for each variable. Note: Not for TT2000 epoch.
    compressVars = {name.label, COMPRESS_LEVEL, name.dce, COMPRESS_LEVEL,...
      name.bitmask, COMPRESS_LEVEL, name.quality, COMPRESS_LEVEL};

    %% Update VariableAttributes
    VATTRIB.CATDESC = {name.epoch, 'Time tags, in TT2000'; ...
      name.label,      'Label'; ...
      name.dce,        'DC E field in SC probe pairs'; ...
      name.bitmask,    'Status bitmask';...
      name.quality,    'Quality indicator';...
      name.sfitsEpoch, 'Time tags, in TT2000'; ...
      name.label2,     'Label'};
    VATTRIB.DEPEND_0 = {name.dce, name.epoch; ...
      name.bitmask, name.epoch;...
      name.quality, name.epoch};
    VATTRIB.DISPLAY_TYPE = {name.dce,     'time_series'};
    VATTRIB.FIELDNAM = {name.epoch, 'Time tags'; ...
      name.label,      'Label'; ...
      name.dce,        'DC E field (probes)'; ...
      name.bitmask,    'Status bitmask'; ...
      name.quality,    'Quality indicator';...
      name.sfitsEpoch, 'Time tags'; ...
      name.label2,     'Label'};
    VATTRIB.FILLVAL = {name.epoch, int64(-9223372036854775808); ...
      name.dce,        single(-1.0E31);...
      name.bitmask,    uint8(255);...
      name.quality,    int16(-32768);...
      name.sfitsEpoch, int64(-9223372036854775808)};
    VATTRIB.FORMAT = {name.label,   'A23'; ...
      name.dce,     'F8.3';...
      name.bitmask, 'I7';...
      name.quality, 'I7';...
      name.label2,  'A23'};
    VATTRIB.LABL_PTR_1 = {name.dce, name.label; name.bitmask, name.label};
    VATTRIB.SI_CONVERSION = {name.dce,     '1.0e3>V/m'};
    VATTRIB.UNITS = {name.dce,     'mV/m';...
      name.bitmask, 'unitless';...
      name.quality, 'unitless'};
    VATTRIB.VALIDMIN = {name.epoch, spdfcomputett2000([1990,01,01,0,0,0,0,0,0]); ...
      name.dce,     -EFIELD_MAX;...
      name.bitmask, uint8(0);...
      name.quality, int16(-32767);...
      name.sfitsEpoch, spdfcomputett2000([1990,01,01,0,0,0,0,0,0])};
    VATTRIB.VALIDMAX = {name.epoch, spdfcomputett2000([2100,01,01,0,0,0,0,0,0]); ...
      name.dce,        EFIELD_MAX;...
      name.bitmask,    uint8(254);...
      name.quality,    QUALITY_MAX;...
      name.sfitsEpoch, spdfcomputett2000([2100,01,01,0,0,0,0,0,0])};
    VATTRIB.VAR_TYPE = {name.epoch, 'support_data'; ...
      name.label,      'metadata'; ...
      name.dce,        'data';...
      name.bitmask,    'support_data';...
      name.quality,    'support_data';...
      name.sfitsEpoch, 'support_data'; ...
      name.label2,     'metadata'};
    VATTRIB.MONOTON = {name.epoch, 'INCREASE';...
      name.sfitsEpoch, 'INCREASE'};

    sdpPair = fields(spinfits.sfit); % Defailt e12 & e34
    switch (size(spinfits.sfit.(sdpPair{1}),2))
      case 3
        % A, B & C, 15 char each.
        label2 = ['Fit y=A        '; 'Fit y=B*sin(wt)'; 'Fit y=C*cos(wt)'];
      case 5
        % A, B, C, D & E, 16 char each.
        label2 = ['Fit y=A         '; 'Fit y=B*sin(wt) ';...
          'Fit y=C*cos(wt) '; 'Fit y=D*sin(2wt)'; 'Fit y=E*cos(2wt)'];
      case 7
        % A, B, C, D, E, F & G, 16 char each.
        label2 = ['Fit y=A         '; 'Fit y=B*sin(wt) ';...
          'Fit y=C*cos(wt) '; 'Fit y=D*sin(2wt)'; 'Fit y=E*cos(2wt)';...
          'Fit y=F*sin(3wt)'; 'Fit y=G*cos(3wt)'];
      otherwise
        errStr = 'Number of terms in spinfit incorrect, must be 3, 5 or 7.';
        irf.log('critical', errStr);
        error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end
    outVars = [outVars {name.label2, label2}];
    varDatatype = [varDatatype {name.label2, 'cdf_char'}];
    compressVars = [compressVars {name.label2, COMPRESS_LEVEL}];

    for iPair=1:numel(sdpPair)
      name.sfits.(sdpPair{iPair}) = [datasetPrefix '_dce_spinfits_' sdpPair{iPair}];
      name.sdev.(sdpPair{iPair}) = [datasetPrefix '_dce_spinfits_sdev_' sdpPair{iPair}];
      % Update output variables
      outVars = [outVars {...
        name.sfits.(sdpPair{iPair}), spinfits.sfit.(sdpPair{iPair}),...
        name.sdev.(sdpPair{iPair}),  spinfits.sdev.(sdpPair{iPair})}];
      recBound = [recBound {...
        name.sfits.(sdpPair{iPair}), name.sdev.(sdpPair{iPair})}];
      varDatatype = [varDatatype {...
        name.sfits.(sdpPair{iPair}), 'cdf_real4', ...
        name.sdev.(sdpPair{iPair}),  'cdf_real4'}];
      compressVars = [compressVars {...
        name.sfits.(sdpPair{iPair}), COMPRESS_LEVEL, ...
        name.sdev.(sdpPair{iPair}),  COMPRESS_LEVEL}];
      % Update VATTRIB
      VATTRIB.CATDESC = [VATTRIB.CATDESC; {...
        name.sfits.(sdpPair{iPair}), ['Spinfit coefficients probes ', sdpPair{iPair}];...
        name.sdev.(sdpPair{iPair}), ['Standard deviation of spinfit coefficients probes ', sdpPair{iPair}]}];
      VATTRIB.DEPEND_0 = [VATTRIB.DEPEND_0; {...
        name.sfits.(sdpPair{iPair}), name.sfitsEpoch;...
        name.sdev.(sdpPair{iPair}),  name.sfitsEpoch}];
      VATTRIB.DISPLAY_TYPE = [VATTRIB.DISPLAY_TYPE; {...
        name.sfits.(sdpPair{iPair}), 'time_series'; ...
        name.sdev.(sdpPair{iPair}),  'time_series'}];
      VATTRIB.FIELDNAM = [VATTRIB.FIELDNAM; {...
        name.sfits.(sdpPair{iPair}), ['Spinfits', sdpPair{iPair}];...
        name.sdev.(sdpPair{iPair}),  ['Sdev spinfits ', sdpPair{iPair}]}];
      VATTRIB.FILLVAL = [VATTRIB.FILLVAL; {...
        name.sfits.(sdpPair{iPair}), single(-1.0E31);...
        name.sdev.(sdpPair{iPair}),  single(-1.0E31)}];
      VATTRIB.FORMAT = [VATTRIB.FORMAT; {...
        name.sfits.(sdpPair{iPair}), 'F8.3'; ...
        name.sdev.(sdpPair{iPair}),  'F8.3'}];
      VATTRIB.LABL_PTR_1 = [VATTRIB.LABL_PTR_1; {...
        name.sfits.(sdpPair{iPair}), name.label2}];
      VATTRIB.SI_CONVERSION = [VATTRIB.SI_CONVERSION; {...
        name.sfits.(sdpPair{iPair}), '1.0e3>V/m';...
        name.sdev.(sdpPair{iPair}),  '1.0e3>V/m'}];
      VATTRIB.UNITS = [VATTRIB.UNITS; {...
        name.sfits.(sdpPair{iPair}), 'mV/m'; ...
        name.sdev.(sdpPair{iPair}),  'mV/m'}];
      VATTRIB.VALIDMIN = [VATTRIB.VALIDMIN; {...
        name.sfits.(sdpPair{iPair}), -EFIELD_MAX;...
        name.sdev.(sdpPair{iPair}),  -EFIELD_MAX}];
      VATTRIB.VALIDMAX = [VATTRIB.VALIDMAX; {...
        name.sfits.(sdpPair{iPair}), EFIELD_MAX;...
        name.sdev.(sdpPair{iPair}),  EFIELD_MAX}];
      VATTRIB.VAR_TYPE = [VATTRIB.VAR_TYPE; {...
        name.sfits.(sdpPair{iPair}), 'data';...
        name.sdev.(sdpPair{iPair}),  'support_data'}];
    end

  case MMS_CONST.SDCProc.scpot
    %% ScPot - get data
    dataType = [tmModeStr '_l2_scpot'];
    dataDesc = sprintf(...
      'MMS %i dual probe %s (%s), Spacecraft potential',...
      scId,procName,tmModeStr);
    
    dcv = mms_sdc_sdp_datamanager('dcv');
    if ~isstruct(dcv) && dcv == MMS_CONST.Error
      errStr = 'Cannot output ''dcv''';
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end
    probe2sc_pot = mms_sdc_sdp_datamanager('probe2sc_pot');
    if ~isstruct(probe2sc_pot) && probe2sc_pot == MMS_CONST.Error
      errStr = 'Cannot output ''probe2sc_pot''';
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end
    sc_pot = mms_sdc_sdp_datamanager('sc_pot');
    if ~isstruct(sc_pot) && sc_pot == MMS_CONST.Error
      errStr = 'Cannot output ''sc_pot''';
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_CDFWRITE:OUT', errStr);
    end
    
    epochTT = dcv.time;
    psp_p = [dcv.v1.data, dcv.v2.data, dcv.v3.data, ...
      dcv.v4.data, dcv.v5.data, dcv.v6.data];
    PSP = probe2sc_pot.data;
    ESCP = sc_pot.data;
    bitmask = uint8(sc_pot.bitmask);
    quality = int16(mms_sdc_sdp_bitmask2quality('e',sc_pot.bitmask));

    name.epoch   = [datasetPrefix '_scpot_epoch']; % Timestamp in TT2000
    name.scpot   = [datasetPrefix '_scpot']; % Estimated Spacecraft potential
    name.psp     = [datasetPrefix '_psp']; % Probe to spacecraft potential, averaged
    name.psp_p   = [datasetPrefix '_dcv']; % Probe to spacecraft potential, indiv probes
    name.bitmask = [datasetPrefix '_scpot_bitmask']; % Bitmask
    name.quality = [datasetPrefix '_scpot_quality']; % Quality
    name.label   = 'LABL_1';
    label = ['PSP_P1'; 'PSP_P2'; 'PSP_P3'; 'PSP_P4'; 'PSP_P5'; 'PSP_P6'];
    
    outVars = {name.epoch, epochTT, ...
      name.label, label, ...
      name.scpot, ESCP, ...
      name.psp, PSP, ...
      name.psp_p, psp_p, ...
      name.bitmask, bitmask,...
      name.quality, quality};

    % RecordBound, ie individual cdf records for each row.
    recBound = {name.epoch, name.scpot, name.psp, name.psp_p, ...
      name.bitmask, name.quality};

    % Variable Datatypes
    varDatatype = {name.epoch, 'cdf_time_tt2000', name.label, 'cdf_char', ...
      name.scpot, 'cdf_real4', name.psp, 'cdf_real4', name.psp_p, 'cdf_real4',...
      name.bitmask, 'cdf_uint1', name.quality, 'cdf_int2'};

    % Compression
    compressVars = {name.label, COMPRESS_LEVEL, name.scpot, COMPRESS_LEVEL, ...
      name.psp, COMPRESS_LEVEL, name.psp_p, COMPRESS_LEVEL, ...
      name.bitmask, COMPRESS_LEVEL, name.quality, COMPRESS_LEVEL};

    %% Update VariableAttributes
    VATTRIB.CATDESC = {name.epoch, 'Time tags, in TT2000'; ...
      name.label,   'Label'; ...
      name.scpot,   'Spacecraft potential';...
      name.psp,     'Probe to spacecraft potential, averaged'; ...
      name.psp_p,   'Probe to spacecraft potential, individual probes'; ...
      name.bitmask, 'Status bitmask';...
      name.quality, 'Quality indicator'};
    VATTRIB.DEPEND_0 = {name.scpot, name.epoch; ...
      name.psp,     name.epoch; ...
      name.psp_p,   name.epoch; ...
      name.bitmask, name.epoch;...
      name.quality, name.epoch};
    VATTRIB.DISPLAY_TYPE = {name.scpot, 'time_series'; ...
      name.psp,     'time_series'; ...
      name.psp_p,   'time_series';...
      name.quality, 'time_series'};
    VATTRIB.FIELDNAM = {name.epoch, 'Time tags'; ...
      name.label,   'Label'; ...
      name.scpot,   'Spacecraft potential'; ...
      name.psp,     'Probe to spacecraft potential'; ...
      name.psp_p,   'Probe to spacecraft potential individual probe'; ...
      name.bitmask, 'Status bitmask';...
      name.quality, 'Quality indicator'};
    VATTRIB.FILLVAL = {name.epoch, int64(-9223372036854775808); ...
      name.scpot,    single(-1.0E31); ...
      name.psp,     single(-1.0E31); ...
      name.psp_p,   single(-1.0E31); ...
      name.bitmask, uint8(255);...
      name.quality, int16(-32768)};
    VATTRIB.FORMAT = {name.scpot, 'F8.3'; ...
      name.label,   'A23'; ...
      name.psp,     'F8.3'; ...
      name.psp_p,   'F8.3'; ...
      name.bitmask, 'I7';...
      name.quality, 'I7'};
    VATTRIB.LABLAXIS = {name.scpot, 'Spacecraft potential';
      name.psp,     'Probe to spacecraft potential';...
      name.quality, 'ScPot quality'};
    VATTRIB.LABL_PTR_1 = {name.psp_p, name.label};
    VATTRIB.SI_CONVERSION = {name.scpot, '1.0>V'; ...
      name.psp,     '1.0>V'; ...
      name.psp_p,   '1.0>V'};
    VATTRIB.UNITS = {name.epoch, 'ns'; ...
      name.scpot,   'V'; ...
      name.psp,     'V'; ...
      name.psp_p,   'V'; ...
      name.bitmask, 'unitless';...
      name.quality, 'unitless'};
    VATTRIB.VALIDMIN = {name.epoch, spdfcomputett2000([1990,01,01,0,0,0,0,0,0]); ...
      name.scpot,   VOLTAGE_MIN; ...
      name.psp,     VOLTAGE_MIN; ...
      name.psp_p,   VOLTAGE_MIN; ...
      name.bitmask, uint8(0);...
      name.quality, int16(-32767)};
    VATTRIB.VALIDMAX = {name.epoch, spdfcomputett2000([2100,01,01,0,0,0,0,0,0]); ...
      name.scpot,   VOLTAGE_MAX; ...
      name.psp,     VOLTAGE_MAX; ...
      name.psp_p,   VOLTAGE_MAX; ...
      name.bitmask, uint8(254);...
      name.quality, QUALITY_MAX};
    VATTRIB.VAR_TYPE = {name.epoch, 'support_data'; ...
      name.label,   'metadata'; ...
      name.scpot,   'data'; ...
      name.psp,     'data'; ...
      name.psp_p,   'data'; ...
      name.bitmask, 'support_data';...
      name.quality, 'support_data'};
    VATTRIB.MONOTON = {name.epoch, 'INCREASE'};

  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

%% Write to file.
% Update GlobalAttributes specific to ql/sitl or l2pre.
GATTRIB.Data_type = {dataType}; % 'fast_l1b_dce2d', 'slow_l1b_dce2d' or 'brst_l1b_dce2d'.
GATTRIB.Logical_source = {[datasetPrefix '_' dataType]}; % ie. mms2_edp_fast_l1b_OptDesc FIXME...
GATTRIB.Logical_source_description = {dataDesc}; % in full words.

% Replace any NaN with correct FILLVAL
replaceNaN()

write_file()

% Return to previous working directory.
cd(oldDir);

%% Help functions
  function [fileName, verStr] = get_file_name
    % Generate output file name incrementing the file version if necessary
    switch procId
      case {MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql}
        subDir = procName; suf = 'dce2d';
      case MMS_CONST.SDCProc.l2pre
        subDir = 'l2pre'; suf = 'dce2d';
      case MMS_CONST.SDCProc.scpot
        subDir = 'l2'; suf = 'scpot';
      case MMS_CONST.SDCProc.l2a
        subDir = 'l2'; suf = 'Adce2d';
      otherwise
        errStr = 'unrecognized procId';
        irf.log('critical', errStr); error(errStr)
    end
    scIdStr = sprintf('mms%d',scId);
    startTime =  HeaderInfo.startTime;
    verStr = sprintf('%d.%d.',MMS_CONST.Version.X,MMS_CONST.Version.Y);
    fileName = [scIdStr '_' INST_NAME, '_' tmModeStr '_' subDir '_' ...
      suf '_' startTime '_v' ];
    
    % Check for preexisting files and increment file version
    dataPathPref = [ENVIR.DATA_PATH_ROOT, filesep,'science',filesep, ...
      scIdStr, filesep, INST_NAME, filesep, tmModeStr, filesep, ...
      subDir, filesep, startTime(1:4), filesep, startTime(5:6), filesep, ...
      startTime(7:8), filesep];
    
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

  function replaceNaN
    % NaN is not ISTP compliant in CDF files and should not be used for the
    % MMS mission, however it can be written/created in CDF files therefor
    % replace any found NaN with the corresponding correct FillValue as
    % specified for that variable.
    for ii=1:2:size(outVars,2)
      % outVars = {'varName', data}, but may include variables without
      % FillVal (such as labels).
      ind = find(strcmp(VATTRIB.FILLVAL(:,1), outVars{1,ii})>0);
      if ind
        % Variable (ind) has a FILLVAL, locate any NaN and replace them
        outVars{1,ii+1}(isnan(outVars{1,ii+1}))=VATTRIB.FILLVAL{ind,2};
      end
    end
  end

  function write_file
   % write file with arguments obtained above, also include md5 checksum.
   spdfcdfwrite(outFileName, outVars, 'Vardatatypes',varDatatype, ...
     'GlobalAttributes', GATTRIB, 'VariableAttributes', VATTRIB, ...
     'RecordBound', recBound, 'VarCompress', compressVars, ...
     'Checksum', 'md5');
  end

  function GATTRIB = getGlobalAttributes
    %% Create GlobalAttribute struct
    % Source: MMS_CDF_Format_Guide, v1.7 draft, issued 2014/10/11.
    GATTRIB=[];
    % Global Attributes REQUIRED:
    GATTRIB.Data_type = cell(0,1); % mode_dataLevel_optionalDescriptor
    GATTRIB.Data_version = cell(0,1); % Same as version number in filename.
    GATTRIB.Descriptor = {'EDP>Electric Double Probe'};
    GATTRIB.Discipline = {'Space Physics>Magnetospheric Science'};
    GATTRIB.Generation_date = {datestr(now,'yyyymmdd')};
    GATTRIB.Instrument_type = {'Electric Fields (space)'};
    GATTRIB.Logical_file_id = cell(0,1); % Same as filename without ".cdf".
    GATTRIB.Logical_source = cell(0,1); % Ex: mms3_edp_fast_ql_swd (mmsSC_instrument_mode_dataLevel_optionalDescriptor)
    GATTRIB.Logical_source_description = cell(0,1); % Description in full words..
    GATTRIB.Mission_group = {'MMS'};
    GATTRIB.PI_affiliation = {'SWRI, LASP, KTH'};
    GATTRIB.PI_name = {'J.Burch, R.Ergun, P.Lindqvist.'};
    GATTRIB.Project = {'STP>Solar-Terrestrial Physics'};
    GATTRIB.Source_name = {sprintf('MMS%i>MMS Satellite Number %i',scId,scId)}; % Or possibly 'MMS>MMS Constellation'.
    GATTRIB.TEXT = {'http://mms.gsfc.nasa.gov/'; ...
      ['The full name of PI affiliations: SWRI - Southwest Research Institute. ',...
      'LASP - Laboratory for Atmospheric and Space Physics. ',...
      'KTH - Kungliga Tekniska Hogskolan (Swedish Royal Institute of Technology). ']};  % FIXME This attribute is an SPDF standard global attribute, which is a text description of the
      %experiment whose data is included in the CDF. A reference to a journal article(s) or to a
      %World Wide Web page describing the experiment is essential, and constitutes the
      %minimum requirement. A written description of the data set is also desirable. This
      %attribute can have as many entries as necessary to contain the desired information.
      %Typically, this attribute is about a paragraph in length and is not shown on CDAWeb.
      % Note: Matlab & spdcdfwrite to cdf files result in issues with ASCII
      % for the charachter "ö", therefor replace ö with o.
    GATTRIB.HTTP_LINK = {'http://mms.gsfc.nasa.gov/'; 'http://mms.space.swri.edu/'}; % FIXME should point to data
    GATTRIB.LINK_TEXT = {'Magnetospheric Multiscale (MMS) mission home page'; 'SMART package home page'}; % FIXME as well
    GATTRIB.LINK_TITLE = {'At NASA GSFC'; 'At SWRI'}; % FIXME as well
    GATTRIB.MODS = MMS_CONST.Version.MODS; % Text describing major version changes, ie. "vX" changes.
    % Global Attributes RECOMMENDED:
    GATTRIB.Acknowledgement = cell(0,1);
    GATTRIB.Generated_by = {['IRFU Matlab', irf('version')]};
    % Global Attributes OPTIONAL:
%    GATTRIB.Parents = cell(0,1); % Req if number of source cdf >= 2.
    GATTRIB.Skeleton_version = {'v0.0.4'};
    GATTRIB.Rules_of_use = cell(0,1);
    GATTRIB.Time_resolution = cell(0,1);
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
