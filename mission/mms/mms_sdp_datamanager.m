function [ varOut ] = mms_sdp_datamanager( param, dataObj )
% mms_sdp_datamanager will store and retrive data for
%  	mms_sdp_datamanager( dataType, dataObj ) will store
%  	appropriate data variables related to dataType in the global variable
%  	DATAC for later retrival.
%  
%   [varOut] = mms_sdp_datamanager( variable ) will return the variable
%   requested to varOut, if no such variable has been created already it
%   will try and calculate it using the stored data.
%   
%  	Example:
%
%  		mms_sdp_datamanager('dce',dceDataObj)
%  		mms_sdp_datamanager('phase')
%  
%   	See also DATAOBJ, MMS_SDP_LOAD.

narginchk(1,2); % One argument is simply retreive variable, Two arguments
% store "dataObj" as "dataType".

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
global DATAC; % Here we store read data.

if nargout, varOut = []; end

if ~ischar(param),
  errStr = 'PARAM must be a string';
  irf.log('critical', errStr); error(errStr);
end

if strcmpi(param, 'init')
  % Initialize
  if nargin==1
    errStr = 'INIT requires second input argument';
    irf.log('critical', errStr); error(errStr);
  elseif ~isstruct(dataObj)
    errStr = 'Second input argument for INIT must be a structure';
    irf.log('critical', errStr); error(errStr);
  elseif ~isfield(dataObj,'scId') || ~isnumeric(dataObj.scId) || ...
      isempty(intersect(dataObj.scId, MMS_CONST.MMSids))
    errStr = 'Invalid input for init_struct.scId';
    irf.log('critical', errStr); error(errStr);
  end
  DATAC.scId = dataObj.scId;
  if ~isfield(dataObj,'tmMode')
    DATAC.tmMode = 1;
    irf.log('warining',['init_struct.tmMode not specified, defaulting to '''...
      MMS_CONST.TmModes{DATAC.tmMode} ''''])
  elseif ~isnumeric(dataObj.tmMode) || ...
      isempty(intersect(dataObj.tmMode, 1:numel(MMS_CONST.TmModes)))
    errStr = 'Invalid input for init_struct.tmMode';
    irf.log('critical', errStr); error(errStr);
  else DATAC.tmMode = dataObj.tmMode;
  end
  if ~isfield(dataObj,'procId')
    DATAC.procId = 1;
    irf.log('warining',['init_struct.procId not specified, defaulting to '''...
      MMS_CONST.SDCProcs{DATAC.procId} ''''])
  elseif ~isnumeric(dataObj.procId) || ...
      isempty(intersect(dataObj.procId, 1:numel(MMS_CONST.SDCProcs)))
    errStr = 'Invalid input for init_struct.procId';
    irf.log('critical', errStr); error(errStr);
  else DATAC.procId = dataObj.procId;
  end
  if( ~isfield(dataObj,'samplerate') )
    % Normal operations, samplerate identified from TmMode
    if(~isfield(MMS_CONST.Samplerate,MMS_CONST.TmModes{DATAC.tmMode}))
      irf.log('warning', ['init_struct.samplerate not specified,'...
        'nor found in MMS_CONST for ',MMS_CONST.TmModes{DATAC.tmMode}]);
      irf.log('warning', ['Defaulting samplerate to ',...
        MMS_CONST.Samplerate.(MMS_CONST.TmModes{1})]);
      DATAC.samplerate = MMS_CONST.Samplerate.(MMS_CONST.TmModes{1});
    else
      % Automatic, determined by tmMode.
      DATAC.samplerate = ...
        MMS_CONST.Samplerate.(MMS_CONST.TmModes{DATAC.tmMode});
    end
  else
    % Commissioning, sample rate identified previously.
    DATAC.samplerate = dataObj.samplerate;
  end
  DATAC.adc_off = [];     % comp ADC offsets
  DATAC.dce = [];         % src DCE file
  DATAC.dce_xyz_dsl = []; % comp E-field xyz DSL-coord
  DATAC.dcv = [];         % src DCV file
  DATAC.defatt = [];      % src DEFATT file
  DATAC.defeph = [];      % src DEFEPH file
  DATAC.hk_101 = [];      % src HK_101 file
  DATAC.hk_105 = [];      % src HK_105 file
  DATAC.hk_10e = [];      % src HK_10E file
  DATAC.l2pre = [];       % src L2Pre file
  DATAC.phase = [];       % comp phase
  DATAC.probe2sc_pot = [];% comp probe to sc potential
  DATAC.sc_pot = [];      % comp sc potential
  DATAC.spinfits = [];    % comp spinfits
  return
end

if ~isfield(DATAC, 'scId')
  errStr = 'Not initialized! Run: mms_sdp_datamanager(''init'',init_struct)';
  irf.log('critical', errStr);
  error('MATLAB:MMS_SDP_DATAMANAGER:INPUT', errStr);
end

% If one input is give, we give the parameter back if it is already there,
% otherwise we compute it
if nargin==1
  if ~isfield(DATAC,param)
    errStr = ['unknown parameter (' param ')'];
    irf.log('critical',errStr), error(errStr)
  end
  
  if ~isempty(DATAC.(param)), varOut = DATAC.(param); return, end
  
  % Check for function fo compute param
  funcName = ['mms_sdp_comp_' param];
  if ~( exist(funcName)==2 ) %#ok<EXIST>
    irf.log('warning',...
      ['Cannot find function (' funcName ') to compute param(' param ')'])
    DATAC.(param) = MMS_CONST.Error; varOut = DATAC.(param); return
  end
  % Run the function
  DATAC.(param) = feval(funcName); varOut = DATAC.(param); return
end

%% nargin==2   
% Make sure first argument is a dataobj class object,
% otherwise a read cdf file.
if isa(dataObj,'dataobj') % do nothing
elseif ischar(dataObj) && exist(dataObj, 'file')
  % If it is not a read cdf file, is it an unread cdf file? Read it.
  irf.log('warning',['Loading ' param ' from file: ', dataObj]);
  dataObj = dataobj(dataObj, 'KeepTT2000');
elseif isstruct(dataObj) && any(strcmp(param, {'defatt', 'defeph'}))
  % Is it the special case of DEFATT/DEFEPH (read as a struct into dataObj).
  % Do nothing..
else
  errStr = 'Unrecognized input argument';
  irf.log('critical', errStr);
  error('MATLAB:MMS_SDP_DATAMANAGER:INPUT', errStr);
end

if( isfield(DATAC, param) ) && ~isempty(DATAC.(param))
  % Error, Warning or Notice for replacing the data variable?
  errStr = ['replacing existing variable (' param ') with new data'];
  irf.log('critical', errStr);
  error('MATLAB:MMS_SDP_DATAMANAGER:INPUT', errStr);
end

varPrefix = sprintf('mms%d_edp_',DATAC.scId);

switch(param)
  case('dce')
    sensors = {'e12','e34','e56'};
    init_param()
    
  case('dcv')
    sensors = {'v1','v2','v3','v4','v5','v6'};
    init_param()
    chk_timeline()
    chk_latched_p()
    %apply_transfer_function()
    v_from_e_and_v()
    chk_bias_guard()
    chk_sweep_on()
    chk_sdp_v_vals()
    apply_nom_amp_corr()
    
  case('hk_101')
    % HK 101, contains sunpulses.
    varPrefix = sprintf('mms%d_101_',DATAC.scId);
    DATAC.(param) = [];
    DATAC.(param).dataObj = dataObj;
    x = getdep(dataObj,[varPrefix 'cmdexec']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    % Add sunpulse times (TT2000) of last recieved sunpulse.
    DATAC.(param).sunpulse = dataObj.data.([varPrefix 'sunpulse']).data;
    % Add sunpulse indicator, real: 0, SC pseudo: 1, CIDP pseudo: 2.
    DATAC.(param).sunssps = dataObj.data.([varPrefix 'sunssps']).data;
    % Add CIDP sun period (in microseconds, 0 if sun pulse not real.
    DATAC.(param).iifsunper = dataObj.data.([varPrefix 'iifsunper']).data;
    
  case('hk_105')
    % HK 101, contains sunpulses.
    varPrefix = sprintf('mms%d_105_',DATAC.scId);
    DATAC.(param) = [];
    DATAC.(param).dataObj = dataObj;
    x = getdep(dataObj,[varPrefix 'sweepstatus']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    % Add sweepstatus which indicates if any of the probes is sweeping
    DATAC.(param).sweepstatus = dataObj.data.([varPrefix 'sweepstatus']).data;

  case('hk_10e')
    % HK 10E, contains bias.
    varPrefix = sprintf('mms%d_10e_',DATAC.scId);
    DATAC.(param) = [];
    DATAC.(param).dataObj = dataObj;
    x = getdep(dataObj,[varPrefix 'seqcnt']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    % Go through each probe and store values for easy access,
    % for instance probe 1 dac values as: "DATAC.hk_10e.beb.dac.v1".
    hk10eParam = {'dac','ig','og','stub'}; % DAC, InnerGuard, OuterGuard & Stub
    for iParam=1:length(hk10eParam)
      for jj=1:6
        % stub only exist if probe is 5 or 6.
        if( ~strcmp(hk10eParam{iParam},'stub') || ...
            (strcmp(hk10eParam{iParam},'stub') && jj>=5))
          tmpStruct = getv(dataObj,...
            [varPrefix 'beb' num2str(jj,'%i') hk10eParam{iParam}]);
          if isempty(tmpStruct)
            errS = ['cannot get ' varPrefix 'beb' num2str(jj,'%i') ...
              hk10eParam{iParam}]; irf.log('warning',errS), warning(errS);
          else
            DATAC.(param).beb.(hk10eParam{iParam}).(sprintf('v%i',jj)) = ...
              tmpStruct.data;
          end
        end
      end % for jj=1:6
    end % for iParam=1:length(hk10eParam)

  case('defatt')
    % DEFATT, contains Def Attitude (Struct with 'time' and 'zphase' etc)
    % As per e-mail discussion of 2015/04/07, duplicated timestamps can
    % occur in Defatt (per design). If any are found, use the last data
    % point and disregard the first duplicate.
    idxBad = diff(dataObj.time)==0; % Identify first duplicate
    fs = fields(dataObj);
    for idxFs=1:length(fs), dataObj.(fs{idxFs})(idxBad) = []; end
    DATAC.(param) = dataObj;
    check_monoton_timeincrease(DATAC.(param).time);

  case('defeph')
    % DEFEPH, contains Def Ephemeris (Struct with 'time', 'Pos_X', 'Pos_Y'
    % and 'Pos_Z')
    DATAC.(param) = dataObj;
    check_monoton_timeincrease(DATAC.(param).time);

  case('l2pre')
    % L2Pre, contain dce data, spinfits, etc. for L2a processing.
    DATAC.(param) = [];
    DATAC.(param).dataObj = dataObj;
    % Split up the various parts (spinfits, dce data [e12, e34, e56],
    % dce bitmask [e12, e34, e56], phase from the l2pre file to their
    % expected locations in DATAC. (so that remaining processing can use
    % same syntax).
    varPre = ['mms', num2str(DATAC.scId), '_edp_dce'];
    varPre2 = '_spinfit_';
    DATAC.spinfits = [];
    sdpPair = {'e12', 'e34'};
    for iPair=1:numel(sdpPair)
      DATAC.spinfits.sfit.(sdpPair{iPair}) = ...
        dataObj.data.([varPre, varPre2, sdpPair{iPair}]).data(:,2:end);
      DATAC.spinfits.sdev.(sdpPair{iPair}) = ...
        dataObj.data.([varPre, varPre2, sdpPair{iPair}]).data(:,1);
    end
    x = getdep(dataObj,[varPre, varPre2, sdpPair{iPair}]);
    DATAC.spinfits.time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.spinfits.time, 'L2Pre spinfits');
    sensors = {'e12', 'e34', 'e56'};
    DATAC.dce = [];
    x = getdep(dataObj,[varPre, '_data']);
    DATAC.dce.time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.dce.time, 'L2Pre dce');
    for iPair=1:numel(sensors);
      DATAC.dce.(sensors{iPair}).data = ...
        dataObj.data.([varPre, '_data']).data(:,iPair);
      DATAC.dce.(sensors{iPair}).bitmask = ...
        dataObj.data.([varPre, '_bitmask']).data(:,iPair);
    end
    DATAC.phase.data = dataObj.data.([varPre, '_phase']).data;

  otherwise
    % Not yet implemented.
    errStr = [' unknown parameter (' param ')'];
    irf.log('critical',errStr);
    error('MATLAB:MMS_SDP_DATAMANAGER:INPUT', errStr);
end

  function chk_latched_p()
    % Check that probe values are varying. If there are 3 identical points,
    % or more, after each other mark this as latched data. If it is latched
    % and the data has a value below MMS_CONST.Limit.LOW_DENSITY_SATURATION
    % it will be Bitmasked with Low density saturation otherwise it will be
    % bitmasked with just Probe saturation.
    
    % For each sensor, check each pair, i.e. V_1 & V_2 and E_12.
    for iSen = 1:2:numel(sensors)
      senA = sensors{iSen};  senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor
      irf.log('notice', ...
        sprintf('Checking for latched probes on %s, %s and %s.', senA, ...
        senB, senE));
      DATAC.dcv.(senA) = latched_mask(DATAC.dcv.(senA));
      DATAC.dcv.(senB) = latched_mask(DATAC.dcv.(senB));
      DATAC.dce.(senE) = latched_mask(DATAC.dce.(senE));
      % TODO: Check overlapping stuck values, if senA stuck but not senB..  
    end
    function sen = latched_mask(sen)
      idx = irf_latched_idx(sen.data);
      if ~isempty(idx)
        sen.bitmask(idx) = bitor(sen.bitmask(idx),...
          MMS_CONST.Bitmask.PROBE_SATURATION);
      end
      idx = sen.data<MMS_CONST.Limit.LOW_DENSITY_SATURATION;
      if any(idx)
        sen.bitmask(idx) = bitor(sen.bitmask(idx), ...
          MMS_CONST.Bitmask.LOW_DENSITY_SATURATION);
      end
    end
  end

  function chk_timeline()
    % Check that DCE time and DCV time overlap and are measured at the same
    % time (within insturument delays). Throw away datapoint which does not
    % overlap between DCE and DCV.
    if isempty(DATAC.dce),
      irf.log('warning','Empty DCE, cannot proceed')
      return
    end
    % 3.8 us per channel and 7 channels between DCV (probe 1) and DCE (12).
    % A total shift of 26600 ns is therefor to be expected, add this then
    % convert to seconds before comparing times.
    %[~, dce_ind, dcv_ind] = intersect(DATAC.dce.time, DATAC.dcv.time+26001);
    % XXX: The above line is faster, but we cannot be sure it works, as it
    % requires exact mathcing of integer numbers
    % XXX: The code below should work in principle, but id does not
    %[dce_ind, dcv_ind] = irf_find_comm_idx(DATAC.dce.time,...
    %  DATAC.dcv.time+26600,int64(40000)); % tolerate 40 us jitter
    % Highest bitrate is to be 8192 samples/s, which correspond to about
    % 122 us between consecutive measurements. The maximum theoretically
    % allowed jitter would be half of this (60us) for the dcv & dce
    % measurements to be completely unambiguous, use 1/3 margin on this.
    
    %This is a HACK. We just take the nearest time, assuming times in DCE
    %and DCV must be identical.
    tE = DATAC.dce.time; tV = DATAC.dcv.time;
     
    if ~(all(median(diff(tE))==diff(tE)) && all(median(diff(tV))==diff(tV)))
      errStr1 = 'Do not know how to handle gaps';
      irf.log('critical',errStr1), error(errStr1)
    end
    
    % Bring together the DCE and DCV time series
    % NOTE: No gaps allowed below this line
    dt = median(diff(tE));
    if tV(1)>tE(1), tStart = tE(1);
    else tStart = tE(1) - ceil((tE(1)-tV(1))/dt)*dt;
    end
    if tE(end)>tV(end), tStop = tE(end);
    else tStop = tE(end) + ceil((tV(end)-tE(end))/dt)*dt;
    end
    nData = (tStop - tStart)/dt + 1;
    newTime = int64((1:nData) - 1)'*dt + tStart;
    [~,idxEonOld,idxEonNew] = intersect(tE,newTime); 
    idxEoffNew = setxor(1:length(newTime),idxEonNew);
    tDiffNew = abs(newTime-tV(1)); tDiffOld = abs(newTime(1)-tV);
    if min(min(tDiffNew),min(tDiffOld))==min(tDiffOld),
      iDcvStartOld = find(tDiffOld==min(tDiffOld));
      tDiffNew = abs(newTime-tV(iDcvStartOld));
      iDcvStartNew = find(tDiffNew==min(tDiffNew));
    else
      iDcvStartNew = find(tDiffNew==min(tDiffNew));
      tDiffOld = abs(newTime(iDcvStartNew)-tV);
      iDcvStartOld = find(tDiffOld==min(tDiffOld));
    end
    idxVonOld = (1:length(tV))'-1 +iDcvStartOld;
    idxVonNew = (1:length(tV))'-1 +iDcvStartNew;
    idxVoffNew = setxor(1:length(tV),idxVonNew);
    
    for iSen = 1:2:numel(sensors)  % Loop over e12, e34, e56
      senA = sensors{iSen};  senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor
      save_restore('dce',senE,idxEonOld,idxEonNew,idxEoffNew)
      save_restore('dcv',senA,idxVonOld,idxVonNew,idxVoffNew)
      save_restore('dcv',senB,idxVonOld,idxVonNew,idxVoffNew)
    end
    DATAC.dce.time = newTime; DATAC.dcv.time = newTime;
    
    function save_restore(sig,sen,idxOnOld,idxOnNew,idxOffNew)
      % Save old values, expand the variables and restore the old values
      SAVE = DATAC.(sig).(sen);
      DATAC.(sig).(sen).data = NaN(size(newTime));
      DATAC.(sig).(sen).bitmask = zeros(size(newTime));
      DATAC.(sig).(sen).data(idxOnNew)    = SAVE.data(idxOnOld);
      DATAC.(sig).(sen).bitmask(idxOnNew) = SAVE.bitmask(idxOnOld);
      DATAC.(sig).(sen).bitmask(idxOffNew) = MMS_CONST.Bitmask.SIGNAL_OFF;
    end
  end % CHK_TIMELINE

  function chk_bias_guard()
    % Check that bias/guard setting, found in HK_10E, are nominal. If any
    % are found to be non nominal set bitmask value in both V and E.


    if(~isempty(DATAC.hk_10e))  % is a hk_10e file loaded?

      NomBias = uint32(32768); % FIXME, read from a table file

      %DATAC.hk_10e.beb.ig.v1(600:630) = uint32(10); % FIXME, forced to non-moninal, for debug.
      %DATAC.hk_10e.beb.og.v2(600) = uint32(0); % FIXME, forced to non-nominal

      irf.log('notice','Checking for non nominal bias settings.');
      for iSen = 1:2:numel(sensors)
        senA = sensors{iSen};  senB = sensors{iSen+1};
        senE = ['e' senA(2) senB(2)]; % E-field sensor

        hk10eParam = {'ig','og'}; % InnerGuard, OuterGuard
        for iiParam = 1:length(hk10eParam);

          % FIXME, proper test of existing fields?
          if( ~isempty(DATAC.hk_10e.beb.(hk10eParam{iiParam}).(senA)) && ...
            ~isempty(DATAC.hk_10e.beb.(hk10eParam{iiParam}).(senB)) );

            % Interpolate HK_10E to match with DCV timestamps, using the
            % previous HK value.
            interp_DCVa = interp1(double(DATAC.hk_10e.time), ...
              double(DATAC.hk_10e.beb.(hk10eParam{iiParam}).(senA)), ...
              double(DATAC.dcv.time), 'previous', 'extrap');

            interp_DCVb = interp1(double(DATAC.hk_10e.time), ...
              double(DATAC.hk_10e.beb.(hk10eParam{iiParam}).(senB)), ...
              double(DATAC.dcv.time), 'previous', 'extrap');

            % Locate Non Nominal values
            indA = interp_DCVa ~= NomBias;
            indB = interp_DCVb ~= NomBias;
            indE = or(indA,indB); % Either senA or senB => senE non nominal.

            if(any(indE))
              irf.log('notice',['Non-nominal bias settings found on pair ',...
                senE,' from ',hk10eParam{iiParam},'. Setting bitmask.']);

              % Add bitmask values to SenA, SenB and SenE for these ind.
              bits = MMS_CONST.Bitmask.BAD_BIAS;
              % Add value to the bitmask, leaving other bits untouched.
              DATAC.dcv.(senA).bitmask(indA) = ...
                bitor(DATAC.dcv.(senA).bitmask(indA), bits);
              DATAC.dcv.(senB).bitmask(indB) = ...
                bitor(DATAC.dcv.(senB).bitmask(indB), bits);
              DATAC.dce.(senE).bitmask(indE) = ...
                bitor(DATAC.dce.(senE).bitmask(indE), bits);
            end
          else
            irf.log('Warning',['HK_10E : no proper values for ',...
              senA,' and ',senB,'.']);
          end % if ~isempty()
        end % for iiParam
      end % for iSen
    else
      irf.log('Warning','No HK_10E file : cannot perform bias/guard check');
    end % if ~isempty(hk_10e)
  end

  function chk_sweep_on()
    % Check if sweep is on for all probes
    % if yes, set bit in both V and E bitmask

    if isempty(DATAC.dce),
      irf.log('warning','Empty DCE, cannot proceed')
      return
    end
    
    varPref = sprintf('mms%d_sweep_', DATAC.scId);
    if ~isfield(DATAC.dce.dataObj.data,[varPref 'start'])
      errS = ['Did not find ',varPref,'start'];
      irf.log('critical',errS); error(errS)
    end
    
    % Get sweep status and sweep timestamp
    sweepStart = DATAC.dce.dataObj.data.([varPref 'start']).data;
    sweepStop = DATAC.dce.dataObj.data.([varPref 'stop']).data;
    sweepSwept = DATAC.dce.dataObj.data.([varPref 'swept']).data;
    
    if isempty(sweepStart)
      irf.log('warning','No sweep status in DCE file');
      % Alternative approach for finding sweep times using hk_105
      if isempty(DATAC.hk_105)
        irf.log('warning','No HK_105 file loaded: cannot identify sweeps.');
        return
      end
      sweepStatus = logical(DATAC.hk_105.sweepstatus);
      sweepStart = DATAC.hk_105.time([diff(sweepStatus)==1; false]);
      if sweepStatus(1) % First point nas sweep ON, start at t(0)-4s
        sweepStart = [DATAC.hk_105.time(1)-int64(4e9) sweepStart];
      end
      sweepStop = DATAC.hk_105.time([false; diff(sweepStatus)==-1]);
      if sweepStatus(end)  % Last point has sweep ON, stop at t(end)+4s
        sweepStop = [sweepStop DATAC.hk_105.time(end)+int64(4e9)];
      end
      if length(sweepStop)~=length(sweepStart)
        % Sanity check, should never be here
        errSt = 'length(sweepStop) != length(sweepStart)!!';
        irf.log('critical',errSt), error(errSt)
      end
      % No info on which probe is swept in hk_105, can be any pair
      sweepSwept = zeros(size(sweepStart));
    end
    
    % For each pair, E_12, E_34, E_56.
    for iSen = 1:2:numel(sensors)
      senA = sensors{iSen};  senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor
      irf.log('notice', ['Checking for sweep status on probe pair ', senE]);
      % Locate probe pair senA and senB, SweepSwept = 1 (for pair 12), etc.
      if all(sweepSwept==0), senN = 0; else senN = str2double(senA(2)); end
      ind = find(sweepSwept==senN);
      sweeping = false(size(DATAC.dce.time)); % First assume no sweeping.
      for ii = 1:length(ind)
        % Each element in ind correspond to a sweep_start and sweep_stop
        % with the requested probe pair. Identify which index these times
        % correspond to in DATAC.dce.time. Each new segment, where
        % (sweep_start<=dce.time<=sweep_stop), is added with 'or' to the
        % previous segments.
        sweeping = or( and(DATAC.dce.time>=sweepStart(ind(ii)), ...
          DATAC.dce.time<=sweepStop(ind(ii))), sweeping);
      end
      if(any(sweeping))
        % Set bitmask on the corresponding pair, leaving the other 16 bits
        % untouched.
        irf.log('notice','Sweeping found, bitmasking it.');
        bits = MMS_CONST.Bitmask.SWEEP_DATA;
        DATAC.dcv.(senA).bitmask(sweeping) = ...
          bitor(DATAC.dcv.(senA).bitmask(sweeping), bits);
        DATAC.dcv.(senB).bitmask(sweeping) = ....
          bitor(DATAC.dcv.(senB).bitmask(sweeping), bits);
        DATAC.dce.(senE).bitmask(sweeping) = ...
          bitor(DATAC.dce.(senE).bitmask(sweeping), bits);
      else
        irf.log('debug',['Did not find any sweep for probe pair ', senE]);
      end % if any(sweeping)
    end % for iSen
  end

  function chk_sdp_v_vals()
    % check if probe-to-spacecraft potentials  averaged over one spin for 
    % all probes are similar (within TBD %, or V). 
    % If not, set bit in both V and E bitmask.
    
    %XXX: Does nothing at the moment
  end

  function apply_nom_amp_corr()
    % Apply a nominal amplitude correction factor to DCE & DCV for p1..4
    % values after cleanup but before any major processing has occured.
    factor = MMS_CONST.NominalAmpCorr;
    for iSen = 1:2:min(numel(sensors),4)
      senA = sensors{iSen}; senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)];
      logStr = sprintf(['Applying nominal amplitude correction factor, '...
        '%d, to probes %s, %s and %s'], factor, senA, senB, senE);
      irf.log('notice',logStr);
      DATAC.dcv.(senA).data = DATAC.dcv.(senA).data * factor;
      DATAC.dcv.(senB).data = DATAC.dcv.(senB).data * factor;
      DATAC.dce.(senE).data = DATAC.dce.(senE).data * factor;
    end
  end

  function v_from_e_and_v
    % Compute V from E and the other V
    % typical situation is V2 off, V1 on
    % E12[mV/m] = ( V1[V] - V2[V] ) / L[km]
    if isempty(DATAC.dce),
      irf.log('warning','Empty DCE, cannot proceed')
      return
    end
    
    % Nominal boom length used in L1b processor
    NOM_BOOM_L = .001; % 1m, XXX the tru value should be 120 m
    
    MSK_OFF = MMS_CONST.Bitmask.SIGNAL_OFF;
    for iSen = 1:2:numel(sensors)
      senA = sensors{iSen}; senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor
      senA_off = bitand(DATAC.dcv.(senA).bitmask, MSK_OFF);
      senB_off = bitand(DATAC.dcv.(senB).bitmask, MSK_OFF);
      senE_off = bitand(DATAC.dce.(senE).bitmask, MSK_OFF);
      idxOneSig = xor(senA_off,senB_off);
%      if ~any(idxOneSig), return, end
      iVA = idxOneSig & ~senA_off;
      if any(iVA),
        irf.log('notice',...
          sprintf('Computing %s from %s and %s for %d data points',...
          senB,senA,senE,sum(iVA)))
        DATAC.dcv.(senB).data(iVA) = DATAC.dcv.(senA).data(iVA) - ...
          NOM_BOOM_L*DATAC.dce.(senE).data(iVA);
      end
      iVB = idxOneSig & ~senB_off;
      if any(iVB),
        irf.log('notice',...
          sprintf('Computing %s from %s and %s for %d data points',...
          senA,senB,senE,sum(iVA)))
        DATAC.dcv.(senA).data(iVB) = DATAC.dcv.(senB).data(iVB) + ...
          NOM_BOOM_L*DATAC.dce.(senE).data(iVB);
      end
      % For comissioning data we will have all DCE/DCV, verify consistency.
      idxBoth = and(and(~senA_off, ~senB_off), ~senE_off);
      if any(idxBoth)
        irf.log('notice',...
          sprintf('Verifying %s = (%s - %s)/NominalLength for %d points',...
          senE, senA, senB, sum(idxBoth)));
        remaining = abs(NOM_BOOM_L*DATAC.dce.(senE).data(idxBoth) - ...
          DATAC.dcv.(senA).data(idxBoth) + ...
          DATAC.dcv.(senB).data(idxBoth));
        if(any(remaining>MMS_CONST.Limit.DCE_DCV_DISCREPANCY))
          irf.log('critical',...
            'Datapoints show a discrepancy between DCE and DCV!');
          % FIXME: Bitmasking them or exit with Error?
        end
      end % if any(idxBoth)
    end % for iSen
  end

  function init_param
    DATAC.(param) = [];
    if ~all(diff(dataObj.data.([varPrefix 'samplerate_' param]).data)==0)
      err_str = ...
        'MMS_SDP_DATAMANAGER changing sampling rate not yet implemented.';
      irf.log('warning', err_str);
      %error('MATLAB:MMS_SDP_DATAMANAGER:INPUT', err_str);
    end
    DATAC.(param).dataObj = dataObj;
    fileVersion = DATAC.(param).dataObj.GlobalAttributes.Data_version{:};
    % Skip the intial "v" and split it into [major minor revision].
    fileVersion = str2double(strsplit(fileVersion(2:end),'.'));
    DATAC.(param).fileVersion = struct('major', fileVersion(1), 'minor',...
        fileVersion(2), 'revision', fileVersion(3));
    % Make sure it is not too old to work properly.
    if DATAC.(param).fileVersion.major < MMS_CONST.MinFileVer
      err_str = sprintf('File too old: major version %d < %d',...
        DATAC.(param).fileVersion.major, MMS_CONST.MinFileVer);
      irf.log('critical',err_str), error(err_str); %#ok<SPERR>
    end
    x = getdep(dataObj,[varPrefix param '_sensor']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    sensorData = dataObj.data.([varPrefix param '_sensor']).data;
    if isempty(sensors), return, end
    probeEnabled = resample_probe_enable(sensors);
    %probeEnabled = are_probes_enabled;
    for iSen=1:numel(sensors)
      DATAC.(param).(sensors{iSen}) = struct(...
        'data',sensorData(:,iSen), ...
        'bitmask',zeros(size(sensorData(:,iSen)),'uint16'));
      %Set disabled bit
      idxDisabled = probeEnabled(:,iSen)==0;
      if(any(idxDisabled>0))
        irf.log('notcie', ['Probe ',sensors{iSen}, ' disabled for ',...
          num2str(sum(idxDisabled)),' points. Bitmask them and set to NaN.']);
        DATAC.(param).(sensors{iSen}).bitmask(idxDisabled) = ...
          bitor(DATAC.(param).(sensors{iSen}).bitmask(idxDisabled), ...
          MMS_CONST.Bitmask.SIGNAL_OFF);
        DATAC.(param).(sensors{iSen}).data(idxDisabled,:) = NaN;
      end
    end
  end

  function res = are_probes_enabled
    % Use FILLVAL of each sensor to determine if probes are enabled or not.
    % Returns logical of size correspondig to sensor.
    sensorData = dataObj.data.([varPrefix param '_sensor']).data;
    FILLVAL = getfillval(dataObj, [varPrefix, param, '_sensor']);
    if( ~ischar(FILLVAL) )
      % Return 'true' for all data not equal to specified FILLVAL
      res = (sensorData ~= FILLVAL);
    else
      errStr = 'Unable to get FILLVAL.';
      irf.log('critical',errStr); error(errStr);
    end
  end

  function res = resample_probe_enable(fields)
  % resample probe_enabled data to E-field cadense
    probe = fields{1};
    flag = get_variable(dataObj,[varPrefix probe '_enable']);
    dtSampling = median(diff(flag.DEPEND_0.data));
    switch DATAC.tmMode
%      case MMS_CONST.TmMode.srvy, error('kaboom')
      case MMS_CONST.TmMode.slow, dtNominal = [20, 160]; % seconds
      case MMS_CONST.TmMode.fast, dtNominal = 5;
      case MMS_CONST.TmMode.brst, dtNominal = [0.625, 0.229 0.0763];
      case MMS_CONST.TmMode.comm, dtNominal = [0.500, 1.250, 2.500, 5.000];
      otherwise
        errS = 'Unrecognized tmMode';
        irf.log('critical',errS), error(errS)
    end
    dtNominal = int64(dtNominal*1e9); % convert to ns
    
    flagOK = false;
    for i=1:numel(dtNominal)
      if dtSampling > dtNominal(i)*.95 && dtSampling < dtNominal(i)*1.05
        dtSampling = dtNominal(i); flagOK = true; break
      end
    end
    if ~flagOK
      errS = ['bad sampling for ' varPrefix probe '_enable'];
      irf.log('critical',errS), error(errS)
    end
    enabled.time = flag.DEPEND_0.data;
    nData = numel(enabled.time);
    enabled.data = zeros(nData,numel(fields));
    enabled.data(:,1) = flag.data;
    for iF=2:numel(fields)
      probe = fields{iF};
      flag = getv(dataObj,[varPrefix probe '_enable']);
      if isempty(flag)
        errS = ['cannot get ' varPrefix probe '_enable'];
        irf.log('critical',errS), error(errS)
      elseif numel(flag.data) ~= nData
        errS = ['bad size for ' varPrefix probe '_enable'];
        irf.log('critical',errS), error(errS)
      end
      enabled.data(:,iF) = flag.data;
    end
    newT = DATAC.(param).time;
    % Default to zero - probe disabled
    res = zeros(numel(newT), numel(fields));
    if all(diff(enabled.data))==0,
      ii = newT>enabled.time(1)-dtSampling & newT<=enabled.time(end);
      for iF=1:numel(fields), 
        res(ii,iF) = enabled.data(1,iF); 
      end
    else
      % TODO: implements some smart logic.
      errS = 'MMS_SDP_DATAMANAGER enabling/disabling probes not yet implemented.';
      irf.log('critical', errS); error(errS);
    end
  end
end


% Short function for verifying Time is increasing.
function check_monoton_timeincrease(time, dataType)
    
if(any(diff(time)<=0))
        err_str = ['Time is NOT increasing for the datatype ', dataType];
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDP_DATAMANAGER:TIME:NONMONOTON', err_str);
end

end
