classdef mms_sdp_dmgr < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    adc_off = [];     % comp ADC offsets
    dce = [];         % src DCE file
    dce_xyz_dsl = []; % comp E-field xyz DSL-coord
    dcv = [];         % src DCV file
    defatt = [];      % src DEFATT file
    defeph = [];      % src DEFEPH file
    hk_101 = [];      % src HK_101 file
    hk_105 = [];      % src HK_105 file
    hk_10e = [];      % src HK_10E file
    l2pre = [];       % src L2Pre file
    phase = [];       % comp phase
    probe2sc_pot = [];% comp probe to sc potential
    sc_pot = [];      % comp sc potential
    spinfits = [];    % comp spinfits
  end
  properties (SetAccess = immutable)
    CONST = [];
    samplerate = [];
    procId = [];
    tmMode = [];
    scId = [];
  end
  
  methods
    function DATAC = mms_sdp_dmgr(scId,procId,tmMode,samplerate)
      DATAC.CONST = mms_constants();
      MMS_CONST = DATAC.CONST;
      if nargin == 0,
        errStr = 'Invalid input for scId';
        irf.log('critical', errStr); error(errStr);
      end
      
      if  ~isnumeric(scId) || ...
          isempty(intersect(scId, MMS_CONST.MMSids))
        errStr = 'Invalid input for scId';
        irf.log('critical', errStr); error(errStr);
      end
      DATAC.scId = scId;
      
      if nargin < 2 || isempty(procId)
        DATAC.procId = 1;
        irf.log('warning',['procId not specified, defaulting to '''...
          MMS_CONST.SDCProcs{DATAC.procId} ''''])
      elseif ~isnumeric(procId) || ...
          isempty(intersect(procId, 1:numel(MMS_CONST.SDCProcs)))
        errStr = 'Invalid input for init_struct.procId';
        irf.log('critical', errStr); error(errStr);
      else DATAC.procId = procId;
      end
      
      if nargin < 3 || isempty(tmMode)
        DATAC.tmMode = 1;
        irf.log('warning',['tmMode not specified, defaulting to '''...
          MMS_CONST.TmModes{DATAC.tmMode} ''''])
      elseif ~isnumeric(tmMode) || ...
          isempty(intersect(tmMode, 1:numel(MMS_CONST.TmModes)))
        errStr = 'Invalid input for init_struct.tmMode';
        irf.log('critical', errStr); error(errStr);
      else DATAC.tmMode = tmMode;
      end
      
      if nargin < 4 || isempty(samplerate)
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
        DATAC.samplerate = samplerate;
      end
    end
  end
  
end

