function off = mms_sdp_get_offset(scId, procId, Tint)
% Function to return offsets to be applied to DSL Ex (sunward) and DSL Ey
% (duskward) in processing of MMS despun electric fields.
% Example: 
%   MMS_CONST = mms_constants;
%   off = mms_sdp_dsl_offset(1, MMS_CONST.SDCProc.ql);
%  returns 
%   off.ex      - Sunward offset (DSL Ex), for QL dce2d, L2 dce2d products
%   off.ey      - Duskward offset (DSL Ey), for QL dce2d, L2 dce2d products
%   off.p2p     - Probe2Plasma potential, for L2 scpot products
%   off.calFile - Calibration filename used (to be stored in GATTRIB in
%                 output cdf file).
%
% Note: For now only a fixed value per s/c and for QL dce only.
%
% NOTE: Keep in mind the signs used here and when applying the offsets and
% scale factors.
%
% See also: MMS_SDP_DMGR

narginchk(2,3);
global MMS_CONST ENVIR
if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
% ENVIR gets loaded from mms_sdc_sdp_init, should not be empty here.
if isempty(ENVIR), errStr='Empty ENVIR'; irf.log(errStr); error(errStr); end

try 
  calFiles = [ENVIR.CAL_PATH_ROOT, filesep 'mms', scId, filesep 'edp', filesep, ...
    'mms', scId, '_edp_sdp_', MMS_CONST.SDCProcs{procId},'_v*.txt' ];
  list = dir(calFiles);
  if(~isempty(list))
    % Found at least one cal file for EDP_SDP_{SDCProcs}. Use the last
    % version, ie. end of list.
    % Load the txt file
    [~,off.calFile, ~] = fileparts(list(end));
    offTmp = load(list(end),'-ascii');
    % Process it, interp1(,,'extrap');
  else
    off = useStaticOffset;
  end
catch ME
  irf.log('warning', ME.message);
  off = useStaticOffset;
end


  function off = useStaticOffset
    % Did not locate any calibration offset file or failed to read it, use
    % static offsets instead.
    irf.log('warning', 'Did not locate any CALIBRATION files, using static values');
    switch procId
      case MMS_CONST.SDCProc.ql
        switch lower(scId)
          % See mms_sdp_dmgr;
          % dE = mms_sdp_despin(...);
          % dE(:,1) = dE(:,1) - offs.ex;
          % dE(:,2) = dE(:,2) - offs.ey;
          case 1
            off.ex = 2.2;
            off.ey = 0.0;
          case 2
            off.ex = 3.0;
            off.ey = 0.0;
          case 3
            off.ex = 2.55;
            off.ey = 0;
          case 4
            off.ex = 1.32;
            off.ey = 0;
          otherwise
            errStr = 'Unexpected scId';
            irf.log('critical',errStr); error(errStr);
        end
      case MMS_CONST.SDCProc.scpot
        switch lower(scId)
          % See mms_sdp_dmgr;
          % SCpot = - Probe2sc_pot.data(:) .* off.shortening(:) + offs.p2p;
          case 1
            off.p2p = 1.3;
            off.shortening = 1.2;
          case 2
            off.p2p = 1.5;
            off.shortening = 1.2;
          case 3
            off.p2p = 1.2;
            off.shortening = 1.2;
          case 4
            off.p2p = 0.0;
            off.shortening = 1.2;
          otherwise
            errStr = 'Unexpected scId';
            irf.log('critical',errStr); error(errStr);
        end
      otherwise
        irf.log('critical', 'Process does not yet have offsets. SETTING ALL TO ZERO');
        off.ex = 0.0;
        off.ey = 0.0;
        off.p2p = 0;
    end
    % Return a static, non empty, string to indicate these values was used.
    off.calFile = '20160204 static';
  end % useStaticOffset

end