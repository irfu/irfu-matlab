function off = mms_sdp_get_offset(scId, procId, time)
% Return offsets etc. to be applied in MMS processing (currently only 
% QL dce2d, L2Pre dce2d or Scpot.
%
% Example: 
%   MMS_CONST = mms_constants;
%   off = mms_sdp_dsl_offset(1, MMS_CONST.SDCProc.ql, dceTime);
%  returns 
%   off.ex         - Sunward offset (DSL Ex), for QL dce2d or L2Pre dce2d
%   off.ey         - Duskward offset (DSL Ey), for QL dce2d or L2Pre dce2d
%   off.p2p        - Probe2Plasma potential, for L2 scpot
%   off.shortening - Shortening factor, for L2 scpot
%   off.calFile    - Calibration filename used (to be stored in GATTRIB in
%                    the resulting cdf file).
%
% If calibration files are found, offsets will be the last values just
% preceeding the start of time. (Start of data interval processing). If no
% files are found a static value is used.
%
% NOTE: Keep in mind the signs used here and when applying the offsets and
% scale factors.
%
% See also: MMS_SDP_DMGR

narginchk(3,3);
global MMS_CONST ENVIR
if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
% ENVIR gets loaded from mms_sdc_sdp_init, should not be empty here.
if isempty(ENVIR), errStr='Empty ENVIR'; irf.log(errStr); error(errStr); end
if isa(time,'GenericTime'), time = time.ttns; end

switch procId
  case {MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2pre}
    calStr = 'dsl';
  case MMS_CONST.SDCProc.scpot
    calStr = 'scpot';
  otherwise
    calStr = ''; errStr = 'Unexpected procId';
    irf.log('critical',errStr); %error(errStr);
end

try
  calPath = [ENVIR.CAL_PATH_ROOT, filesep 'mms', num2str(scId), filesep, ...
    'edp'];
  calFiles = [calPath, filesep, 'mms', num2str(scId), ...
    '_edp_sdp_', calStr, '_v*.txt' ];
  list = dir(calFiles);
  if(~isempty(list))
    % Found at least one cal file for EDP_SDP_{SDCProcs}. Use the last
    % version, ie. end of list.
    % Load the txt file
    list = list(end);
    [~, off.calFile, ~] = fileparts(list.name);
    fmt = '%s\t%f\t%f';
    fileID = fopen([calPath, filesep, list.name]);
    C = textscan(fileID, fmt, 'HeaderLines', 1);
    fclose(fileID);
    offTime = EpochTT(cell2mat(C{1}));
    ind = offTime.ttns <= time(1); % Preceeding offsets only.
    if(~any(ind)), ind = true; end % Sanity check if time was incorrect...
    % Use only the last of the preceeding offset values for the entire data
    % interval.
    offIntrp = [C{2}(ind), C{3}(ind)];
    switch procId
      case {MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2pre}
        off.ex = offIntrp(end,1);
        off.ey = offIntrp(end,2);
      case MMS_CONST.SDCProc.scpot
        off.p2p = offIntrp(end,1);
        off.shortening = offIntrp(end,2);
   end
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
    % Default output if nothing matching the requested process.
    off = struct('ex', 0.0, 'ey', 0.0, 'p2p', 0.0, 'shortening', 1.0, ...
      'calFile', '20160204 static');
    switch procId
      case {MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2pre}
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
    end
  end % useStaticOffset

end
