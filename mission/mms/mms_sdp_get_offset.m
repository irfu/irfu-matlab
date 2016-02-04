function off = mms_sdp_get_offset(scId, procId, time)
% Return offsets etc. to be applied in MMS processing (currently only 
% QL dce2d or Scpot.
%
% Example: 
%   MMS_CONST = mms_constants;
%   off = mms_sdp_dsl_offset(1, MMS_CONST.SDCProc.ql, dceTime);
%  returns 
%   off.ex         - Sunward offset (DSL Ex), for QL dce2d
%   off.ey         - Duskward offset (DSL Ey), for QL dce2d
%   off.p2p        - Probe2Plasma potential, for L2 scpot
%   off.shortening - Shortening factor, for L2 scpot
%   off.calFile    - Calibration filename used (to be stored in GATTRIB in
%                    the resulting cdf file).
%
% If calibration files are found, offsets will be interpolated with
% 'previous' and 'extrap'. If failing to read latest calibration file a
% static value is returned.
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

try
  calPath = [ENVIR.CAL_PATH_ROOT, filesep 'mms', num2str(scId), filesep, ...
    'edp'];
  calFiles = [calPath, filesep, 'mms', num2str(scId), ...
    '_edp_sdp_', MMS_CONST.SDCProcs{procId},'_v*.txt' ];
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
    offIntrp = interp1(double(offTime.ttns-offTime.start.ttns), ...
      [C{2}, C{3}], double(time-offTime.start.ttns), 'previous', 'extrap');
    switch procId
      case MMS_CONST.SDCProc.ql
        off.ex = offIntrp(:,1);
        off.ey = offIntrp(:,2);
      case MMS_CONST.SDCProc.scpot
        off.p2p = offIntrp(:,1);
        off.shortening = offIntrp(:,2);
      otherwise
        errStr = 'Unexpected procId';
        irf.log('critical',errStr); error(errStr);
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
    end
  end % useStaticOffset

end
