function off = mms_sdp_get_offset(scId, procId, time, TMmode, Vpsp)
% Return offsets etc. to be applied in MMS processing (currently only
% QL dce2d, L2Pre dce2d or Scpot.
%
% Example:
%   MMS_CONST = mms_constants;
%   TMmode = MMS_CONST.TmMode.brst;
%   off = mms_sdp_get_offset(1, MMS_CONST.SDCProc.ql, dceTime, TMmode);
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
% NOTE: if processing Slow L2pre than also include the otherwise optional
%      "Vpsp" - a TSeries of "probe2scpot" for corresponding interval.
%
% NOTE: Keep in mind the signs used here and when applying the offsets and
% scale factors.
%
% See also: MMS_SDP_DMGR

narginchk(3,5);
global MMS_CONST ENVIR
if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
if ~exist('TMmode','var'), TMmode=[]; end
% ENVIR gets loaded from mms_sdc_sdp_init, should not be empty here.
if isempty(ENVIR)
  errStr='Empty ENVIR'; irf.log('critical',errStr); error(errStr);
end
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
  if procId == MMS_CONST.SDCProc.l2pre && TMmode == MMS_CONST.TmMode.slow
    useVpspdependentoff
  else
    calPath = [ENVIR.CAL_PATH_ROOT, filesep 'mms', num2str(scId), filesep, ...
      'edp', filesep, 'sdp', filesep, calStr];
    calFiles = [calPath, filesep, 'mms', num2str(scId), ...
      '_edp_sdp_', calStr, '_*v*.txt' ];
    list = dir(calFiles);
    if(~isempty(list))
      % Found at least one cal file for EDP_SDP_{SDCProcs}. Use the last
      % version, ie. end of list.
      % Load the txt file
      list = list(end);
      [~, off.calFile, ~] = fileparts(list.name);

      verStr = regexp(off.calFile, ['mms[1-4]_edp_sdp_', calStr, '_\d{8,8}_v(?<VERSION>\d{1,}.\d{1,}.\d{1,})'],'tokens');
      if is_version_geq(verStr{1}{1}, '1.0.0')
        % NEW format of Calibration files
        % As discussed during meeting at KTH, 2017/10/?? the new calibration
        % files have two extra columns providing the time margins to be
        % applied before and after each timestamp.
        %         __________
        %        /
        %      ./
        %_____/
        %
        %      | <-- Times given by first column, but the value does not come
        %      into full effect until after the margin of the last column.
        %      But it starts to come into effect before the margin of the
        %      second to last column. If margins are zero, then the old value
        %      is used until the timestamp in the first column and after this
        %      timestamp the new offsets are used. Otherwise a linear
        %      interpolation is applied in the intermediate interval, before
        %      which the previous offset is used and correspondingly after
        %      the margin the new offset is fully used.
        fmt = '%s\t%f\t%f\t%f\t%f';
        fileID = fopen([calPath, filesep, list.name]);
        C = textscan(fileID, fmt, 'HeaderLines', 1);
        fclose(fileID);
        % Interpret the time strings and convert it to int64 (tt2000).
        offTime = EpochTT(cell2mat(C{1}));
        % Offset start time (based on ROI).
        time1 = offTime.ttns;
        % Add margins to start time (ROI),
        % positive for slow to ensure same offset is used for continous
        % segments (each ROI start close to switch from slow to fast tmmode).
        % ROI region 2015-09-20T06:27:14Z/2015-09-20T20:26:54 UTC
        % MMS1 edp fast l1b start 2015-09-20T06:13:45.097915313 UTC
        % ie. continous fast mode start some 13 min 29 seconds before ROI.
        if TMmode == MMS_CONST.TmMode.slow
          % +15 minutes margin (undo the built in margin written into calibration files (ROI start time has been adjusted).
          time1 = time1 + int64(15*60*10^9);
        end
        % Ensure we have at least a margin of 5 us, even if zero was entered
        % for both margins.
        % NOTE: While it could be as small as 1ns with int64 representation it
        % can in reallity not be that small as "interp1" in Matlab R2013b req.
        % double representation. And must be kept as distinct value or offset
        % values will not have the same number of elements as our offset time
        % (Highest MMS burst rate used as of 2017/05/05 is 16384Hz => dt=61 us).
        ind = (abs(C{4}) + abs(C{5}) <= 5e-6);
        if any(ind)
          irf.log('notice', 'Increasing time margins read from file to at least 5e-6 seconds.');
          C{5}(ind) = 5e-6;
        end
        if any(C{4}>0)
          errStr = 'Expected negative or zero values for the negative time margin "-dt". Aborting and using static values';
          irf.log('warning', errStr); error(errStr);
        end
        time_end   = time1 + int64(10^9 * C{4}); % The previous offsets end here (and interpolation start)
        time_start = time1 + int64(10^9 * C{5}); % The new offset fully start here (and interpolation ends)
        ind = repmat(1:length(C{2}), 2, 1); % 1 1 2 2 3 3 4 4 ...
        dataOff = [C{2}(ind(:)), C{3}(ind(:))]; % Repeated data
        % Let the last offset be valid "forever" after (interp1 with "linear"
        % and "extrap" require one extra datapoint to ensure it is interpolated
        % as a static value and not a linear trend between the penultimate and
        % ultimate offset value.) Add one extra point one year after the last.
        time_end(1) = time_start(end) + int64(365*86400e9);
        timeComb = unique(sort([time_start, time_end]));

      else
        % OLD format of Calibration files
        fmt = '%s\t%f\t%f';
        fileID = fopen([calPath, filesep, list.name]);
        C = textscan(fileID, fmt, 'HeaderLines', 1);
        fclose(fileID);
        % Interpret the time strings and convert it to int64 (tt2000).
        offTime = EpochTT(cell2mat(C{1}));
        % Offset start time (based on ROI).
        time1 = offTime.ttns;
        % Add margins to start time (ROI),
        % positive for slow to ensure same offset is used for continous
        % segments (each ROI start close to switch from slow to fast tmmode).
        % ROI region 2015-09-20T06:27:14Z/2015-09-20T20:26:54 UTC
        % MMS1 edp fast l1b start 2015-09-20T06:13:45.097915313 UTC
        % ie. continous fast mode start some 13 min 29 seconds before ROI.
        switch TMmode
          case MMS_CONST.TmMode.slow
            Tmargin = +int64(15*60*10^9); % +15 minutes margin (undo the built in margin written into calibration files (ROI start time has been adjusted).
          otherwise % fast, brst have a built in correction to times in the offset files compared with official ROI start time.
            Tmargin = int64(0);
        end
        time1 = time1 + Tmargin;
        % Let each offset be valid until 5 us before next offset begin
        % NOTE: While it could be as small as 1ns with int64 representation it
        % can in reallity not be that small as "interp1" in Matlab R2013b req.
        % double representation.
        % (Highest MMS burst rate used as of 2017/05/05 is 16384Hz => dt=61 us).
        time2 = time1(2:end) - int64(5000);
        % Let the last offset be valid "forever" after (interp1 with "linear"
        % and "extrap" require one extra datapoint to ensure it is interpolated
        % as a static value and not a linear trend between the penultimate and
        % ultimate offset value.) Add one extra point one year after the last.
        time2(end+1) = time2(end) + int64(365*86400e9);

        data3 = [C{2}, C{3}; ...
          C{2}, C{3}]; % Repeated data
        [timeSort, indSort] = sort([time1(:); time2(:)]); % Almost repeated time (5 us diff), then sorted
        [timeComb, indUniq] = unique(timeSort); % Ensure no duplicated values
        dataSort = data3(indSort, :); % Sorted data (based on time)
        dataOff = dataSort(indUniq, :); % Ensure no duplicated values (based on time)

      end

      % Covert time and timeComb to double after subtracting the start time
      % (2015/01/01) from both to ensure interp1 works as expected.
      timeReq = double(time - timeComb(1));
      timeOff = double(timeComb - timeComb(1));
      % Interpolate offsets with linear and extrapolation to the requested
      % time
      offIntrp = interp1(timeOff, dataOff, timeReq, 'linear', 'extrap');


      switch procId
        case {MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2pre}
          off.ex = offIntrp(:,1);
          off.ey = offIntrp(:,2);
        case MMS_CONST.SDCProc.scpot
          off.p2p = offIntrp(:,1);
          off.shortening = offIntrp(:,2);
      end
    else
      off = useStaticOffset;
    end
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

  function useVpspdependentoff
    timeTS = EpochTT(time);
    Tint = irf.tint(timeTS);
    epoch1min = ceil(Tint.start.epochUnix/60)*60:20:fix(Tint.stop.epochUnix/60)*60;
    Epoch20s = EpochUnix(epoch1min); % Define 20 baseline
    if ~exist('Vpsp','var') || isempty(Vpsp) || ~isa(Vpsp, 'TSeries')
      errStr='NO SCPOT TSeries loaded/given as argument. Fall back to static offset!';
      irf.log('critical', errStr);
      error(errStr);
    end
    Vpsp20s = Vpsp.resample(Epoch20s);
    switch lower(scId)
      case 1
        V0 = -0.2; % Define constants, IF changed remember to update the date in "off.calFile" below
        a = 0.9;
        b = 2.5;
        c = 0.5;
        d = 3.0;
        e = 1.2;
      case 2
        V0 = -0.3;
        a = 1.7;
        b = 3.7;
        c = 0.5;
        d = 2.5;
        e = 1.2;
      case 3
        V0 = -0.25;
        a = 1.9;
        b = 5.3;
        c = 1.7;
        d = 4.2;
        e = 1.5;
      case 4
        V0 = -0.35;
        a = 2.2;
        b = 15;
        c = 5;
        d = 0.12;
        e = 1.0;
      otherwise
        errStr = 'Unexpected scId';
        irf.log('critical',errStr); error(errStr);
    end

    f = a - b/abs(V0 - c) + d/abs(V0 + e)^4;

    dEx = zeros(size(Epoch20s));
    dEy = zeros(size(Epoch20s));
    idx1 = Vpsp20s.data <= V0;
    idx2 = Vpsp20s.data > V0;

    dEx(idx1) = a-b./abs(Vpsp20s.data(idx1) - c);
    dEx(idx2) = f-d./abs(Vpsp20s.data(idx2) + e).^4;

    dE = irf.ts_vec_xy(Epoch20s,[dEx dEy]);
    dEslow = dE.resample(timeTS);

    off.ex = dEslow.data(:,1);
    off.ey = dEslow.data(:,2);
    off.calFile = 'SCPOT dependent with 20200302 static coef.'; % Date, YYYYMMDD, of when coef. above was last changed
  end


end
