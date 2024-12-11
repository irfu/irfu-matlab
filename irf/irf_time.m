function tOutput = irf_time(tInput,flag)
%IRF_TIME  Convert time between different formats
%
%   t_out=IRF_TIME(t_in,'in>out');
%   t_out=IRF_TIME(t_in,'in2out');  % deprecated
%
%   Input:
%       t_in: column vector of input time in format 'in'
%   Output:
%       t_out: column vector of output time in format 'out'
%
%   Formats 'in' and 'out' can be (default 'in' is 'epoch'):%
%       epoch: seconds since the epoch 1 Jan 1970, default, used by the ISDAT system.
%      vector: [year month date hour min sec] last five columns optional
%     vector6: [year month date hour min sec]
%     vector9: [year month date hour min sec msec micros nanosec]
%         iso: deprecated, same as 'utc'
%        date: MATLAB datenum format
%     datenum: same as 'date'
%         doy: [year, doy]
%        doy8: only input, [year, doy, hour, min, sec, msec, micros, nanosec]
%          tt: Terrestrial Time, seconds past  January 1, 2000, 11:58:55.816 (UTC)
%        ttns: Terrestrial Time, nanoseconds past  January 1, 2000, 11:58:55.816 (UTC)
%         utc: UTC string (see help spdfparsett2000 for supported formats)
%  utc_format: only output, where 'format' can be any string where 'yyyy' is
%              changed to year, 'mm' month, 'dd' day, 'HH' hour, 'MM'
%              minute, 'SS' second, 'mmm' milisceonds, 'uuu' microseconds,
%              'nnn' nanoseconds. E.g. 'utc_yyyy-mm-dd HH:MM:SS.mmm'
%              Values exceeding the requested precision are truncated,
%              e.g. 10:40.99 is returned as "10:40" using format "utc_HH:MM".
%    cdfepoch: miliseconds since 1-Jan-0000
%  cdfepoch16: [seconds since 1-Jan-0000, picoseconds within the second]
%     epochtt: return class EPOCHTT
%
%  t_out=IRF_TIME(t_in,'out') equivalent to t_out=IRF_TIME(t_in,'epoch>out');
%  t_out=IRF_TIME(t_in) equivalent to t_out=IRF_TIME(t_in,'vector>epoch');
%
%  Example: t_out=irf_time('2011-09-13T01:23:19.000000Z','utc>epoch');
%
%  There are also commands to convert time intervals
%   time_int=IRF_TIME(tint,'tint>utc')
%           convert time interval to utc (first column epoch start time, 2nd epoch end time)


persistent tlastcall strlastcall
if nargin==0 % return string with current time (second precision)
  % datestr is slow function therefore if second has not passed use the old
  % value of string as output without calling datestr function. Increases speed!
  if isempty(tlastcall)
    tlastcall=0; % initialize
  end
  if 24*3600*(now-tlastcall)>1
    tlastcall=now;
    strlastcall=irf_time(now,'date>utc_yyyy-mm-dd HH:MM:SS');
  end
  tOutput=strlastcall;
  return
elseif nargin==1
  flag='vector>epoch';
end
if isempty(tInput),tOutput=[];return;end
flagTint=strfind(flag,'tint'); % check if we work with time interval (special case)
if isempty(flagTint)   %#ok<STREMP>       % default transformation
  flag2=strfind(flag,'2');   % see if flag has number 2 in it
  if isempty(flag2)
    flag2=strfind(flag,'>');   % see if flag has '>' in it
  else
    irf.log('warning',['irf_time(..,''' flag '''): irf_time(..,''in2out'') is deprecated, please use irf_time(..,''in>out'')']);
    flag(flag2)='>';
  end
  if isempty(flag2)         % if no '2' convert from epoch
    format_in='epoch';
    format_out=flag;
    flag=[format_in '>' format_out];
  else
    format_in=flag(1:flag2-1);
    format_out=flag(flag2+1:end);
  end
  if strcmp(format_in,format_out) % if in and out equal return
    tOutput=tInput;
    return;
  elseif ~strcmp(format_in,'ttns') && ~strcmp(format_out,'ttns')
    % if there is no epoch in the format then
    % first convert from 'in' to 'epoch'
    % and then from 'epoch' to 'out'
    t_temp=irf_time(tInput,[format_in '>ttns']);
    tOutput=irf_time(t_temp,['ttns>' format_out]);
    return
  end
else
  flag2=strfind(flag,'2');   % see if flag has number 2 in it
  if ~isempty(flag2)
    irf.log('warning',['irf_time(..,''' flag '''): irf_time(..,''in2out'') is deprecated, please use irf_time(..,''in>out'')']);
    flag(flag2)='>';
  end
end


%
% flag should include 'tint' or 'epoch'
% no other flags are allowed below
%
switch lower(flag)
  case 'vector>tt'
    % Convert a [YYYY MM DD hh mm ss ms micros ns] to TT2000 in double s
    % the last columns can be omitted, default values MM=1,DD=1,other=0
    tOutput = double(irf_time(tInput,'vector>ttns')/1e9);
  case 'vector>ttns'
    % Convert a [YYYY MM DD hh mm ss ms micros ns] to TT2000 in int64 ns
    % the last columns can be omitted, default values MM=1,DD=1,other=0
    nCol = size(tInput,2);
    if nCol > 9
      error('irf_time:badInputs',...
        'input should be column vector [YYYY MM DD hh mm ss ms micros ns], last columns can be omitted.')
    elseif nCol == 6
      tOutput = irf_time(tInput,'vector6>ttns');
    elseif nCol < 9
      defValues = [2000 1 1 0 0 0 0 0 0];
      tInput = [tInput repmat(defValues(nCol+1:9),size(tInput,1),1)];
      tOutput = spdfcomputett2000(tInput);
    end
  case 'vector6>ttns'
    % Convert a [YYYY MM DD hh mm ss.xxxx] to TT2000 in int64 ns
    nCol = size(tInput,2);
    if nCol ~= 6
      error('irf_time:badInputs',...
        'input should be column vector with 6 columns [YYYY MM DD hh mm ss.xxx].')
    end
    tSecRound = floor(tInput(:,6));
    tmSec = 1e3*(tInput(:,6)-tSecRound);
    tmSecRound = floor(tmSec);
    tmicroSec = 1e3*(tmSec - tmSecRound);
    tmicroSecRound = floor(tmicroSec);
    tnSecRound = floor(1e3*(tmicroSec - tmicroSecRound));
    tInput(:,6:9) = [tSecRound tmSecRound tmicroSecRound tnSecRound];
    tOutput = spdfcomputett2000(tInput);
  case {'ttns>utc','ttns>iso'}
    if any(strfind(flag,'iso'))
      irf.log('warning','irf_time: ''iso'' is deprecated and will be removed, please use ''utc'', see help.');
    end
    tOutput = GenericTimeArray.ttns2utc(tInput);
  case 'tt>ttns'
    tOutput = int64(tInput)*1e9;
  case 'ttns>tt'
    tOutput = double(tInput)/1e9;
  case 'ttns>vector9'
    tOutput = spdfbreakdowntt2000(tInput);
  case 'ttns>vector'
    tVec9 = spdfbreakdowntt2000(tInput);
    tOutput = tVec9(:,1:6);
    tOutput(:,6) = tVec9(:,6)+1e-3*tVec9(:,7)+1e-6*tVec9(:,8)+1e-9*tVec9(:,9);
  case {'utc>ttns','iso>ttns'}
    if any(strfind(flag,'iso'))
      irf.log('warning','irf_time: ''iso'' is deprecated and will be removed, please use ''utc'', see help.');
    end
    if any(strfind(tInput(1,:),'T'))
      tOutput = GenericTimeArray.utc2ttns(tInput);
    else
      mask = '%4d-%2d-%2d %2d:%2d:%f%*c';
      s=tInput';
      s(end+1,:)=sprintf(' ');
      a = sscanf(s,mask);
      N = numel(a)/6;
      if N~=fix(N) || N~=size(s,2)
        irf.log('warning','something is wrong with iso input format, returning empty!'),
        tOutput=[];
        return;
      end
      a = reshape(a,6,N);
      a = a';
      tOutput = irf_time(a,'vector6>ttns');
    end
  case 'ttns>epoch'
    tOutput = toepoch(irf_time(tInput,'ttns>vector'));
  case 'epoch>ttns'
    tOutput = EpochUnix.to_ttns(tInput);
  case {'ttns>date','ttns>datenum'} % matlab date
    tOutput = spdftt2000todatenum(tInput);
  case {'date>ttns','datenum>ttns'}
    tOutput = spdfdatenumtott2000(tInput);
  case 'ttns>doy'
    ttBreak = irf_time(tInput,'ttns>vector');
    ttBreak(:,2:end) = 1;
    tOutput = [ttBreak(:,1) ...
      floor(irf_time(tInput,'ttns>datenum'))-floor(irf_time(ttBreak,'vector>datenum'))+1];
  case 'doy>ttns'
    tOutput=irf_time([tInput(:,1) tInput(:,1).*0+1 tInput(:,1).*0+1 ...
      tInput(:,2).*24-12 tInput(:,1).*0 tInput(:,1).*0],'vector6>ttns');
  case 'doy8>ttns'
    tOutput = spdfcomputett2000([tInput(:,1) tInput(:,1).*0+1 tInput(:,2) ...
      tInput(:,3) tInput(:,4) tInput(:,5) tInput(:,6) tInput(:,7) tInput(:,8)]);
  case 'cdfepoch>ttns'
    sz = size(tInput);
    if sz(2)>1
      irf.log('warning','irf_time(cdfepoch>ttns: input is not column vector! output is column vector!');
      tInput = tInput(:);
    end
    ttBreak = spdfbreakdownepoch(tInput); % cdfepoch, not other "epoch"
    % Assume 0 microsec and 0 nanosec.
    tOutput = spdfcomputett2000([ttBreak zeros(size(ttBreak,1),2)]);
  case 'ttns>cdfepoch'
    ttBreak = spdfbreakdowntt2000(tInput);
    % Should round up/down column 7, depending on column 8 to 9
    ttBreak(:,7) = round(ttBreak(:,7)+ttBreak(:,8)/10^3+ttBreak(:,9)/10^6);
    tOutput = spdfcomputeepoch(ttBreak(:,1:7)); % cdfepoch, not other "epoch"
  case 'cdfepoch16>ttns'
    ttBreak = spdfbreakdownepoch16(tInput);
    % Should round up/down column 9, depending on column 10
    ttBreak(:,9) = round(ttBreak(:,9)+ttBreak(:,10)/10^3);
    tOutput = spdfcomputett2000(ttBreak(:,1:9));
  case 'ttns>cdfepoch16'
    ttBreak = spdfbreakdowntt2000(tInput);
    % Assume 0 picoseconds.
    tOutput = spdfcomputeepoch16([ttBreak zeros(size(ttBreak,1),1)]);

  case 'ttns>epochtt'
    tOutput = EpochTT(tInput);
  case 'epochtt>ttns'
    tOutput = tInput.ttns;
    %
    % Time interval conversions
    %

  case {'tint>utc','tint>iso'}
    if any(strfind(flag,'iso'))
      irf.log('warning','irf_time(..,''tint>iso'') is deprecated and will be removed, please use irf_time(..,''tint>utc'').');
    end
    tOutput = irf_time(tInput,'tint>utc_'); % fmt = []

  case {'utc>tint','iso>tint'}
    if any(strfind(flag,'iso'))
      irf.log('warning','irf_time(..,''iso>tint'') is deprecated and will be removed, please use irf_time(..,''utc>tint'').');
    end
    % assume column array where each row is interval in iso format
    ii=strfind(tInput(1,:),'/');
    t1=irf_time(tInput(:,1:ii-1),'utc>epoch');
    t2=irf_time(tInput(:,ii+1:end),'utc>epoch');
    if isempty(t1) || isempty(t2)
      tOutput=[];
    else
      tOutput=[t1 t2];
    end

  case 'tint>isoshort'
    tOutput = irf_time(tInput,'tint>utc_yyyy-mm-ddTHH:MM:SS.mmmZ');

  otherwise
    if numel(flag)>=9 && strcmp(flag(1:9),'tint>utc_')
      fmt = flag(10:end);
      if isa(tInput,'GenericTimeArray') && length(tInput) == 2
        t1iso=tInput(1).utc(fmt);
        t2iso=tInput(2).utc(fmt);
      elseif isa(tInput,'int64')
        t1iso = irf_time(tInput(:,1),['ttns>utc_' fmt]);
        t2iso = irf_time(tInput(:,2),['ttns>utc_' fmt]);
      else
        t1iso = irf_time(tInput(:,1),['epoch>utc_' fmt]);
        t2iso = irf_time(tInput(:,2),['epoch>utc_' fmt]);
      end
      t1iso(:,end+1)='/';
      tOutput = [t1iso t2iso];
    elseif numel(flag)>=9 && strcmp(flag(1:9),'ttns>utc_')
      fmt = flag(10:end);
      tOutput = GenericTimeArray.ttns2utc(tInput,fmt);
    else
      disp(['!!! irf_time: unknown flag ''' lower(flag) ''', not converting.'])
      tOutput=tInput;
    end
end
