function t_out = irf_time(t_in,flag)
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
if nargin==0, % return string with current time (second precision)
	% datestr is slow function therefore if second has not passed use the old
	% value of string as output without calling datestr function. Increases speed!
	if isempty(tlastcall),
		tlastcall=0; % initialize
	end
	if 24*3600*(now-tlastcall)>1,
		tlastcall=now;
		strlastcall=irf_time(now,'date>utc_yyyy-mm-dd HH:MM:SS');
	end
	t_out=strlastcall;
	return
elseif nargin==1,
	flag='vector>epoch';
end
if isempty(t_in),t_out=[];return;end
flagTint=strfind(flag,'tint'); % check if we work with time interval (special case)
if isempty(flagTint),          % default transformation
	flag2=strfind(flag,'2');   % see if flag has number 2 in it
	if isempty(flag2)
		flag2=strfind(flag,'>');   % see if flag has '>' in it
	else
		irf.log('warning',['irf_time(..,''' flag '''): irf_time(..,''in2out'') is deprecated, please use irf_time(..,''in>out'')']);
		flag(flag2)='>';
	end
	if isempty(flag2),         % if no '2' convert from epoch
		format_in='epoch';
		format_out=flag;
		flag=[format_in '>' format_out];
	else
		format_in=flag(1:flag2-1);
		format_out=flag(flag2+1:end);
	end
	if strcmp(format_in,format_out) % if in and out equal return
		t_out=t_in;
		return;
	elseif ~strcmp(format_in,'ttns') && ~strcmp(format_out,'ttns')
		% if there is no epoch in the format then
		% first convert from 'in' to 'epoch'
		% and then from 'epoch' to 'out'
		t_temp=irf_time(t_in,[format_in '>ttns']);
		t_out=irf_time(t_temp,['ttns>' format_out]);
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
		% the last columns can be ommitted, default values MM=1,DD=1,other=0
		t_out = double(irf_time(t_in,'vector>ttns')/1e9);
	case 'vector>ttns'
		% Convert a [YYYY MM DD hh mm ss ms micros ns] to TT2000 in int64 ns
		% the last columns can be ommitted, default values MM=1,DD=1,other=0
		nCol = size(t_in,2);
		if nCol > 9,
			error('irf_time:vector>tt:badInputs',...
				'input should be column vector [YYYY MM DD hh mm ss ms micros ns], last columns can be ommitted.')
		elseif nCol < 9
			defValues = [2000 1 1 0 0 0 0 0 0];
			t_in = [t_in repmat(defValues(nCol+1:9),size(t_in,1),1)];
		end
		t_out = spdfcomputett2000(t_in);
	case 'vector6>ttns'
		% Convert a [YYYY MM DD hh mm ss.xxxx] to TT2000 in int64 ns
		nCol = size(t_in,2);
		if nCol ~= 6,
			error('irf_time:vector6>ttns:badInputs',...
				'input should be column vector with 6 columns [YYYY MM DD hh mm ss.xxx].')
		end
		tSecRound = floor(t_in(:,6));
		tmSec = 1e3*(t_in(:,6)-tSecRound);
		tmSecRound = floor(tmSec);
		tmicroSec = 1e3*(tmSec - tmSecRound);
		tmicroSecRound = floor(tmicroSec);
		tnSecRound = floor(1e3*(tmicroSec - tmicroSecRound));
		t_in(:,6:9) = [tSecRound tmSecRound tmicroSecRound tnSecRound];
		t_out = spdfcomputett2000(t_in);
	case {'ttns>utc','ttns>iso'}
		if any(strfind(flag,'iso')),
			irf.log('warning','irf_time: ''iso'' is deprecated and will be removed, please use ''utc'', see help.');
		end
		t_out = GenericTimeArray.ttns2utc(t_in);
	case 'tt>ttns'
		t_out = int64(t_in)*1e9;
	case 'ttns>tt'
		t_out = double(t_in)/1e9;
	case 'ttns>vector9'
		t_out = spdfbreakdowntt2000(t_in);
	case 'ttns>vector'
		tVec9 = spdfbreakdowntt2000(t_in);
		t_out = tVec9(:,1:6);
		t_out(:,6) = tVec9(:,6)+1e-3*tVec9(:,7)+1e-6*tVec9(:,8)+1e-9*tVec9(:,9);
	case {'utc>ttns','iso>ttns'}
		if any(strfind(flag,'iso')),
			irf.log('warning','irf_time: ''iso'' is deprecated and will be removed, please use ''utc'', see help.');
		end
		if any(strfind(t_in(1,:),'T'))
			t_out = GenericTimeArray.utc2ttns(t_in);
		else
			mask = '%4d-%2d-%2d %2d:%2d:%f%*c';
			s=t_in';
			s(end+1,:)=sprintf(' ');
			a = sscanf(s,mask);
			N = numel(a)/6;
			if N~=fix(N) || N~=size(s,2),
				irf.log('warning','something is wrong with iso input format, returning empty!'),
				t_out=[];
				return;
			end
			a = reshape(a,6,N);
			a = a';
			t_out = irf_time(a,'vector6>ttns');
		end
	case 'ttns>epoch'
		t_out = toepoch(irf_time(t_in,'ttns>vector'));
	case 'epoch>ttns'
		t_out = EpochUnix.to_ttns(t_in);
	case {'ttns>date','ttns>datenum'} % matlab date
		t_out = spdftt2000todatenum(t_in);
	case {'date>ttns','datenum>ttns'}
		t_out = spdfdatenumtott2000(t_in);
	case 'ttns>doy'
		ttBreak = irf_time(t_in,'ttns>vector');
		ttBreak(:,2:end) = 1;
		t_out = [ttBreak(:,1) ...
			floor(irf_time(t_in,'ttns>datenum'))-floor(irf_time(ttBreak,'vector>datenum'))+1];
		
	case 'doy>ttns'
		t_out=irf_time([t_in(:,1) t_in(:,1).*0+1 t_in(:,1).*0+1 ...
			t_in(:,2).*24-12 t_in(:,1).*0 t_in(:,1).*0],'vector6>ttns');
		
	case 'cdfepoch>ttns'
        sz = size(t_in);
        if sz(2)>1
            irf.log('warning','irf_time(cdfepoch>ttns: input is not column vector! output is column vector!');
            t_in = t_in(:);
        end
        ttBreak = spdfbreakdownepoch(t_in); % cdfepoch, not other "epoch"
        % Assume 0 microsec and 0 nanosec.
        t_out = spdfcomputett2000([ttBreak zeros(size(ttBreak,1),2)]);
	case 'ttns>cdfepoch'
        ttBreak = spdfbreakdowntt2000(t_in);
        % Should round up/down column 7, depending on column 8 to 9
        ttBreak(:,7) = round(ttBreak(:,7)+ttBreak(:,8)/10^3+ttBreak(:,9)/10^6);
        t_out = spdfcomputeepoch(ttBreak(:,1:7)); % cdfepoch, not other "epoch"
	case 'cdfepoch16>ttns'
        ttBreak = spdfbreakdownepoch16(t_in);
        % Should round up/down column 9, depending on column 10
        ttBreak(:,9) = round(ttBreak(:,9)+ttBreak(:,10)/10^3);
        t_out = spdfcomputett2000(ttBreak(:,1:9));
	case 'ttns>cdfepoch16'
        ttBreak = spdfbreakdowntt2000(t_in);
        % Assume 0 picoseconds.
        t_out = spdfcomputeepoch16([ttBreak zeros(size(ttBreak,1),1)]);

	case 'ttns>epochtt'
		t_out = EpochTT(t_in);
	case 'epochtt>ttns'
		t_out = t_in.ttns;
		%
		% Time interval conversions
		%
		
	case {'tint>utc','tint>iso'}
		if any(strfind(flag,'iso')),
			irf.log('warning','irf_time(..,''tint>iso'') is deprecated and will be removed, please use irf_time(..,''tint>utc'').');
		end
		t1iso=irf_time(t_in(:,1),'epoch>utc');
		t2iso=irf_time(t_in(:,2),'epoch>utc');
		t_out=[t1iso repmat('/',size(t1iso,1),1) t2iso];
	case {'utc>tint','iso>tint'}
		if any(strfind(flag,'iso')),
			irf.log('warning','irf_time(..,''iso>tint'') is deprecated and will be removed, please use irf_time(..,''utc>tint'').');
		end
		% assume column array where each row is interval in iso format
		ii=strfind(t_in(1,:),'/');
		t1=irf_time(t_in(:,1:ii-1),'utc>epoch');
		t2=irf_time(t_in(:,ii+1:end),'utc>epoch');
		if isempty(t1) || isempty(t2)
			t_out=[];
		else
			t_out=[t1 t2];
		end
		
	case 'tint>isoshort'
		t1iso=irf_time(t_in(:,1),'epoch>utc_yyyy-mm-ddTHH:MM:SS.mmmZ');
		t2iso=irf_time(t_in(:,2),'epoch>utc_yyyy-mm-ddTHH:MM:SS.mmmZ');
		t_out=[t1iso repmat('/',size(t1iso,1),1) t2iso];
		
	otherwise
		if numel(flag)>9 && strcmp(flag(1:9),'tint>utc_')
			fmt = flag(10:end);
			t1iso = irf_time(t_in(:,1),['epoch>utc_' fmt]);
			t2iso = irf_time(t_in(:,2),['epoch>utc_' fmt]);
			t1iso(:,end+1)='/';
			t_out = [t1iso t2iso];
		elseif numel(flag)>9 && strcmp(flag(1:9),'ttns>utc_')
			fmt = flag(10:end);
			t_out = GenericTimeArray.ttns2utc(t_in,fmt);
		else
		disp(['!!! irf_time: unknown flag ''' lower(flag) ''', not converting.'])
		t_out=t_in;
		end
end
