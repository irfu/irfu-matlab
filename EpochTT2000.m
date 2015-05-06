classdef EpochTT2000 < GenericTimeArray
  %EpochTT2000 Class representing T2000 epoch, nanoseconds since 2000.
  %
	% EpochTT2000(t) - initialize class, where t can be:
	%                   - vector of integer number (int64) of nanoseconds as TT2000
	%                   - UTC string array
	
% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

  methods
    function obj = EpochTT2000(inp)
      if nargin==0, return, end
      if isa(inp,'int64'),
        if min(size(inp))>1
          error('irf:EpochTT2000:EpochTT2000:badInputs',...
            'int64 input (nanoseconds since 2000) must be a columt or row vector')
        end
        obj.epoch = inp(:); % column vector
      elseif isa(inp,'char')
        if ~GenericTimeArray.validate_utc_time_str(inp)
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
        tmpStr = pad_utc(inp);
        obj.epoch = spdfparsett2000(tmpStr);
        if obj.epoch==int64(-9223372036854775805)
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
      else
        error('irf:EpochUnix:EpochUnix:badInputs',...
          'Expected inputs: int64 (nanoseconds since 2000), double (seconds since 1970) or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
      end
      function utcNew = pad_utc(utc)
        MAX_NUM_IDX = 29; IDX_DOT = 20;
        utcNew = utc; lUtc = size(utc,2); flagAddZ = true;
        if all(all(utc(:,end)=='Z')), lUtc = lUtc - 1; flagAddZ = false; end
        if lUtc == MAX_NUM_IDX && ~flagAddZ, return, end
        if lUtc == IDX_DOT-1, utcNew(:,IDX_DOT) = '.'; lUtc = lUtc + 1; end
        if lUtc < MAX_NUM_IDX, utcNew(:,(lUtc+1):MAX_NUM_IDX) = '0'; end % Pad with zeros
        utcNew(:,MAX_NUM_IDX+1) = 'Z';
      end
    end
    function s = toUtc(obj,format)
      % s = toUtc(obj,format)
      if nargin<2, format = 2; end
      s_tmp = char(spdfencodett2000(obj.epoch));
      switch format
        case 0, s_tmp = s_tmp(:,1:26);
        case 1, s_tmp = s_tmp(:,1:23);
        case 2,
        otherwise
          error('irf:EpochUnix:toUtc:badInputs',...
            'wrong format value')
      end
      s_tmp(:,end+1) = 'Z';
      s = s_tmp;
    end
    
    function res = toEpochUnix(obj)
      s_tmp = spdfencodett2000(obj.epoch(1)); epoch0 = iso2epoch(s_tmp{:});
      epoch = double(obj.epoch - obj.epoch(1))*1e-9 + epoch0;
      if numel(epoch) == 1, res = EpochUnix(epoch); return; end
      
      % Check for leap seconds during the time interval of interest
      lSecs = GenericTimeArray.leap_seconds();
      yyyy = str2double(s_tmp{:}(1:4));
      lSecs = lSecs(lSecs(:,1)>=yyyy,:); lSecs(:,4:6) = 0;
      lSecsEpoch = [];
      if ~isempty(lSecs), lSecsEpoch = toepoch(lSecs); end
      if ~isempty(lSecsEpoch) && any( lSecsEpoch>=epoch(1)-1 & lSecsEpoch< epoch(end))
        % Found: convert via UTC string
        res = EpochUnix(toUtc(obj));
      else res = EpochUnix(epoch);
      end
    end
  end
end