classdef EpochUnix < GenericTimeArray
  %EpochUnix Class representing UNIX epoch, double seconds since 1970.
  %   Detailed explanation goes here
  %
  % See also GenericTimeArray
  
  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------
  
  methods
    function obj = EpochUnix(inp)
      if nargin==0, return, end
      if isa(inp,'double')
        if min(size(inp))>1
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'Double input (seconds since 1970) must be a column or row vector')
        end
        if size(inp,2)~=1, inp = inp'; end % to column
        obj.epoch = inp;
      elseif isa(inp,'int64')
        if min(size(inp))>1
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'int64 input (nanoseconds since 2000) must be a column or row vector')
        end
        obj.epoch = EpochUnix.from_ttns(inp(:)); % column vector
      elseif isa(inp,'char')
        if GenericTimeArray.validate_utc(inp)
          obj.epoch = iso2epoch(inp);
        else
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
      elseif isa(inp,'GenericTimeArray')
        if isa(inp,'EpochUnix')
          obj = inp;
        else
          obj = EpochUnix(inp.utc);
        end
      else
        error('irf:EpochUnix:EpochUnix:badInputs',...
          'Expected inputs: double (seconds since 1970), int64 (nanoseconds since 2000), or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
      end
    end
    function objOut = plus(obj,arg)
      if isnumeric(arg)
        if isa(arg,'double'), inp = arg;
        elseif isa(arg,'int64'), inp = double(arg)*1e-9;
        else
          error('Input type not defined');
        end
        objOut = obj;
        objOut.epoch = obj.epoch + inp(:);
      end
    end
    function outObj = colon(obj,varargin)
      if nargin == 2 && isa(varargin{1},'GenericTimeArray')
        epoch = obj.start.epochUnix:1:varargin{1}.stop.epochUnix;
      elseif nargin == 3 && isa(varargin{2},'GenericTimeArray') && isnumeric(varargin{1})
        epoch = obj.start.epochUnix:double(varargin{1}):varargin{2}.stop.epochUnix;
      else, error('Invalid input(s)')
      end
      outObj = EpochUnix(epoch);
    end
  end
  
  methods (Static)
    function ttns = to_ttns(epoch)
      % Convert unix epoch to TT2000 (ns)
      vector6        = fromepoch(epoch);
      tSecRound      = floor(vector6(:,6));
      tmSec          = 1e3*(vector6(:,6)-tSecRound);
      tmSecRound     = floor(tmSec);
      tmicroSec      = 1e3*(tmSec - tmSecRound);
      tmicroSecRound = floor(tmicroSec);
      tnSecRound     = floor(1e3*(tmicroSec - tmicroSecRound));
      vector6(:,6:9) = [tSecRound tmSecRound tmicroSecRound tnSecRound];
      ttns           = spdfcomputett2000(vector6);
    end
    function epoch = from_ttns(ttns)
      % Convert TT2000 (ns) to unix epoch
      % Empirically: Only seems to work if ttns is int64.
      s_tmp = spdfencodett2000(ttns(1));
      epoch0 = iso2epoch(s_tmp{:});
      epoch = double(ttns - ttns(1))*1e-9 + epoch0;
      if numel(epoch) == 1 % WHY THIS CEHCK??
        return;
      end
      % Check for leap seconds during the time interval of interest
      lSecs = GenericTimeArray.leap_seconds();
      yyyy = str2double(s_tmp{:}(1:4));
      lSecs = lSecs(lSecs(:,1)>=yyyy,:); lSecs(:,4:6) = 0;
      lSecsEpoch = [];
      if ~isempty(lSecs), lSecsEpoch = toepoch(lSecs); end
      if ~isempty(lSecsEpoch) && any( lSecsEpoch>=epoch(1)-1 & lSecsEpoch< epoch(end))
        % Found: convert via UTC string
        epoch = EpochUnix(GenericTimeArray.ttns2utc(ttns)).epoch;
      end
    end
  end
end

