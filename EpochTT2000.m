classdef EpochTT2000 < GenericTimeArray
  %EpochTT2000 Class representing T2000 epoch, nanoseconds since 2000.
  %   Detailed explanation goes here
  methods
    function obj = EpochTT2000(inp)
      if nargin==0, return, end
      if isa(inp,'int64'), 
        if min(size(inp))>1
          error('irf:EpochTT2000:EpochTT2000:badInputs',...
            'int64 input (nanoseconds since 2000) must be a columt or row vector')
        end
        if size(inp,2)~=1, inp = inp'; end % to column
        obj.epoch = inp;
      elseif isa(inp,'char')
        if GenericTimeArray.validate_iso_time_str(inp)
          obj.epoch = parsett2000(inp);
        else
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
      else
        error('irf:EpochUnix:EpochUnix:badInputs',...
          'Expected inputs: double (seconds since 1970) or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
      end
    end
    function s = toUtc(obj,format)
      % s = toUtc(obj,format)
      if nargin<2, format = 2; end
      s_tmp = char(encodett2000(obj.epoch));
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
  end
end