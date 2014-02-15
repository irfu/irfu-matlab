classdef EpochUnix < GenericTimeArray
  %EpochUnix Class representing UNIX epoch, double seconds since 1970.
  %   Detailed explanation goes here
  
  methods
    function obj = EpochUnix(inp)
      if nargin==0, return, end
      if isa(inp,'double'), 
        if min(size(inp))>1
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'Double input (seconds since 1970) must be a columt or row vector')
        end
        if size(inp,2)~=1, inp = inp'; end % to column
        obj.epoch = inp;
      elseif isa(inp,'char')
        if GenericTimeArray.validate_iso_time_str(inp)
          obj.epoch = iso2epoch(inp);
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
      if nargin<2, format = 0; end
      s = epoch2iso(obj.epoch, format);
    end
  end
  
end

