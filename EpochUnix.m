classdef EpochUnix < GenericTimeArray
%EpochUnix Class representing UNIX epoch, double seconds since 1970.
%   Detailed explanation goes here
%
% See also GenericTimeArray
  
% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

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
        if GenericTimeArray.validate_utc_time_str(inp)
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
      %toUtc  convert to UTC string
      if nargin<2, format = 0; end
      s = epoch2iso(obj.epoch, format);
    end
  end
end

