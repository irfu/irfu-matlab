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
        tmpStr = GenericTimeArray.pad_utc(inp);
        obj.epoch = spdfparsett2000(tmpStr);
				if obj.epoch==int64(-9223372036854775805)
					error('irf:EpochUnix:EpochUnix:badInputs',...
						'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
				end
			elseif isa(inp,'GenericTimeArray')
				if isa(inp,'EpochTT2000'),
					obj = inp;
				else
					obj = EpochTT2000(inp.ttns);
				end
      else
        error('irf:EpochUnix:EpochUnix:badInputs',...
          'Expected inputs: int64 (nanoseconds since 2000), double (seconds since 1970) or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
      end
		end  
	end
	
	methods (Static)
		function output = from_ttns(input,index) % for consistency with other EpochXX routines
			if nargin == 1,
				output = input;
			else
				output = input(index);
			end
		end
		function output = to_ttns(input,index)
			if nargin == 1,
				output = input;
			else
				output = input(index);
			end
		end
	end
	
end