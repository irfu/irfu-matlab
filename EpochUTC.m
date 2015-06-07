classdef EpochUTC < GenericTimeArray
	%EpochUTC Class representing Terrestrial Time epoch
	%
	% EpochTT(t) - initialize class whith time vector being char array 
	%               - UTC string array
	
	methods
		function obj = EpochTT(inp)
			if nargin==0, return, end
			if isa(inp,'double'),
				if min(size(inp))>1
					error('irf:EpochTT:EpochTT:badInputs',...
						'input must be a column or row vector')
				end
				obj.epoch = EpochUTC.from_ttns(int64(inp(:)*1e9)); 
			elseif isa(inp,'int64'),
				if min(size(inp))>1
					error('irf:EpochTT:EpochTT:badInputs',...
						'input must be a column or row vector')
				end
				obj.epoch = EpochUTC.from_ttns(inp(:)); % column vector
			elseif isa(inp,'char')
				if ~GenericTimeArray.validate_utc_time_str(inp)
					error('irf:EpochUnix:EpochUnix:badInputs',...
						'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
				end
				obj.epoch = inp;
			elseif isa(inp,'GenericTimeArray')
				if isa(inp,'EpochUTC'),
					obj = inp;
				else
					obj = EpochUTC(inp.ttns);
				end
			else
				error('irf:EpochUnix:EpochUnix:badInputs',...
					'Expected inputs: int64 (nanoseconds since 2000), double (seconds since 1970) or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
			end
		end
		function out = epochUnix(obj)
			out = iso2epoch(obj.epoch);
		end
		function s = toUtc(obj,varargin)
			s = utc(obj,varargin{:});
		end
		function s = utc(obj,format)
			% s = utc(obj,format)
			if nargin<2,
				format = '';
			else
				format = ['_' format];
			end
			s = irf_time(obj.epoch,['ttns>utc' format]);
		end
		function s = tts(obj,index)
			% s = tt(obj,index)
			% return index points, if not given return all
			if nargin == 1,
				s = double(obj.epoch)/1e9;
			elseif nargin == 2 && isnumeric(index),
				s = double(obj.epoch(index))/1e9;
			end
		end
		function s = ttns(obj,index)
			% s = ttns(obj,index)
			if nargin == 1,
				s = obj.epoch;
			elseif nargin == 2 && isnumeric(index),
				s = obj.epoch(index);
			end
		end
		
		function objOut = plus(obj,arg)
			if isnumeric(arg)
				if isa(arg,'double'),
					inp = int64(arg*1e9);
				elseif isa(arg,'int64'),
					inp = arg;
				else
					error('Input type not defined');
				end
				objOut = obj;
				objOut.epoch = obj.epoch + inp(:);
			end
		end
		function outObj = colon(obj,varargin)
			if nargin == 2 && isa(varargin{1},'EpochTT')
				tns = obj.start.ttns:int64(1e9):varargin{1}.stop.ttns;
				outObj = EpochTT(tns);
			elseif nargin == 3 && isa(varargin{2},'EpochTT') && isnumeric(varargin{1})
				tns = obj.start.ttns:int64(varargin{1}*1e9):varargin{2}.stop.ttns;
				outObj = EpochTT(tns);
			end
		end
	end
	methods (Static)
		function ttns = to_ttns(utc)
			ttns = spdfparsett2000(irf.utc_validate_and_pad(utc));
			if ttns==int64(-9223372036854775805)
				error('EpochUTC:to_ttns',...
					'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
			end
		end
		function utc = from_ttns(ttns,format)
      if nargin<2, format = 2; end
			utc =  char(spdfencodett2000(int64(ttns)));
			if isnumeric(format)
				switch format
					case 0, utc = utc(:,1:26);
					case 1, utc = utc(:,1:23);
					case 2,
					otherwise
						error('irf:EpochUnix:toUtc:badInputs',...
							'wrong format value')
				end
				utc(:,end+1)='Z';
			elseif ischar(format)
				fmt = format;
				iyyyy = strfind(fmt,'yyyy');
				imm   = strfind(fmt,'mm');
				idd   = strfind(fmt,'dd');
				iHH   = strfind(fmt,'HH');
				iMM   = strfind(fmt,'MM');
				iSS   = strfind(fmt,'SS');
				immm  = strfind(fmt,'mmm');
				iuuu  = strfind(fmt,'uuu');
				innn  = strfind(fmt,'nnn');
				tVec9 = irf_time(ttns,'ttns>vector9');
				utc = repmat(fmt,size(ttns,1),1);
				if iyyyy, utc(:,iyyyy:iyyyy+3)= num2str(tVec9(:,1),'%04u'); end
				if imm,   utc(:,imm:imm+1)    = num2str(tVec9(:,2),'%02u'); end
				if idd,   utc(:,idd:idd+1)    = num2str(tVec9(:,3),'%02u'); end
				if iHH,   utc(:,iHH:iHH+1)    = num2str(tVec9(:,4),'%02u'); end
				if iMM,   utc(:,iMM:iMM+1)    = num2str(tVec9(:,5),'%02u'); end
				if iSS,   utc(:,iSS:iSS+1)    = num2str(tVec9(:,6),'%02u'); end
				if immm,  utc(:,immm:immm+2)  = num2str(tVec9(:,7),'%03u'); end
				if iuuu,  utc(:,iuuu:iuuu+2)  = num2str(tVec9(:,8),'%03u'); end
				if innn,  utc(:,innn:innn+2)  = num2str(tVec9(:,9),'%03u'); end
			end
		end
	end
end