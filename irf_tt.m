function out=irf_tt(varargin)
% IRF_TT work with time tables IN DEVELOPMENT!!!!
%   time table tt is a structure with tt.start being vector with start times in epoch and
%   tt.end being vector with end times in epoch. tt.description includes general description
%   of the time table. There can be other fields to structure, those are added after #
%   when exporting time table. Any comments should be put in cell array
%   tt.comment that should have the same size.
%
%   IRF_TT(tt) display time table on screen
%
%   ttascii=IRF_TT(tt,'export_ASCII');	% export time table in ascii format
%	IRF_TT(tt,'export_ASCII',filename); % save time table in ascii format to	file
%
%   tt=IRF_TT('import_ASCII',filename);% import time table from file
%
%   tt=IRF_TT('read_www',http_address) read time table from http_address
%
%   tt=IRF_TT('read_IRF',tt_id) read time table tt_id from IRF
%
%   IRF_TT(tt,'write_IRF',tt_id) write time table tt_id to IRF
%
%   IRF_TT('list_IRF') list time tables in IRF repository
%
%  =========== BELOW NOT IMPLEMENTED =========
%
%   [tt]=IRF_TT('read_AMDA',tt_id) read time table tt_id from AMDA
%
%   hashid=IRF_TT(tt,'save_AMDA',tt_id,'user',username,'passw',password)
%     save time table tt_id. Return hashid, if saving did not succeed return NaN.
%
%   out=IRF_TT(tt,'save_AMDA',hashid) replace time table with hashid with a new one
%    return NaN if operation did not succeed.
%
%   out=IRF_TT(tt,'add_AMDA',hashid) add time intervals to the table hashid
%    return NaN if operation did not succeed.
%
%   Example:
%       TT=struct('start',irf_time([2000 1 1 0 0 0]) + [0 10],...
%        'end',irf_time([2000 1 1 0 0 0]) + [1 11],'description','two second intervals');
%        irf_tt(TT); % display time table
%		tt=irf_tt('read_IRF','Example');
%		irf_tt(tt,'write_IRF','new_id'); % write time table tt to brain.irfu.se with new id 'new_id'
%		% read Cluster Master Science Plan
%		msp=irf_tt('read_www','http://jsoc1.bnsc.rl.ac.uk/msp/full_msp_ascii.lst');

% $Id$
% $Revision$  $Date$

%% Check inputs

if nargin==0,
	help irf_tt;
	return;
end

if isTT(varargin{1}) % first argument is TT in matlab format (structure)
	tt=varargin{1};
	if nargin==1 && nargout == 0, % IRF_TT(tt) display time table on screen
		action='display';
	elseif nargin>1 && ischar(varargin{2}), % second argument is action
		action=varargin{2};
		if strcmpi(action,'export_ASCII') && nargin >2 % IRF_TT(tt,'export_ASCII',filename);
			filename=varargin{3};
			if ~ischar(filename),
				disp('ERROR! dont understand irf_tt arguments. See help!');
				return
			end
		elseif strcmp(action,'write_IRF'), % [out]=IRF_TT(tt,'write_IRF',tt_id)
			if nargin==2,
				disp('ERROR! dont understand irf_tt arguments. See help!');
				return
			else
				tt_id=varargin{3};
				if ~ischar(varargin{3}), % filename should be character array
					disp('ERROR! dont understand irf_tt arguments. See help!');
					return
				end
			end
		end
	end
elseif ischar(varargin{1}) % first argument is action
	action=varargin{1};
	if nargin>1 && ischar(varargin{2}) && strcmpi(varargin{1},'import_ASCII')
		filename=varargin{2};
	end
	if nargin>1 && ischar(varargin{2}) && ...
			(strcmpi(varargin{1},'read_IRF') || strcmpi(varargin{1},'read_www'))
		tt_id=varargin{2};
	end
end

%% Act
switch action
	case 'display'
		display(asciiTT(tt));
	case 'export_ASCII'
		out=asciiTT(tt);
		if exist('filename','var') % read from file
			fid = fopen(filename,'w');
			fwrite(fid,out);
			fclose(fid);
		end
	case 'read_IRF'
		remoteFile=['brain.irfu.se:/share/Cluster/TT/' tt_id];
		tempFile=tempname;
		eval(['!scp ' remoteFile ' ' tempFile]);
		out=irf_tt('import_ASCII',tempFile);
		delete(tempFile);
	case 'read_www'
		tempTT=urlread(tt_id);
		out=tt_from_ascii(tempTT);
	case 'write_IRF'
		remoteFile=['brain.irfu.se:/share/Cluster/TT/' tt_id];
		tempFile=tempname;
		irf_tt(tt,'export_ASCII',tempFile);
		eval(['!scp ' tempFile ' ' remoteFile]);
		delete(tempFile);
	case 'list_IRF'
		eval('!ssh brain.irfu.se ls /share/Cluster/TT');
	case 'import_ASCII'
		if exist('filename','var') % read from file
			out=readasciiTT(filename);
		else
			disp('ERROR! dont understand irf_tt arguments. See help!');
			return
		end
	otherwise
		disp('ERROR! dont understand irf_tt arguments. See help!');
		return
end

if nargout==0,
	clear out;
end

end


%% functions
function out=isTT(x) % check if input is time table structure
if isstruct(x) && ... % is structure
		isfield(x,'start') && isfield(x,'end') &&... % start and end field should be defined
		isnumeric(x.start) && isnumeric(x.end) && ...	% start and end should be numeric
		any(size(x.start) == size(x.end)) && ...		% start and end fields should be the same size
		min(size(x.start)) == 1 && min(size(x.end))==1 % start end fields should be vectors
	out = 1;
else
	out=0;
end
end
function out=asciiTT(tt) % tt in ascii format
out=[];
if isfield(tt,'description'),
	if ischar(tt.description),
		description={tt.description}; % make description cell array
	else
		description=tt.description;
	end
end
if isfield(tt,'comment'),
	isComment=1;
	comment=tt.comment;
else
	isComment=0;
end
if numel(description)<numel(tt.start) % fill with empty
	for j=(numel(description)+1):numel(tt.start),
		description{j}=[];
	end
end
description=regexprep(description,'\n','\n# '); % new lines in description should be commented

for iTT=1:numel(tt.start)
	if ~isempty(description{iTT}),
		out=[out sprintf('# %s\n',description{iTT})];
	end
	if isComment && ~isempty(comment{iTT}),
		if numel(comment{iTT})>1 && any(strcmp(comment{iTT}(1),'#')),
			fmt='%s %s %s\n';
		else
			fmt='%s %s #%s\n';
		end
		out=[out sprintf(fmt,irf_time(tt.start(iTT),'iso'),...
			irf_time(tt.end(iTT),'iso'),comment{iTT})];
	else
		out=[out sprintf('%s %s\n',irf_time(tt.start(iTT),'iso'),...
			irf_time(tt.end(iTT),'iso'))];
	end
end
end
function out=readasciiTT(filename) % read in tt from ascii format
if ischar(filename) && exist(filename,'file') % is ok
else % not ok
	out = NaN;
end
fid = fopen(filename);
fileContents=fscanf(fid,'%c',inf);
fclose(fid);
ttAscii=textscan(fileContents,'%s','delimiter',sprintf('\n'));
out=tt_from_ascii(ttAscii);
end

function out=tt_from_ascii(ttascii) % convert ascii tt to matlab format

if ~iscell(ttascii) && ischar(ttascii), % character array
	ttAscii=textscan(ttascii,'%s','delimiter',sprintf('\n'));
else
	ttAscii=ttascii; % input is already cell array with one line per cell
end

nTimeInterval			= 0;
isFirstDescriptionLine	= 1;
description				= '';

for iTtAscii=1:numel(ttAscii{1})
	str=ttAscii{1}{iTtAscii};
	if ~isempty(regexp(str,'^\s*[#,!]','match')) % commented line
		if isFirstDescriptionLine,
			isFirstDescriptionLine=0;
		end
		linetext=regexp(str,'^\s*[#,!]\s?(.*)','tokens');
		if isempty(description)
			description=[description linetext{1}{1}];
		else
			description=[description sprintf('\n') linetext{1}{1}];
		end
	else % assume number in iso format
		timeInterval=regexp(str,'^\s*(?<start>[\d-]*T[\d:\.]*Z+)\s*(?<end>[\d-]*T[\d:\.]*Z+)\s?(.*)','names');
		if ~isempty(timeInterval) && isfield(timeInterval,'start') && isfield(timeInterval,'end')
			nTimeInterval=nTimeInterval+1;
			out.start(nTimeInterval)		= irf_time(timeInterval.start,'iso2epoch');
			out.end(nTimeInterval)			= irf_time(timeInterval.end,'iso2epoch');
			if isfield(timeInterval,'comment'),
				out.comment{nTimeInterval}		= timeInterval.comment;
			else
				out.comment{nTimeInterval}		= [];
			end
			out.description{nTimeInterval}	= description;
			isFirstDescriptionLine			= 1;
			description						= '';
		else
			description=[description sprintf('\n') str];
		end
	end
end
end