function res = c_load(vs,cl_id)
%C_LOAD load cluster data
%
% C_LOAD(VS) will attempth to load a varible given by string VS from a MAT
% file in the current directory. Returns is 1 in case of sucess.
%
% C_LOAD(VS,CL_ID) will replace '?' in VS by CL_ID and perform load.
%
% Input:
%	VS - variable string
%
% Example:
%	c_load('diE3')
% c_load('diE?',3) 
%   % Loads variable diE3 from the current directory into workspace.
%
% See also C_DESC, AV_SSUB
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%
error(nargchk(1,2,nargin))

if ~isstr(vs)
	error('Input argument must be a string')
end

if nargin==2 
	if ~isnumeric(cl_id), error('Second input argument must be a number 1..4'),end
end

if cl_id<0 | cl_id>4
	error('Second input argument must be a number 1..4')
end

if nargin==2, vs = av_ssub(vs,cl_id); end

d = c_desc(vs);

% Try to load from file
if exist([d.file '.mat'],'file')
	warning off
	eval(['load -mat ' d.file ' ' vs])
	warning on
end

% Return the result
if exist(vs,'var')
	assignin('caller',vs,eval(vs));
	if nargout>0, res = 1; end
else
	if nargout>0, res = 0; 
	else, c_log('load',['cannot load ' vs])
	end
end
