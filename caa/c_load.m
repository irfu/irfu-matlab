function [res,v] = c_load(vs,cl_id,mode_s)
%C_LOAD load cluster data
%
% C_LOAD(V_S) will attempth to load a varible given by string VS from a MAT
%   file in the current directory.
%
% C_LOAD(V_S,CL_ID) will replace '?' in VS by CL_ID and perform load.
%
% [RES,VAR] = C_LOAD(V_S,[CL_ID]) in case of success returns RES=1 and
%   loaded variable in VAR. Otherwise RES=0 and VAR=[].
%
% RES = C_LOAD(V_S[,CL_ID,MODE_S]) returns variable in RES if MODE_S 
%   is set to 'var';
%
% Input:
%	V_S - variable string
% CL_ID - Cluster id
% MODE_S - string 'var' or 'res' (default) defining the output in case of
%   one input argument.
%
% Examples:
%	c_load('diE3')
% c_load('diE?',3)
% ok = c_load('diE?',3)
%   % Loads variable diE3 from the current directory into workspace. 
%   % ok is 1 if load was sucessfull. 
% B_tmp = c_load('B2','var');
% B_tmp = c_load('B?',2,'var');
% [ok,B_tmp] = c_load('B2');
% [ok,B_tmp] = c_load('B?',2);
%   % Loads variable B2 from the current directory into variable B_tmp. 
%   % ok is 1 if load was sucessfull. 
%
% See also C_DESC, AV_SSUB
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%
error(nargchk(1,3,nargin))
if nargout==2 & nargin==3
	error('Invalid number of input and output arguments. See HELP C_LOAD')
end

if ~isstr(vs), error('V_S must be a string'), end

switch nargin
case 3
	if ~isnumeric(cl_id), error('Second input argument must be a number 1..4'),end
	if cl_id<0 | cl_id>4
		error('CL_ID must be in a range 1..4')
	end
	vs = av_ssub(vs,cl_id);
case 2
	if isnumeric(cl_id)
		if cl_id<0 | cl_id>4
			error('CL_ID must be in a range 1..4')
		end
		vs = av_ssub(vs,cl_id);
		mode_s = 'res';
	elseif isstr(cl_id), mode_s = cl_id;
	else, error('Second input argument must be eather a number 1..4 or a string.')
	end
case 1
	mode_s = 'res';
end

if strcmp(mode_s,'var'), ret_var = 1;
elseif strcmp(mode_s,'res'), ret_var = 0;
else, c_log('fcal','Invalid value of MODE_S. Defaulting to ''res''')
end

d = c_desc(vs);

% Try to load from file
if exist([d.file '.mat'],'file')
	warning off
	eval(['load -mat ' d.file ' ' vs])
	warning on
end

% Return the result
if exist(vs,'var')
	switch nargout
	case 2
		res = 1;
		v = eval(vs);
	case 1
		if ret_var, res = eval(vs);
		else, res = 1; assignin('caller',vs,eval(vs));
		end
	case 0
		assignin('caller',vs,eval(vs));
	end
else
	switch nargout
	case 2
		res = 0;
		v = [];
	case 1
		if ret_var, res = [];
		else, res = 0;
		end
	case 0
		c_log('load',['cannot load ' vs])
	end
end
