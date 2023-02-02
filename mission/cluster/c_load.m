function [res,v,msg] = c_load(vs,cl_id,mode_s,probe_id)
%C_LOAD load cluster data
%
% C_LOAD(V_S) will attempth to load a varible given by string VS from a MAT
%   file in the current directory.
%
% C_LOAD(V_S,SC_LIST) will replace '?' in VS by SC_LIST and perform load.
%
% [RES,VAR,MSG] = C_LOAD(V_S,[SC_LIST]) in case of success returns RES=1 and
%	loaded variable in VAR. Otherwise RES=0 and VAR=[]. If SC_LIST has more
%	than one entry, RES will be an array and VAR a cell array.
%	MSG contains an error message
%
% RES = C_LOAD(V_S[,SC_LIST,MODE_S])
%     returns variable in RES if MODE_S	is set to 'var';
%
% [RES,V,MSG] = C_LOAD(V_S,CL_ID,'res',PROBE_LIST])
%     Expand '?' in V_S using PROBE_LIST, or if both '?' and '!' are
%     present in V_S, expand '?' using CL_ID; and '!' using PROBE_LIST
%
% Input:
%	V_S - variable string
% CL_ID - Cluster id
% MODE_S - string 'var' or 'res' (default) defining the output in case of
%	one input argument.
%
% Examples:
%	c_load('diE3')
%	c_load('diE?',3)
%	ok = c_load('diE?',3)
%	% Loads variable diE3 from the current directory into workspace.
%	% ok is 1 if load was sucessfull.
%	c_load('diE?')
%	c_load('diE?',1:4)
%	ok = c_load('diE?')
%	% Loads variables diE1..4 from the current directory into workspace.
%	% ok is [1 1 1 1] if load was sucessfull.
%	B_tmp = c_load('B2','var');
%	B_tmp = c_load('B?',2,'var');
%	[ok,B_tmp] = c_load('B2');
%	[ok,B_tmp] = c_load('B?',2);
%	% Loads variable B2 from the current directory into variable B_tmp.
%	% ok is 1 if load was sucessfull.
%
% See also C_DESC, IRF_SSUB

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

%error(nargchk(1,3,nargin))
narginchk(1,4)
if nargout==2 && nargin==3
  error('Invalid number of input and output arguments. See HELP C_LOAD')
end

ERR_RET = -157e8;

if ~ischar(vs), error('V_S must be a string'), end

switch nargin
  case 4
    if ~isnumeric(cl_id) || length(cl_id)>1
      error('Second input argument must be a (scalar) number 1..4')
    end
    if cl_id<0 || cl_id>4
      error('CL_ID must be in a range 1..4')
    end
    if (strcmp(mode_s, 'res') == 0)
      error('mode_s must be ''res''.')
    end
    if regexp(vs,'?')
      if regexp(vs, '!')
        vs = irf_ssub(vs, cl_id);     % Change '?' --> cl_id and '!' --> '?' in vs ...
        vs = strrep(vs, '!', '?');    % ...and '!' --> '?' in vs
      end
      cl_id = probe_id;    % Fool rest of program into thinking probe_id is cl_id.
    end
  case 3
    if ~isnumeric(cl_id), error('Second input argument must be a number 1..4'),end
    if cl_id<0 || cl_id>4
      error('CL_ID must be in a range 1..4')
    end
  case 2
    if isnumeric(cl_id)
      if any(find(cl_id<0)) || any(find(cl_id>4))
        error('CL_ID must be in a range 1..4')
      end
      mode_s = 'res';
    elseif ischar(cl_id)
      mode_s = cl_id;
      if regexp(vs,'?'), cl_id = 1:4;
      else, cl_id = 1;
      end
    else, error('Second input argument must be either a number 1..4 or a string.')
    end
  case 1
    mode_s = 'res';
    if regexp(vs,'?'), cl_id = 1:4;
    else, cl_id = 1;
    end
end

if strcmp(mode_s,'var'), ret_var = 1;
elseif strcmp(mode_s,'res'), ret_var = 0;
else, irf.log('warning','Invalid value of MODE_S. Defaulting to ''res''')
end

kk = 1;
for cli=cl_id
  vs_tmp = irf_ssub(vs,cli);
  d = c_desc(vs_tmp);
  if isempty(d), error(['Unknown variable: ' vs_tmp]), end
  % Try to load from file
  if exist(['./' d.file '.mat'],'file')
    warning off %#ok<WNOFF>
    eval(['load -mat ./' d.file '.mat ' vs_tmp])
    warning on %#ok<WNON>
  elseif exist(['./' d.file_old '.mat'],'file')
    warning off %#ok<WNOFF>
    eval(['load -mat ./' d.file_old '.mat ' vs_tmp])
    warning on %#ok<WNON>
  end
  
  % Return the result
  if exist(vs_tmp,'var')
    switch nargout
      case {2,3}
        res(kk) = 1; %#ok<AGROW>
        if length(cl_id)>1, v(kk) = {eval(vs_tmp)};  %#ok<AGROW>
        else, v = eval(vs_tmp);
        end
        msg = ['sucessfully loaded ' vs_tmp];
      case 1
        if ret_var
          if length(cl_id)>1, res(kk) = {eval(vs_tmp)}; %#ok<AGROW>
          else, res = eval(vs_tmp);
          end
        else, res(kk) = 1; assignin('caller',vs_tmp,eval(vs_tmp)); %#ok<AGROW>
        end
      case 0
        assignin('caller',vs_tmp,eval(vs_tmp));
    end
  else
    switch nargout
      case {2,3}
        res(kk) = 0;  %#ok<AGROW>
        if length(cl_id)>1, v(kk) = {[]}; %#ok<AGROW>
        else, v = [];
        end
        if d.lev, OBJ = 'ClusterProc';
        else, OBJ = 'ClusterDB';
        end
        msg = ['No ' vs_tmp ' in ' d.file '. Use getData(' OBJ...
          ',...,cl_id,''' d.quant ''')'];
      case 1
        if ret_var
          if length(cl_id)>1, res(kk) = {ERR_RET};  %#ok<AGROW>
          else, res = ERR_RET;
          end
        else, res(kk) = 0;  %#ok<AGROW>
        end
      case 0
        irf.log('warning',['cannot load ' vs_tmp])
    end
  end
  kk = kk + 1;
end
