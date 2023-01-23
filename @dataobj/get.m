function res = get(dobj,var_s)
%GET(dobj, var_s)  get VariableAttributes, GlobalAttributes or a variable
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

if ~ischar(var_s), error('VAR_S must be a string'), end

switch var_s
  case {'va','VariableAttributes'}
    res = dobj.VariableAttributes;
  case {'ga','GlobalAttributes'}
    res = dobj.GlobalAttributes;
  otherwise
    res = getv(dobj,var_s);
end