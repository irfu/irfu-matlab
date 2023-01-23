function res = findva(dobj,field,var_s)
%FINDVA(dobj, var_s) find attribute for a variable in a list
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(3,3)

if ~ischar(var_s), error('VAR_S must be a string'), end

res = '';

nvars = size(dobj.vars,1);
if nvars>0
  for v=1:nvars
    if strcmpi(dobj.vars{v,1},var_s) || strcmpi(dobj.vars{v,2},var_s)
      if isfield( dobj.VariableAttributes ,field)
        for vv=1:size( dobj.VariableAttributes.(field) ,1)
          if strcmp( dobj.VariableAttributes.(field){vv,1} , dobj.vars{v,2})
            res = dobj.VariableAttributes.(field){vv,2};
            return
          end
        end
      end
      return
    end
  end
end

error(['No such variable : ' var_s])