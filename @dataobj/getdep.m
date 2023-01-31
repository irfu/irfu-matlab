function res = getdep(dobj,var_s)
%GETDEP(dobj, var_s)  get dependencies for a variable

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

if ~ischar(var_s), error('VAR_S must be a string'), end

nvars = size(dobj.vars,1);
if nvars>0
  for v=1:nvars
    if strcmpi(dobj.vars{v,1},var_s) || strcmpi(dobj.vars{v,2},var_s)
      v1_s = dobj.vars{v,2};
      maxdep = length(dobj.data.(dobj.vars{v,1}).variance(3:end));
      
      dep_x = cell(maxdep,2);
      found_any = 0;
      for d=1:maxdep
        found = 0;
        
        field = sprintf('LABEL_%d',d);
        tt = findva(dobj,field,v1_s);
        if ~isempty(tt)
          dep_x(d,:) = {tt,field};
          found = 1;
          found_any = 1;
        end
        
        field = sprintf('DEPEND_%d',d);
        tt = findva(dobj,field,v1_s);
        if ~isempty(tt)
          if found
            irf.log('warning','found both LABEL_X and DEPEND_X'),
          else
            dep_x(d,:) = {tt,field};
            found_any = 1;
          end
        end
      end
      if dobj.data.(dobj.vars{v,1}).variance(1) == 'T'
        if isfield(dobj.VariableAttributes,'DEPEND_0')
          va = dobj.VariableAttributes.DEPEND_0;
          tvar = [];
          for vv=1:size(va,1)
            if strcmp(va{vv,1},v1_s)
              tvar = getv(dobj,va{vv,2});
            end
          end
          if isempty(tvar)  % can be time variable itself because there is no DEPEND_O for a T variable
            res=[];
            irf.log('warning','No DEPEND_O for a T variable');
            return;
          end
        end
        if found_any
          res = struct('DEPEND_O',tvar,'DEPEND_X',{dep_x});
        else
          res = struct('DEPEND_O',tvar,'DEPEND_X',[]);
        end
      else
        if found_any, res = struct('DEPEND_X',{dep_x});
        else, res = struct('DEPEND_X',[]);
        end
      end
      return
    end
  end
end

errStr = ['No such variable : ' var_s];
irf.log('critical',errStr), error(errStr)

