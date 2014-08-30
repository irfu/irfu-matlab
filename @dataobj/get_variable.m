function res = get_variable(dobj,varName)
%GETV(dobj, varName)  get a variable
%

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)

if ~ischar(varName)
  error('variable name must be a string')
end

% fix variable in matlab format, remove minuses, dots ...
varName = variable_mat_name(varName);

nVars = size(dobj.vars,1);
if nVars>0
  iVar=find(sum(strcmp(varName,dobj.vars(:,1:2)),2)>0); % find variable
  if isempty(iVar)
    irf.log('notice',['No such variable : ' varName])
    res = [];
    return;
  else
    res = dobj.data.(dobj.vars{iVar,1});
    res.name = varName;
    res.GlobalAttributes = dobj.GlobalAttributes;
    % Add Variable attributes to the returned variable
    variableAttributeNames=fieldnames(dobj.VariableAttributes);
    for j=1:length(variableAttributeNames),
      iattr=find(strcmpi(dobj.vars{iVar,2},...
        dobj.VariableAttributes.(variableAttributeNames{j})(:,1))==1);
      if iattr,
        attr = dobj.VariableAttributes.(variableAttributeNames{j}){iattr,2};
        if ischar(attr), 
			if strcmp(attr,varName),
        irf.log('critical',['Variable depends on itself : ' varName])
				error('Recursive variable dependence');
			end
			varTmp = get_variable(dobj,attr);
        else varTmp = [];
        end
        if isempty(varTmp), res.(variableAttributeNames{j})= attr;
        else res.(variableAttributeNames{j}) = varTmp;
        end
      end
    end
    return
  end
end

disp(['No such variable : ' varName])
res = [];
