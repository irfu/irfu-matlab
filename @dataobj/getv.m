function res = getv(dobj,varName)
%GETV(dobj, varName)  get a variable
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
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
    disp(['No such variable : ' varName])
    res = [];
    return;
  end
  res = dobj.data.(dobj.vars{iVar,1});
  res.name = varName;
  % Add Variable attributes to the returned variable
  variableAttributeNames=fieldnames(dobj.VariableAttributes);
  for j=1:length(variableAttributeNames)
    iattr=find(strcmpi(dobj.vars{iVar,2},...
      dobj.VariableAttributes.(variableAttributeNames{j})(:,1))==1);
    if iattr
      res.(variableAttributeNames{j})=dobj.VariableAttributes.(variableAttributeNames{j}){iattr,2};
    end
  end
  return
end

disp(['No such variable : ' varName])
res = [];
