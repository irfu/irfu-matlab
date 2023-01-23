function res = get_variable(dobj,varName,parent)
%GET_VARIABLE  extract a variable from dataObj with all dependencies
%
%  var = get_variable(dobj,varName)
%
%  Exctracts a varibale from dataobj by recurcively pulling in all
%  dependent variables so that the resulting VAR is self consistent.
%
% See also: dataobj/getv

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,3)
if nargin==2, parent = ''; end
if ~ischar(varName)
  errStr = 'variable name must be a string';
  irf.log('critical',errStr), error(errStr)
end

% Fix variable in matlab format, remove minuses, dots ...
varName = variable_mat_name(varName);
nVars = size(dobj.vars,1);
if ~nVars
  irf.log('warning',['No such variable : ' varName])
  res = []; return
end

iVar=find(sum(strcmp(varName,dobj.vars(:,1:2)),2)>0); % find variable
if isempty(iVar), res = []; return; end % Variale not found

res = dobj.data.(dobj.vars{iVar,1});
res.name = varName;
res.GlobalAttributes = dobj.GlobalAttributes;

% Add all revelant variable attributes to the returned variable
varAtts = dobj.VariableAttributes; varAttNames = fieldnames(varAtts);
allVars = dobj.Variables(:,1);
for iName=1:length(varAttNames)
  iattr = find(strcmpi(dobj.vars{iVar,2},varAtts.(varAttNames{iName})(:,1))==1);
  if isempty(iattr), continue, end
  
  % Copy attribute
  attr = varAtts.(varAttNames{iName}){iattr,2};
  res.(varAttNames{iName}) = attr;
  
  % Char attr may be a dependence on another variable, so we pull it in
  % Some buggy files may have identical key/value pairs, we ignore those
  if ~ischar(attr) || strcmp(attr,varName) || ...
      isempty(intersect(allVars,attr))
    continue
  end
  attrMatName = variable_mat_name(attr);
  % Detect recursive deps
  if strcmp(parent,attrMatName)
    irf.log('warning',['dataobj:get_variable:cyclic_dep',...
      ['Cyclic dependency: ' parent ' <-> ' varName]]);
    continue
  end
  varTmp = get_variable(dobj,attr,varName);
  if ~isempty(varTmp), res.(varAttNames{iName}) = varTmp; end
end % for

end % function


