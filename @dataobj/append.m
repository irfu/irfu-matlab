function dataObject = append(dataObject1,dataObject2)
% APPEND(dobj1,dobj2)  appends dataobject dobj1 to dobj2
%
% dobj = APPEND(dobj1,dobj2) returns resulting dataobject
%
% Routine is checking if the meta information is the same
% in which case it just appends the values of dobj2 to dobj1
%

% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(nargin,2,2))

if nargout==1, dataObject=[];end % default empty output

if ~isa(dataObject1,'dataobj') || ~isa(dataObject2,'dataobj')
	irf_log('fcal','Both input parameters should be dataobj');
	return;
end

% compare if dataobjects have the same variables
ok=compare_cell_arrays(dataObject1.Variables(:,[1 2 4 5 6]),dataObject2.Variables(:,[1 2 4 5 6]));
if ~ok, 
	irf_log('fcal','Databojects do not have the same variables');
	return;
end

% should compare also attributes, but thats for later TODO

% append variable values

dataObject = dataObject1;
variableNameArray=fieldnames(dataObject.data);
for iVariable=1:numel(variableNameArray)
	variableName=variableNameArray{iVariable};
	if strcmp(dataObject.data.(variableName).variance(1),'T') % variable changes with time
		dataObject.data.(variableName).nrec=dataObject1.data.(variableName).nrec + dataObject2.data.(variableName).nrec;
		dataObject.data.(variableName).data=[dataObject1.data.(variableName).nrec ; dataObject2.data.(variableName).nrec];		
		dataObject.Variables{iVariable,3}=dataObject.data.(variableName).nrec;
	end
end

function ok=compare_cell_arrays(c1,c2)
ok=false;
if ~iscell(c1),			return; end
if ~iscell(c2),			return; end
if size(c1)~=size(c2),	return;end
ok=any(any(cellfun(@(x,y) isequal(x,y),c1,c2)));

