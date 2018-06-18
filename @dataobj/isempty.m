function isemp = isempty(dobj)
%ISEMPTY True for empty dataobject
%
% Empty dataobject means that there is no data, there still can be meta-data
% information.
%
%    ISEMPTY(X) returns 1 if X is an empty data object and 0 otherwise. An

isemp = false; % default

% check nrec of time variable
numRecords = dobj.Variables{1,3};
if isempty(numRecords) || (numRecords == 0)
	isemp = true;
	return;
end

numTimeData =  numel(dobj.data.(variable_mat_name(dobj.Variables{1,1})).data);
if isempty(numTimeData) || (numTimeData == 0)
	isemp = true;
	irf_log('dsrc',['WARNING!!!! Number of records is ' num2str(numRecords) ...
		' but there is no data!']);
	disp(dobj);
	return;
end

