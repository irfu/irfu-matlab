function [datasetName,varName]=caa_get_dataset_name(fullVarName,flag)
% CAA_GET_DATASET_NAME construct dataset name form variable name
%
% [datasetName,varName]=CAA_GET_DATASET_NAME(fullVarName) where fullVarName
% is string of syntax 'varName__datasetName'. fullVarName can also be a
% cell array of strings (TODO).
%
% CAA_GET_DATASET_NAME(fullVarName,'_') substitutes all minus signs '-' in
% the output datasetName with underscores '_'
%
% Example:
%	[d,v] = caa_get_dataset_name('pressure__C3_CP_CIS-CODIF_HS_H1_MOMENTS','_');

%% Check input
if nargin == 0, help caa_get_dataset_name; return; end
if nargin == 1, flag = '';end
assert(ischar(fullVarName) | iscell(fullVarName),'fullVarName should be string or cell array of strings');
if nargin == 2 && strcmpi(flag,'_')
	doSubstituteUnderscoreToMinus = true;
else
	doSubstituteUnderscoreToMinus = false;
end

%% Obtain dataset names
if ischar(fullVarName)
	[datasetName,varName] = dataset_name(fullVarName);
elseif iscell(fullVarName)
	datasetName = cell(size(fullVarName));
	varName = datasetName;
	for iVar = 1:numel(fullVarName)
		[datasetName{iVar},varName{iVar}] ...
			= dataset_name(fullVarName{iVar});
	end
end

%% Output
if nargout == 0
elseif nargout == 1
	varName = [];
elseif nargout == 2
else
	error('caa_get_dataset_name: Too many output parameters');
end

return;

% Nested functions
	function [dName,vName] = dataset_name(fullVar)
		dd=regexp(fullVar, '__', 'split');
		vName = dd{1};
		if length(dd)==2 % data object can be properly identified
			dName=dd{2};
			if strcmp(dName,'C3_CP_PEA_') % the bad case of PEACE
				dName='C3_CP_PEA_MOMENTS';
			elseif strcmp(dName,'C2_CP_PEA_') % the bad case of PEACE
				dName='C2_CP_PEA_MOMENTS';
			elseif strcmp(dName,'C1_CP_PEA_') % the bad case of PEACE
				dName='C1_CP_PEA_MOMENTS';
			elseif strcmp(dName,'C4_CP_PEA_') % the bad case of PEACE
				dName='C4_CP_PEA_MOMENTS';
			end
		elseif length(dd)==3 % the case of PEACE moments
			if strcmp(dd{3},'MOMENTS')
				dName=[dd{2}(1:2) '_CP_PEA_' dd{3}];
			end
		end
		if doSubstituteUnderscoreToMinus
			dName(strfind(dName,'-'))='_'; % substitute '-' with '_'
		end
	end
end
