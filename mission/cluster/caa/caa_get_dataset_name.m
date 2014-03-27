function [datasetName,varName]=caa_get_dataset_name(fullVarName,flag)
% CAA_GET_DATASET_NAME construct dataset name form variable name
%
% [datasetName,varName]=CAA_GET_DATASET_NAME(fullVarName) where fullVarName
% is of syntax 'varName__datasetName'. 
%
% CAA_GET_DATASET_NAME(fullVarName,'_') substitutes all minus signs '-' in
% the output datasetName with underscores '_'
%
% Example:
%	[d,v] = caa_get_dataset_name('pressure__C3_CP_CIS-CODIF_HS_H1_MOMENTS','_');

if nargin == 0; help caa_get_dataset_name;return; end

assert(ischar(fullVarName),'fullVarName should be string');

dd=regexp(fullVarName, '__', 'split');
varName = dd{1};
if length(dd)==2, % data object can be properly identified
	datasetName=dd{2};
	if strcmp(datasetName,'C3_CP_PEA_'), % the bad case of PEACE
		datasetName='C3_CP_PEA_MOMENTS';
	elseif strcmp(datasetName,'C2_CP_PEA_'), % the bad case of PEACE
		datasetName='C2_CP_PEA_MOMENTS';
	elseif strcmp(datasetName,'C1_CP_PEA_'), % the bad case of PEACE
		datasetName='C1_CP_PEA_MOMENTS';
	elseif strcmp(datasetName,'C4_CP_PEA_'), % the bad case of PEACE
		datasetName='C4_CP_PEA_MOMENTS';
	end
elseif length(dd)==3, % the case of PEACE moments
	if strcmp(dd{3},'MOMENTS'),
		datasetName=[dd{2}(1:2) '_CP_PEA_' dd{3}];
	end
end

if nargin == 2 && strcmpi(flag,'_')
	datasetName(strfind(datasetName,'-'))='_'; % substitute '-' with '_'
end
