function out=c_caa_meta(varargin)
% LOCAL.C_CAA_META return meta data structure
%	IRF_TEMPLATE(AX,...) applies function on axis with handle AX
%
%
%	Example:
%		LOCAL.C_CAA_META({'C1' 'PEA'})
%		IRF_TEMPLATE(gca,[1 2 3])
%
% 	See also IRF_FUNCTION_A, IRF_FUNCTION_B.
%
% $Id$

%	LOCAL.C_CAA_META('create') create file with all structures


%- (The line above is needed for CVS, keep it! REMOVE THIS LINE!!)
%-----------------------------------------
%- This is a template function that shows how the help should be formatted
%- When writing the code please follow guidelines in http://www.datatool.com/prod02.htm
%- Remove all lines from this help starting with '-' in your code
%-----------------------------------------

if nargin==0, 
	help local.c_caa_meta;
	return;
end

if nargin==1 && ischar(varargin{1}) && strcmp(varargin{1},'create')
		irf_log('fcal','Creating structures');
	    d = dir('./*.XML');
		isub = [d(:).isdir]; %# returns logical vector
		d(isub)=[];
		nameFiles = {d.name}';
		delme=[];
		save -v7 index delme;
		for iFile=1:numel(nameFiles)
			fileName=nameFiles{iFile};
			if fileName(1)==' ', disp(['Starts with space:' fileName]);continue;end % in case space in beginning of filename
			nameDataset=fileName(1:find(fileName=='.',1)-1);
			nameDataset(strfind(nameDataset,'-')) = '_';
			disp(nameDataset);
			try
				s=xml2struct(fileName);
			catch 
				irf_log('dsrc','Error reading!');
				continue;
			end
			if isfield(s,'DATASETS'),
				s=s.DATASETS;
			end
			if isfield(s,'DATASET_METADATA'),
				s=s.DATASET_METADATA;
			end
			eval(['meta_' nameDataset '=s;']);
%			save(['meta_' nameDataset],'-v7',['meta_' nameDataset]);
			save('index',['meta_' nameDataset],'-append');
		end
else
  load index s;
  metaNames={s(:).name};
  metaNames(1)=[];
  datasetNames=cellfun(@(x) x(6:end),metaNames,'uniformoutput',false);
  datasetNames=datasetNames(:);
  iSelected=true(numel(datasetNames),1);
  for jInp=1:numel(varargin)
	filter=varargin{jInp};
	if ischar(filter)
	  ii=cellfun(@(x) any(strfind(x,filter)),datasetNames);
	  iSelected=iSelected & ii(:);
	end
  end
  if sum(iSelected)>0
	if sum(iSelected)==1 && nargout == 1,
	  load('index',metaNames{iSelected});
	  eval(['out = ' metaNames{iSelected} ';']);
	else
		disp(vertcat(datasetNames(iSelected)));
  	end
  end
end
