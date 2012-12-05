function dobj = dataobj(varargin)
%DATAOBJ  constructor function for DATAOBJ object
%
% DATAOBJ(FILENAME)
%    Construct dataobj form file FILENAME. FILENAME can also contain
%    wildcards ('*').
%
% DATAOBJ(FILENAME,'tint',tint)
%       tint - limit dataobject to time interval (good for large files)
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
persistent usingNasaPatchCdf

if isempty(usingNasaPatchCdf), % check only once if using NASA cdf
	usingNasaPatchCdf=check_if_using_nasa_patch_cdf;   % stop with error if using NASA cdf
end
shouldReadAllData= 1; % default read all data
noDataReturned   = 0; % default expects data to be returned
if nargin==0, action='create_default_object'; end
if nargin==1, action='read_data_from_file'; end
if nargin==3 && ...
		ischar(varargin{2}) && strcmp(varargin{2},'tint') && ...
		isnumeric(varargin{3}) && (length(varargin{3})==2),
	tint=varargin{3};
	action='read_data_from_file';
end
switch action
	case 'create_default_object'
		% if no input arguments, create a default object
		dobj.FileModDate = datestr(now);
		dobj.VariableAttributes = {};
		dobj.GlobalAttributes = {};
		dobj.Variables = {};
		dobj.vars = {};
		dobj = class(dobj,'dataobj');
	case 'read_data_from_file'
		% if single argument of class ClusterDB, return it
		if (isa(varargin{1},'dataobj'))
			dobj = varargin{1};
			
		elseif ischar(varargin{1})
			
			if strfind(varargin{1},'*')
				cdf_files = dir(varargin{1});
				cdf_files([cdf_files(:).isdir])=[]; % remove directories from the list
				if strfind(varargin{1},filesep) % if there is directory in file name
					filesep_indexes=strfind(varargin{1},filesep);
					directory_name=varargin{1}(1:filesep_indexes(end));
				else
					directory_name='';
				end
				switch numel(cdf_files)
					case 0
						error('no cdf files specified')
					case 1
						cdf_file = [directory_name cdf_files.name];
					otherwise
						for j=1:numel(cdf_files),
							disp([num2str(j) '. ' cdf_files(j).name]);
						end
						j = irf_ask('Choose cdf-file? [%]>','cdf_file',1);
						cdf_file = [directory_name cdf_files(j).name];
				end
				clear cdf_files
			else cdf_file = varargin{1};
			end
			
			if ~exist(cdf_file,'file')
				error(['file ' cdf_file ' does not exist'])
			end
			
			irf_log('dsrc',['Reading: ' cdf_file]);
			%% read in file
			if usingNasaPatchCdf
				[data,info] = cdfread(cdf_file,'CombineRecords',true);
			else
				cdfid   = cdflib.open(cdf_file);
				
				% check if epoch 16
				usingCdfepoch16=strcmpi('cdf_epoch16',...
					getfield(cdflib.inquireVar(cdfid,0),'datatype'));
				
				% read in file
				if usingCdfepoch16,
					irf_log('dsrc',['EPOCH16 time in cdf file:' cdf_file]);
					shouldReadAllData=1; % read all data
					info = cdflib.inquire(cdfid);
					vars=cell(info.numVars-1,1);
					vars_i16 = ones(size(1:info.numVars)); % array indicating which of the variables are EPOCH16
					for jj=1:info.numVars
						vars{jj}=cdflib.getVarName(cdfid,jj-1);
						inq=cdflib.inquireVar(cdfid,jj-1);
						if strcmpi(inq.datatype,'cdf_epoch16')
							vars_i16(jj)=0;
						end
					end
					data=cell(1,info.numVars);
					data(vars_i16==1) = cdfread(cdf_file,'variables',vars(vars_i16==1),'CombineRecords',true);
					info=cdfinfo(cdf_file);
					ii = find(vars_i16==0);
					numrecs = cdflib.getVarAllocRecords(cdfid,0);
					for i=1:length(ii)
						% get time axis
						tc=zeros(2,numrecs);
						for jj=1:numrecs,
							tc(:,jj) = cdflib.getVarRecordData(cdfid,ii(i)-1,jj-1);
						end
						data(ii(i))={tc'};
					end
				else
					[data,info] = cdfread(cdf_file,...
						'ConvertEpochToDatenum',true,...
						'CombineRecords',true);
				end
				cdflib.close(cdfid);
				if ~shouldReadAllData, % check which records to return later
					info=cdfinfo(cdf_file);
					timevar=info.Variables{strcmpi(info.Variables(:,4),'epoch')==1,1};
					timeline = irf_time(cdfread(cdf_file,'Variable',{timevar},'ConvertEpochToDatenum',true,'CombineRecords',true),'date2epoch');
					records_within_interval=find((timeline > tint(1)) & (timeline < tint(2)));
				end
			end
			%% construct data object
			dobj.FileModDate		= info.FileModDate;
			dobj.VariableAttributes = info.VariableAttributes;
			dobj.GlobalAttributes	= info.GlobalAttributes;
			dobj.Variables			= info.Variables;
			% test if there are some data
			if ~(any(strcmpi(info.Variables(:,4),'epoch')==1) || ...
                    any(strcmpi(info.Variables(:,4),'epoch16')==1)),
				nVariables=0; % no time variable, return nothing
				irf_log('dsrc','CDF FILE IS EMPTY!')
			else
				nVariables = size(info.Variables,1);
			end
			dobj.vars = cell(nVariables,2);
			if nVariables>0
				dobj.vars(:,1) = info.Variables(:,1);
				dobj.vars(:,2) = info.Variables(:,1);
				for v=1:nVariables
					% Replace minuses with underscores
					dobj.vars{v,1}(strfind(dobj.vars{v,1},'-')) = '_';
					% Remove training dots
					while (dobj.vars{v,1}(end) == '.')
						dobj.vars{v,1}(end) = [];
					end
					% Take care of '...'
					d3 = strfind(dobj.vars{v,1},'...');
					if d3, dobj.vars{v,1}( d3 + (1:2) ) = []; end
					% Replace dots with underscores
					dobj.vars{v,1}(strfind(dobj.vars{v,1},'.')) = '_';
					% Add "x" if the varible name starts with a number
					if ~isletter(dobj.vars{v,1}(1)),
						dobj.vars{v,1}=['x' dobj.vars{v,1}];
					end
					% Take care of names longer than 63 symbols (Matlab limit)
					if length(dobj.vars{v,1})>63
						dobj.vars{v,1} = dobj.vars{v,1}(1:63);
						disp(['orig var : ' dobj.vars{v,2}])
						disp(['new var  : ' dobj.vars{v,1}])
					end
					if shouldReadAllData, % return all data
						if usingNasaPatchCdf
							dobj.data.(dobj.vars{v,1}).data = data{v};
							dobj.data.(dobj.vars{v,1}).nrec = info.Variables{v,3};
						else
							dobj.data.(dobj.vars{v,1}).data = [data{:,v}];
							dobj.data.(dobj.vars{v,1}).nrec = info.Variables{v,3};
						end
					else
						data_all_records=[data{:,v}];
						if numel(size(data_all_records))==2,
							data_records_within_interval=data_all_records(records_within_interval,:);
						elseif numel(size(data_all_records))==3,
							data_records_within_interval=data_all_records(records_within_interval,:,:);
						elseif numel(size(data_all_records))==4,
							data_records_within_interval=data_all_records(records_within_interval,:,:,:);
						elseif numel(size(data_all_records))==5,
							data_records_within_interval=data_all_records(records_within_interval,:,:,:,:);
						end
						dobj.data.(dobj.vars{v,1}).data = data_records_within_interval;
						dobj.data.(dobj.vars{v,1}).nrec = numel(records_within_interval);
						if numel(records_within_interval)==0,
							irf_log('dsrc','No data within specified time interval');
							noDataReturned=1;
							break;
						end
					end
					dobj.data.(dobj.vars{v,1}).dim = info.Variables{v,2};
					dobj.data.(dobj.vars{v,1}).type = info.Variables{v,4};
					dobj.data.(dobj.vars{v,1}).variance = info.Variables{v,5};
					dobj.data.(dobj.vars{v,1}).sparsity = info.Variables{v,6};
					%Convert to isdat epoch
					if strcmp(dobj.data.(dobj.vars{v,1}).type,'epoch')
						if numel(dobj.data.(dobj.vars{v,1}).data)==1 && dobj.data.(dobj.vars{v,1}).data == 1,
							irf_log('dsrc','CDF FILE IS EMPTY!')
							noDataReturned=1;
						else
							convert_cdfepoch_to_isdat_epoch % update all the structure fields necessary when doing conversion
						end
					elseif strcmp(dobj.data.(dobj.vars{v,1}).type,'epoch16')
						if usingNasaPatchCdf
							convert_cdfepoch16_to_isdat_epoch % update all the structure fields necessary when doing conversion							
						else
							dobj.data.(dobj.vars{v,1}).data = irf_time(dobj.data.(dobj.vars{v,1}).data,'cdfepoch162epoch');
						end
					end
				end
			end
			if noDataReturned
				dobj.data = [];
				irf_log('dsrc','No data returned!')
			end
			dobj = class(dobj,'dataobj');
		else
			error('Wrong argument type')
		end
	otherwise
		error('Wrong number of input arguments')
end

	function convert_cdfepoch_to_isdat_epoch() % nested function
		dobj.data.(dobj.vars{v,1}).data = irf_time(dobj.data.(dobj.vars{v,1}).data,'date2epoch');
		isFieldUnits        = isfield(dobj.VariableAttributes,'UNITS');
		isFieldSIConversion = isfield(dobj.VariableAttributes,'SI_CONVERSION');
		isFieldDeltaPlus    = isfield(dobj.VariableAttributes,'DELTA_PLUS');
		isFieldDeltaMinus   = isfield(dobj.VariableAttributes,'DELTA_MINUS');
		if isFieldUnits,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.UNITS(:,1))==1);
			if iattr, dobj.VariableAttributes.UNITS(iattr,2)={'s'};end % change from ms to s UNITS of epoch if present
		end
		if isFieldSIConversion,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.SI_CONVERSION(:,1))==1);
			if iattr, dobj.VariableAttributes.SI_CONVERSION(iattr,2)={'1.0>s'};end % change from ms to s SI_CONVERSION of epoch if present
		end
		if isFieldDeltaPlus,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.DELTA_PLUS(:,1))==1); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_PLUS{iattr,2}),
				dobj.VariableAttributes.DELTA_PLUS{iattr,2}=dobj.VariableAttributes.DELTA_PLUS{iattr,2}/1000;
			end
		end
		if isFieldDeltaMinus,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.DELTA_MINUS(:,1))==1); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_MINUS{iattr,2}),
				dobj.VariableAttributes.DELTA_MINUS{iattr,2}=dobj.VariableAttributes.DELTA_MINUS{iattr,2}/1000;
			end
		end
	end
	function convert_cdfepoch16_to_isdat_epoch() % nested function
		epochData=cell2mat(dobj.data.(dobj.vars{v,1}).data);
		DIn=epochData(:,1:20);
		Dout=zeros(size(DIn,1),6);
		Dx   = double(DIn - '0');        % For faster conversion of numbers
		Dout(:,1)=Dx(:,8) * 1000 + Dx(:,9) * 100 + Dx(:,10) * 10 + Dx(:,11); % Year
		monthN=sum(Dx(:,4:6),2)-143; % Apr,Sep ok
		monthN(monthN==-3)=1; % jan
		monthN(monthN==-18)=2; % feb
		monthN(monthN==1)=3; % mar
		monthN(monthN==8)=5; % may
		monthN(monthN==14)=6; % jun
		monthN(monthN==12)=7; % jul
		monthN(monthN==-2)=8; % aug
		monthN(monthN==7)=10; % oct
		monthN(monthN==20)=11; % nov
		monthN(monthN==-19)=12; % dec
		Dout(:,2)=monthN;
		Dout(:,3) = Dx(:,1)  * 10 + Dx(:,2);   % Day
		Dout(:,4) = Dx(:,13) * 10 + Dx(:,14);  % Hour
		Dout(:,5) = Dx(:,16) * 10 + Dx(:,17);  % Minute
		Dout(:,6) = Dx(:,19) * 10 + Dx(:,20);     % Second
		psStr=epochData(:,[22:24 26:28 30:32 34:36]); % ps
		psVec=double(psStr-'0');
		power=repmat(10.^(-1:-1:-12),size(psStr,1),1);
		ps=sum(psVec.*power,2);
		Dout(:,6)=Dout(:,6)+ps;
		t=irf_time(Dout,'vector2epoch');
		dobj.data.(dobj.vars{v,1}).data = t;
		isFieldUnits        = isfield(dobj.VariableAttributes,'UNITS');
		isFieldSIConversion = isfield(dobj.VariableAttributes,'SI_CONVERSION');
		isFieldDeltaPlus    = isfield(dobj.VariableAttributes,'DELTA_PLUS');
		isFieldDeltaMinus   = isfield(dobj.VariableAttributes,'DELTA_MINUS');
		if isFieldUnits,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.UNITS(:,1))==1);
			if iattr, dobj.VariableAttributes.UNITS(iattr,2)={'s'};end % change from ms to s UNITS of epoch if present
		end
		if isFieldSIConversion,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.SI_CONVERSION(:,1))==1);
			if iattr, dobj.VariableAttributes.SI_CONVERSION(iattr,2)={'1.0>s'};end % change from ms to s SI_CONVERSION of epoch if present
		end
		if isFieldDeltaPlus,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.DELTA_PLUS(:,1))==1); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_PLUS{iattr,2}),
				dobj.VariableAttributes.DELTA_PLUS{iattr,2}=dobj.VariableAttributes.DELTA_PLUS{iattr,2}/1e12;
			end
		end
		if isFieldDeltaMinus,
			iattr=find(strcmpi(dobj.vars{v,2},dobj.VariableAttributes.DELTA_MINUS(:,1))==1); % to convert DELTA_PLUS
			if iattr && isnumeric(dobj.VariableAttributes.DELTA_MINUS{iattr,2}),
				dobj.VariableAttributes.DELTA_MINUS{iattr,2}=dobj.VariableAttributes.DELTA_MINUS{iattr,2}/1e12;
			end
		end
	end
	function ok=check_if_using_nasa_patch_cdf
		fid=fopen(which('cdfread'));
		while 1
			tline = fgetl(fid);
			if ~ischar(tline), break, end
			if strfind(tline,'Mike Liu')
				fprintf('\n\n\n');
				disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
				disp(' You are using NASA cdfread patch which we have bad support!')
				disp(' This may give errors reading in multidimensional data sets!')
				disp(' Also option ''tint'' in routine databoj is disabled.');
				disp(' We suggest you to use the MATLAB cdfread!');
				disp(' To use MATLAB cdfread please remove path to NASA cdfread patch.');
				disp(' You can execute and then continue:');
				a=which('cdfread');
				ai=strfind(a,'/');
				disp(['> rmpath ' a(1:ai(end))]);
				disp('> clear databoj');
				disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
				fprintf('\n\n\n');
				irf_log('fcal','Using NASA cdf is not supported!');
				ok=true;
				break;
			else
				ok=false;
			end
		end
		fclose(fid);
	end

end % end of main functions

