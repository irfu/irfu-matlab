function out=c_read(varargin)
% LOCAL.C_READ read local cluster aux information
%	[out]=LOCAL.C_READ(variable,tint) 
%		read variable for given time interval in matlab format (matrix)
%
% Variable can be CAA variable or shortcuts
%  'R1'  - Cluster 1 position
%  'dR1' - Cluster 1 relative position wrt center
%   'R'   - Cluster center posiiton
%
% Possible variables: (question mark needs to be sustituted to Cluster number)
% status__CL_SP_AUX
% sc_orbit_num__CL_SP_AUX
% sc_r_xyz_gse__CL_SP_AUX
% sc_v_xyz_gse__CL_SP_AUX
% sc_dr?_xyz_gse__CL_SP_AUX
% sc_at?_lat__CL_SP_AUX
% sc_at?_long__CL_SP_AUX
% sc_config_QG__CL_SP_AUX
% sc_config_QR__CL_SP_AUX
% sc_dr_min__CL_SP_AUX
% sc_dr_max__CL_SP_AUX
% gse_gsm__CL_SP_AUX
% dipole_tilt__CL_SP_AUX
% sc_geom_size__CL_SP_AUX
% sc_geom_elong__CL_SP_AUX
% sc_geom_planarity__CL_SP_AUX
% sc_geom_E_dir_gse__CL_SP_AUX
% sc_geom_P_nor_gse__CL_SP_AUX
%
% Invar_Lat__C1_JP_PMP
% Mag_Local_time__C1_JP_PMP
% L_value__C1_JP_PMP
% Pred_B_mag__C1_JP_PMP
% 
%	Example:
%		tint='2005-01-01T05:00:00.000Z/2005-01-05T05:10:00.000Z';
%		R1=local.c_read('R1',tint);
%		DipoleTilt=local.c_read('dipole_tilt__CL_SP_AUX',tint);
%
% to update file index run LOCAL.C_UPDATE
%
% $Id$

persistent index % to make fast access read only once

caaDir='/data/caa/CAA/';
if isempty(index), index=struct('dummy',[]);end
% sc_r_xyz_gse__CL_SP_AUX
% sc_v_xyz_gse__CL_SP_AUX
% sc_dr1_xyz_gse__CL_SP_AUX
% sc_at1_lat__CL_SP_AUX
% sc_at1_long__CL_SP_AUX

if nargin==2,
	varName=varargin{1};
	tint=varargin{2};
	if ischar(tint),
		tint=irf_time(tint,'iso2tint');
	end
else
	irf_log('fcal','Only 2 arguments supported');
	return
end

out=[];
specialCaseCis=0;
switch lower(varName)
	case {'r'}
		varToRead={'sc_r_xyz_gse__CL_SP_AUX'};
		ok=readdata;
		if ok, out=[data{1} double(data{2})]; end
	case {'r1','r2','r3','r4'}
		varToRead={'sc_r_xyz_gse__CL_SP_AUX',['sc_dr' varName(2) '_xyz_gse__CL_SP_AUX']};
		ok=readdata;
		if ok, out=[data{1} double(data{2}+data{3})];end
	case {'dr1','dr2','dr3','dr4'}
		varToRead={['sc_dr' varName(3) '_xyz_gse__CL_SP_AUX']};
		ok=readdata;
		if ok, out=[data{1} double(data{2})];end
	otherwise
		irf_log('fcal',['Reading variable (assume to exist): ' varName]);
		if strfind(varName,'CIS'),specialCaseCis=1;end
		varToRead={varName};
		ok=readdata;
		if ok, out=[data{1} double(data{2})];end
end

	function status=readdata
		%% find index
		ii=strfind(varToRead{1},'__');
		if ii,
			dataset=varToRead{1}(ii+2:end);
			if ~isfield(index,dataset)
				s=load([caaDir 'caa'],['index_' dataset]);
				index.(dataset)=s.(['index_' dataset]);
			end
			index=index.(dataset);
		else
			irf_log('dsrc',['Do not now how to read variable: ' varToRead{1}]);
			status=0;
			return
		end
		%% find records within time interval
		istart=find(index.tend>tint(1),1);
		iend=find(index.tstart<tint(2),1,'last');
		
		if isempty(istart) || isempty(iend),
			status=0;
			return
		end
		%% read in records
		for iFile=istart:iend
			cdf_file=index.filename(iFile,:);
			if specialCaseCis,
				dataset=strrep(dataset,'CIS_','CIS-');
				varToRead=strrep(varToRead,'CIS_','CIS-');
			end
			irf_log('dsrc',['Reading: ' cdf_file]);
			%% check if epoch16
			cdfid=cdflib.open([caaDir cdf_file]);
			useCdfepoch16=strcmpi(getfield(cdflib.inquireVar(cdfid,0),'datatype'),'cdf_epoch16');
			if useCdfepoch16,
				irf_log('dsrc',['EPOCH16 time in cdf file:' cdf_file]);
				tmptime=readCdfepoch16(cdfid,0); % read time which has variable number 0
				tt=irf_time(tmptime','cdfepoch162epoch');
				tmpdata=cell(1,numel(varToRead));
				for iVar=1:numel(varToRead),
					tmp=readCdfepoch16(cdfid,varToRead{iVar}); % currently only first variable read
					tmpdata{iVar}=[tt tmp'];
				end
			else
				[tmpdata,~] = cdfread([caaDir cdf_file],'ConvertEpochToDatenum',true,'CombineRecords',true,'Variables', [{['time_tags__' dataset]},varToRead{:}]);
				tmpdata{1}=irf_time(tmpdata{1},'date2epoch');
			end
			if iFile==istart, data=cell(size(tmpdata));end
			iist=1;iien=numel(tmpdata{1});
			if iFile==istart
				iist=find(tmpdata{1}>tint(1),1);
			end
			if iFile==iend
				iien=find(tmpdata{1}<tint(2),1,'last');
			end
			%% check for NaNs
			for iVar=1:numel(varToRead),
				fillVal=value_of_variable_attribute(cdfid,varToRead{iVar},'FILLVAL');
				tmpdata{iVar+1}(tmpdata{iVar+1}==fillVal)=NaN;
			end
			%% attach to result
			for j=1:numel(data),
				data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:));
			end
			cdflib.close(cdfid);
		end
		status=1;
	end
end
function data = readCdfepoch16(cdfid,varName)
if isnumeric(varName),
	varnum=varName;
elseif ischar(varName)
	varnum  = cdflib.getVarNum(cdfid,varName);
else 
	error('varName should be variable number or name');
end
numrecs = cdflib.getVarNumRecsWritten(cdfid,varnum);
numElements = getfield(cdflib.inquireVar(cdfid,varnum),'numElements');
data=zeros(2+numElements,numrecs);

for j = 0:numrecs-1
	% This reads in the data in raw epoch 16 format
	% That implies each Epoch16 value is a 2-element double precision
	% value in MATLAB
	data(:,1+j) = cdflib.getVarRecordData(cdfid,varnum,j);
end
data=data';
end
function value=value_of_variable_attribute(cdfid,varName,attrName)
attrnum = cdflib.getAttrNum(cdfid,attrName);
varnum = cdflib.getVarNum(cdfid,varName);
value = cdflib.getAttrEntry(cdfid,attrnum,varnum);
end
