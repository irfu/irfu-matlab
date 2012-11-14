function out=c_read(varargin)
% LOCAL.C_READ read local cluster aux information
%	[out]=LOCAL.C_READ(variable,tint) read variable for given time interval
%
% Variable can be CAA variable or shortcuts
%  'R1'  - Cluster 1 position
%  'dR1' - Cluster 1 relative position wrt center
%   'R'   - Cluster center posiiton
% status__CL_SP_AUX
% sc_orbit_num__CL_SP_AUX
% sc_r_xyz_gse__CL_SP_AUX
% sc_v_xyz_gse__CL_SP_AUX
% sc_dr1_xyz_gse__CL_SP_AUX
% sc_dr2_xyz_gse__CL_SP_AUX
% sc_dr3_xyz_gse__CL_SP_AUX
% sc_dr4_xyz_gse__CL_SP_AUX
% sc_at1_lat__CL_SP_AUX
% sc_at1_long__CL_SP_AUX
% sc_at2_lat__CL_SP_AUX
% sc_at2_long__CL_SP_AUX
% sc_at3_lat__CL_SP_AUX
% sc_at3_long__CL_SP_AUX
% sc_at4_lat__CL_SP_AUX
% sc_at4_long__CL_SP_AUX
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
% to update file index run LOCA.C_UPDATE
%
% $Id$

persistent index % to make fast access read only once

cd('/data/caa/CAA');
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
				s=load('caa',['index_' dataset]);
				eval(['index.' dataset '=s.index_' dataset ';']);
			end
			eval(['index=index.' dataset ';']);
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
			[tmpdata,~] = cdfread(cdf_file,'ConvertEpochToDatenum',true,'CombineRecords',true,'Variables', [{['time_tags__' dataset]},varToRead{:}]);
			if iFile==istart, data=cell(size(tmpdata));end
			iist=1;iien=numel(tmpdata{1});
			tmpdata{1}=irf_time(tmpdata{1},'date2epoch');
			if iFile==istart
				iist=find(tmpdata{1}>tint(1),1);
			end
			if iFile==iend
				iien=find(tmpdata{1}<tint(2),1,'last');
			end
			%% check for NaNs
			cdfInfo=cdfinfo(cdf_file);
			fillValList=cdfInfo.VariableAttributes.FILLVAL;
			for iVar=1:numel(varToRead),
				ii=strcmp(fillValList,varToRead{iVar});
				fillVal=fillValList{ii,2};
				tmpdata{iVar+1}(tmpdata{iVar+1}==fillVal)=NaN;
			end
			%% attach to result
			for j=1:numel(data),
				data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:));
			end
		end
		status=1;
	end
end


