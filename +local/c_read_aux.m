function out=c_read_aux(varargin)
% LOCAL.C_READ_AUX read local cluster aux information
%	[out]=LOCAL.C_READ_AUX(variable,tint) read variable for given time interval
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
%	Example:
%		tint='2005-01-01T05:00:00.000Z/2005-01-05T05:10:00.000Z';
%		R1=LOCAL.C_READ_AUX('R1',tint);
%		DipoleTilt=LOCAL.C_READ_AUX('dipole_tilt__CL_SP_AUX',tint);
%
% to update file index run LOCA.C_UPDATE
%
persistent aux % to make fast access read only once

cd('~/data/CAA');
if isempty(aux)
	load('caa','aux');
end
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

istart=find(aux.tend>tint(1),1);
iend=find(aux.tstart<tint(2),1,'last');

if isempty(istart) || isempty(iend),
	out=[];
	return
end

switch lower(varName)
	case {'r'}
		varToRead={'sc_r_xyz_gse__CL_SP_AUX'};
		readdata;
		out=[data{1} double(data{2})];
	case {'r1','r2','r3','r4'}
		varToRead={'sc_r_xyz_gse__CL_SP_AUX',['sc_dr' varName(2) '_xyz_gse__CL_SP_AUX']};
		readdata;
		out=[data{1} double(data{2}+data{3})];
	case {'dr1','dr2','dr3','dr4'}
		varToRead={['sc_dr' varName(3) '_xyz_gse__CL_SP_AUX']};
		readdata;
		out=[data{1} double(data{2})];
	otherwise
		irf_log('fcal',['Reading variable (assume to exist): ' varName]);
		varToRead={varName};
		readdata;
		out=[data{1} double(data{2})];		
end

	function readdata
		for iFile=istart:iend
			cdf_file=aux.filename(iFile,:);
			[tmpdata,~] = cdfread(cdf_file,'ConvertEpochToDatenum',true,'CombineRecords',true,'Variables', {'time_tags__CL_SP_AUX',varToRead{:}});
			if iFile==istart, data=cell(size(tmpdata));end
			iist=1;iien=numel(tmpdata{1});
			tmpdata{1}=irf_time(tmpdata{1},'date2epoch');
			if iFile==istart
				iist=find(tmpdata{1}>tint(1),1);
			end
			if iFile==iend
				iien=find(tmpdata{1}<tint(2),1,'last');
			end
			for j=1:numel(data),
				data{j}=vertcat(data{j},tmpdata{j}(iist:iien,:));
			end
		end
	end

end




% $Id$
% $Revision$  $Date$

