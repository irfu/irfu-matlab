function data=readCSDS(data_path,start_time,dt,cl_id,quantity)
%readCSDS read CSDS data from ISDAT or disk
%
% data = readCSDS(data_path,start_time,dt,cl_id,quantity)
%
% Input:
%	data_path - ISDAT database strings and directories containing CSDS 
%   subdirectory with data separated by '|'. These data sources are	being tried
%	in the same order. A typical example would be: 
%		'ice:10|disco:10|/data/cluster'
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	cl_id - SC#
%	quantity - one of the following:
%		'b' : B FGM CSDS PP
%		'slat' : spin axis latitude
%		'slong' : spin axis longitude
%
% $Revision$  $Date$
%
% see also C_GET, TOEPOCH

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(5,5,nargin))

cl_id_s = num2str(cl_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define request parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (quantity)
case 'b'
	r.file	= ['PP/FGM/C' cl_id_s '/C' cl_id_s '_PP_FGM_'];
	r.var	= ['B_xyz_gse__C' cl_id_s '_PP_FGM'];
	r.pr	= 'CSDS_PP';
	r.mem	= ['C' cl_id_s];
	r.inst  = 'FGM';
	r.sen	= r.var;
case 'edi'
	r.file	= ['PP/EDI/C' cl_id_s '/C' cl_id_s '_PP_EDI_'];
	r.var	= ['E_xyz_gse__C' cl_id_s '_PP_EDI'];
	r.pr	= 'CSDS_PP';
	r.mem	= ['C' cl_id_s];
	r.inst  = 'EDI';
	r.sen	= r.var;
case 'slat'
	r.file	= 'SP/AUX/CL_SP_AUX_';
	r.var	= ['sc_at' cl_id_s '_lat__CL_SP_AUX'];
	r.pr	= 'CSDS_SP';
	r.mem	= 'CL';
	r.inst  = 'AUX';
	r.sen	= r.var;
case 'slong'
	r.file	= 'SP/AUX/CL_SP_AUX_';
	r.var	= ['sc_at' cl_id_s '_long__CL_SP_AUX'];
	r.pr	= 'CSDS_SP';
	r.mem	= 'CL';
	r.inst  = 'AUX';
	r.sen	= r.var;
otherwise
	error('caa:noSuchQuantity','Quantity ''%s'' is not recongized',quantity)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
start_date_s = strrep(datestr(fromepoch(start_time),29),'-','');

p = tokenize(data_path,'|');

for i=1:length(p)

	if ~isempty(regexp(p{i},':\d{1,2}\>'))
		% use ISDAT
		useISDAT = 1;
	else
		if exist(p{i},'dir')
			% read from directory
			useISDAT = 0;
		else
			warning('caa:noSuchDir','Directory %s does not exist',p{i})
			continue
		end
	end


	if useISDAT
		warning('caa:dataSource','Using ISDAT')
		lasterr('')
		try
			dbase = Mat_DbOpen(p{i});
		catch
			warning('ISDAT:dbOpen','Cannot open ISDAT database')
			continue
		end
		if exist('dbase','var')
			[t, dat] = isGetDataLite( dbase, start_time, dt, ...
			r.pr, r.mem, r.inst, r.sen);
			Mat_DbClose(dbase)
			clear dbase

			if ~isempty(dat)
				% if dat has more the one column, we need to transpose it
				sz = size(dat);
				i_s = find(sz~=length(t));
				if sz(i_s)>1, dat = dat'; end

				data = [double(t) double(dat)];
				return
			else
				% warning('caa:noData','No data')
			end
		end

	else
		warning('caa:dataSource','Using FILE')
		disp([p{i} '/CSDS/' r.file start_date_s '*']);
		data = av_read_cdf([p{i} '/CSDS/' r.file start_date_s '*'], r.var);
		if ~isempty(data)
			data = av_t_lim(data,start_time + [0 dt]);
			return
		else 
			% warning('caa:noData','No data')
		end
	end
end

% if we are here means there was no data
warning('caa:noData','No data')
