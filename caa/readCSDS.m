function data=readCSDS(data_path,start_time,dt,cl_id,quantity)
%readCSDS read CSDS data from ISDAT or disk
%
% data = readCSDS(data_path,start_time,dt,cl_id,quantity)
%
% Input:
%	data_path - eather directory containing CSDS subdirectory with data, or 
%	ISDAT database string host:NN
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	cl_id - SC#
%	quantity - one of the following:
%		'b' : B FGM CSDS PP
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
	r.sen	= ['B_xyz_gse__C' cl_id_s '_PP_FGM'];
otherwise
	error('caa:noSuchQuantity','Quantity ''%s'' is not recongized',quantity)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];

if ~isempty(regexp(data_path,':\d{1,2}\>'))
	% use ISDAT
	useISDAT = 1;
	lasterr('')
	try
		dbase = Mat_DbOpen(data_path);
	catch
		error('ISDAT:dbOpen','Cannot open ISDAT database')
	end
else
	if exist(data_path,'dir')
		% read from directory
		useISDAT = 0;
	else
		error('caa:noSuchDir','Directory %s does not exist',data_path)
	end
end

start_date_s = strrep(datestr(fromepoch(start_time),29),'-','');

if useISDAT
	warning('caa:dataSource','Using ISDAT')
	[t, dat] = isGetDataLite( dbase, start_time, dt, ...
	r.pr, r.mem, r.inst, r.sen);
	if ~isempty(dat)
		data = [double(t) double(dat)'];
		clear t dat
	else, warning('caa:noData')
	end

else
	warning('caa:dataSource','Using FILE')
	data = av_read_cdf([data_path '/CSDS/' r.file start_date_s '*'], r.var);
	if ~isempty(data)
		data = av_t_lim(data,start_time + [0 dt]);
	else, warning('caa:noData')
	end
end
