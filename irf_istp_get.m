function data=irf_istp_get(data_path,start_time,dt_int,sc_id,quantity,lev)
%IRF_ISTP_GET read ISTP data from CDF files
%
% data = irf_istp_get(data_path,start_time,dt,sc_id,quantity,[lev])
%
% Input:
%	data_path - directories containing ISTP data 
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	sc_id - SC name:
%		'ace' : ACE
%	quantity - one of the following:
%		'b', 'bgsm' : B GSM
%		'bgse' : B GSE
%		'n'    : density
%		'v'    : velocity
%	lev - data level:
%		'h2'   : final data [DEFAULT]
%       'k1', 'k2' : preliminary data
%
% See also C_GET, TOEPOCH
%
% $Id$

% Copyright 2006 Yuri Khotyaintsev

error(nargchk(5,6,nargin))

switch (sc_id)
case 'ace'
	sc_pref = 'ac';
otherwise
	error('caa:noSuchQuantity','SC_ID ''%s'' is not recongized',sc_id)
end

if nargin>5
	if ~strcmp(lev,'h2') || ~strcmp(lev,'k1') || ~strcmp(lev,'k2')
		error('caa:noSuchQuantity','LEV ''%s'' is not recongized',lev)
	end
else lev = 'h2';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define request parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (quantity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSDS PP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'b', 'bgsm'}
	inst  = 'mfi';
	var_s = 'BGSM';
case 'bgse'
	inst  = 'mfi';
	var_s = 'BGSEc';
case 'n'
	inst  = 'swe';
	var_s = 'Np';
case 'v'
	inst  = 'swe';
	var_s = 'Vp';
otherwise
	error('caa:noSuchQuantity','Quantity ''%s'' is not recongized',quantity)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = [];
st =[]; dt =[];
start_time_a = fromepoch(start_time);

if start_time+dt_int > toepoch([start_time_a(1:3) 0 0 0]) +86400
	%interval spans multiple days
	st = start_time; dt = toepoch([start_time_a(1:3) 0 0 0]) +86400 - st;
	while 1
		st(end+1) = st(end) + dt(end);
		if st(end) +86400 > start_time+dt_int
			dt(end+1) = start_time+dt_int -st(end);
			break
		else dt(end+1) = 86400;
		end
	end
else st = start_time; dt = dt_int;
end

for k=1:length(st)
	st_a = fromepoch(st(k));
	
	data_subdir = [inst '_' lev '/' num2str(st_a(1))];
	file_mask = [sc_pref '_' lev '_' inst '_' ...
			strrep(datestr(st_a,29),'-','') ];
	
	data_tmp = [];
	data_tmp = irf_cdf_read([data_path '/' sc_id '/' data_subdir '/' file_mask '*'], ...
			var_s,'latest');
	
	if data_tmp
		if dt(k)~=86400, data_tmp = irf_tlim(data_tmp,st(k) + [0 dt(k)]); end
		if data, data = [data; data_tmp];
		else data = data_tmp;
		end
	end
end
if isempty(data), irf_log('dsrc','No data'), end
