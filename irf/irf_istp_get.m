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
%		'ace'  : ACE
%		'omni2' : OMNI2 merged data
%	quantity - one of the following:
%		'b', 'bgsm' : B GSM
%		'bgse' : B GSE
%		'n'    : density
%		'v'    : velocity
%	lev - data level:
%		'h2'   : final data [DEFAULT]
%		'h0'   : combined, definitive data [DEFAULT for OMNI2]
%       'k1', 'k2' : preliminary data
%
% See also C_GET, TOEPOCH
%

% Copyright 2006 Yuri Khotyaintsev

narginchk(5,6)

switch (sc_id)
  case 'ace'
    sc_pref = 'ac';
    % ACE data is split into daily files
    month_step=0;
    day_step=1;
  case 'omni2'
    sc_pref = 'omni2';
    % OMNI2 data is split into 6-month files
    month_step=6;
    day_step=0;
  otherwise
    error('caa:noSuchQuantity','SC_ID ''%s'' is not recongized',sc_id)
end

if nargin>5
  if ~strcmp(lev,'h2') || ~strcmp(lev,'h0') || ~strcmp(lev,'k1') || ~strcmp(lev,'k2')
    error('caa:noSuchQuantity','LEV ''%s'' is not recongized',lev)
  end
else
  switch (sc_id)
    case 'ace'
      lev = 'h2';
    case 'omni2'
      lev = 'h0';
    otherwise
      error('Unknown default level.')
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define request parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (sc_id)
  case 'ace'
    switch (quantity)
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
        error('caa:noSuchQuantity','Quantity ''%s'' is not recongized for spacecraft %s',quantity,sc_id)
    end
  case 'omni2'
    inst  = 'mrg1hr';
    switch (quantity)
      case {'b', 'bgsm'}
        var_s = {'BX_GSE','BY_GSM','BZ_GSM'};
      case 'bgse'
        var_s = {'BX_GSE','BY_GSE','BZ_GSE'};
      case 'n'
        var_s = 'N';
      case 'v'
        var_s = 'V';
      otherwise
        error('caa:noSuchQuantity','Quantity ''%s'' is not recongized for spacecraft %s',quantity,sc_id)
    end
  otherwise
    error('caa:noSuchQuantity','Quantity ''%s'' is not recongized for spacecraft %s',quantity,sc_id)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate intervals and construct filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

st = start_time;
start_time_a = fromepoch(start_time);
if month_step ~= 0
  start_time_a(2)=fix((start_time_a(2)-1)/month_step)*month_step + 1;
  start_time_a(3)=1;
end
file_mask = {[sc_pref '_' lev '_' inst '_' ...
  strrep(datestr(start_time_a,29),'-','') ]};
eof_a=date_increment(start_time_a,month_step,day_step);

if start_time+dt_int > toepoch(eof_a)
  %interval spans multiple files
  dt = toepoch(eof_a) - st;
  while 1
    st(end+1) = st(end) + dt(end);
    st_a = fromepoch(st(end));
    file_mask(end+1) = {[sc_pref '_' lev '_' inst '_' ...
      strrep(datestr(st_a,29),'-','') ]};
    eof_a=date_increment(st_a,month_step,day_step);
    if toepoch(eof_a) > start_time+dt_int
      dt(end+1) = start_time+dt_int -st(end);
      break
    else, dt(end+1) = toepoch(eof_a)-toepoch(st_a);
    end
  end
else
  dt = dt_int;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
for k=1:length(st)
  st_a = fromepoch(st(k));
  
  data_subdir = [inst '_' lev '/' num2str(st_a(1))];
  if iscell(var_s)
    data_tmp = irf_cdf_read([data_path '/' sc_id '/' data_subdir '/' file_mask{k} '*'], ...
      var_s(1),'latest');
    for var_s_indx=2:length(var_s)
      data_tmp_2 = irf_cdf_read([data_path '/' sc_id '/' data_subdir '/' file_mask{k} '*'], ...
        var_s(var_s_indx),'latest');
      if any(data_tmp_2(:,1) ~= data_tmp_2(:,1)),error('Unable to concatenate variables.'),end
      data_tmp=[data_tmp data_tmp_2(:,2)];
    end
  else
    data_tmp = irf_cdf_read([data_path '/' sc_id '/' data_subdir '/' file_mask{k} '*'], ...
      var_s,'latest');
  end
  
  if ~isempty(data_tmp)
    if dt(k) > 0, if k==1 || k==length(st), data_tmp = irf_tlim(data_tmp,st(k) + [0 dt(k)]); end, end
    if ~isempty(data), data = [data; data_tmp];
    else, data = data_tmp;
    end
  end
end
if isempty(data), irf_log('dsrc','No data'), end
end

function date_out=date_increment(date_in,month_step,day_step)
date_out=date_in;
date_out(2)=date_out(2)+month_step;
date_out(3)=date_out(3)+day_step;
if date_out(2) > 12
  date_out(1)=date_out(1)+1;
  date_out(2)=date_out(2)-12;
end
if date_out(3) > eomday(date_out(1),date_out(2))
  date_out(3)=date_out(3)-eomday(date_out(1),date_out(2));
  date_out(2)=date_out(2)+1;
  if date_out(2) > 12
    date_out(1)=date_out(1)+1;
    date_out(2)=date_out(2)-12;
  end
end
end
