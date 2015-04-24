function res = get_ts(dobj,var_s)
%GETMAT(dobj, var_s)  get a variable in the matlab format
%
% Output:
%			empty		if variable does not exist
%			matrix		if variable depends only on time (e.g. fields, position,..)
%						first column time, other columns variable values
%			structure 	if variable has additional dependencies (e.g. spectra)
%						res.t - time 
%						res.dt - time step (can be also defined res.dt.plus, res.dt.minus)
%						res.data - data
%						res.dep_x - additional dependencies

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

data = get_variable(dobj,var_s);
if isempty(data), % no such variable, return empty
	res=[];
	return;
end
fillv = getfillval(dobj,var_s);
if ~ischar(fillv), data.data(data.data==fillv) = NaN;
else irf.log('warning','fill value is character: discarding')
end

if strcmpi(data.DEPEND_0.type,'tt2000')
  Time = EpochTT2000(data.DEPEND_0.data);
else
  Time = EpochUnix(data.DEPEND_0.data);
end

tensorOrder = length(data.variance(3:end));

switch tensorOrder
  case 0 % scalar
  case 1 % vector
  case 2 % tensor
  otherwise
    error('TensorOrder>2 not supported')
end

res = TSeries(Time,data.data,'TensorOrder',tensorOrder);
res.name = data.name;
ud = data; ud = rmfield(ud,'DEPEND_0'); ud = rmfield(ud,'data');
ud = rmfield(ud,'nrec'); ud = rmfield(ud,'dim'); ud = rmfield(ud,'name');
ud = rmfield(ud,'variance');
res.userData = ud;
