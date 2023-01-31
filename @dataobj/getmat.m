function res = getmat(dobj,var_s)
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
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

data = getv(dobj,var_s);
if isempty(data) % no such variable, return empty
  res=[];
  return;
end
fillv = getfillval(dobj,var_s);
if ~ischar(fillv)
  data.data(data.data==fillv) = NaN;
else
  irf_log('dsrc','fill value is character: discarding')
end

dep = getdep(dobj,var_s);
dim = length(data.variance(3:end));

if dim <=1
  plot_data = double(data.data)';
  if isfield(dep,'DEPEND_O') % add time column as first column
    if isfield(dep.DEPEND_O,'type') && strcmpi(dep.DEPEND_O.type,'tt2000')
      timeLine = EpochTT(dep.DEPEND_O.data).epochUnix;
    else, timeLine = dep.DEPEND_O.data;
    end
    res = [timeLine plot_data'];
  else                        % time variable requested, return column vector
    res = plot_data(:);
  end
else
  dep_x = cell(size(dep.DEPEND_X,1),1);
  for d = 1:size(dep.DEPEND_X,1)
    dep_x{d} = getv(dobj,dep.DEPEND_X{d,1});
    dep_x{d}.s = dep.DEPEND_X{d,1};
    dep_x{d}.fillv = getfillval(dobj,dep_x{d}.s);
    if isnumeric(dep_x{d}.fillv) % only implemented for numeric data
      dep_x{d}.data(dep_x{d}.data==dep_x{d}.fillv) = NaN;
    end
    dep_x{d}.units = getunits(dobj,dep_x{d}.s);
    dep_x{d}.lab = getlablaxis(dobj,dep_x{d}.s);
    % check if DELTA_PLUS and  DELTA_MINUS are given
    if isfield(dep_x{d},'DELTA_PLUS') && isfield(dep_x{d},'DELTA_MINUS')
      dep_x{d}.df=struct('plus',dep_x{d}.DELTA_PLUS,'minus',dep_x{d}.DELTA_MINUS);
      if ischar(dep_x{d}.DELTA_PLUS)
        deltaplus  = getv(dobj,dep_x{d}.DELTA_PLUS);
        deltaminus = getv(dobj,dep_x{d}.DELTA_MINUS);
        dep_x{d}.df.plus = double(deltaplus.data);
        dep_x{d}.df.minus= double(deltaminus.data);
      end
    else
      dep_x{d}.df=[];
    end
  end
  % Obtain time DELTA_PLUS and  DELTA_MINUS if given
  % Also do necessary tome conversion if needed
  if strcmpi(dep.DEPEND_O.type,'tt2000')
    timeLine = EpochTT(dep.DEPEND_O.data).epochUnix; factor = 1e9;
  else, timeLine = dep.DEPEND_O.data; factor = 1;
  end
  timevar=getv(dobj,dobj.VariableAttributes.DEPEND_0{1,2});
  if isfield(timevar,'DELTA_PLUS') && isfield(timevar,'DELTA_MINUS')
    dt=struct('plus',timevar.DELTA_PLUS,'minus',timevar.DELTA_MINUS);
    if ischar(timevar.DELTA_PLUS)
      deltaplus = getv(dobj,timevar.DELTA_PLUS);
      dt.plus = double(deltaplus.data(1,:))/factor; % TODO check if "(1,:)" is correct order after latest NASA CDF update, see bug #23
    elseif isnumeric(timevar.DELTA_PLUS)
      dt.plus = double(timevar.DELTA_PLUS)/factor;
    end
    if ischar(timevar.DELTA_MINUS)
      deltaminus = getv(dobj,timevar.DELTA_MINUS);
      dt.minus = double(deltaminus.data(1,:))/factor;
    elseif isnumeric(timevar.DELTA_MINUS)
      dt.minus = double(timevar.DELTA_MINUS)/factor;
    end
  else
    dt=[];
  end
  if isnumeric(data.FILLVAL) % put fillvalues to NaN
    data.data(data.data==data.FILLVAL) = NaN;
  end
  if isempty(dt)
    res = struct('t',timeLine,'dep_x',[],'data',data.data);
  else
    res = struct('t',timeLine,'dep_x',[],'data',data.data,'dt',dt);
  end
  res.dep_x = dep_x;
end
