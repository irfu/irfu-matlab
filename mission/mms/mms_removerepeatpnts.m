function newdata = mms_removerepeatpnts(tsdata)
%MMS_REMOVEREPEATPNTS  Remove repeated elements in TSERIES or
%structure data. Important when using DEFATT products. Must have a time
%variable.
%
% newtsdata = mms_removerepeatpnts(tsdata);
%
% Written by D. B. Graham


threshold = 100; % Points separated in time by less than 100ns are treated as repeats

if isa(tsdata,'TSeries')
  diffs = diff(tsdata.time.epoch);

  norepeat = ones(length(tsdata.time),1);
  norepeat(diffs < threshold) = 0;

  newtstime = tsdata.time(norepeat == 1);
  newtsdata = tsdata.data(norepeat == 1,:);

  newdata = TSeries(newtstime,newtsdata,'to',tsdata.tensorOrder);

elseif isstruct(tsdata) && isfield(tsdata,'time')
  if isa(tsdata.time,'EpochTT')
    diffs = diff(tsdata.time.epoch);
  else
    diffs = diff(tsdata.time);
  end

  norepeat = ones(length(tsdata.time),1);
  norepeat(diffs < threshold) = 0;

  varnames = fieldnames(tsdata);
  for ii=1:length(varnames)
    tsdata.(varnames{ii}) = tsdata.(varnames{ii})(norepeat == 1,:);
  end
  newdata = tsdata;
else
  newdata = tsdata; % no change to input if it's not a TSERIES or structure
end

end

