function TS = ts_scalar(time,data)
%IRF.TS_SCALAR  Factory for scalar (TSeries)
%
% TsScalar = irf.ts_scalar(time,data)
%
% Create TSeries object - scalar
% time  - Epoch object or TT2000 epoch, see EpochTT. 
% data  - vector of length equal to time or it can be a NxM matrix with the 
%         number of rows N being equal to the length of time.

if ~isa(time,'GenericTimeArray')
	epoch = EpochTT(time);
else
	epoch = time;
end

if size(data,1) == epoch.length
	TS = TSeries(epoch,data,'TensorOrder',0);
elseif (ndims(data)==2) && (size(data,2) == epoch.length)
	data = data';
	TS = TSeries(epoch,data,'TensorOrder',0);
end
