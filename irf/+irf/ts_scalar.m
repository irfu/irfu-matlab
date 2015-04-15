function TS = ts_scalar(time,data)
%TS_SCALAR  Factory for scalar (TSeries)
%
% TsScalar = irf.ts_scalar(time,data)
%
% Create TSeries object - scalar

if ~isa(time,'GenericTimeArray'), epoch = EpochTT2000(time);
else epoch = time;
end

TS = TSeries(epoch,data,'TensorOrder',0);
