function TS = ts_scalar(time,data)
%IRF.TS_SCALAR  Factory for scalar (TSeries)
%
% TsScalar = irf.ts_scalar(time,data)
%
% Create TSeries object - scalar

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else epoch = time;
end

TS = TSeries(epoch,data,'TensorOrder',0);
