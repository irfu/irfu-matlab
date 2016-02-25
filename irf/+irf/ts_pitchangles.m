function TS = ts_skymap(time,data,energy,pitchangles,varargin)
%IRF.TS_SCALAR  Factory for scalar (TSeries)
%
% TsSkymap = irf.ts_skymap(time,data,energy,pitchangles,'optionalArgs1',Val1,...)
%
% Create TSeries object - skymap

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else epoch = time;
end

TS = PDist(epoch,data,energy,pitchangles,varargin{:});