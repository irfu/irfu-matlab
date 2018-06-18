function TS = ts_vec_xy(time,data)
%TS_VEC_XYZ  Factory for 2D vector [X,Y] (TSeries)
%
% TsVec3DXYZ = ts_vec_xy(time,data)
%
% Create TSeries object - 2D vector [X,Y]

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else, epoch = time;
end

TS = TSeries(epoch,data,'vec_xy');
