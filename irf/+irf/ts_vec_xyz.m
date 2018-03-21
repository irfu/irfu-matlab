function TS = ts_vec_xyz(time,data)
%TS_VEC_XYZ  Factory for 3D vector [X,Y,Z] (TSeries)
%
% TsVec3DXYZ = ts_vec_xyz(time,data)
%
% Create TSeries object - 3D vector [X,Y,Z]

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else, epoch = time;
end

TS = TSeries(epoch,data,'vec_xyz');
