function TS = ts_tensor_xyz(time,data)
%TS_VEC_XYZ  Factory for Tensor [XX,XY,XZ; YX,YY,YZ; ZX,ZY,ZZ] (TSeries)
%
% TsTensorXYZ = ts_tensor_xyz(time,data)
%
% Create TSeries object - 3D vector [X,Y,Z]

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else, epoch = time;
end

TS = TSeries(epoch,data,'tensor_xyz');
