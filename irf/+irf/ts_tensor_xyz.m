function TS = ts_tensor_xyz(time,data)
%TS_VEC_XYZ  Factory for Tensor [XX,XY,XZ; YX,YY,YZ; ZX,ZY,ZZ] (TSeries)
%
% TsTensorXYZ = ts_tensor_xyz(time,data)
%
% Create TSeries object - 3D vector [X,Y,Z]

if ~isa(time,'GenericTimeArray'), epoch = EpochTT(time);
else, epoch = time;
end

% if only one entry (3x3) add singleton dimension to make first dim time (1x3x3)
if ndims(data)==2; data_temp(1,:,:)=data; data=data_temp; end

TS = TSeries(epoch,data,'tensor_xyz');
