function [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = rotate_tensor_fac(PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ,Bback)
%
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor_fac(PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ,Bback)
%
% Function to coordinate transform the pressure tensor term into
% field-aligned coordinates. PeZZp is aligned with background magnetic
% field, PeXXp is closest to the x-direction and PeYYp completes the
% system.
% Written by D. B. Graham.
%
% Input: (All data must be in TSeries format)
%       PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ - Pressure terms or temperature   
%       Bback - Background magnetic field
% Output: 
%       PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp - Pressure or temperature terms in
%       field-aligned coordinates
% 

Ptensor = zeros(length(PeXX.time),3,3);
Ptensor(:,1,1) = PeXX.data;
Ptensor(:,2,1) = PeXY.data;
Ptensor(:,3,1) = PeXZ.data;
Ptensor(:,1,2) = PeXY.data;
Ptensor(:,2,2) = PeYY.data;
Ptensor(:,3,2) = PeYZ.data;
Ptensor(:,1,3) = PeXZ.data;
Ptensor(:,2,3) = PeYZ.data;
Ptensor(:,3,3) = PeZZ.data;

% Make transformation matrix
Bback = Bback.resample(PeXX);
Bmag = Bback.abs.data;
Rz = Bback.data./([Bmag Bmag Bmag]);
Rx = [1 0 0];
Ry = irf_cross(Rz,Rx);
Rmag = irf_abs(Ry,1);
Ry = Ry./[Rmag Rmag Rmag];
Rx = irf_cross(Ry, Rz);
Rmag = irf_abs(Rx,1);
Rx = Rx./[Rmag Rmag Rmag];
Rotmat = zeros(length(PeXX.time),3,3);
Rotmat(:,1,:) = Rx;
Rotmat(:,2,:) = Ry;
Rotmat(:,3,:) = Rz;


Ptensorp = zeros(length(PeXX.time),3,3);
for ii = 1:length(PeXX.time);
    rottemp = squeeze(Rotmat(ii,:,:));
    Ptensorp(ii,:,:) = rottemp*(squeeze(Ptensor(ii,:,:))*transpose(rottemp));
end

PeXXp = TSeries(PeXX.time,Ptensorp(:,1,1));
PeXYp = TSeries(PeXX.time,Ptensorp(:,1,2));
PeXZp = TSeries(PeXX.time,Ptensorp(:,1,3));
PeYYp = TSeries(PeXX.time,Ptensorp(:,2,2));
PeYZp = TSeries(PeXX.time,Ptensorp(:,2,3));
PeZZp = TSeries(PeXX.time,Ptensorp(:,3,3));


