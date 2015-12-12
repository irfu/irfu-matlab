function [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = rotate_tensor_fac(varargin)
%
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor_fac(PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ,Bback)
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor_fac(Peall,Bback)
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
%       Peall - TSERIES of all tensor terms with column order PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ
% Output: 
%       PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp - Pressure or temperature terms in
%       field-aligned coordinates
% 

if nargin == 2;
    Peall = varargin{1};
    times = Peall.time;
    Bback = varargin{2};
    Ptensor = zeros(length(Peall.time),3,3);
    Ptensor(:,1,1) = Peall.data(:,1);
    Ptensor(:,2,1) = Peall.data(:,2);
    Ptensor(:,3,1) = Peall.data(:,3);
    Ptensor(:,1,2) = Peall.data(:,2);
    Ptensor(:,2,2) = Peall.data(:,4);
    Ptensor(:,3,2) = Peall.data(:,5);
    Ptensor(:,1,3) = Peall.data(:,3);
    Ptensor(:,2,3) = Peall.data(:,5);
    Ptensor(:,3,3) = Peall.data(:,6);
elseif nargin == 7;
    times = varargin{1}.time;
    Ptensor = zeros(length(varargin{1}.time),3,3);
    Ptensor(:,1,1) = varargin{1}.data;
    Ptensor(:,2,1) = varargin{2}.data;
    Ptensor(:,3,1) = varargin{3}.data;
    Ptensor(:,1,2) = varargin{2}.data;
    Ptensor(:,2,2) = varargin{4}.data;
    Ptensor(:,3,2) = varargin{5}.data;
    Ptensor(:,1,3) = varargin{3}.data;
    Ptensor(:,2,3) = varargin{5}.data;
    Ptensor(:,3,3) = varargin{6}.data;
    Bback = varargin{7};
else
    help mms.rotate_tensor_fac;
    return;
end

% Make transformation matrix
Bback = Bback.resample(times);
Bmag = Bback.abs.data;
Rz = Bback.data./([Bmag Bmag Bmag]);
Rx = [1 0 0];
Ry = irf_cross(Rz,Rx);
Rmag = irf_abs(Ry,1);
Ry = Ry./[Rmag Rmag Rmag];
Rx = irf_cross(Ry, Rz);
Rmag = irf_abs(Rx,1);
Rx = Rx./[Rmag Rmag Rmag];
Rotmat = zeros(length(times),3,3);
Rotmat(:,1,:) = Rx;
Rotmat(:,2,:) = Ry;
Rotmat(:,3,:) = Rz;


Ptensorp = zeros(length(times),3,3);
for ii = 1:length(times);
    rottemp = squeeze(Rotmat(ii,:,:));
    Ptensorp(ii,:,:) = rottemp*(squeeze(Ptensor(ii,:,:))*transpose(rottemp));
end

PeXXp = irf.ts_scalar(times,Ptensorp(:,1,1));
PeXYp = irf.ts_scalar(times,Ptensorp(:,1,2));
PeXZp = irf.ts_scalar(times,Ptensorp(:,1,3));
PeYYp = irf.ts_scalar(times,Ptensorp(:,2,2));
PeYZp = irf.ts_scalar(times,Ptensorp(:,2,3));
PeZZp = irf.ts_scalar(times,Ptensorp(:,3,3));

end


