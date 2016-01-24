function [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = rotate_tensor(varargin)
% MMS.ROTATE_TENSOR rotate pressure or temperature tensor into another
% coordinate system
%
% Examples:
% Rotate tensor into field-aligned coordinates
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor(PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ,'fac',Bback)
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor(Peall,'fac',Bback)
% 
% Rotate tensor into user-defined coordinate system
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor(Peall,'rot',xnew,[ynew,znew])
%
% Rotate tensor from spacecraft coordinates into GSE coordinates
% [PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp] = mms.rotate_tensor(Peall,'gse',MMSnum)
%
%
% Function to rotate the pressure tensor term into field-aligned coordinates or 
% another coordinate system. 
% For rotate into field-aligned coordinates PeZZp is aligned with background magnetic
% field, PeXXp is closest to the x-direction and PeYYp completes the system.
% 
% Written by D. B. Graham.
%
% Input: 
%       PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ - Pressure terms or temperature in TSeries 
%       Peall - TSERIES of all tensor terms with column order PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ
%       'fac' - Transform tensor into field-aligned coordinates. PeZZp is
%           field-aligned. PeXXp and PeYYp are closest orthogonal components to
%           the x and y directions
%           * Bback - Background magnetic field (TSERIES format)
%       'rot' - Transform tensor into an arbitrary coordinate system
%           * xnew - new x-direction (required after 'rot')
%           * ynew, znew - new y and z directions (if not included y and 
%               z directions are orthogonal to xnew and closest to the orginal 
%               y and z directions).
%       'gse' - Transform tensor into GSE coordinates (not yet
%           implemented).
%           * MMSnum - MMS spacecraft number 1--4.
% Output: 
%       PeXXp,PeXYp,PeXZp,PeYYp,PeYZp,PeZZp - Pressure or temperature terms in
%       field-aligned, user-defined, or GSE coordinates
% 

% Check input and load pressure/temperature terms
if isa(varargin{2},'char'),
    rotflag = varargin{2};
    rotflagpos = 2;
    Peall = varargin{1};
    Petimes = Peall.time;
    Ptensor = zeros(length(Petimes),3,3);
    Ptensor(:,1,1) = Peall.data(:,1);
    Ptensor(:,2,1) = Peall.data(:,2);
    Ptensor(:,3,1) = Peall.data(:,3);
    Ptensor(:,1,2) = Peall.data(:,2);
    Ptensor(:,2,2) = Peall.data(:,4);
    Ptensor(:,3,2) = Peall.data(:,5);
    Ptensor(:,1,3) = Peall.data(:,3);
    Ptensor(:,2,3) = Peall.data(:,5);
    Ptensor(:,3,3) = Peall.data(:,6);
elseif isa(varargin{7},'char'),
    rotflag = varargin{7};
    rotflagpos = 7;
    Petimes = varargin{1}.time;
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
else
    irf_log('proc','Something is wrong with the input.')
    help mms.rotate_tensor;
    return;
end    

Rotmat = zeros(length(Petimes),3,3);

if (rotflag(1) == 'f'),
    irf_log('proc','Transforming tensor into field-aligned coordinates.')
    if (nargin == rotflagpos)
        irf_log('proc','B TSeries is missing.')
        return;
    end
    Bback = varargin{rotflagpos+1};
    Bback = Bback.resample(Petimes);
    Bmag = Bback.abs.data;
    Rz = Bback.data./([Bmag Bmag Bmag]);
    Rx = [1 0 0];
    Ry = irf_cross(Rz,Rx);
    Rmag = irf_abs(Ry,1);
    Ry = Ry./[Rmag Rmag Rmag];
    Rx = irf_cross(Ry, Rz);
    Rmag = irf_abs(Rx,1);
    Rx = Rx./[Rmag Rmag Rmag];
    Rotmat(:,1,:) = Rx;
    Rotmat(:,2,:) = Ry;
    Rotmat(:,3,:) = Rz;
elseif (rotflag(1) == 'r'),
    irf_log('proc','Transforming tensor into user defined coordinate system.')
    if (nargin == rotflagpos),
        irf_log('proc','Vector(s) is(are) missing.')
        return;
    end
    vectors = varargin(rotflagpos+1:end);
    if (numel(vectors) == 1),
        Rx = vectors{1};
        if(length(Rx) ~= 3);
            irf_log('proc','Vector format not recognized.')
            return;
        end
        Rx = Rx/norm(Rx);
        Ry = [0 1 0];
        Rz = cross(Rx,Ry);
        Rz = Rz/norm(Rz);
        Ry = cross(Rz,Rx);
        Ry = Ry/norm(Ry);
    elseif (numel(vectors) == 3)
        Rx = vectors{1};
        Ry = vectors{2};
        Rz = vectors{3};
        Rx = Rx/norm(Rx);
        Ry = Ry/norm(Ry);
        Rz = Rz/norm(Rz);
        % TO DO: add check that vectors are orthogonal
    else
        irf_log('proc','Vector format not recognized.')
        return; 
    end
    Rotmat(:,1,:) = ones(length(Petimes),1)*Rx;
    Rotmat(:,2,:) = ones(length(Petimes),1)*Ry;
    Rotmat(:,3,:) = ones(length(Petimes),1)*Rz;    
elseif (rotflag(1) == 'g'),
    irf_log('proc','Transforming tensor into GSE coordinates.')
    SCnum = varargin{rotflagpos+1};
	Tint = irf.tint(Petimes.start.utc,Petimes.stop.utc);
	c_eval('defatt = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',Tint);',SCnum);
	c_eval('defatt.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',Tint).zdec;',SCnum);
    defatt = mms_removerepeatpnts(defatt);
    
    % Development of transformation matrix follows modified version of mms_dsl2gse.m
    ttDefatt = EpochTT(defatt.time);
    zra = irf.ts_scalar(ttDefatt,defatt.zra);
    zdec = irf.ts_scalar(ttDefatt,defatt.zdec);
    zra = zra.resample(Petimes);
    zdec = zdec.resample(Petimes);
    [x,y,z] = sph2cart(zra.data*pi/180,zdec.data*pi/180,1);
    saxTmp = irf.geocentric_coordinate_transformation([Petimes.epochUnix x y z],'gei>gse');
    spin_axis = saxTmp(:,2:4);
    Rx = spin_axis(:,1); 
    Ry = spin_axis(:,2); 
    Rz = spin_axis(:,3);
    a = 1./sqrt(Ry.^2+Rz.^2);
    Rotmat(:,1,:) = [a.*(Ry.^2+Rz.^2) -a.*Rx.*Ry -a.*Rx.*Rz];
    Rotmat(:,2,:) = [0*a a.*Rz	-a.*Ry];
    Rotmat(:,3,:) = [Rx	Ry Rz];
else
    irf_log('proc','Flag is not recognized.')
    return;
end

Ptensorp = zeros(length(Petimes),3,3);
for ii = 1:length(Petimes);
    rottemp = squeeze(Rotmat(ii,:,:));
    Ptensorp(ii,:,:) = rottemp*(squeeze(Ptensor(ii,:,:))*transpose(rottemp));
end

PeXXp = irf.ts_scalar(Petimes,Ptensorp(:,1,1));
PeXYp = irf.ts_scalar(Petimes,Ptensorp(:,1,2));
PeXZp = irf.ts_scalar(Petimes,Ptensorp(:,1,3));
PeYYp = irf.ts_scalar(Petimes,Ptensorp(:,2,2));
PeYZp = irf.ts_scalar(Petimes,Ptensorp(:,2,3));
PeZZp = irf.ts_scalar(Petimes,Ptensorp(:,3,3));

end