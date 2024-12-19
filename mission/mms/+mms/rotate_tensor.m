function Pe = rotate_tensor(varargin)
% MMS.ROTATE_TENSOR rotate pressure or temperature tensor into another
% coordinate system
%
% Examples:
% Rotate tensor into field-aligned coordinates
% Pe = mms.rotate_tensor(PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ,'fac',Bback)
% Pe = mms.rotate_tensor(Peall,'fac',Bback)
% Pe = mms.rotate_tensor(Peall,'fac',Bback,'pp' or 'qq')
%
% Rotate tensor into user-defined coordinate system
% Pe = mms.rotate_tensor(Peall,'rot',xnew)
% Pe = mms.rotate_tensor(Peall,'rot',xnew,ynew,znew)
%
% Rotate tensor from spacecraft coordinates into GSE coordinates
% Pe = mms.rotate_tensor(Peall,'gse',MMSnum)
%
% Function to rotate the pressure tensor term into field-aligned coordinates or
% another coordinate system.
%
% Written by D. B. Graham.
%
% Input:
%       PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ - Pressure terms or temperature in TSeries
%       Peall - TSERIES of all tensor terms with column order PeXX,PeXY,PeXZ,PeYY,PeYZ,PeZZ
%               OR 3x3 data (can have tensorOrder=2) with [PeXX,PeXY,PeXZ;PeYX,PeYY,PeYZ;PeZX,PeZY,PeZZ]
%       'fac' - Transform tensor into field-aligned coordinates.
%           * Bback - Background magnetic field (TSERIES format)
%           * 'pp' - optional flag to rotate perpendicular components so
%             P_perp1 = P_perp2.
%           * 'qq' - optional flag to rotate perpendicular components so
%             P_perp1 and P_perp2 are most unequal, sets P23 to zero
%       'rot' - Transform tensor into an arbitrary coordinate system
%           * xnew - new x-direction (required after 'rot')
%           * ynew, znew - new y and z directions (if not included y and
%               z directions are orthogonal to xnew and closest to the orginal
%               y and z directions).
%       'gse' - Transform tensor into GSE coordinates
%           * MMSnum - MMS spacecraft number 1--4.
%
% Output:
%       Pe - Pressure or temperature terms in field-aligned, user-defined,
%       or GSE coordinates. Tseries with 3*3 in data.
%       For 'fac' Pe = [Ppar P12 P13; P12 Pperp1 P23; P13 P23 Pperp2]
%       For 'rot' and 'gse' Pe = [Pxx Pxy Pxz; Pxy Pyy Pyz; Pxz Pyz Pzz]
%

rtntensor = 0;
% Check input and load pressure/temperature terms
if isa(varargin{2},'char')
  rotflag = varargin{2};
  rotflagpos = 2;
  Peall = varargin{1};
  Petimes = Peall.time;
  if ndims(Peall.data) == 3
    Ptensor = Peall.data;
    rtntensor = 1;
  else
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
  end
elseif isa(varargin{7},'char')
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
  irf.log('critical','Something is wrong with the input.');
  help mms.rotate_tensor;
  return;
end

ppeq = 0;
qqeq = 0;
Rotmat = zeros(length(Petimes),3,3);

if (rotflag(1) == 'f')
  irf.log('notice','Transforming tensor into field-aligned coordinates.');
  if (nargin == rotflagpos)
    irf.log('critical','B TSeries is missing.');
    return;
  end
  Bback = varargin{rotflagpos+1};
  Bback = Bback.resample(Petimes);
  if (nargin == 4)
    if (isa(varargin{4},'char') && varargin{4}(1) == 'p')
      ppeq = 1;
    elseif (isa(varargin{4},'char') && varargin{4}(1) == 'q')
      qqeq = 1;
    else
      irf.log('critical','Flag not recognized no additional rotations applied.');
    end
  end
  if (nargin == 9)
    if (isa(varargin{9},'char') && varargin{9}(1) == 'p')
      ppeq = 1;
    elseif (isa(varargin{9},'char') && varargin{9}(1) == 'q')
      qqeq = 1;
    else
      irf.log('critical','Flag not recognized no additional rotations applied.');
    end
  end
  Bvec = Bback/Bback.abs;
  Rx = Bvec.data;
  Ry = [1 0 0];
  Rz = irf_cross(Rx,Ry);
  Rmag = irf_abs(Rz,1);
  Rz = Rz./[Rmag Rmag Rmag];
  Ry = irf_cross(Rz, Rx);
  Rmag = irf_abs(Ry,1);
  Ry = Ry./[Rmag Rmag Rmag];
  Rotmat(:,1,:) = Rx;
  Rotmat(:,2,:) = Ry;
  Rotmat(:,3,:) = Rz;
elseif (rotflag(1) == 'r')
  irf.log('notice','Transforming tensor into user defined coordinate system.');
  if (nargin == rotflagpos)
    irf.log('critical','Vector(s) is(are) missing.');
    return;
  end
  vectors = varargin(rotflagpos+1:end);
  if isscalar(vectors)
    Rx = vectors{1};
    if(length(Rx) ~= 3)
      irf.log('critical','Vector format not recognized.');
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
    irf.log('critical','Vector format not recognized.');
    return;
  end
  Rotmat(:,1,:) = ones(length(Petimes),1)*Rx;
  Rotmat(:,2,:) = ones(length(Petimes),1)*Ry;
  Rotmat(:,3,:) = ones(length(Petimes),1)*Rz;
elseif (rotflag(1) == 'g')
  irf.log('notice','Transforming tensor into GSE coordinates.');
  SCnum = varargin{rotflagpos+1};
  Tint = irf.tint(Petimes.start.utc,Petimes.stop.utc); %#ok<NASGU>
  c_eval('defatt = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',Tint);',SCnum);
  c_eval('defatt.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',Tint).zdec;',SCnum);
  defatt = mms_removerepeatpnts(defatt); %#ok<NODEF>

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
  irf.log('critical','Flag is not recognized.');
  return;
end

Ptensorp = zeros(length(Petimes),3,3);
for ii = 1:length(Petimes)
  rottemp = squeeze(Rotmat(ii,:,:));
  Ptensorp(ii,:,:) = rottemp*(squeeze(Ptensor(ii,:,:))*transpose(rottemp));
end

if ppeq
  irf.log('notice','Rotating tensor so perpendicular diagonal components are equal.');
  theta = 0.5*atan((Ptensorp(:,3,3)-Ptensorp(:,2,2))./(2*Ptensorp(:,2,3)));
  for ii = 1:length(Petimes)
    if isnan(theta(ii));    theta(ii) = 0;  end
    rottemp = [1 0 0; 0 cos(theta(ii)) sin(theta(ii)); 0 -sin(theta(ii)) cos(theta(ii))];
    Ptensorp(ii,:,:) = rottemp*(squeeze(Ptensorp(ii,:,:))*transpose(rottemp));
  end
end

if qqeq
  irf.log('notice','Rotating tensor so perpendicular diagonal components are most unequal.');
  theta = 0.5*atan((2*Ptensorp(:,2,3))./(Ptensorp(:,3,3)-Ptensorp(:,2,2)));
  for ii = 1:length(Petimes)
    rottemp = [1 0 0; 0 cos(theta(ii)) -sin(theta(ii)); 0 sin(theta(ii)) cos(theta(ii))];
    Ptensorp(ii,:,:) = rottemp*(squeeze(Ptensorp(ii,:,:))*transpose(rottemp));
  end
end

% Construct output
if rtntensor
  Pe = irf.ts_tensor_xyz(Petimes,Ptensorp);
else % for backwards compability
  Pe = TSeries(Petimes,Ptensorp);
end
Pe.units = varargin{1}.units;
end
