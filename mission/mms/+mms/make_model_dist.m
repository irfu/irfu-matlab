function modelPDist = make_model_dist(varargin)
%MAKE_MODEL_DIST Make a general bi-Maxwellian distribution function based
%on particle moment data in the same format as PDist.
%
% Example:
% modelPDist = mms.make_model_dist(PDist,Bxyz,SCpot,ne,Ve,Te)
%
%   Input:
%       PDist - particle distribution (skymap). Must be in PDist format.
%       Bxyz - Background magnetic field (TSeries)
%       SCpot - Spacecraft potential (TSeries)
%       n - number density (TSeries)
%       V - bulk velocity (TSeries)
%       T - temperature (TSeries, must be full tensor)
%
%   Options:
%       'isotropic' - set to 1 to make an isotropic model distribution
%       funciton.
%
%   Output:
%       modelPDist - Distribution function in the same format as PDist.
%       Ouput in units of s^3 km^{-6}
%
%   Written by D. B. Graham
%
% See also mms.calculate epsilon

% First input check
if (nargin < 6)
  modelPDist = NaN;
  help mms.calculate_epsilon;
  return;
end

PDist = varargin{1};
Bxyz = varargin{2};
SCpot = varargin{3};
n = varargin{4};
V = varargin{5};
T = varargin{6};

isisotropic = 0;

args=varargin(7:end);
if numel(args)>0
  options=1;
else
  options=0;
end

tic;
% Check that PDist and moments have the same times
if abs(median(diff(PDist.time-n.time))) > 0
  modelPDist = NaN;
  irf.log('critical','PDist and moments have different times.')
  return;
end

while options
  l = 2;
  switch(lower(args{1}))
    case 'isotropic'
      if numel(args)>1 && isnumeric(args{2})
        if args{2} == 1
          isisotropic = 1;
          irf.log('notice','Using isotropic Maxwellian distribution.');
        end
      end
    otherwise
      irf.log('critical',['Unknown flag: ' args{1}]);
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), options=0; end
end


% Convert to SI units
PDist = PDist.convertto('s^3/m^6');

% Resample Bxyz and SCpot to particle data resolution
Bxyz = Bxyz.resample(n);
SCpot = SCpot.resample(n);


% Define directions based on Bxyz and V, calculate relevant temperatures
Tfac = mms.rotate_tensor(T,'fac',Bxyz,'pp'); % N.B makes final distribution gyrotropic
if isisotropic
  Tpar = irf.ts_scalar(Tfac.time,Tfac.trace.data/3);
  Trat = irf.ts_scalar(Tfac.time,ones(size(Tfac.time)));
else
  Tpar = irf.ts_scalar(Tfac.time,Tfac.data(:,1,1));
  Trat = irf.ts_scalar(Tfac.time,Tfac.data(:,1,1)./Tfac.data(:,2,2));
end

[Vpar,Vperp]=irf_dec_parperp(Bxyz,V);
Vpmag = Vperp.abs;
Vpdir = Vperp/Vpmag;
Bdir = Bxyz/Bxyz.abs;

% Define constants
Units = irf_units;
qe = Units.e;

% Check whether particles are electrons or ions
if (PDist.species(1) == 'e')
  pmass = Units.me;
  irf.log('notice','Particles are electrons');
elseif (PDist.species(1) == 'i')
  pmass = Units.mp;
  SCpot.data = -SCpot.data;
  irf.log('notice','Particles are Ions');
else
  modelPDist = NaN;
  irf.log('critical','Could not identify the particle type');
  return;
end

% Convert moments to SI units
vthpar = sqrt(2*Tpar.data*qe/pmass);
vthpar = irf.ts_scalar(n.time,vthpar);
Vpmag.data = Vpmag.data*1e3;
Vpar.data = Vpar.data*1e3;
n.data = n.data*1e6;

% Defines dimensions of array below
lengthphi = length(PDist.depend{2}(1,:));
lengththeta = length(PDist.depend{3});
lengthenergy = length(PDist.ancillary.energy0);

% Get phi
phi = PDist.depend{2};
if size(phi,1) < 2
  phi = ones(size(PDist.time))*phi;
end

% Get energy array
energyarr = PDist.depend{1,1};

% Define Cartesian coordinates
x = zeros(length(PDist.time),lengthphi,lengththeta);
y = zeros(length(PDist.time),lengthphi,lengththeta);
z = zeros(length(PDist.time),lengthphi,lengththeta);
r = zeros(length(PDist.time),lengthenergy);

for ii = 1:length(PDist.time)
  x(ii,:,:) = -cosd(phi(ii,:)')*sind(PDist.depend{3});
  y(ii,:,:) = -sind(phi(ii,:)')*sind(PDist.depend{3});
  z(ii,:,:) = -ones(lengthphi,1)*cosd(PDist.depend{3});
  r(ii,:) = real(sqrt(2*(energyarr(ii,:)-SCpot.data(ii))*qe/pmass));
end
r(r == 0) = 0.0;

% Define rotation vectors based on B and Ve directions
Rz = Bdir.data;
Rx = Vpdir.data;
Ry = cross(Bdir,Vpdir);
Ry = Ry.data;

% Rotated coordinate system for computing bi-Maxwellian distribution
xp = zeros(length(PDist.time),lengthphi,lengththeta);
yp = zeros(length(PDist.time),lengthphi,lengththeta);
zp = zeros(length(PDist.time),lengthphi,lengththeta);

for ii = 1:length(PDist.time)
  xp(ii,:,:) = x(ii,:,:)*Rx(ii,1)+y(ii,:,:)*Rx(ii,2)+z(ii,:,:)*Rx(ii,3);
  yp(ii,:,:) = x(ii,:,:)*Ry(ii,1)+y(ii,:,:)*Ry(ii,2)+z(ii,:,:)*Ry(ii,3);
  zp(ii,:,:) = x(ii,:,:)*Rz(ii,1)+y(ii,:,:)*Rz(ii,2)+z(ii,:,:)*Rz(ii,3);
end

% Make 4D position matrix
xp = repmat(xp,1,1,1,lengthenergy);
%xp = squeeze(permute(xp,[1 4 2 3]));
xp = permute(xp,[1 4 2 3]);
yp = repmat(yp,1,1,1,lengthenergy);
%yp = squeeze(permute(yp,[1 4 2 3]));
yp = permute(yp,[1 4 2 3]);
zp = repmat(zp,1,1,1,lengthenergy);
%zp = squeeze(permute(zp,[1 4 2 3]));
zp = permute(zp,[1 4 2 3]);
rmat = repmat(r,1,1,lengthphi,lengththeta);

% Can use 4D matrices. Too much memory is used for my computer; too
% slow for large files.
% Make 4D matrices required for distribution calculation
%nmat = repmat(n.data,1,lengthenergy,lengthphi,lengththeta);
%Tratmat = repmat(Trat.data,1,lengthenergy,lengthphi,lengththeta);
%Vpmagmat = repmat(Vpmag.data,1,lengthenergy,lengthphi,lengththeta);
%Vparmat = repmat(Vpar.data,1,lengthenergy,lengthphi,lengththeta);
%vthparmat = repmat(vthpar.data,1,lengthenergy,lengthphi,lengththeta);

% Calculate bi-Maxwellian distribution function
%bimaxdist = nmat.*Tratmat./(sqrt(pi^3)*vthparmat.^3);
%bimaxdist = bimaxdist.*exp(-(xp.*rmat-Vpmagmat).^2./(vthparmat.^2).*Tratmat);
%bimaxdist = bimaxdist.*exp(-(yp.*rmat).^2./(vthparmat.^2).*Tratmat);
%bimaxdist = bimaxdist.*exp(-(zp.*rmat-Vparmat).^2./(vthparmat.^2));

% Construct biMaxwellian distribution function
bimaxdist = zeros(size(rmat));
for ii = 1:length(PDist.time)
  coeff = n.data(ii)*Trat.data(ii)/(sqrt(pi^3)*vthpar.data(ii)^3);
  bimaxtemp = coeff*exp(-(xp(ii,:,:,:).*rmat(ii,:,:,:)-Vpmag.data(ii)).^2./(vthpar.data(ii).^2).*Trat.data(ii));
  bimaxtemp = bimaxtemp.*exp(-(yp(ii,:,:,:).*rmat(ii,:,:,:)).^2./(vthpar.data(ii).^2).*Trat.data(ii));
  bimaxtemp = bimaxtemp.*exp(-(zp(ii,:,:,:).*rmat(ii,:,:,:)-Vpar.data(ii)).^2./(vthpar.data(ii).^2));
  bimaxdist(ii,:,:,:) = bimaxtemp;
end

% Make modelPDist file for output
modelPDist = PDist;
modelPDist.data = bimaxdist;
modelPDist = modelPDist.convertto('s^3/km^6');

toc;

end

