function epsilon = calculate_epsilon(varargin)
%CALCULATE_EPSILON Calculate epsilon parameter using model distribution
%
%   Examples:
%     epsilon = mms.calculate_epsilon(PDist,modelPDist,n,SCpot,'enchannels',4:32)
%
%   Input:
%     PDist - observed particle distribution (skymap). Must be in PDist format.
%     modelPDist - model particle distribution (skymap). Must be in PDist format.
%     ne - number density (TSeries).
%     SCpot - Spacecraft potential (TSeries).
%   Options:
%     'enchannels' - set energy channels to integrate over [min max]; min and max
%       between must be between 1 and num.
%     'energyrange' - set energy range in eV to integrate over [E_min E_max].
%
%   Output:
%     epsilon - epsilon parameter (TSeries).
%
% Written by D. B. Graham
%
% See also mms.make_model_dist

% First input check
if (nargin < 4)
  epsilon = NaN;
  help mms.calculate_epsilon;
  return;
end

tic;

flag_dE = 0;
flag_same_e = 0;
flag_dphi = 0;
flag_dtheta = 0;

PDist = varargin{1};
modelPDist = varargin{2};
n = varargin{3};
SCpot = varargin{4};
args=varargin(5:end);

energy = PDist.depend{1};
phi = PDist.depend{2};
theta = PDist.depend{3};
particletype = PDist.species;

if numel(args)>0
  options=1;
else
  options=0;
end

if abs(median(diff(PDist.time-n.time))) > 0
  epsilon = NaN;
  irf.log('critical','PDist and moments have different times.')
  return;
end

% Default energy channels used to compute epsilon.
intenergies = 1:size(energy,2);

if isfield(PDist.ancillary,'energy0') && isfield(PDist.ancillary,'energy1') && isfield(PDist.ancillary,'esteptable')
  if sum(abs(PDist.ancillary.energy0-PDist.ancillary.energy1)) < 0.0001
    flag_same_e = 1;
  end
end

%check theta dimensions
thetasize = size(theta);
if thetasize(1) > thetasize(2)
  theta = theta';
end

while options
  l = 2;
  switch(lower(args{1}))
    case 'energyrange'
      if numel(args)>1 && isnumeric(args{2})
        Eminmax = args{2};
        starte = find(energy(1,:) > Eminmax(1));
        starte = starte(1);
        ende = find(energy(1,:) < Eminmax(2));
        ende = ende(end);
        intenergies = starte:ende;
        irf.log('notice','Using partial energy range');
      end
    case 'enchannels'
      if numel(args)>1 && isnumeric(args{2})
        intenergies = args{2}(1):args{2}(end);
      end
    otherwise
      irf.log('critical',['Unknown flag: ' args{1}]);
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), options=0; end
end

% Resample SCpot
SCpot = SCpot.resample(n);

% Remove zero count points from final calculation
%modelPDist.data(PDist.data <= 0) = 0;

modelPDist = modelPDist.convertto('s^3/m^6');
PDist = PDist.convertto('s^3/m^6');

Pdiff = abs(PDist.data-modelPDist.data);

% Define constants
Units = irf_units;
qe = Units.e;

% Check whether particles are electrons or ions
if (particletype(1) == 'e')
  pmass = Units.me;
  irf.log('notice','Particles are electrons');
elseif (particletype(1) == 'i')
  pmass = Units.mp;
  SCpot.data = -SCpot.data;
  irf.log('notice','Particles are Ions');
else
  epsilon = NaN;
  irf.log('critical','Could not identify the particle type');
  return;
end

% Calculate angle differences
if isfield(PDist.ancillary,'delta_phi_minus') && isfield(PDist.ancillary,'delta_phi_plus')
  deltaphi = PDist.ancillary.delta_phi_plus+PDist.ancillary.delta_phi_minus;
  flag_dphi = 1;
else
  deltaphi = median(diff(PDist.depend{2}(1,:)));
end
if isfield(PDist.ancillary,'delta_theta_minus') && isfield(PDist.ancillary,'delta_theta_plus')
  deltatheta = PDist.ancillary.delta_theta_plus+PDist.ancillary.delta_theta_minus;
  flag_dtheta = 1;
else
  deltatheta = median(diff(PDist.depend{3}(1,:)));
end
deltaphi = deltaphi*pi/180;
deltatheta = deltatheta*pi/180;
deltaang = deltaphi*deltatheta';

if size(PDist.depend{2},1) > 1
  phitr = phi;
else
  phitr = ones(size(PDist.time))*phi;
end

if size(PDist.depend{3},1) > 1
  thetatr = theta;
else
  thetatr = ones(size(PDist.time))*theta;
end

if isfield(PDist.ancillary,'delta_energy_minus') && isfield(PDist.ancillary,'delta_energy_plus')
  flag_dE = 1;
  energy_minus = PDist.ancillary.delta_energy_minus;
  energy_plus = PDist.ancillary.delta_energy_plus;
end

% Calculate speed widths associated with each energy channel.
energycorr = energy - SCpot.data*ones(size(energy(1,:)));
v = real(sqrt(2*qe*(energycorr)/pmass));
if flag_same_e && flag_dE
  energyupper = energy + energy_plus;
  energylower = energy - energy_minus;
  vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
  vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
elseif flag_same_e && ~flag_dE
  temp0 = 2*energy(:,1)-energy(:,2);
  tempend = 2*energy(:,end)-energy(:,end-1);
  energyall = [temp0 energy tempend];
  diffenall = diff(energyall,1,2);
  energyupper = 10.^(log10(energy+diffenall(:,2:end)/2));
  energylower = 10.^(log10(energy-diffenall(:,1:end-1)/2));
  vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
  vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
elseif ~flag_same_e && flag_dE
  energyupper = energy + energy_plus;
  energylower = energy - energy_minus;
  vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
  vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
elseif ~flag_same_e && ~flag_dE
  energy0 = PDist.ancillary.energy0;
  energy1 = PDist.ancillary.energy1;
  if size(energy1,2) == 1
    energy1 = energy1';
  end
  if size(energy0,2) == 1
    energy0 = energy0';
  end
  esteptable = PDist.ancillary.esteptable;
  temp0 = 2*energy0(1)-energy0(2);
  tempend = 2*energy0(end)-energy0(end-1);
  energyall = [temp0 energy0 tempend];
  diffenall = diff(energyall);
  energyupper0 = 10.^(log10(energy0+diffenall(2:1:end)/2));
  energylower0 = 10.^(log10(energy0-diffenall(1:1:end-1)/2));
  temp0 = 2*energy1(1)-energy1(2);
  tempend = 2*energy1(end)-energy1(end-1);
  energyall = [temp0 energy1 tempend];
  diffenall = diff(energyall);
  energyupper1 = 10.^(log10(energy1+diffenall(2:1:end)/2));
  energylower1 = 10.^(log10(energy1-diffenall(1:1:end-1)/2));
  esteptablemat = single(esteptable)*ones(size(energy0));
  energyupper = esteptablemat.*(ones(size(PDist.time))*energyupper1)+abs(esteptablemat-1).*(ones(size(PDist.time))*energyupper0);
  energylower = esteptablemat.*(ones(size(PDist.time))*energylower1)+abs(esteptablemat-1).*(ones(size(PDist.time))*energylower0);
  vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
  vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
end

vupper(vupper < 0) = 0;
vlower(vlower < 0) = 0;
vupper = real(vupper);
vlower = real(vlower);

deltav = (vupper-vlower);
vmat = repmat(v,1,1,length(phitr(1,:)),length(theta));
deltavmat = repmat(deltav,1,1,length(phitr(1,:)),length(theta));
vmat = vmat(:,intenergies,:,:);
deltavmat = deltavmat(:,intenergies,:,:);

thetamat = permute(repmat(thetatr,1,1,length(intenergies),length(phitr(1,:))),[1 3 4 2]);
%phimat = permute(repmat(phitr,1,1,length(intenergies),length(thetatr(1,:))),[1 3 2 4]);

if flag_dphi && flag_dtheta
  deltaang = repmat(deltaang,1,1,length(PDist.time),length(intenergies));
  deltaang = permute(deltaang,[3 4 1 2]);
end


M = ones(size(thetamat)).*sind(thetamat).*deltaang;
epsilon = irf.nansum(irf.nansum(irf.nansum(M.*Pdiff(:,intenergies,:,:).*vmat.^2.*deltavmat,4),3),2);

epsilon = epsilon/1e6./(n.data*2);
epsilon = irf.ts_scalar(PDist.time,epsilon);

toc;

end
