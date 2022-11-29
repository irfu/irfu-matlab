function particlemoments = psd_moments(varargin)
% PSD_MOMENTS compute particle moments from a 3D distribution function
%
% Input:
%   pdist - TSeries of the full particle distribution of electrons or ions
%   in PDist data format (burst and fast)
%   SCpot - TSeries of spacecraft potential.
%
% Optional Inputs:
%   'energyrange' - set energy range in eV to integrate over [E_min E_max].
%   energy range is applied to energy0 (if applicable) and the same elements are used for energy1 to
%   ensure that the same number of points are integrated over.
%   'noscpot' - set to 1 to set spacecraft potential to zero. Calculates moments without
%   correcting for spacecraft potential.
%   'enchannels' - set energy channels to integrate over [min max]; min and max
%   between must be between 1 and total number of energy channels.
%   'partialmoms' - use a binary array (or TSeries) (pmomsarr) to select which psd points are used
%   in the moments calculation. pmomsarr must be a binary array (1s and 0s, 1s correspond to points used).
%   Array (or data of TSeries) must be the same size as pdist.data. For
%   examples see Example_MMS_partialmoments.
%   'innerelec' ('on') - innerelectron potential for electron moments
%
% Output:
%   psd_moments - structure containing the particle moments: density, bulk
%   velocity, pressure, temperature, particle heat flux vector as TSeries'.
%
% Tint = irf.tint('2015-10-30T05:15:20.00Z/2015-10-30T05:16:20.00Z');
% ePDist = mms.get_data('PDe_fpi_brst_l2',Tint,1);
% SCpot = mms.get_data('V_edp_brst_l2',Tint,1);
% particlemoments = mms.psd_moments(ePDist,SCpot,'energyrange',[1 1000]);

% 1. basic
flag_dE = 0;
flag_dE_SameDim = 0;            % flag for dimension check of energy, energy_plus & energy_minus
flag_same_e = 0;
flag_innerelec = 0;
flag_dphi = 0;
flag_dtheta = 0;
%W_innerelec = 3.5;         % [eV] scpot + W_innerelec for electron moments calculation; 2018-01-26, wy;

% First input check
if (nargin < 2)
  nargin
  help psd_moments;
  return;
end

pdist = varargin{1};
SCpot = varargin{2};
args=varargin(3:end);
pdist = pdist.convertto('s^3/m^6');

energy = pdist.depend{1};
phi = pdist.depend{2};
theta = pdist.depend{3};
particletype = pdist.species;

intenergies = 1:size(energy,2);

if isfield(pdist.ancillary,'energy0') && isfield(pdist.ancillary,'energy1') && isfield(pdist.ancillary,'esteptable')
	if sum(abs(pdist.ancillary.energy0-pdist.ancillary.energy1)) < 0.0001
    flag_same_e = 1;
  end
end

%resample SCpot to same resolution as particle distributions
SCpot = SCpot.resample(pdist);

%check theta dimensions
thetasize = size(theta);
if thetasize(1) > thetasize(2)
  theta = theta';
end

if numel(args)>0
  options=1;
else
  options=0;
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
    case 'noscpot'
      if numel(args)>1 && args{2}
        SCpot.data = zeros(size(SCpot.data));
        irf.log('notice','Setting spacecraft potential to zero.');
      end
    case 'enchannels'
      if numel(args)>1 && isnumeric(args{2})
        intenergies = args{2}(1):args{2}(2);
      end
    case 'partialmoms'
      if numel(args)>1
        partialmoms = args{2};
        if isa(partialmoms,'TSeries')
          partialmoms = partialmoms.data;
        end
        % Check size of partialmoms
        if (size(partialmoms) == size(pdist.data))
          sumones = sum(sum(sum(sum(partialmoms))));
          sumzeros = sum(sum(sum(sum(-partialmoms+1))));
          if ((sumones+sumzeros) == numel(pdist.data))
            irf.log('notice','partialmoms is correct. Partial moments will be calculated');
            pdist.data = pdist.data.*partialmoms;
          else
            irf.log('notice','All values are not ones and zeros in partialmoms. Full moments will be calculated.');
          end
        else
          irf.log('notice','Size of partialmoms is wrong. Full moments will be calculated.');
        end
      end
    case 'innerelec'
      if numel(args)>1
        innerelec_tmp = args{2};
        if (strcmp(innerelec_tmp, 'on') && (particletype(1) == 'e'))
          flag_innerelec = 1;
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

irf.log('notice',strcat('Integrating over energy levels ',num2str(min(intenergies)),' to ',num2str(max(intenergies))));

% Define constants
Units = irf_units; % Use IAU and CODATA values for fundamental constants.
qe = Units.e;
kb = Units.kB;

if (particletype(1) == 'e')
  pmass = Units.me;
  irf.log('notice','Particles are electrons');
elseif (particletype(1) == 'i')
  pmass = Units.mp;
  SCpot.data = -SCpot.data;
  irf.log('notice','Particles are Ions');
else
  particlemoments = NaN;
  irf.log('critical','Could not identify the particle type');
  return;
end

% Define arrays for output
%sizedist = size(pdist.data);
%n_psd_e32 = zeros(length(pdist.time), 32);
%n_psd_e32_phi_theta = zeros(sizedist(1), sizedist(2), sizedist(3), sizedist(4));
V_psd = zeros(length(pdist.time), 3);
P_psd = zeros(length(pdist.time), 3, 3);
P2_psd = zeros(length(pdist.time), 3, 3);
H_psd = zeros(length(pdist.time), 3);

tic

% Calculate angle differences
if isfield(pdist.ancillary,'delta_phi_minus') && isfield(pdist.ancillary,'delta_phi_plus') 
  deltaphi = pdist.ancillary.delta_phi_plus+pdist.ancillary.delta_phi_minus;
  flag_dphi = 1;
else
  deltaphi = median(diff(pdist.depend{2}(1,:)));
end
if isfield(pdist.ancillary,'delta_theta_minus') && isfield(pdist.ancillary,'delta_theta_plus') 
  deltatheta = pdist.ancillary.delta_theta_plus+pdist.ancillary.delta_theta_minus;
  flag_dtheta = 1;
else
  deltatheta = median(diff(pdist.depend{3}(1,:)));
end
deltaphi = deltaphi*pi/180;
deltatheta = deltatheta*pi/180;
deltaang = deltaphi*deltatheta';

if size(pdist.depend{2},1) > 1
  phitr = phi;
else
  phitr = ones(size(pdist.time))*phi;
end

if size(pdist.depend{3},1) > 1
  thetatr = theta;
else
  thetatr = ones(size(pdist.time))*theta;
end

if isfield(pdist.ancillary,'delta_energy_minus') && isfield(pdist.ancillary,'delta_energy_plus')
  flag_dE = 1;
  energy_minus = pdist.ancillary.delta_energy_minus;
  energy_plus = pdist.ancillary.delta_energy_plus;
  % check energy & energy_minus & energy_plus dimensions
  energy_size = size(energy);
  energy_minus_size = size(energy_minus);
  energy_plus_size = size(energy_plus);
  if isequal(energy_size, energy_minus_size) && isequal(energy_size, energy_plus_size)
      flag_dE_SameDim = 1;      % same dimensions goto Line 230 section. 
  end
end

% Calculate speed widths associated with each energy channel.
energycorr = energy - SCpot.data*ones(size(energy(1,:)));
v = real(sqrt(2*qe*(energycorr)/pmass));
if flag_same_e && flag_dE && ~flag_dE_SameDim
  energyupper = energy + ones(size(pdist.time)).*energy_plus;
  energylower = energy - ones(size(pdist.time)).*energy_minus;
  vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
  vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
elseif flag_same_e && ~flag_dE && ~flag_dE_SameDim
  temp0 = 2*energy(:,1)-energy(:,2);
  tempend = 2*energy(:,end)-energy(:,end-1);
  energyall = [temp0 energy tempend];
  diffenall = diff(energyall,1,2);
  energyupper = 10.^(log10(energy+diffenall(:,2:end)/2));
  energylower = 10.^(log10(energy-diffenall(:,1:end-1)/2));
  vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
  vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
elseif flag_same_e && flag_dE && flag_dE_SameDim
    energyupper = energy + energy_plus;
    energylower = energy - energy_minus;
    vupper = sqrt(2*qe*(energyupper - SCpot.data*ones(size(energy(1,:))))/pmass);
    vlower = sqrt(2*qe*(energylower - SCpot.data*ones(size(energy(1,:))))/pmass);
elseif ~flag_same_e % && ~flag_dE && ~flag_dE_SameDim
  energy0 = pdist.ancillary.energy0;
  energy1 = pdist.ancillary.energy1;
  if size(energy1,2) == 1
    energy1 = energy1';
  end
  if size(energy0,2) == 1
    energy0 = energy0';
  end
  esteptable = pdist.ancillary.esteptable;
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
  energyupper = esteptablemat.*(ones(size(pdist.time))*energyupper1)+abs(esteptablemat-1).*(ones(size(pdist.time))*energyupper0);
  energylower = esteptablemat.*(ones(size(pdist.time))*energylower1)+abs(esteptablemat-1).*(ones(size(pdist.time))*energylower0);
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
phimat = permute(repmat(phitr,1,1,length(intenergies),length(thetatr(1,:))),[1 3 2 4]);

if flag_dphi && flag_dtheta
  deltaang = repmat(deltaang,1,1,length(pdist.time),length(intenergies));
  deltaang = permute(deltaang,[3 4 1 2]);
end

% Compute moments
% Density
M = ones(size(thetamat)).*sind(thetamat).*deltaang;
n_psd = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^2.*deltavmat,4),3),2);
irf.log('notice','Density Calculated');

% Particle flux and heat flux vector
M = -cosd(phimat).*sind(thetamat).*sind(thetamat).*deltaang;
V_psd(:,1) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^3.*deltavmat,4),3),2);
H_psd(:,1) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^5.*deltavmat,4),3),2);
M = -sind(phimat).*sind(thetamat).*sind(thetamat).*deltaang;
V_psd(:,2) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^3.*deltavmat,4),3),2);
H_psd(:,2) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^5.*deltavmat,4),3),2);
M = -ones(size(phimat)).*sind(thetamat).*cosd(thetamat).*deltaang;
V_psd(:,3) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^3.*deltavmat,4),3),2);
H_psd(:,3) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^5.*deltavmat,4),3),2);
irf.log('notice','Velocity Calculated');

% Particle thermal pressure
M = cosd(phimat).^2.*(sind(thetamat).^3).*deltaang; % xx
P_psd(:,1,1) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^4.*deltavmat,4),3),2);
M = sind(phimat).^2.*(sind(thetamat).^3).*deltaang; % yy
P_psd(:,2,2) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^4.*deltavmat,4),3),2);
M = ones(size(phimat)).*sind(thetamat).*cosd(thetamat).^2.*deltaang; % zz
P_psd(:,3,3) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^4.*deltavmat,4),3),2);
M = cosd(phimat).*sind(phimat).*sind(thetamat).^3.*deltaang; % xy
P_psd(:,1,2) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^4.*deltavmat,4),3),2);
M = cosd(phimat).*sind(thetamat).^2.*cosd(thetamat).*deltaang; % xz
P_psd(:,1,3) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^4.*deltavmat,4),3),2);
M = sind(phimat).*sind(thetamat).^2.*cosd(thetamat).*deltaang; % yz
P_psd(:,2,3) = irf.nansum(irf.nansum(irf.nansum(M.*pdist.data(:,intenergies,:,:).*vmat.^4.*deltavmat,4),3),2);
irf.log('notice','Pressure Calculated');

% Compute moments in SI units
P_psd = pmass*P_psd;
V_psd = V_psd./[n_psd n_psd n_psd];
P2_psd(:,1,1) = P_psd(:,1,1);
P2_psd(:,1,2) = P_psd(:,1,2);
P2_psd(:,1,3) = P_psd(:,1,3);
P2_psd(:,2,2) = P_psd(:,2,2);
P2_psd(:,2,3) = P_psd(:,2,3);
P2_psd(:,3,3) = P_psd(:,3,3);
P2_psd(:,2,1) = P2_psd(:,1,2); P2_psd(:,3,1) = P2_psd(:,1,3); P2_psd(:,3,2) = P2_psd(:,2,3);

P_psd(:,1,1) = P_psd(:,1,1)-pmass*n_psd.*V_psd(:,1).*V_psd(:,1);
P_psd(:,1,2) = P_psd(:,1,2)-pmass*n_psd.*V_psd(:,1).*V_psd(:,2);
P_psd(:,1,3) = P_psd(:,1,3)-pmass*n_psd.*V_psd(:,1).*V_psd(:,3);
P_psd(:,2,2) = P_psd(:,2,2)-pmass*n_psd.*V_psd(:,2).*V_psd(:,2);
P_psd(:,2,3) = P_psd(:,2,3)-pmass*n_psd.*V_psd(:,2).*V_psd(:,3);
P_psd(:,3,3) = P_psd(:,3,3)-pmass*n_psd.*V_psd(:,3).*V_psd(:,3);
P_psd(:,2,1) = P_psd(:,1,2); P_psd(:,3,1) = P_psd(:,1,3); P_psd(:,3,2) = P_psd(:,2,3);

ntemp = reshape([n_psd n_psd n_psd;n_psd n_psd n_psd;n_psd n_psd n_psd],length(n_psd),3,3);

Ptrace = (P_psd(:,1,1)+P_psd(:,2,2)+P_psd(:,3,3));
T_psd = P_psd./ntemp/kb;
T_psd(:,2,1) = T_psd(:,1,2); T_psd(:,3,1) = T_psd(:,1,3); T_psd(:,3,2) = T_psd(:,2,3);

H_psd = pmass/2*H_psd;
Vabs2 = V_psd(:,1).^2+V_psd(:,2).^2+V_psd(:,3).^2;
H_psd(:,1) = H_psd(:,1)-(V_psd(:,1).*P_psd(:,1,1)+V_psd(:,2).*P_psd(:,1,2)+V_psd(:,3).*P_psd(:,1,3))-0.5*V_psd(:,1).*Ptrace-0.5*pmass*n_psd.*Vabs2.*V_psd(:,1);
H_psd(:,2) = H_psd(:,2)-(V_psd(:,1).*P_psd(:,1,2)+V_psd(:,2).*P_psd(:,2,2)+V_psd(:,3).*P_psd(:,2,3))-0.5*V_psd(:,2).*Ptrace-0.5*pmass*n_psd.*Vabs2.*V_psd(:,2);
H_psd(:,3) = H_psd(:,3)-(V_psd(:,1).*P_psd(:,1,3)+V_psd(:,2).*P_psd(:,2,3)+V_psd(:,3).*P_psd(:,3,3))-0.5*V_psd(:,3).*Ptrace-0.5*pmass*n_psd.*Vabs2.*V_psd(:,3);

% Convert to typical units (/cc, km/s, nP, eV, and ergs/s/cm^2).
n_psd = n_psd/1e6;
%n_psd_e32 = n_psd_e32/1e6;
%n_psd_e32_phi_theta = n_psd_e32_phi_theta / 1e6;
V_psd = V_psd/1e3;
P_psd = P_psd*1e9;
P2_psd = P2_psd*1e9;
T_psd = T_psd*kb/qe;
H_psd = H_psd*1e3;

% Construct TSeries'
n_psd = irf.ts_scalar(pdist.time,n_psd);
%n_psd_e32 = irf.ts_scalar(pdist.time, n_psd_e32);
%if isstruct(phi)
%  n_psd_skymap = PDist(pdist.time, n_psd_e32_phi_theta, 'skymap', energy, phi.data, thetak);
%else
%  n_psd_skymap = PDist(pdist.time, n_psd_e32_phi_theta, 'skymap', energy, phi, thetak);
%end
%n_psd_skymap.userData = pdist.userData;
%n_psd_skymap.units = 'cm^{-3}';
V_psd = irf.ts_vec_xyz(pdist.time,V_psd);
P_psd = irf.ts_tensor_xyz(pdist.time,P_psd);
P2_psd = irf.ts_tensor_xyz(pdist.time,P2_psd);
T_psd = irf.ts_tensor_xyz(pdist.time,T_psd);
H_psd = irf.ts_vec_xyz(pdist.time,H_psd);

toc;

% make structure for output
particlemoments =struct('n_psd',n_psd,'V_psd',V_psd,'P_psd',P_psd,'P2_psd',P2_psd, ...
  'T_psd',T_psd,'H_psd',H_psd); %, 'n_psd_e32', n_psd_e32, 'n_psd_skymap', n_psd_skymap);

end
