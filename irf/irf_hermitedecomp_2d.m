function hermitestruct = irf_hermitedecomp_2d(varargin)
% IRF_HERMITEDECOMP_2D - decompose a series of 2D distributions in the
% hermite basis of functions
%
%   Detailed explanation goes here
%
% Input:
%       PDist2d - Structure or TSeries of PDist of 2D particle distributions.
% Must contain the distributions and speeds. TSeries should be 2D reduced
% distributions (Cartesian coordinates).
%       numhermites - number of hermite functions used in the decomposition.
%
% Options:
%       species - set particle species. 'e', 'p', or 'a' for electrons,
%       protons, or alpha particles.
%
% Output:
%       hermitestruct - structure containing hermite function weights
%
% Notes: Velocities in the structure are converted to km s^{-1}
%
% Example:
%   hermitestruct = irf_hermitedecomp_2d(PDist2D,64,'species','e');
%
% Written by D. B. Graham
tic

if nargin < 2
  help irf_hermitedecomp_2d;
  return;
end

PDist2d = varargin{1};
numhermites = varargin{2};
args=varargin(3:end);
if numel(args)>0
  flag_have_options=true;
else
  flag_have_options=false;
end

species = 'e'; % set for electrons as default
massnotset = true;

while flag_have_options
  l = 2;
  switch(lower(args{1}))
    case 'species'
      if numel(args)>1 && ischar(args{2})
        species = args{2}(1);
        massnotset = false;
      end
    otherwise
      irf_log('fcal',['Unknown flag: ' args{1}])
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), flag_have_options=0; end
end

isTSeries = false;

if isa(PDist2d,'TSeries')
  if massnotset
    species = PDist2d.species(1);
  end
  distributions = PDist2d.data; % Should be units of s m^{-4}
  times = PDist2d.time;
  vxmat = PDist2d.depend{1,1}*1e3; % Convert to m s^{-1}
  vymat = PDist2d.depend{1,2}*1e3; % Convert to m s^{-1}
  dvx = median(diff(PDist2d.depend{1,1}(1,:)))*1e3; % Convert to m s^{-1}
  dvy = median(diff(PDist2d.depend{1,2}(1,:)))*1e3; % Convert to m s^{-1}
  lengthvx = length(vxmat(1,:));
  lengthvy = length(vymat(1,:));
  vxmat = repmat(vxmat,1,1,lengthvy);
  vymat = repmat(vymat,1,1,lengthvx);
  vymat = permute(vymat,[1 3 2]);
  isTSeries = true;
end

if isTSeries == false
  times = PDist2d.t;
  distributions = PDist2d.dist;
  if isfield(PDist2d,'species') && massnotset
    if ischar(PDist2d.species)
      species = PDist2d.species(1);
    end
  end

  timesformat = size(times);
  if timesformat(1) < timesformat(2)
    times = times';
  end

  vxvec = PDist2d.vx*1e3; % convert from km s^-1 to m s^-1
  lengthvx = length(vxvec);
  dvx = median(diff(vxvec));

  vyvec = PDist2d.vy*1e3; % convert from km s^-1 to m s^-1
  lengthvy = length(vyvec);
  dvy = median(diff(vyvec));

  if min(size(vxvec)) == 1
    vxmat = ones(size(times))*vxvec;
  end
  if min(size(vyvec)) == 1
    vymat = ones(size(times))*vyvec;
  end
  vxmat = repmat(vxmat,1,1,lengthvy);
  vymat = repmat(vymat,1,1,lengthvx);
  vymat = permute(vymat,[1 3 2]);
end

% Define constants
Units = irf_units;
mp = Units.mp;
me = Units.me;
qe = Units.e;

if isequal(species,'e')
  ms = me;
elseif isequal(species,'p')
  ms = mp;
elseif isequal(species,'i')
  ms = mp;
elseif isequal(species,'a')
  ms = mp*4;
else
  irf_log('warning','Particle species unknown; defaulting to electrons.')
  ms = me;
end

% Compute moments
nmoms = sum(sum(distributions,3),2)*dvx*dvy;
Vxmoms = sum(sum(vxmat.*distributions,3),2)*dvx*dvy;
Vxmoms = Vxmoms./nmoms;
Vymoms = sum(sum(vymat.*distributions,3),2)*dvx*dvy;
Vymoms = Vymoms./nmoms;
Vxmomsmat = repmat(Vxmoms,1,lengthvx,lengthvy);
Vymomsmat = repmat(Vymoms,1,lengthvx,lengthvy);
Pxxmoms = ms*sum(sum(distributions.*(vxmat-Vxmomsmat).^2,3),2)*dvx*dvy;
Pyymoms = ms*sum(sum(distributions.*(vymat-Vymomsmat).^2,3),2)*dvx*dvy;
Pxymoms = ms*sum(sum(distributions.*(vxmat-Vxmomsmat).*(vymat-Vymomsmat),3),2)*dvx*dvy;
Txxmoms = Pxxmoms./nmoms/qe;
Tyymoms = Pyymoms./nmoms/qe;
Txymoms = Pxymoms./nmoms/qe;
vxth = sqrt(qe*Txxmoms/ms); % Use RMS v_th rather than mean v_th for Hermite decomposition
vyth = sqrt(qe*Tyymoms/ms);
vxthmat = repmat(vxth,1,lengthvx,lengthvy);
vythmat = repmat(vyth,1,lengthvx,lengthvy);

sqrtvxthmat = sqrt(vxthmat);
sqrtvythmat = sqrt(vythmat);

HermitePs = zeros([length(times) numhermites+1 numhermites+1]);
vxnormmat = (vxmat-Vxmomsmat)./vxthmat;
vynormmat = (vymat-Vymomsmat)./vythmat;

c_eval('hermitejj? = irf_hermitefunction(vynormmat,?)./sqrtvythmat;',0:numhermites);

for ii = 0:numhermites
  hermiteii = irf_hermitefunction(vxnormmat,ii)./sqrtvxthmat;
  for jj = 0:numhermites
    c_eval('HermitePs(:,ii+1,jj+1) = sum(sum(hermiteii.*hermitejj?.*distributions,3),2)*dvx*dvy;',jj);
  end
end

hermitestruct = struct;
hermitestruct.hermitespectra = HermitePs;
hermitestruct.orders = 0:numhermites;
hermitestruct.times = times;
hermitestruct.nmoms = nmoms/1e6; % in cm^{-3};
hermitestruct.Vxmoms = Vxmoms/1e3; % in km s^{-1};
hermitestruct.Vymoms = Vymoms/1e3; % in km s^{-1};
hermitestruct.Pxxmoms = Pxxmoms*1e9; % in nPa;
hermitestruct.Pyymoms = Pyymoms*1e9; % in nPa;
hermitestruct.Pxymoms = Pxymoms*1e9; % in nPa;
hermitestruct.Txxmoms = Txxmoms; % in eV
hermitestruct.Tyymoms = Tyymoms; % in eV
hermitestruct.Txymoms = Txymoms; % in eV
hermitestruct.vxth = vxth/1e3; % in km s^{-1};
hermitestruct.vyth = vyth/1e3; % in km s^{-1};
hermitestruct.vxmat = vxmat/1e3; % in km s^{-1};
hermitestruct.vymat = vymat/1e3; % in km s^{-1};
hermitestruct.species = species;

toc

end