function hermitestruct = irf_hermitedecomp_1d(varargin)
% IRF_HERMITEDECOMP_1D - decompose a series of 1D distributions in the
% hermite basis of functions
%
%   Detailed explanation goes here
%
% Input:
%       PDist1d - Structure or TSeries of PDist of 1D particle distributions.
% Must contain the distributions and speeds. TSeries should be 1D reduced
% distributions.
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
%   hermitestruct = irf_hermitedecomp_1d(PDist1D,64,'species','e');
%
% Written by D. B. Graham

if nargin < 2
  help irf_hermitedecomp_1d;
  return;
end

PDist1d = varargin{1};
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

if isa(PDist1d,'TSeries')
  if massnotset
    species = PDist1d.species(1);
  end
  distributions = PDist1d.data; % Should be units of s m^{-4}
  times = PDist1d.time;
  vmat = PDist1d.depend{1,1}*1e3; % Convert to m s^{-1}
  dv = median(diff(PDist1d.depend{1,1}(1,:)))*1e3; % Convert to m s^{-1}
  lengthv = length(vmat(1,:));
  isTSeries = true;
end

if isTSeries == false
  times = PDist1d.t;
  distributions = PDist1d.dist;
   if isfield(PDist1d,'species') && massnotset
     if ischar(PDist1d.species)
        species = PDist1d.species(1);
     end
   end

  timesformat = size(times);
  if timesformat(1) < timesformat(2)
    times = times';
  end

  vvec = PDist1d.v*1e3; % convert from km s^-1 to m s^-1
  lengthv = length(vvec);
  dv = median(diff(vvec));

  if min(size(vvec)) == 1
    vmat = ones(size(times))*vvec;
  end
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
nmoms = sum(distributions,2)*dv;
Vmoms = sum(vmat.*distributions,2)*dv;
Vmoms = Vmoms./nmoms;
Vmomsmat = repmat(Vmoms,1,lengthv);
Pmoms = ms*sum(distributions.*(vmat-Vmomsmat).^2*dv,2);
Tmoms = Pmoms./nmoms/qe;
vth = sqrt(qe*Tmoms/ms); % Use RMS v_th rather than mean v_th for Hermite decomposition
vthmat = repmat(vth,1,lengthv);
sqrtvthmat = sqrt(vthmat);

HermitePs = zeros([length(times) numhermites+1]);
vnormmat = (vmat-Vmomsmat)./vthmat;

for ii = 0:numhermites
  hermiteii = irf_hermitefunction(vnormmat,ii)./sqrtvthmat;
  HermitePs(:,ii+1) = sum(hermiteii.*distributions,2)*dv;
end

hermitestruct = struct;
hermitestruct.hermitespectra = HermitePs;
hermitestruct.orders = 0:numhermites;
hermitestruct.times = times;
hermitestruct.vmat = vmat/1e3; % in km s^{-1}
hermitestruct.nmoms = nmoms/1e6; % in cm^{-3}
hermitestruct.Vmoms = Vmoms/1e3; % in km s^{-1}
hermitestruct.Pmoms = Pmoms*1e9; % in nPa
hermitestruct.Tmoms = Tmoms; % in eV
hermitestruct.vth = vth/1e3; % in km s^{-1}
hermitestruct.species = species;

end
