function Distpart = get_hiresdistributions(varargin)
% GET_PARTIALDISTRIBUTIONS Obtains particle distributions at 2 or 4 times
% the standard resolution of FPI. 
% Written by D. B. Graham
%
% Distpart = mms.get_partialdistributions(Dist,4);
%
% Input:
%   Dist - Burst mode distribution function. Must be in PDist format and
%   include usec_offsets in ancillary data.
%   pf - Optional parameter. Determines how many distributions are
%   constructed out of the original distributions. Must be 2 or 4. If not included 4 is
%   assumed. 
%
% Output: 
%   Distpart - Higher temporal resolution distribution function. Parts of
%   distribution not measured are given by NaN. Same energy and angular
%   format as the original distribution.
%
% See also: mms.hiresmoments

pf = 4;

if (nargin < 1)
    nargin
    help get_partialdistributions;
    return;
end

Dist = varargin{1};

if isa(Dist,'PDist')
  if isempty('Dist.ancillary.usec_offsets')
    Distpart = NaN;
    help get_partialdistributions;
  return;
  end
else
  irf.log('critical','Distribution is not PDist format');
  Distpart = NaN;
  return;
end

if nargin > 1
  pf = varargin{2};
  if (pf ~= 2) && (pf ~= 4)
    pf = 4;
    irf.log('warning','pf value is not valid. pf set to 4.');
  end
end

% Define new timeline and new metadata
dt = Dist.ancillary.dt_minus+Dist.ancillary.dt_plus;
dtpart = dt/pf;
halfpartdt = dtpart/2;
OriginalTS = Dist.time+-Dist.ancillary.dt_minus;

PartialTS = int64(zeros(length(OriginalTS)*pf,1));
energyp = single(zeros(length(OriginalTS)*pf,length(Dist.depend{1,1}(1,:))));
phip = single(zeros(length(OriginalTS)*pf,length(Dist.depend{1,2}(1,:))));
Distp = double(ones(length(OriginalTS)*pf,length(Dist.depend{1,1}(1,:)),length(Dist.depend{1,2}(1,:)),length(Dist.depend{1,3})))*NaN;
esteptablep = uint8(zeros(length(OriginalTS)*pf,1));
delta_energy_minusp = single(zeros(length(OriginalTS)*pf,length(Dist.depend{1,1}(1,:))));
delta_energy_plusp = single(zeros(length(OriginalTS)*pf,length(Dist.depend{1,1}(1,:))));
dt_minusp = halfpartdt;
dt_plusp = halfpartdt;

for ii = 1:length(OriginalTS)
  PartialTS(pf*(ii-1)+1:pf*(ii-1)+pf) = OriginalTS(ii).ttns + int64((-halfpartdt + [1:pf]*dtpart)*1e9);
  delta_energy_minusp(pf*(ii-1)+1:pf*(ii-1)+pf,:) = ones(pf,1)*Dist.ancillary.delta_energy_minus(ii,:);
  delta_energy_plusp(pf*(ii-1)+1:pf*(ii-1)+pf,:) = ones(pf,1)*Dist.ancillary.delta_energy_plus(ii,:);
  energyp(pf*(ii-1)+1:pf*(ii-1)+pf,:) = ones(pf,1)*Dist.depend{1,1}(ii,:);
  esteptablep(pf*(ii-1)+1:pf*(ii-1)+pf) = Dist.ancillary.esteptable(ii);
  phip(pf*(ii-1)+1:pf*(ii-1)+pf,:) = ones(pf,1)*Dist.depend{1,2}(ii,:);
end

if pf == 2
  for ii = 1:length(OriginalTS)
    Pdatatemp = squeeze(Dist.data(ii,:,:,:));
    Toffsetstemp = squeeze(Dist.ancillary.usec_offsets(ii,:,:,:));
    idx1 = Toffsetstemp < dt/2*1e6;
    idx2 = Toffsetstemp >= dt/2*1e6;
    Distptemp1 = ones(size(Pdatatemp))*NaN;
    Distptemp2 = ones(size(Pdatatemp))*NaN;
    Distptemp1(idx1) = Pdatatemp(idx1);
    Distptemp2(idx2) = Pdatatemp(idx2);
    Distp(pf*(ii-1)+1,:,:,:) = Distptemp1;
    Distp(pf*(ii-1)+2,:,:,:) = Distptemp2;
  end
end

if pf == 4
  for ii = 1:length(OriginalTS)
    Pdatatemp = squeeze(Dist.data(ii,:,:,:));
    Toffsetstemp = squeeze(Dist.ancillary.usec_offsets(ii,:,:,:));
    idx1 = Toffsetstemp < dt/4*1e6;
    idx2 = Toffsetstemp < dt/2*1e6 & Toffsetstemp >= dt/4*1e6;
    idx3 = Toffsetstemp < 3*dt/4*1e6 & Toffsetstemp >= dt/2*1e6;
    idx4 = Toffsetstemp >= 3*dt/4*1e6;
    Distptemp1 = ones(size(Pdatatemp))*NaN;
    Distptemp2 = ones(size(Pdatatemp))*NaN;
    Distptemp3 = ones(size(Pdatatemp))*NaN;
    Distptemp4 = ones(size(Pdatatemp))*NaN;
    Distptemp1(idx1) = Pdatatemp(idx1);
    Distptemp2(idx2) = Pdatatemp(idx2);
    Distptemp3(idx3) = Pdatatemp(idx3);
    Distptemp4(idx4) = Pdatatemp(idx4);
    Distp(pf*(ii-1)+1,:,:,:) = Distptemp1;
    Distp(pf*(ii-1)+2,:,:,:) = Distptemp2;
    Distp(pf*(ii-1)+3,:,:,:) = Distptemp3;
    Distp(pf*(ii-1)+4,:,:,:) = Distptemp4;
  end
end

Distpart = PDist(EpochTT(PartialTS),Distp,'skymap',energyp,phip,Dist.depend{1,3});
Distpart.userData = Dist.userData;
Distpart.name = Dist.name;
Distpart.units = Dist.units;
Distpart.siConversion = Dist.siConversion;
Distpart.species = Dist.species;
Distpart.ancillary.dt_minus = dt_minusp;
Distpart.ancillary.dt_plus = dt_plusp;
Distpart.ancillary.energy = energyp;
Distpart.ancillary.energy0 = Dist.ancillary.energy0;
Distpart.ancillary.energy1 = Dist.ancillary.energy1;
Distpart.ancillary.esteptable = esteptablep;
Distpart.ancillary.delta_energy_minus = delta_energy_minusp;
Distpart.ancillary.delta_energy_plus = delta_energy_plusp;
Distpart.ancillary.usec_offsets = [];
Distpart.ancillary.pf = pf;

end