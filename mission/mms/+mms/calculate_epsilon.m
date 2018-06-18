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
%       between must be between 1 and 32.
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

PDist = varargin{1};
modelPDist = varargin{2};
n = varargin{3};
SCpot = varargin{4};
args=varargin(5:end);

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

% Default energy channels used to compute epsilon, lowest energy channel
% should not be used.
intenergies = 2:32; 

while options
    l = 2;
    switch(lower(args{1}))
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
if (PDist.species(1) == 'e')
    pmass = Units.me;
    irf.log('notice','Particles are electrons');
elseif (PDist.species(1) == 'i')
    pmass = Units.mp;
    SCpot.data = -SCpot.data;
    irf.log('notice','Particles are Ions');
else
    epsilon = NaN;
    irf.log('critical','Could not identify the particle type');
    return;
end

% Define lengths of variables
lengthphi = length(PDist.depend{1,2}(1,:));

% Define corrected energy levels using SCpot
energyarr = PDist.ancillary.energy;
v = zeros(size(energyarr));
deltav = zeros(size(energyarr));

for nt = 1:length(PDist.time)
    energyvec = energyarr(nt,:);
    v(nt,:) = real(sqrt(2*(energyvec-SCpot.data(nt))*qe/pmass));
    energylog = log10(energyvec);
    temp0 = 2*energylog(1)-energylog(2);
    temp33 = 2*energylog(end)-energylog(end-1);
    energyall = [temp0 energylog temp33];
    diffenall = diff(energyall);    
    energyupper = 10.^(energylog+diffenall(2:1:33)/2);
    energylower = 10.^(energylog-diffenall(1:1:32)/2);
    vupper = sqrt(2*qe*energyupper/pmass);
    vlower = sqrt(2*qe*energylower/pmass);
    vlower(isnan(vlower)) = 0;
    vupper(isnan(vupper)) = 0;
    deltav(nt,:) = (vupper-vlower);
end
v(v < 0) = 0;

% Calculate density of Pdiff
deltaang = (11.25*pi/180)^2;
thetak = PDist.depend{1,3};

epsilon = zeros(length(PDist.time), 1);

Mpsd2n = ones(lengthphi,1) * sind(thetak);

for nt = 1:length(PDist.time)
    for ii = intenergies 
        tmp = squeeze(Pdiff(nt, ii, :, :));
        epsilon(nt) = epsilon(nt) + irf.nansum(irf.nansum(tmp .* Mpsd2n, 1), 2) * v(nt,ii)^2 * deltav(nt,ii) * deltaang;
    end
end

epsilon = epsilon/1e6./(n.data*2);
epsilon = irf.ts_scalar(PDist.time,epsilon);

toc;

end
