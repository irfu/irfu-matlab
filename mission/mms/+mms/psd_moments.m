function particlemoments = psd_moments(varargin)
% PSD_MOMENTS compute moments from the FPI particle phase-space densities 
%
% particlemoments = mms.psd_moments(pdist,phi,theta,stepTable,energy0,energy1,SCpot,particle,option,option_value)
%
% Input:
%   pdist - TSeries of the full particle distribution of electrons or ions (must be in s^3/cm^6)
%   phi - TSeries of all phi angles of distribution
%   theta - 1D array or structure of theta angles
%   stepTable - TSeries of stepping table between energies
%   energy0 - 1D array or structure of energy table 0
%   energy1 - 1D array or structure of energy table 1
%   SCpot - TSeries of spacecraft potential. (Make sure sign is correct, should be typically positive)
%   particle - indicate particle type: 'electron' or 'ion'
% Optional Inputs:
%   'energyrange' - set energy range in eV to integrate over [E_min E_max].
%   Not yet implemented because of energy channel switching. 
%   'noscpot' - set to 1 to set spacecraft potential to zero. Calculates moments without
%   correcting for spacecraft potential. 
%   'enchannels' - set energy channels to integrate over [min max]; min and max
%   between must be between 1 and 32.
%
% Output: 
%   psd_moments - structure containing the particle moments: density, bulk
%   velocity, pressure and temperature (n_psd, V_psd, P_psd, T_psd, and H_psd,
%   respectively) as TSeries'. For temperature and
%   pressure tensors the order of the columns is XX, XY, XZ, YY, YZ, ZZ.
%
% Still testing but seems to work okej.
%
% Notes: 
% Regarding the spacecraft potential, the best estimate of is -1.2*(probe
% to spacecraft voltage)+2. Note that in most plasmas the spacecraft
% potential is positive. E.g.
% do = dataobj('data/mms1_edp_brst_l2_scpot_20151202011414_v1.0.0.cdf');
% SCpot = mms.variable2ts(get_variable(tmpDataObj,'mms?_edp_psp'));
% SCpot.data = -SCpot.data*1.2+2; 
% Apply correction for input. Correction is not applied in this script. 
%
% Currently the heat flux vector does not match with the FPI ion moments. Currently
% using Eq. (6.8) of Analysis Methods for Multi-Spacecraft Data. This needs
% to be investigated further. 
%

% Check input
if nargin < 8,
    nargin
    help psd_moments;
    return;
end

pdist = varargin{1};
phi = varargin{2};
thetak = varargin{3};
stepTable = varargin{4};
energy0 = varargin{5};
energy1 = varargin{6};
SCpot = varargin{7};
particletype = varargin{8};

%resample SCpot to same resolution as particle distributions
SCpot = SCpot.resample(pdist);

if isstruct(thetak),
    thetak = thetak.data;
end
if isstruct(energy0),
    energy0 = energy0.data;
end
if isstruct(energy1),
    energy1 = energy1.data;
end

%check theta dimensions
thetasize = size(thetak);
if thetasize(1) > thetasize(2)
    thetak = thetak';
end

args=varargin(9:end);

if numel(args)>0,
    options=1;
else
    options=0;
end

intenergies = [1:32];
while options
    l = 2;
    switch(lower(args{1}))
        %case 'energyrange'
            %if numel(args)>1 && isnumeric(args{2})
            %    Eminmax = args{2};
            %    irf_log('proc','Using partial energy range')
            %end
      	case 'noscpot'
            if numel(args)>1 && args{2};
                SCpot.data = zeros(size(SCpot.data));
                irf_log('proc','Setting spacecraft potential to zero.')
            end
        case 'enchannels'
            if numel(args)>1 && isnumeric(args{2});
                intenergies = args{2}(1):args{2}(2);
            end
        otherwise
            irf_log('fcal',['Unknown flag: ' args{1}])
            l=1;
            break
    end
    args = args(l+1:end);
    if isempty(args), options=0; end
end

irf_log('proc',strcat('Integrating over energy levels ',num2str(min(intenergies)),' to ',num2str(max(intenergies))));

% Define constants
Units = irf_units; % Use IAU and CODATA values for fundamental constants.
qe = Units.e;
kb = Units.kB;

if (particletype(1) == 'e')
    pmass = Units.me;
    irf_log('proc','Particles are electrons')
elseif (particletype(1) == 'i')
    pmass = Units.mp;
    SCpot.data = -SCpot.data;
    irf_log('proc','Particles are Ions')
else
    particlemoments = NaN;
    irf_log('proc','Could not identify the particle')
    return;
end

pdist.data = pdist.data*1e12; % convert to SI units

% Define arrays for output
n_psd = zeros(length(pdist.time), 1);
V_psd = zeros(length(pdist.time), 3);
P_psd = zeros(length(pdist.time), 6);
T_psd = zeros(length(pdist.time), 6);
H_psd = zeros(length(pdist.time), 3);

tic

% angle between theta and phi points is 360/32 = 11.25
deltaang = (11.25*pi/180)^2;

phitr = phi.data';

for nt = 1:length(pdist.time)
	energy = energy0;
	if stepTable.data(nt), 
        energy = energy1;
    end 
    
    v = real(sqrt(2*qe*(energy-SCpot.data(nt))/pmass));    
    v(energy-SCpot.data(nt)<0) = 0;
        
    
    % estimate of dv - incorperate upper and lower energies when available
    energy = [2*energy(1)-energy(2) energy 2*energy(32)-energy(31)];
    v2 = sqrt(2*qe*energy/pmass);
    
    phij = phitr(:,nt);
    Mpsd2n = ones(length(phij),1) * sind(thetak);
    Mpsd2Vx = -cosd(phij) * (sind(thetak) .* sind(thetak));
    Mpsd2Vy = -sind(phij) * (sind(thetak) .* sind(thetak));
    Mpsd2Vz = -ones(length(phij),1) * (sind(thetak) .* cosd(thetak));   
    Mpsdmfxx = cosd(phij).^2*(sind(thetak).^3);
    Mpsdmfyy = sind(phij).^2*(sind(thetak).^3);
    Mpsdmfzz = ones(length(phij),1) * (sind(thetak).* cosd(thetak).^2);
    Mpsdmfxy = cosd(phij).*sind(phij) * (sind(thetak).^3);
    Mpsdmfxz = cosd(phij) * (sind(thetak).^2.*cosd(thetak));
    Mpsdmfyz= sind(phij) * (sind(thetak).^2.*cosd(thetak));
   
    for ii = intenergies; 
        deltav = (v2(ii+2) - v2(ii))/2;
        tmp = squeeze(pdist.data(nt, ii, :, :));
        n_psd(nt) = n_psd(nt) + irf.nansum(irf.nansum(tmp .* Mpsd2n, 1), 2) * v(ii)^2 * deltav * deltaang;
        Vxtemp = irf.nansum(irf.nansum(tmp .* Mpsd2Vx, 1), 2) * v(ii)^3 * deltav * deltaang;
        Vytemp = irf.nansum(irf.nansum(tmp .* Mpsd2Vy, 1), 2) * v(ii)^3 * deltav * deltaang;
        Vztemp = irf.nansum(irf.nansum(tmp .* Mpsd2Vz, 1), 2) * v(ii)^3 * deltav * deltaang;
        V_psd(nt, 1) = V_psd(nt, 1) + Vxtemp;
        V_psd(nt, 2) = V_psd(nt, 2) + Vytemp;
        V_psd(nt, 3) = V_psd(nt, 3) + Vztemp;
        P_psd(nt, 1) = P_psd(nt, 1) + irf.nansum(irf.nansum(tmp .* Mpsdmfxx, 1), 2) * v(ii)^4 * deltav * deltaang;
        P_psd(nt, 2) = P_psd(nt, 2) + irf.nansum(irf.nansum(tmp .* Mpsdmfxy, 1), 2) * v(ii)^4 * deltav * deltaang;
        P_psd(nt, 3) = P_psd(nt, 3) + irf.nansum(irf.nansum(tmp .* Mpsdmfxz, 1), 2) * v(ii)^4 * deltav * deltaang;
        P_psd(nt, 4) = P_psd(nt, 4) + irf.nansum(irf.nansum(tmp .* Mpsdmfyy, 1), 2) * v(ii)^4 * deltav * deltaang;
        P_psd(nt, 5) = P_psd(nt, 5) + irf.nansum(irf.nansum(tmp .* Mpsdmfyz, 1), 2) * v(ii)^4 * deltav * deltaang;
        P_psd(nt, 6) = P_psd(nt, 6) + irf.nansum(irf.nansum(tmp .* Mpsdmfzz, 1), 2) * v(ii)^4 * deltav * deltaang;
        H_psd(nt, 1) = H_psd(nt, 1) + Vxtemp * v(ii)^2;
        H_psd(nt, 2) = H_psd(nt, 2) + Vytemp * v(ii)^2;
        H_psd(nt, 3) = H_psd(nt, 3) + Vztemp * v(ii)^2;
    end
end
toc

    
% Compute moments in SI units    
P_psd = pmass*P_psd;
V_psd = V_psd./[n_psd n_psd n_psd];
P_psd(:,1) = P_psd(:,1)-pmass*n_psd.*V_psd(:,1).*V_psd(:,1);
P_psd(:,2) = P_psd(:,2)-pmass*n_psd.*V_psd(:,1).*V_psd(:,2);
P_psd(:,3) = P_psd(:,3)-pmass*n_psd.*V_psd(:,1).*V_psd(:,3);
P_psd(:,4) = P_psd(:,4)-pmass*n_psd.*V_psd(:,2).*V_psd(:,2);
P_psd(:,5) = P_psd(:,5)-pmass*n_psd.*V_psd(:,2).*V_psd(:,3);
P_psd(:,6) = P_psd(:,6)-pmass*n_psd.*V_psd(:,3).*V_psd(:,3);
Ptrace = (P_psd(:,1)+P_psd(:,4)+P_psd(:,6));
T_psd = P_psd./([n_psd n_psd n_psd n_psd n_psd n_psd])/kb;
H_psd = pmass/2*H_psd;
H_psd(:,1) = H_psd(:,1)-(V_psd(:,1).*P_psd(:,1)+V_psd(:,2).*P_psd(:,2)+V_psd(:,3).*P_psd(:,3))-0.5*V_psd(:,1).*Ptrace;
H_psd(:,2) = H_psd(:,2)-(V_psd(:,1).*P_psd(:,2)+V_psd(:,2).*P_psd(:,4)+V_psd(:,3).*P_psd(:,5))-0.5*V_psd(:,2).*Ptrace;
H_psd(:,3) = H_psd(:,3)-(V_psd(:,1).*P_psd(:,3)+V_psd(:,2).*P_psd(:,5)+V_psd(:,3).*P_psd(:,6))-0.5*V_psd(:,3).*Ptrace;


% Convert to typical units (/cc, km/s, nP, eV, and ergs/s/cm^2).
n_psd = n_psd/1e6;
V_psd = V_psd/1e3;
P_psd = P_psd*1e9;
T_psd = T_psd*kb/qe;
H_psd = H_psd*1e3;

% Construct TSeries'
n_psd = irf.ts_scalar(pdist.time,n_psd);
V_psd = irf.ts_vec_xyz(pdist.time,V_psd);
P_psd = TSeries(pdist.time,P_psd);
T_psd = TSeries(pdist.time,T_psd);
H_psd = irf.ts_vec_xyz(pdist.time,H_psd);

% make structure for output
particlemoments =struct('n_psd',n_psd,'V_psd',V_psd,'P_psd',P_psd,'T_psd',T_psd,'H_psd',H_psd);

end