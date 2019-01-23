function particlemoments = psd_moments(varargin)
% PSD_MOMENTS compute moments from the FPI particle phase-space densities 
%
% For brst mode data
% particlemoments = mms.psd_moments(pdist,phi,theta,stepTable,energy0,energy1,SCpot,particle,option,option_value)
%
% For fast mode data
% particlemoments = mms.psd_moments(pdist,phi,theta,energy,SCpot,particle,'fast',option,option_value)
%
% Input:
%   pdist - TSeries of the full particle distribution of electrons or ions
%   (must be in s^3/cm^6) (burst and fast)
%   phi - TSeries of all phi angles of distribution for burst data. 1D array or
%   structure for fast data.
%   theta - 1D array or structure of theta angles (burst and fast)
%   stepTable - TSeries of stepping table between energies (burst)
%   energy0 - 1D array or structure of energy table 0 (burst)
%   energy1 - 1D array or structure of energy table 1 (burst)
%   energy - 1D array or structure of energy table (fast)
%   SCpot - TSeries of spacecraft potential (burst and fast). 
%   (Make sure sign is correct, should be typically positive)
%   particle - indicate particle type: 'electron' or 'ion'
%
%   See Example_MMS_EDRsignatures for example of loading the necessary data 
%   and running the function.
%
% Optional Inputs:
%   'energyrange' - set energy range in eV to integrate over [E_min E_max].
%   energy range is applied to energy0 and the same elements are used for energy1 to 
%   ensure that the same number of points are integrated over. 
%   'noscpot' - set to 1 to set spacecraft potential to zero. Calculates moments without
%   correcting for spacecraft potential. 
%   'enchannels' - set energy channels to integrate over [min max]; min and max
%   between must be between 1 and 32.
%   'partialmoms' - use a binary array (or TSeries) (pmomsarr) to select which psd points are used
%   in the moments calculation. pmomsarr must be a binary array (1s and 0s, 1s correspond to points used).
%   Array (or data of TSeries) must be the same size as pdist.data. For
%   examples see Example_MMS_partialmoments.
%   'energy_minus' / 'energy_plus' - energy range: [energy - energy_minus,
%   energy - energy_plus]; 
%   'innerelec' ('on') - innerelectron potential for electron moments
% 
%  how to prepare data [example for electron]:
%       c_eval('eenergy = mms.db_get_variable(''mms?_fpi_brst_l2_des-dist'', ''mms?_des_energy_brst'', Tint);',ic);    
%       energy_minus = squeeze(eenergy.DELTA_MINUS_VAR.data(index, :));
%       energy_plus = squeeze(eenergy.DELTA_PLUS_VAR.data(index, :));
%
% Output: 
%   psd_moments - structure containing the particle moments: density, bulk
%   velocity, pressure, temperature, particle heat flux, density contribution
%   from 32 energy bins, and density contribution skymap (n_psd, V_psd, P_psd, 
%   T_psd, and H_psd, n_psd_e32, n_psd_skymap, respectively) as TSeries'.
%
% Notes: 
% Regarding the spacecraft potential, the best estimate of is -1.2*(probe
% to spacecraft voltage)+MMSoffset. Note that in most plasmas the spacecraft
% potential is positive. E.g.
% ic = 1,2,3, or 4;
% c_eval('do = dataobj(''data/mms?_edp_brst_l2_scpot_20151202011414_v1.0.0.cdf'');',ic);
% c_eval('SCpot = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_psp''));',ic);
% offset1 = 1.3; offset2 = 1.5; offset3 = 1.2; offset4 = 0.0; %For v1 data
% c_eval('SCpot.data = -SCpot.data*1.2+offset?;',ic);
% Apply correction for input. Correction is not applied in this script. 
% This correction is applied to v2 spacecraft potential so use 
% c_eval('SCpot = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_scpot_fast_l2''));',ic);
%
% For Phase 1B electron density:
%   tmpNe = irf.ts_scalar(epsd2mom.n_psd_e32.time,
%   irf.nansum(epsd2mom.n_psd_e32.data(:, 2:32), 2)); close to the official data;           % 2017-05-08, irfu.
%   or use 'innerelec' to remove some of the inner photoelectrons;
%
% Currently the heat flux vector does not match with the FPI ion moments. Currently
% using Eq. (6.8) of Analysis Methods for Multi-Spacecraft Data. This needs
% to be investigated further. 
%

% 1. basic
    flag_eminus = 0;
    flag_eplus = 0;
    flag_same_e = 0;
    flag_innerelec = 0;
    W_innerelec = 3.5;         % [eV] scpot + W_innerelec for electron moments calculation; 2018-01-26, wy;
    
% First input check
if (nargin < 7)
    nargin
    help psd_moments;
    return;
end

% Check if data is fast or burst resolution
if isa(varargin{7},'TSeries')
    isbrstdata = 1;
    irf.log('notice','Burst resolution data is used.');
elseif (varargin{7}(1) == 'f')
    isbrstdata = 0;
    irf.log('notice','Fast resolution data is used.');
else
    irf.log('critical','Something is wrong with the input.');
    nargin
    help psd_moments;
    return;
end

% Second input check
if (nargin < 8 && isbrstdata)
    nargin
    help psd_moments;
    return;
end

pdist = varargin{1};
phi = varargin{2};
thetak = varargin{3};
if isstruct(thetak)
    thetak = thetak.data;
end

if isbrstdata
    stepTable = varargin{4};
    energy0 = varargin{5};
    energy1 = varargin{6};
    SCpot = varargin{7};
    particletype = varargin{8};
    args=varargin(9:end);
    if isstruct(energy0)
        energy0 = energy0.data;
    end
    if isstruct(energy1)
        energy1 = energy1.data;
    end
    etmp = energy1 - energy0;
    if (all(etmp) == 0); flag_same_e = 1; end
else
    energy = varargin{4};
    SCpot = varargin{5};
    particletype = varargin{6};
    args=varargin(8:end);
    if isstruct(energy)
        energy = energy.data;
    end
end

%resample SCpot to same resolution as particle distributions
SCpot = SCpot.resample(pdist);

%check theta dimensions
thetasize = size(thetak);
if thetasize(1) > thetasize(2)
    thetak = thetak';
end

if isbrstdata
    args=varargin(9:end);
else
    args=varargin(8:end);
end
    
if numel(args)>0
    options=1;
else
    options=0;
end

intenergies = 1:32;
while options
    l = 2;
    switch(lower(args{1}))
        case 'energyrange'
            if numel(args)>1 && isnumeric(args{2})
                if isbrstdata==0
                    energy0 = energy;
                end
                Eminmax = args{2};
                starte = find(energy0 > Eminmax(1));
                starte = starte(1);
                ende = find(energy0 < Eminmax(2));
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
        case 'energy_plus'
            if numel(args)>1
                energy_plus = args{2};
                flag_eplus = 1;      
            end         
        case 'energy_minus'
            if numel(args)>1
                energy_minus = args{2};
                
                flag_eminus = 1;        
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

pdist.data = pdist.data*1e12; % convert to SI units

% Define arrays for output
n_psd = zeros(length(pdist.time), 1);
    sizedist = size(pdist.data);
    n_psd_e32 = zeros(length(pdist.time), 32);
    n_psd_e32_phi_theta = zeros(sizedist(1), sizedist(2), sizedist(3), sizedist(4));    
V_psd = zeros(length(pdist.time), 3);
P_psd = zeros(length(pdist.time), 3, 3);
P2_psd = zeros(length(pdist.time), 3, 3);
H_psd = zeros(length(pdist.time), 3);

tic

% angle between theta and phi points is 360/32 = 11.25 degrees
deltaang = (11.25*pi/180)^2;

if isbrstdata
    phitr = phi.data';
else
    phitr = phi;
    if isstruct(phitr)
        phitr = phitr.data;
    end
    phisize = size(phitr);
    if phisize(2) > phisize(1)
        phitr = phitr';
    end
end

% Calculate speed widths associated with each energy channel. 
% Lowest energy channels extended to reproduce fpi moments
if isbrstdata % Burst mode energy/speed widths
    if flag_same_e && flag_eminus && flag_eplus
        energy = energy0;
        energyupper = energy + energy_plus;
        energylower = energy - energy_minus;
        vupper = sqrt(2*qe*energyupper/pmass);
        vlower = sqrt(2*qe*energylower/pmass);
        deltav = (vupper-vlower);
    else
    energyall = [energy0 energy1];
    energyall = log10(sort(energyall));
    if abs(energyall(2) - energyall(1)) > 0.0001
        temp0 = 2*energyall(1)-energyall(2);
    else
        temp0 = 2*energyall(2)-energyall(3); 
    end
    if abs(energyall(64) - energyall(63)) > 0.0001
        temp65 = 2*energyall(64)-energyall(63);
    else
        temp65 = 2*energyall(64)-energyall(62); 
    end    
    energyall = [temp0 energyall temp65];
    diffenall = diff(energyall);
    energy0upper = 10.^(log10(energy0)+diffenall(2:2:64)/2);
    energy0lower = 10.^(log10(energy0)-diffenall(1:2:63)/2);
    energy1upper = 10.^(log10(energy1)+diffenall(3:2:65)/2);
    energy1lower = 10.^(log10(energy1)-diffenall(2:2:64)/2);
    v0upper = sqrt(2*qe*energy0upper/pmass);
    v0lower = sqrt(2*qe*energy0lower/pmass);
    v1upper = sqrt(2*qe*energy1upper/pmass);
    v1lower = sqrt(2*qe*energy1lower/pmass);
    deltav0 = (v0upper-v0lower)*2.0; %factor of two is applied because half the energy channels are used in a single sweep
    deltav1 = (v1upper-v1lower)*2.0;
    %deltav0(1) = deltav0(1)*2.7;
    %deltav1(1) = deltav1(1)*2.7;
    end
else % Fast mode energy/speed widths
    energyall = log10((energy));
    temp0 = 2*energyall(1)-energyall(2);
    temp33 = 2*energyall(32)-energyall(31); 
    energyall = [temp0 energyall temp33];
    diffenall = diff(energyall);    
    energyupper = 10.^(log10(energy)+diffenall(2:1:33)/4);
    energylower = 10.^(log10(energy)-diffenall(1:1:32)/4);
    vupper = sqrt(2*qe*energyupper/pmass);
    vlower = sqrt(2*qe*energylower/pmass);
    deltav = (vupper-vlower)*2.0;
    deltav(1) = deltav(1)*2.7;
end
    
for nt = 1:length(pdist.time)
    if isbrstdata
        if flag_same_e && flag_eminus && flag_eplus
        else
            energy = energy0;
            deltav = deltav0;
            if stepTable.data(nt) 
                energy = energy1;
                deltav = deltav1;
            end            
        end
    end
        
    v = real(sqrt(2*qe*(energy-SCpot.data(nt))/pmass));    
    v(energy-SCpot.data(nt) - flag_innerelec * W_innerelec < 0) = 0;

    if isbrstdata
        phij = phitr(:,nt);
    else
        phij = phitr;
    end
    
    Mpsd2n = ones(length(phij),1) * sind(thetak);
    Mpsd2Vx = -cosd(phij) * (sind(thetak) .* sind(thetak));
    Mpsd2Vy = -sind(phij) * (sind(thetak) .* sind(thetak));
    Mpsd2Vz = -ones(length(phij),1) * (sind(thetak) .* cosd(thetak));   
    Mpsdmfxx = cosd(phij).^2*(sind(thetak).^3);
    Mpsdmfyy = sind(phij).^2*(sind(thetak).^3);
    Mpsdmfzz = ones(length(phij),1) * (sind(thetak).* cosd(thetak).^2);
    Mpsdmfxy = cosd(phij).*sind(phij) * (sind(thetak).^3);
    Mpsdmfxz = cosd(phij) * (sind(thetak).^2.*cosd(thetak));
    Mpsdmfyz = sind(phij) * (sind(thetak).^2.*cosd(thetak));
   
    for ii = intenergies
        tmp = squeeze(pdist.data(nt, ii, :, :));
        n_psd_tmp1 = tmp .* Mpsd2n * v(ii)^2 * deltav(ii) * deltaang;
        n_psd_e32_phi_theta(nt, ii, :, :) = n_psd_tmp1;        
        n_psd_tmp = irf.nansum(irf.nansum(tmp .* Mpsd2n, 1), 2) * v(ii)^2 * deltav(ii) * deltaang;
        n_psd_e32(nt, ii) = n_psd_tmp;
        n_psd(nt) = n_psd(nt) + n_psd_tmp;
        %n_psd(nt) = n_psd(nt) + irf.nansum(irf.nansum(tmp .* Mpsd2n, 1), 2) * v(ii)^2 * deltav(ii) * deltaang;
        Vxtemp = irf.nansum(irf.nansum(tmp .* Mpsd2Vx, 1), 2) * v(ii)^3 * deltav(ii) * deltaang;
        Vytemp = irf.nansum(irf.nansum(tmp .* Mpsd2Vy, 1), 2) * v(ii)^3 * deltav(ii) * deltaang;
        Vztemp = irf.nansum(irf.nansum(tmp .* Mpsd2Vz, 1), 2) * v(ii)^3 * deltav(ii) * deltaang;
        V_psd(nt, 1) = V_psd(nt, 1) + Vxtemp;
        V_psd(nt, 2) = V_psd(nt, 2) + Vytemp;
        V_psd(nt, 3) = V_psd(nt, 3) + Vztemp;
        P_psd(nt, 1, 1) = P_psd(nt, 1, 1) + irf.nansum(irf.nansum(tmp .* Mpsdmfxx, 1), 2) * v(ii)^4 * deltav(ii) * deltaang;
        P_psd(nt, 1, 2) = P_psd(nt, 1, 2) + irf.nansum(irf.nansum(tmp .* Mpsdmfxy, 1), 2) * v(ii)^4 * deltav(ii) * deltaang;
        P_psd(nt, 1, 3) = P_psd(nt, 1, 3) + irf.nansum(irf.nansum(tmp .* Mpsdmfxz, 1), 2) * v(ii)^4 * deltav(ii) * deltaang;
        P_psd(nt, 2, 2) = P_psd(nt, 2, 2) + irf.nansum(irf.nansum(tmp .* Mpsdmfyy, 1), 2) * v(ii)^4 * deltav(ii) * deltaang;
        P_psd(nt, 2, 3) = P_psd(nt, 2, 3) + irf.nansum(irf.nansum(tmp .* Mpsdmfyz, 1), 2) * v(ii)^4 * deltav(ii) * deltaang;
        P_psd(nt, 3, 3) = P_psd(nt, 3, 3) + irf.nansum(irf.nansum(tmp .* Mpsdmfzz, 1), 2) * v(ii)^4 * deltav(ii) * deltaang;
        H_psd(nt, 1) = H_psd(nt, 1) + Vxtemp * v(ii)^2;
        H_psd(nt, 2) = H_psd(nt, 2) + Vytemp * v(ii)^2;
        H_psd(nt, 3) = H_psd(nt, 3) + Vztemp * v(ii)^2;
    end
end
toc
   
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
n_psd_e32 = n_psd_e32/1e6;
n_psd_e32_phi_theta = n_psd_e32_phi_theta / 1e6;    
V_psd = V_psd/1e3;
P_psd = P_psd*1e9;
P2_psd = P2_psd*1e9;
T_psd = T_psd*kb/qe;
H_psd = H_psd*1e3;

% Construct TSeries'
n_psd = irf.ts_scalar(pdist.time,n_psd);
n_psd_e32 = irf.ts_scalar(pdist.time, n_psd_e32);
if isstruct(phi)
    n_psd_skymap = PDist(pdist.time, n_psd_e32_phi_theta, 'skymap', energy, phi.data, thetak);
else
    n_psd_skymap = PDist(pdist.time, n_psd_e32_phi_theta, 'skymap', energy, phi, thetak);    
end
    n_psd_skymap.userData = pdist.userData; 
n_psd_skymap.units = 'cm^{-3}';  
V_psd = irf.ts_vec_xyz(pdist.time,V_psd);
P_psd = irf.ts_tensor_xyz(pdist.time,P_psd);
P2_psd = irf.ts_tensor_xyz(pdist.time,P2_psd);
T_psd = irf.ts_tensor_xyz(pdist.time,T_psd);
H_psd = irf.ts_vec_xyz(pdist.time,H_psd);

% make structure for output
particlemoments =struct('n_psd',n_psd,'V_psd',V_psd,'P_psd',P_psd,'P2_psd',P2_psd, ...
    'T_psd',T_psd,'H_psd',H_psd, 'n_psd_e32', n_psd_e32, 'n_psd_skymap', n_psd_skymap);

end
