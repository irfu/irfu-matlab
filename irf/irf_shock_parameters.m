function dspec = irf_shock_parameters(spec)
% IRF_SHOCK_PARAMETERS Calculate shock related plasma parameters.
%
% dspec = IRF_SHOCK_PARAMETERS(spec) Returns struct dspec with derived plasma
% parameters from input struct spec with measured plasma parameters.
%
% Input struct has parameters input with fixed names. After the parameter
% name, a name of the region can be given, e.g. "u" and "d". All parameters
% except B are optional.
%
% INPUT PARAMETERS:
% Magnetic field vector [nT]        -   B
% Plasma velocity [km/s]            -   V
% Plasma number density [cm^-3]     -   n
% Ion/electron temperature [eV]     -   Ti/Te
% Reference system for Mach numbers -   ref_sys
% Normal vector                     -   nvec
% Shock speed along nvec            -   Vsh
%
% OUTPUT PARAMETERS - symbol    - required input parameters
% Speeds [km/s]
% Alfven speed      - Va        - B,n
% fast speed        - Vf        - B,V,n,Ti,Te
% sound speed       - Vts       - Ti,Te
%
% Frame velocities [km/s]
% Normal incidence frame  - Vnif  - V
% deHoffmann-Teller frame - Vhtf  - B,V
%
% Frequencies [s^-1] (not radians)
% proton gyrofreq   - Fcp       - B
%
% Lengths [km]
% ion inert.len     - Li        - n
% proton gyroradius - Rcp       - B,V     (not thermal motion)
%
% Dimensionless
% Alfven Mach #     - Ma        - B,V,n
% fast Mach #       - Mf        - B,V,n,Ti,Te
% ion beta          - bi        - n,Ti
% electron beta     - be        - n,Te
%
% It is possible to calculate the Mach numbers in the Normal Incidence
% Frame (NIF). Set spec.ref_sys = 'NIF' and specify normal vector and
% optionally shock speed.
%
% To transform a velocity to e.g., the NI frame do the following
% transformation:
%   V_in_nif = V - Vnif
%
% ----------------------
% Example 1:
% sp = [];
% sp.B = [10,0,0];
% sp.n = 1;
% sp.Te = 100; sp.Ti=1000;
% dsp = irf_shock_parameters(sp)
%
% dsp =
%
%      Va: 218.1204
%     Vts: 544.9255
%     Fcp: 0.1525
%      Li: 227.7108
%      bi: 4.0267
%      be: 0.4027
% ----------------------
% Example 2:
% sp = [];
% sp.Bu = 5; sp.Bd = 20;
% sp.nu = 3; sp.nd = 12;
% dsp = irf_shock_parameters(sp)
%
% dsp =
%
%      Vau: 62.9659
%      Vad: 125.9318
%     Fcpd: 0.3049
%     Fcpd: 0.3049
%      Liu: 131.4689
%      Lid: 65.7344
% ----------------------
% Example 3:
% spec = [];
% spec.B = [0,5,0];
% spec.V = [-400,0,0];
% spec.n = 5;
% spec.Ti = 10;
% spec.Te = 10;
% spec.ref_sys = 'nif';
% spec.nvec = [1,1,0]/sqrt(2);
% dspec = irf_shock_parameters(spec)
%   [warning: irf_shock_parameters(119)] Setting shock speed, Vsh, to 0.
%
% dspec =
%
%   struct with fields:
%
%      Va: 48.7732
%     Vts: 61.8994
%      Vf: 78.8058
%     Fcp: 0.0762
%      Li: 101.8354
%     Rcp: 835.2014
%      Ma: 5.7991
%      Mf: 3.8633
%      Ms: 4.5694
%      bi: 0.8053
%      be: 0.8053
% ----------------------
%
% See also: IRF_PLASMA_CALC, IRF_SHOCK_NORMAL, IRF_SHOCK_GUI
%


%% Handle input
fn = fieldnames(spec);

% find Bs
idB = ~cellfun(@isempty,strfind(fn,'B')); %#ok<STRCLFH>
% regions
rgsB = fn(idB);
nR = numel(rgsB);
rgs = cell(1,nR);
for k = 1:nR
  rgs{k} = rgsB{k}(end);
  if strcmp(rgs{k},'B')
    rgs{k} = '';
  end
end

% reference system for Mach numbers
if ~ismember(fn,'ref_sys')
  spec.ref_sys = 'sc';
end

if strcmpi(spec.ref_sys,'nif')
  if ~ismember(fn,'nvec')
    irf.log('c','Calculation in the NIF requires a normal vector, nvec.');
  end
  if ~ismember(fn,'Vsh')% set shock velocity to 0 if not entered
    irf.log('w','Setting shock speed, Vsh, to 0.');
    spec.Vsh = 0;
  end
end

if find(ismember(fn,['V',rgs{1}])); hasV = 1; else, hasV = 0; end
if find(ismember(fn,['n',rgs{1}])); hasN = 1; else, hasN = 0; end
if find(ismember(fn,['Ti',rgs{1}])); hasTi = 1; else, hasTi = 0; end
if find(ismember(fn,['Te',rgs{1}])); hasTe = 1; else, hasTe = 0; end

% Frame velocities need to know wether there is a normal vector
if find(ismember(fn,'nvec')); hasNvec = 1; else, hasNvec = 0; end


%% Calculate parameters
dspec = [];


%% Speeds
% Alfven
if hasN
  for k = 1:nR
    dspec.(['Va',rgs{k}]) = v_alfv(spec.(['B',rgs{k}]),spec.(['n',rgs{k}]));
  end
end
% Sound speed
if hasTi && hasTe
  for k = 1:nR
    dspec.(['Vts',rgs{k}]) = v_sound(spec.(['Ti',rgs{k}]),spec.(['Te',rgs{k}]));
  end
end
% Fast
if hasN && hasTi && hasV
  for k = 1:nR
    dspec.(['Vf',rgs{k}]) = v_fast(spec.(['B',rgs{k}]),spec.(['V',rgs{k}]),spec.(['n',rgs{k}]),spec.(['Ti',rgs{k}]),spec.(['Te',rgs{k}]));
  end
end

%% Frame velocities

if hasV && hasNvec
  for k = 1:nR
    dspec.(['Vnif',rgs{k}]) = nif_speed(spec.(['V',rgs{k}]),spec);
  end
end

if hasV && hasNvec
  for k = 1:nR
    dspec.(['Vhtf',rgs{k}]) = htf_speed(spec.(['B',rgs{k}]),spec.(['V',rgs{k}]),spec);
  end
end


%% Frequencies
% Ion gyrofrequency
for k = 1:nR
  dspec.(['Fcp',rgs{k}]) = ion_gyro_freq(spec.(['B',rgs{k}]));
end


%% Lenghts
% Ion inertial length
if hasN
  for k = 1:nR
    dspec.(['Li',rgs{k}]) = ion_in_len(spec.(['n',rgs{k}]));
  end
end
% Ion gyroradius
if hasV
  for k = 1:nR
    dspec.(['Rcp',rgs{k}]) = ion_gyro_rad(spec.(['B',rgs{k}]),spec.(['V',rgs{k}]));
  end
end


%% Dimensionless
% Alfven Mach
if hasV && hasN
  for k = 1:nR
    dspec.(['Ma',rgs{k}]) = alfv_mach(spec.(['B',rgs{k}]),spec.(['V',rgs{k}]),spec.(['n',rgs{k}]),spec);
  end
end
% fast Mach
if hasV && hasN && hasTi && hasTe
  for k = 1:nR
    dspec.(['Mf',rgs{k}]) = fast_mach(spec.(['B',rgs{k}]),spec.(['V',rgs{k}]),spec.(['n',rgs{k}]),spec.(['Ti',rgs{k}]),spec.(['Te',rgs{k}]),spec);
  end
end
% sound Mach
if hasV && hasTi && hasTe
  for k = 1:nR
    dspec.(['Ms',rgs{k}]) = sonic_mach(spec.(['V',rgs{k}]),spec.(['Ti',rgs{k}]),spec.(['Te',rgs{k}]),spec);
  end
end
% ion beta
if hasN && hasTi
  for k = 1:nR
    dspec.(['bi',rgs{k}]) = beta_i(spec.(['B',rgs{k}]),spec.(['n',rgs{k}]),spec.(['Ti',rgs{k}]));
  end
end
% electron beta
if hasN && hasTe
  for k = 1:nR
    dspec.(['be',rgs{k}]) = beta_e(spec.(['B',rgs{k}]),spec.(['n',rgs{k}]),spec.(['Te',rgs{k}]));
  end
end


end


%% Help functions

function Va =  v_alfv(B,n)
Va = irf_plasma_calc(norm(B),norm(n),0,0,0,'Va')*1e-3;
end

function Vts = v_sound(Ti,Te) %
Vts = irf_plasma_calc(0,0,0,Te,Ti,'Vts')*1e-3;
end

function Vf =  v_fast(B,V,n,Ti,Te,th)
if nargin == 5
  th = acosd(dot(B,V,2)./(irf_abs(B,1)*irf_abs(V,1)));
end
Va = v_alfv(B,n);
cs = v_sound(Ti,Te);
cms0 = sqrt(Va.^2+cs.^2);

Vf = sqrt(cms0.^2/2+sqrt(cms0.^4/4-Va.^2.*cs.^2*cosd(th).^2));
end

function Fcp = ion_gyro_freq(B)
Fcp = irf_plasma_calc(norm(B),0,0,0,0,'Fcp');
end

function Li =  ion_in_len(n)
Li = irf_plasma_calc(0,n,0,0,0,'Li')*1e-3; % returns in m
end

function Rcp =  ion_gyro_rad(B,V)
u = irf_units;
Ep_kin = 1/2*u.mp*norm(V)^2/u.e; % kinetic energy of proton (eV)
Rcp = irf_plasma_calc(norm(B),0,0,0,Ep_kin,'Rop'); % returns in km
end

% frame velocities
function vNIF = nif_speed(V,spec)
vNIF = V-(dot(V,spec.nvec)-spec.Vsh)*spec.nvec;
end

function vHTF = htf_speed(B,V,spec) %
% first get the velocity in a shock rest frame
V_in_shock_rest_frame = V-spec.Vsh*spec.nvec;
% then get the dHT frame speed in the shock rest frame
vHTF_srf = cross(spec.nvec,cross(V_in_shock_rest_frame,B))/(dot(B,spec.nvec));
% then the dHT frame speed in the sc frame is the shock speed plus
% the dHT frame speed in the shock rest frame (I think)
vHTF = spec.Vsh*spec.nvec+vHTF_srf;
end

function Ma = alfv_mach(B,V,n,spec)
switch lower(spec.ref_sys)
  case 'sc'
    Ma = norm(V)/v_alfv(B,n);
  case 'nif'
    Ma = norm(dot(V,spec.nvec)-spec.Vsh)/v_alfv(B,n);
end
end

function Mf = fast_mach(B,V,n,Ti,Te,spec)
switch lower(spec.ref_sys)
  case 'sc'
    Mf = norm(V)/v_fast(B,V,n,Ti,Te);
  case 'nif'
    thBn = acosd(dot(B,spec.nvec)/norm(B));
    Mf = norm(dot(V,spec.nvec)-spec.Vsh)/v_fast(B,V,n,Ti,Te,thBn);
end
end

function Ms = sonic_mach(V,Ti,Te,spec)
switch lower(spec.ref_sys)
  case 'sc'
    Ms = norm(V)/v_sound(Ti,Te);
  case 'nif'
    Ms = norm(dot(V,spec.nvec)-spec.Vsh)/v_sound(Ti,Te);
end
end

function bi =  beta_i(B,n,Ti)
u = irf_units;
bi = n*1e6*Ti*u.e/((norm(B)*1e-9)^2/(2*u.mu0));
end

function be =  beta_e(B,n,Te)
u = irf_units;
be = n*1e6*Te*u.e/((norm(B)*1e-9)^2/(2*u.mu0));
end

