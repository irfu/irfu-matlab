% IRF_UNITS create global variable Units with key units and constants
%
% uses units.m from web matlab_exchange depository 
% adds also space relevant stuff
% see help units.m for more details
%
% Examples: 
%   Units.R_Earth                               % display Earth radius in SI
%   Units.R_Earth/Units.R_Sun                   % display Earth radius in solar radia
%   T_in_eV = Units.kB*T_in_MK*1e6 / Units.e    % to convert from MK to eV
%   Units                                       % to see all defined units
%
% Suggestion: execute irf_units in your startup and you can everywhere
% access the units and constants

global Units
clear u;

u=units;

% added some more common units in space missing in units.m

%-------- UNITS ------------------------------
%------- length ----
u.AU = 1.496e11*u.m;
u.R_Earth = 6370e3*u.m;                     % Earth radius
u.R_Sun = 6.96e8*u.m;                       % Solar radius 
u.pc = 3.0857e16*u.m;                        % parsec

%------- mass ----
u.M_Earth = 5.9742e24*u.kg;                 % Mass of the Earth
u.M_Sun = 1.98892e30*u.kg;                  % Mass of the Sun
u.me = 9.1094e-31*u.kg;                      % electron mass
u.mp = 1.6726e-27*u.kg;                      % proton mass

%---- frequency ---- 
u.mHz = 1e-3*u.Hz;

%----- energy -----
u.keV = 1e3*u.eV;
u.MeV = 1e6*u.eV;
u.erg = 1e-7*u.J;

%---- temperature ---
u.MK = 1e6*u.K;

%----magnetic field -----
u.nT = 1e-9*u.T;

%----fundamental constants ----
u.G = 6.6726e-11*u.N*u.m^2/u.kg^2;           % gravitational constant
u.kB = 1.38e-23*u.J/u.K;                     % Boltzmann constant


Units=u; 
clear u;


