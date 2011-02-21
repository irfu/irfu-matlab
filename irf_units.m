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

%%%%%%%%%%%%%%%%%%%%%%%%% including original units.m code %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% from matlab_exchange %%%%%%%%%%%%%%%%%%%%


% The units.m function returns a struct containing the SI values of 
% many common units. The example below demonstrates how to use 
% this struct to cleanly and effeciently perform matlab calculations 
% using commonly encountered units and physical constants. The code 
% is easily modified to include any non-standard unit or constant desired.
%
% To determine the exact syntax of a particular unit you would like to use,
% simply run the function with no semicolon and all the units will be
% displayed.
%
%  Using the units struct:
%  --------------------------------------------------------
%    First create a units struct by including in your code the following
%    line
%          u = units;
%    Then:
%    
%          %To enter a number in a given unit, MULTIPLY by the unit:
%                L = 5*Units.in   % matlab automatically displays L in SI units
%  
%          % To display in a desired unit, simply divide by that unit
%                L/Units.ft       % displays L in ft.
%  
%          % To convert between units, MULTIPLY by the starting unit, and
%          % DIVIDE by the ending unit:
%                Units.mi^2/Units.ft^2  %displays the number of feet^2 in one mile^2
%  
%          %More complicated units can be obtained through arithmatic
%                mach1 = 340.29*Units.m/Units.s;  %speed of sound 
%  
%             %Note... to make the speed of sound available wherever your
%             %units struct is defined, simply write:
%                Units.mach1 = 340.29*Units.m/Units.s;   %mach1 now part of units struct
% 
%
%
% %------  BEGIN EXAMPLE CODE --------------------------------
% %This is an example calculation that uses the units mfile to calculate the
% %pressure at the bottom of a long vertically oriented pipe that is capped 
% %at the bottom and filled with oil.
% 
% u = units;
% pipeInnerDiameter = 4*Units.in;     %4 inch inner diameter
% pipeHeight = 30*Units.ft;           %pipe sticks 30 feet up into the air
% densityOfOil = 0.926*Units.gm/Units.cc; %density of oil as found on some random web site = .926 gm/cc
% pipeCrossSectionArea = pi*(pipeInnerDiameter/2)^2;  %cross sectional area of pipe bore
% volumeOfOil = pipeCrossSectionArea * pipeHeight;    %volume of oil that the pipe can hold
% pressurePipeBottom = densityOfOil * Units.g * pipeHeight;  %pressure formula from physics: P = rho*g*h.
% forceOnPipeBottom = pressurePipeBottom * pipeCrossSectionArea;  %force exterted on bottom cap of the pipe.
% 
% %Note that each variable holds its value as expressed in SI units.  To
% %express these values in different units, simply divide by the desired
% %unit as shown below.
% line1 = sprintf('A %2.3g inch diameter pipe sticking %3.3g meters into the air and filled',pipeInnerDiameter/Units.in, pipeHeight/Units.m);
% line2 = sprintf('with %3.3g fluid ounces of oil will have a pressure at the bottom of %4.4g psi.',volumeOfOil/Units.floz, pressurePipeBottom/Units.psi);
% line3 = sprintf('This will cause a total force of %5.5g lbs to press on the bottom cap of the pipe.',forceOnPipeBottom/Units.lbf);
% 
% textVal = sprintf('\n\n%s\n%s\n%s\n',line1,line2,line3);
% disp(textVal);
%%------  END EXAMPLE CODE --------------------------------



%============ START THE ACTUAL CODE TO DEFINE THE UNITS STRUCT =========
%-------- UNITS ------------------------------
%------- length ----
Units.m = 1;
Units.km = 1e3*Units.m;
Units.cm = 1e-2*Units.m;
Units.mm = 1e-3*Units.m;
Units.um = 1e-6*Units.m;
Units.nm = 1e-9*Units.m;
Units.ang = 1e-10*Units.m;
Units.in = 2.54*Units.cm;
Units.mil = 1e-3*Units.in;
Units.ft = 12*Units.in;
Units.yd = 3*Units.ft;
Units.mi = 5280*Units.ft;
Units.a0 = .529e-10*Units.m;

%------- Volume -------
Units.cc = (Units.cm)^3;
Units.L = 1000*Units.cc;
Units.mL = Units.cc;
Units.floz = 29.5735297*Units.cc;
Units.pint = 473.176475*Units.cc;
Units.quart = 946.35295*Units.cc;
Units.gal = 3.78541197*Units.L;

%----- mass ---------
Units.kg = 1;
Units.gm = 1e-3*Units.kg;
Units.mg = 1e-3*Units.gm;
Units.lb = 0.45359237*Units.kg;
Units.oz = (1/16)*Units.lb;
Units.amu = 1.66e-27*Units.kg;

%---- time -------
Units.s = 1;
Units.ms = 1e-3*Units.s;
Units.us = 1e-6*Units.s;
Units.ns = 1e-9*Units.s;
Units.ps = 1e-12*Units.s;
Units.min = 60*Units.s;
Units.hr = 60*Units.min;
Units.day = 24*Units.hr;
Units.yr = 365.242199*Units.day; 

%---- frequency ---- 
Units.Hz = 1/Units.s;
Units.kHz = 1e3 *Units.Hz;
Units.MHz = 1e6 *Units.Hz;
Units.GHz = 1e9 *Units.Hz;

%---- force -------
Units.N = 1;
Units.dyne = 1e-5*Units.N;
Units.lbf = 4.44822*Units.N;


%----- energy -----
Units.J = 1;
Units.MJ = 1e6*Units.J;
Units.kJ = 1e3*Units.J;
Units.mJ = 1e-3*Units.J;
Units.uJ = 1e-6*Units.J;
Units.nJ = 1e-9*Units.J;
Units.eV = 1.6022e-19*Units.J;
Units.BTU = 1.0550559e3*Units.J;
Units.kWh = 3.6e6*Units.J;
Units.cal = 4.1868*Units.J;
Units.kcal = 1e3*Units.cal;

%---- temperature ---
Units.K = 1;
Units.mK = 1e-3*Units.K;
Units.uK = 1e-6*Units.K;
Units.nK = 1e-9*Units.K;

%---- pressure -----
Units.Pa = 1;
Units.torr = 133.322*Units.Pa;
Units.mtorr = 1e-3*Units.torr;
Units.bar = 1e5*Units.Pa;
Units.mbar = 1e-3*Units.bar;
Units.atm = 1.013e5*Units.Pa;
Units.psi = 6.895e3*Units.Pa;



%----- power --- ---
Units.W = 1;
Units.MW = 1e6*Units.W;
Units.kW = 1e3*Units.W;
Units.mW = 1e-3*Units.W;
Units.uW = 1e-6*Units.W;
Units.nW = 1e-9*Units.W;
Units.pW = 1e-12*Units.W;
Units.hp = 745.69987*Units.W;

%------ charge ------
Units.coul = 1;
Units.e = 1.6022e-19*Units.coul;                        % elementary charge


%------ Voltage -----
Units.V = 1;
Units.kV = 1e3*Units.V;
Units.mV = 1e-3*Units.V;
Units.uV = 1e-6*Units.V;

%----- Current ------
Units.A = 1;
Units.mA = 1e-3*Units.A;
Units.uA = 1e-6*Units.A;
Units.nA = 1e-9*Units.A;

%----magnetic field -----
Units.T = 1;
Units.gauss = 1e-4*Units.T;



%----fundamental constants ----
Units.g = 9.80665*Units.m/Units.s^2;
Units.kB = 1.38e-23*Units.J/Units.K;                        % Boltzman constant
Units.sigma_SB = 5.670e-8 * Units.W/(Units.m^2 * Units.K^4);
Units.h = 6.626e-34 * Units.J*Units.s;
Units.hbar = Units.h/(2*pi);
Units.mu_B = 9.274e-24 * Units.J/Units.T;
Units.mu_N = 5.0507866e-27 * Units.J/Units.T;
Units.c = 2.99792458e8*Units.m/Units.s;
Units.eps0 = 8.8541878176204e-12* Units.coul/(Units.V*Units.m);
Units.mu0 = 1.2566370614359e-6 * Units.J/(Units.m*Units.A^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% added some more common units in space missing in units.m

%-------- UNITS ------------------------------
%------- length ----
Units.AU = 1.496e11*Units.m;
Units.R_Earth = 6370e3*Units.m;                     % Earth radius
Units.RE = Units.R_Earth;                          % Earth radius
Units.R_Sun = 6.96e8*Units.m;                       % Solar radius 
Units.pc = 3.0857e16*Units.m;                        % parsec

%------- mass ----
Units.M_Earth = 5.9742e24*Units.kg;                 % Mass of the Earth
Units.M_Sun = 1.98892e30*Units.kg;                  % Mass of the Sun
Units.me = 9.1094e-31*Units.kg;                      % electron mass
Units.mp = 1.6726e-27*Units.kg;                      % proton mass

%---- frequency ---- 
Units.mHz = 1e-3*Units.Hz;

%----- energy -----
Units.keV = 1e3*Units.eV;
Units.MeV = 1e6*Units.eV;
Units.erg = 1e-7*Units.J;

%---- temperature ---
Units.MK = 1e6*Units.K;

%----magnetic field -----
Units.nT = 1e-9*Units.T;

%----fundamental constants ----
Units.G = 6.6726e-11*Units.N*Units.m^2/Units.kg^2;           % gravitational constant
Units.kB = 1.38e-23*Units.J/Units.K;                     % Boltzmann constant




