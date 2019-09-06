function Units=irf_units(varargin)
% IRF_UNITS returns structure containing key units and constants
%  uses units.m from web matlab_exchange depository
%  adds also space relevant stuff
%
%  Units = IRF_UNITS returns structure with all the values
%
%  IRF_UNITS('pattern')
%  IRF_UNITS pattern
%    displays all the units matching pattern on command line
%
% Examples:
%  Units=irf_units;
%  Units                                     % display all values
%  Units.Earth                               % display all values for Earth
%   Units.R_Earth/Units.R_Sun                % display Earth radius in solar radia
%   T_in_eV = Units.kB*T_in_MK*1e6 / Units.e % to convert from MK to eV

% Sources used:
%   CODATA 2018 - https://physics.nist.gov/cuu/Constants/Table/allascii.txt (using "BIPM 2019")
%   IAU 2012    - https://www.iau.org/static/resolutions/IAU2012_English.pdf
%   IAU 2015    - https://www.iau.org/static/resolutions/IAU2015_English.pdf
%   CGPM 1903   - https://www.bipm.org/en/CGPM/db/3/2/
%   CGPM 1954   - https://www.bipm.org/en/CGPM/db/10/4/
%   BIPM 2019   - https://www.bipm.org/utils/common/pdf/si-brochure/SI-Brochure-9.pdf
%   IAU WG 2009 - DOI 10.1007/s10569-010-9320-4

if nargin == 1 % display matching units
  fid = fopen(which('irf_units'));
  while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline) && ~strcmp(tline(1),'%') % check that not comment line
      if strfind(lower(tline),lower(varargin{:})) %#ok<STRIFCND>
        disp(tline);
      end
    end
  end
  fclose(fid);
  return
end

if nargout==0 && nargin==0
  disp('!!!WARNING!!!!');
  disp('IRF_UNITS usage has been changed, from script it has become function.');
  disp('please check the help!');
  disp('Instead of :')
  disp('> irf_units');
  disp('please run :');
  disp('> Units=irf_units;');
  clear Units;
  return;
end

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
Units.amu = 1.66053906660e-27*Units.kg; % Src: CODATA 2018

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
Units.eV = 1.602176634e-19*Units.J; % Exact, Src: BIPM 2019 - Table 8
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
Units.bar = 1e5*Units.Pa;       % Exact
Units.mbar = 1e-3*Units.bar;
Units.atm = 101325*Units.Pa;    % Exact, Src: CGPM 1954
Units.psi = 6.895e3*Units.Pa;
Units.torr = (1/760)*Units.atm; % Exact
Units.mtorr = 1e-3*Units.torr;

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
Units.e = 1.602176634e-19*Units.coul; % elementary charge, exact. Src: BIPM 2019 - Table 1

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
Units.nT = 1e-9*Units.T;
Units.gauss = 1e-4*Units.T;

%----fundamental constants ----
Units.g = 9.80665*Units.m/Units.s^2;                            % gravitational acceleration, Src: CGPM 1903
Units.G = 6.67430e-11*Units.m^3/Units.kg/Units.s^2;             % Newtonian graviational constant, Src: CODATA 2018
Units.kB = 1.380649e-23*Units.J/Units.K;                        % Boltzman constant, exact. Src: BIPM 2019 - Table 1
Units.sigma_SB = 5.670374419e-8 * Units.W/(Units.m^2 * Units.K^4); % Stefan-Boltzmann constant, exact but here given only values as in Src: CODATA 2018
Units.h = 6.62607015e-34 * Units.J*Units.s;                     % Planck constant, exact. Src: BIPM 2019 - Table 1
Units.hbar = Units.h/(2*pi);                                    % Planck constant
Units.mu_B = 9.2740100783e-24 * Units.J/Units.T;                % Bohr magneton, Src: CODATA 2018
Units.mu_N = 5.0507837461e-27 * Units.J/Units.T;                % Nuclear magneton, Src: CODATA 2018
Units.c = 2.99792458e8*Units.m/Units.s;                         % Speed of light, exact. Src: BIPM 2019 - Table 1
Units.mu0 = 1.25663706212e-6 * Units.J/(Units.m*Units.A^2);     % Magnetic constant, Src: CODATA 2018
Units.eps0 = 8.8541878128e-12; %*Units.F/Units.m;               % Vacuum permittivity, Src: CODATA 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% added some more common units in space missing in units.m

%-------- UNITS ------------------------------
%------- length ----
Units.AU = 149597870700*Units.m;                % astronomical unit, exact, Src: IAU 2012
Units.R_Sun = 6.96e8*Units.m;                   % Solar radius, Src: IAU WG 2009 - Table 4
Units.Sun.radius = Units.R_Sun;
Units.pc = (648000/pi)*Units.AU;                % parsec, exact, Src: IAU 2015
Units.Uranus.distanceTSun = 19.1914*Units.AU;   % Uranus orbit, semimajor axis
Units.Neptune.distanceToSun = 30.0611*Units.AU; % Neptune orbit, semimajor axis

%------- mass ----
Units.M_Earth = 5.9742e24*Units.kg;    % Mass of the Earth
Units.Earth.mass = Units.M_Earth;
Units.M_Sun = 1.98892e30*Units.kg;     % Mass of the Sun
Units.Sun.mass = Units.M_Sun;
Units.me = 9.1093837015e-31*Units.kg;  % electron mass, Src: CODATA 2018
Units.mp = 1.67262192369e-27*Units.kg; % proton mass, Src: CODATA 2018

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
Units.gauss = 1e-4*Units.T;

%---- MERCURY -----
Units.Mercury.semiMajorAxis = 57909100*Units.km;  % Mercury semimajor axis
Units.Merucry.distanceToSun  = 57909100*Units.km; % Mercury orbit, semimajor axis
Units.Mercury.aphelion = 69816900*Units.km;
Units.Mercury.perihelion  = 46001200*Units.km;
Units.Mercury.radius = 2439.7*Units.km;           % Mercury radius (mean), Src: IAU WG 2009 - Table 4

%---- MARS -----
Units.Mars.distanceToSun  = 1.5273*Units.AU;     % Mars orbit, semimajor axis
Units.Mars.radius    = 3396.19*Units.km;         % Mars radius (equatorial), Src: IAU WG 2009 - Table 4
Units.Mars.radiusEquatorial  = 3396.19*Units.km; % Mars radius (equatorial), Src: IAU WG 2009 - Table 4
Units.Mars.radiusPolar    = 3376.2*Units.km;     % Mars radius (average polar), Src: IAU WG 2009 - Table 4

%---- EARTH -----
Units.Earth.semiMajorAxis  = 149598261*Units.km; % Earth semimajor axis
Units.Earth.distanceToSun  = 149598261*Units.km; % Earth distance to Sun (=semimajor axis)
Units.Earth.radius = 6371.2e3*Units.m;           % Earth radius
Units.R_Earth  = 6371.2e3*Units.m;               % Earth radius
Units.RE = Units.R_Earth;                        % Earth radius

%---- VENUS -----
Units.Venus.distanceToSun  = 0.7233*Units.AU; % Venus orbit, semimajor axis
Units.Venus.radius    = 6051.8*Units.km;      % Src: IAU WG 2009 - Table 4

%---- SATURN -----
Units.Saturn.distanceToSun  = 9.5388*Units.AU; % Saturn orbit, semimajor axis
Units.Saturn.radius    = 60268*Units.km;       % Saturn equatorial radius (one-bar surface), Src: IAU WG 2009 - Table 4

%---- JUPITER -----
Units.Jupiter.radius    = 69911*Units.km;         % Jupiter mean radius (one-bar surface), Src: IAU WG 2009 - Table 4
Units.Jupiter.radiusEquatorial = 71492*Units.km;  % Jupiter equatorial radius (one-bar surface), Src: IAU WG 2009 - Table 4
Units.Jupiter.radiusPolar   = 66854*Units.km;     % Jupiter polar radius (one-bar surface), Src: IAU WG 2009 - Table 4
Units.Jupiter.semiMajorAxis = 778547200*Units.km;
Units.Jupiter.distanceToSun  = 778547200*Units.km; % Jupiter semimajor axis
Units.Jupiter.aphelion     = 816520800*Units.km;
Units.Jupiter.perihelion   = 740573600*Units.km;
Units.Jupiter.mass      = 1.8986e27*Units.kg;

%---- GANYMEDE -----
Units.Ganymede.semiMajorAxis = 1070400*Units.km;
Units.Ganymede.radius = 2631.2*Units.km;         % Src: IAU WG 2009 - Table 5

%---- EUROPA -----
Units.Europa.semiMajorAxis  = 670900*Units.km; %
Units.Europa.radius  = 1560.8*Units.km;        % Europa radius (mean), Src: IAU WG 2009 - Table 5

%---- CALLISTO -----
Units.Callisto.semiMajorAxis = 1882700*Units.km; %
Units.Callisto.radius  = 2410.3*Units.km;        % Src: IAU WG 2009 - Table 5

end
