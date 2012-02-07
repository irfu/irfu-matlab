function [Fpe_out,Fce,Fuh,Fpp,Fcp,FpO,FcO,Va,Vte,Le] = irf_plasma_calc(B_inp,n_inp,no_inp,Te_inp,Ti_inp,noshow)
%IRF_PLASMA_CALC   Calculate basic plasma quantities
% Calculate many important plasma parameters (also interactively)
% irf_plasma_calc(B,n,no,Te,Ti)
% irf_plasma_calc
% [Fpe,Fce,Fuh,Fpp,Fcp,FpO,FcO,Va,Vte,Le] = irf_plasma_calc(B,n,no,Te,Ti,flag);
%
%	B - magnetic field [nT]                 (can be time series with first column time)
%	n - density [cc]                        (can be time series with first column time only if B also is time series)
%	no - content of O+ [% of number density](can be time series with first column time only if B also is time series)
%	Te - electron temperature [eV]          (can be time series with first column time only if B also is time series)
%	Ti - ion temperature [eV]               (can be time series with first column time only if B also is time series)
% flag - 1: do not display the output
% flag - string : return only specified output
%        string can be 'Fpe','Fce','Fcp','Flh','Fuh','Fpp'
%             length - 'Li' (c/wpi),'Le' (c/wpe),'Ld' (Debye length)
%             speeds - 'Va' (Alfven), 'Vte','Vtp','VtO','Vts'
%         gyroradius - 'Roe','Rop','RoO'
%
% returned variables are in SI units, frequencies in Hz, lengths in m
%
% Here is also given a small table for simple/fast estimates of
% space plasma paremeters with precision <5%
%  units T[eV] B[nT] n[cc] E[mV/m]	
% 
%  Debye length [m]       = sqrt(200 T/n)
%  e- plasma f. [kHz]     = 9 sqrt(n)
%  e- gyrof. [Hz]         = 28 B
%  e- gyrorad. [km]       = sqrt(10 T) / B
%  e- inert. l [km]       = 5.3 / sqrt(n)
%  e- veloc. [km/s]       = 600 sqrt(T)
%  H+ plasma f. [Hz]      = 0.033 sqrt(n)
%  H+ gyrof. [Hz]         = 0.015 B
%  H+ veloc. [km/s]       = sqrt(200 T)
%  H+ gyror. [km]         = 100 sqrt(2 T) / B
%  H+ inert. l [km]       = 230 / sqrt(n)
%  H+ veloc. [km/s]       = 13.8 sqrt(T)
%  O+ gyrof. [mHz]        = B
%  O+ veloc. [km/s]       = sqrt(10 T)
%  O+ gyror. [km]         = 100 sqrt(33 T) / B
%  low. hybr. f. [Hz]     = sqrt(0.427 (B)^2 /(1+ 9.7e-6 (B)^2)/n);
%  Alfven vel. V_A[km/s]  = 22 B /sqrt(n)
%  Poynting fl. S[uW/m2]  = 0.8 E B
%  Plasma beta = 0.4 n T / B^2 = ([gyroradius]/[inertial length])^2 
%  Magnetic pressure[nPa] = (B/50)^2
%  E_corrotation  [mV/m]  = 0.6e-4 R[RE] B
 

% $Id$
% Copyright 1997-2005 Yuri Khotyaintsev

persistent B np_cc no_rel Te Ti

flag_display_values=1;  % dispaly all values on screen
flag_ask_parameters_interactively=1; % default do ask parameters on command line

if nargin < 6, noshow = 0; end
if nargin >= 6, flag_ask_parameters_interactively=0; end
if nargout,        flag_display_values=0;end

if nargin >= 1, B=B_inp; end
if nargin >= 2, np_cc=n_inp; end
if nargin >= 3, no_rel=no_inp; end
if nargin >= 4, Te=Te_inp; end
if nargin >= 5, Ti=Ti_inp; To=Ti; end % O+ temperature the same as for H+

if flag_ask_parameters_interactively
  if nargin < 1, 
    if numel(B)>1; B=B(1); end             % use only the first element in persistent variable
    B=irf_ask('Magnetic field in nT [%] >','B',10);
  end
  
  if nargin < 2, 
    if numel(np_cc)>1; np_cc=np_cc(1); end % use only the first element in persistent variable
    np_cc=irf_ask('H+ desity in cc [%] >','np_cc',1);
  end
  
  if nargin < 3, 
    if numel(no_rel)>1; no_rel=no_rel(1); end % use only the first element in persistent variable
    no_rel=irf_ask('Oxygen density in percent from H+ density [%] >','no_rel',0);
  end
  
  if nargin < 4, 
    if numel(Te)>1; Te=Te(1); end % use only the first element in persistent variable
    Te=irf_ask('Electron  temperature in eV [%] >','Te',100);
  end
  
  if nargin < 5, 
    if numel(Ti)>1; Ti=Ti(1); end % use only the first element in persistent variable
    Ti=irf_ask('Ion  temperature in eV [%] >','Ti',1000); To=Ti; 
  end
end

%if time series are supplied then time series shoud be returned
if size(B,2)>1,     % we have time series of B
    t=B(:,1);       % time axis
    B(:,1)=[];      % delete time column
    if size(B,2)>3; %asssume that column 4 is amplitude, delete other coolumns
        B(:,[1:3 5:end])=[];
    elseif size(B,2)==3, % assume there are three columns Bx,By,Bz
        B(:,4)=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);
        B(:,1:3)=[];  %leave only amplitude
    else              % do not know what to do if several B columns, take only first
        B(:,2:end)=[];
    end
    flag_time_series='yes';
else
    flag_time_series='no';
end
if strcmp(flag_time_series,'yes'), % check that other variables are time series, if not interpoplate
    variables_to_check={'np_cc','no_rel','Te','Ti','To'};
    for j=1:length(variables_to_check)
        if eval(['size(' variables_to_check{j} ',2)'])>1, % we have time series
            c_eval('?=irf_resamp(?,t);',variables_to_check(j)) % resample to new time axis
            c_eval('?(:,[1,3:end])=[];',variables_to_check(j)) % delete time and  other columns
        elseif eval(['prod(size(' variables_to_check{j} '))==1']), % only one number
            eval([variables_to_check{j} '=repmat(' variables_to_check{j} ',size(t));']);
        else
            c_eval('irf_log(''proc'',''do not understand input <?>'');',variables_to_check{j});
            return;
        end
    end
end

np=np_cc*1e6; % proton density m^-3
no=no_rel/100.*np; % proton density m^-3
n=np+no; % total plasma density  m^-3

irf_units; % read in standard units

Me=Units.me; % electron mass
Mp=Units.mp; % proton mass
c=Units.c; % speed of light
e=Units.e; % elementary charge
epso=Units.eps0; % vacuum dielectric constant 
mu0=Units.mu0; % Mu_0
Mp_Me = Mp/Me; % ratio of proton and electron mass 1836.15;

% in formulas it is more convenient to use variables expresesd in SI units 
B_SI=B*1e-9; % [T]

Wpe = sqrt(n*e^2/Me/epso); % rad/s
Wce = e*B_SI/Me;   % rad/s
Wpp = sqrt(np*e^2/Mp/epso);
WpO = sqrt(no*e^2/Mp/16/epso);
Va = B_SI./sqrt(mu0*(np+16*no)*Mp);
Vae = B_SI./sqrt(mu0*n*Me);
Vte = c*sqrt(1-1./(Te.*e./(Me*c^2)+1).^2);              % m/s (relativ. correct)
Vtp = c*sqrt(1-1./(Ti.*e./(Mp*c^2)+1).^2);              % m/s
Vts = Vtp.*sqrt(Te./Ti); Vts=Vts(:);                    % ? what is relativistic formula???
VtO = c*sqrt(1-1./(To.*e./(16*Mp*c^2)+1).^2);           % m/s
gamma_e=1./sqrt(1-(Vte/c).^2);
gamma_p=1./sqrt(1-(Vtp/c).^2);
gamma_O=1./sqrt(1-(VtO/c).^2);
Le = c./Wpe;
Li = c./Wpp;
Ld = Vte./Wpe ; % Debye length scale


Fpe = Wpe/2/pi; % Hz
Fce = Wce/2/pi;
Fuh = sqrt(Fce.^2+Fpe.^2);
Fpp = Wpp/2/pi;
Fcp = Fce/Mp_Me;
FpO = WpO/2/pi;
FcO = Fce/Mp_Me/16;
Flh = sqrt(Fcp.*Fce./(1+Fce.^2./Fpe.^2)+Fcp.^2);

Roe = Me*c./(e*B_SI).*sqrt(gamma_e.^2-1); % m, relativistically correct
Rop = Mp*c./(e*B_SI).*sqrt(gamma_p.^2-1); % m, relativistically correct
RoO = Mp*16*c./(e*B_SI).*sqrt(gamma_O.^2-1); % m, relativistically correct
Ros = Vts./Fcp/2/pi; % m


% variables to return in case time series input
if strcmp(flag_time_series,'yes')
    var={...
        {'Fpe','Hz'},{'Fce','Hz'}...
        {'Fpp','Hz'},{'Fcp','Hz'}...
        {'Fuh','Hz'},{'Flh','Hz'}...
        {'Le'},{'Li'},{'Ld'}...
        {'Va'},{'Vte'},{'Vtp'},{'VtO'},{'Vts'}...
        {'Roe'},{'Rop'},{'RoO'}...
        };
    for j=1:length(var),
        eval([var{j}{1} '= [t ' var{j}{1} '];']);
    end
end

if ischar(noshow) % return only particular value
    try
        eval(['Fpe_out=' noshow ';']);
        return;
    catch
        disp(['variable ''' noshow ''' not recognized']);
        Fpe_out=[];
        flag_display_values=0;
    end
end

if ~flag_display_values,    return; end 

if strcmp(flag_time_series,'yes'), 
    disp('!!! TIME SERIES. Showing only the values for the first point!');
end

disp('===============================================================')
disp('IRFU plasma calculator, relativistic effects not fully included')
disp('velocities, gyroradia are relativistically correct')
disp('can somebody fix relativstically correct frequencies Fpe, Fce,.. ?')
disp('===============================================================')
disp(['B=' num2str(B(1)) ' [nT]; n_H=' num2str(np(1)*1e-6) ' [cc]; n_O=' num2str(no(1)*1e-6) ' [cc]; ' ...
  'T_e=' num2str(Te(1)) ' [eV]; T_i=' num2str(Ti(1)) ' [eV];']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequencies

freq = [Fpe(:,end) Fce(:,end) Fuh(:,end) Flh(:,end) Fpp(:,end) Fcp(:,end) FpO(:,end) FcO(:,end)];
freqs = {'F_pe'; 'F_ce'; 'F_uh'; 'F_lh'; 'F_pp'; 'F_cp'; 'F_pO'; 'F_cO'};
fprintf('\nPlasma frequencies\n');
for ii = 1:size(freq,2)
    val = freq(1,ii);
    if val > 1e6
        units = '[MHz]';
        koef = 1e-6;
    elseif val < 1e6 && val > 1e3
        units = '[KHz]';
        koef = 1e-3;
    elseif val <1e3 && val >.1
        units = '[Hz]';
        koef = 1;
    else
        units = '[mHz]';
        koef = 1e3;
    end
    fprintf('\n%s = %5.2f %s',freqs{ii},val*koef, units);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocities

v = [Va(:,end) Vae(:,end) Vte(:,end) Vtp(:,end) Vts(:,end) VtO(:,end)];
vs = {'V_a'; 'V_ae'; 'V_Te'; 'V_Tp'; 'C_s'; 'V_TO'};

fprintf('\n\nPlasma velocities\n');
for ii = 1:size(v,2)
    val = v(1,ii);
    if val > 1e6
        units = '[km/s 1e3]';
        koef = 1e-6;
    elseif val < 1e6 && val > 1e3
        units = '[km/s]';
        koef = 1e-3;
    else
        units = '[m/s]';
        koef = 1;
    end
    fprintf('\n%s = %5.2f %s',vs{ii},val*koef, units);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scales

l = [Le(:,end) Li(:,end) Ld(:,end) Roe(:,end) Rop(:,end) RoO(:,end) Ros(:,end)];
ls = {'L_e'; 'L_i'; 'L_d' ; 'Ro_e' ; 'Ro_p'; 'Ro_O'; 'Ro_s'};
fprintf('\n\nPlasma scales\n');
for ii = 1:size(l,2)
    val = l(1,ii);
    if val > 1e3
        units = '[km]';
        koef = 1e-3;
    elseif val < 1e3 && val > 1
        units = '[m]';
        koef = 1;
    else
        units = '[cm]';
        koef = 1e2;
    end
    fprintf('\n%s = %5.2f %s',ls{ii},val*koef, units);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% etc

fprintf('\nPlasma dimensionless parameters\n');

beta  = Vtp.^2./Va.^2;
fprintf('\nbeta  = %1.5f',beta(1));
fprintf('\nbeta*sqrt(Mp/Me)  = %1.5f',beta(1)*sqrt(Mp_Me));
fprintf('\nbeta*(Mp/Me) = %1.5f',beta(1)*Mp_Me);
fprintf('\nGamma_e  = %1.2f\n',gamma_e);

if nargout>0, Fpe_out = Fpe; end
