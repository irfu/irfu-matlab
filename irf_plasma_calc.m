function [Fpe,Fce,Fuh,Fpp,Fcp,FpO,FcO,Va,Vte,Le] = irf_plasma_calc(B_inp,n_inp,no_inp,Te_inp,Ti_inp,noshow)
%IRF_PLASMA_CALC   Calculate basic plasma quantities
%
% irf_plasma_calc(B,n,no,Te,Ti)
%
% [Fpe,Fce,Fuh,Fpp,Fcp,FpO,FcO,Va,Vte,Le] = irf_plasma_calc(B,n,no,Te,Ti,flag);
%
%	B - magnetic field [nT]
%	n - density [cm^-3]
%	no - content of O+ [%]
%	Te - electron temperature [eV]
%	Ti - ion temperature [eV]
% flag - 1: do not display the output
% flag - string : return only specified output
%        string canbe 'Wpe','Wce','Wuh',etc.
%
% returned variables are in SI units, frequencies in Hz, lengths in m
%
% $Id$

persistent B n no Te Ti

if nargin < 6
    noshow = 0;
end
if nargin >= 1, B=B_inp; end
if nargin >= 2, n=n_inp; end
if nargin >= 3, no=no_inp; end
if nargin >= 4, Te=Te_inp; end
if nargin >= 5, To=To_inp; end

% Copyright 1997-2005 Yuri Khotyaintsev
if nargin < 1, B=irf_ask('Magnetic field in nT [%] >','B',10);end
if nargin < 2, n=irf_ask('H+ desity in cc [%] >','n',1);end
if nargin < 3, no=irf_ask('Oxygen density in percent [%] >','no',0);end
if nargin < 4, Te=irf_ask('Electron  temperature in eV [%] >','Te',100);end
if nargin < 5, Ti=irf_ask('Ion  temperature in eV [%] >','Ti',1000);end

%if time series are supplied then time series shoud be returned
if size(B,2)>1, % we have time series of density
    t=B(:,1); % time axis
    if length(t)>1 & ~isstr(noshow), noshow=1;end % if more than  one time point do not print result
    B(:,1)=[]; % delete time column
    if size(B,2)>3; %asssume that column 4 is amplitude, delete other coolumns
        B(:,[1:3 5:end])=[];
    elseif size(B,2)==3, % assume there are three columns Bx,By,Bz
        B(:,4)=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);
        B(:,1:3)=[]; %leave only amplitude
    else % do not know what to do if several B columns, take only first
        B(:,2:end)=[];
    end
    flag_time_series='yes';
else
    flag_time_series='no';
end
if strcmp(flag_time_series,'yes'), % check that other variables are time series, if not interpoplate
    variables_to_check={'n','no','Te','Ti'};
    for j=1,length(variables_to_check)
        if eval(['size(' variables_to_check{j} ',2)'])>1, % we have time series
            c_eval('?=irf_resamp(?,t);',variables_to_check{j}) % resample to new time axis
            c_eval('?(:,[1,3:end])=[];',variables_to_check{j}) % delete time and  other columns
        elseif eval(['prod(size(' variables_to_check{j} '))==1']), % only one number
            c_eval('?=repmat(?,size(t));',variables_to_check{j})
        else
            c_eval('irf_log(''proc'',''do not understand input <?>'');',variables_to_check{j});
            return;
        end
    end
end

Me=9.1094e-31; % electron mass
Mp=1.6726e-27; % proton mass
c=2.9979e8; % speed of light
e=1.6022e-19; % elementary charge
Mp_Me = Mp/Me; % ratio of proton and electron mass 1836.15;

% in formulas it is more convenient to use variables expresesd in SI units 
B_SI=B*1e-9; % [T]

Wpe = 8.973*sqrt(n)*1e3*2*pi; % rad/s
Wce = e*B_SI/Me;   % rad/s
Wpp = 8.973*sqrt(n*(1-no/100)/Mp_Me);
WpO = 8.973*sqrt(n*(no/100)/Mp_Me/16);
Va = 21.8*B./sqrt(n.*(1+15*no/100));
%Vte = 4.19*1e2*sqrt(Te);
Vte = c*sqrt(1-1/(Te*e/(Me*c^2)+1)^2);  % m/s
%Vtp = 9.79*sqrt(Ti);
Vtp = c*sqrt(1-1/(Ti*e/(Mp*c^2)+1)^2);   % m/s
Vts = 9.79*sqrt(Te);  % ? relativistic formula???
%VtO = Vtp/4;
VtO = c*sqrt(1-1/(Ti*e/(16*Mp*c^2)+1)^2);   % m/s
gamma_e=1/sqrt(1-(Vte/c).^2);
gamma_p=1/sqrt(1-(Vtp/c).^2);
gamma_O=1/sqrt(1-(VtO/c).^2);
Le = 5.3e3/sqrt(n);
Li = 3e8/(2*pi*Wpp*1e3);

Fpe = Wpe/2/pi; % Hz
Fce = Wce/2/pi;
Fuh = sqrt(Fce.^2+Fpe.^2);
Fpp = Wpp/2/pi;
Fcp = Fce/Mp_Me;
FpO = WpO/2/pi;
FcO = Fce/Mp_Me/16;
Flh = sqrt((Fpp.^2).*Fce.^2./(Fce.^2+Fpe.^2)+Fcp.^2);

Roe = Me*c/(e*B_SI)*sqrt(gamma_e.^2-1); % m, relativistically correct
Rop = Mp*c/(e*B_SI)*sqrt(gamma_p.^2-1); % m, relativistically correct
RoO = Mp*16*c/(e*B_SI)*sqrt(gamma_O.^2-1); % m, relativistically correct
%Roe = Vte./Fce/2/pi; % m
%Rop = Vtp./Fcp/2/pi; % m
%RoO = VtO./FcO/2/pi; % m
Ros = Vts./Fcp/2/pi; % m


% variables to return in case time series input
if strcmp(flag_time_series,'yes')
    var={...
        {'Fpe','Hz'},{'Fce','Hz'}...
        {'Fpp','Hz'},{'Fcp','Hz'}...
        };
    for j=1:length(var),
        eval([var{j}{1} '= [t ' var{j}{1} '];']);
        if isstr(noshow)
            if strcmp(noshow,var{j}{1}),
                eval(['Wpe=' var{j}{1} ';']);
                return;
            end
        end
    end
end


if noshow > 0
    return
end

disp('===============================================================')
disp('IRFU plasma calculator, relativistic effects not fully included')
disp('velocities, gyroradia are relativistically correct')
disp('can somebody fix relativstically correct frequencies Fpe, Fce,.. ?')
disp('===============================================================')
disp(['B=' num2str(B) ' [nT]; n_H=' num2str(n) ' [cc]; n_O=' num2str(no) ' [cc]; ' ...
  'T_e=' num2str(Te) ' [eV]; T_i=' num2str(Ti) ' [eV];']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequencies

freq = [Fpe Fce Fuh Flh Fpp Fcp FpO FcO];
freqs = {'F_pe'; 'F_ce'; 'F_uh'; 'F_lh'; 'F_pp'; 'F_cp'; 'F_pO'; 'F_cO'};
disp(sprintf('\nPlasma frequencies\n'))
for ii = 1:length(freq)
    val = freq(ii);
    if val > 1e6
        units = '[MHz]';
        koef = 1e-6;
    elseif val < 1e6 & val > 1e3
        units = '[KHz]';
        koef = 1e-3;
    elseif val <1e3 & val >.1
        units = '[Hz]';
        koef = 1;
    else
        units = '[mHz]';
        koef = 1e3;
    end
    disp(sprintf('%s = %5.2f %s',freqs{ii},val*koef, units))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velosities

v = [Va Vte Vtp Vts VtO ];
vs = {'V_a'; 'V_Te'; 'V_Tp'; 'C_s'; 'V_TO'};

disp(sprintf('\nPlasma velocities\n'))
for ii = 1:length(v)
    val = v(ii);
    if val > 1e3
        units = '[km/s 1e3]';
        koef = 1e-3;
    elseif val < 1e3 & val > 1
        units = '[km/s]';
        koef = 1;
    else
        units = '[m/s]';
        koef = 1e3;
    end
    disp(sprintf('%s = %5.2f %s',vs{ii},val*koef, units))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scales

l = [Le Li Roe Rop RoO Ros];
ls = {'L_e'; 'L_i'; 'Ro_e' ; 'Ro_p'; 'Ro_O'; 'Ro_s'};
disp(sprintf('\nPlasma scales\n'))
for ii = 1:length(l)
    val = l(ii);
    if val > 1e3
        units = '[km]';
        koef = 1e-3;
    elseif val < 1e3 & val > 1
        units = '[m]';
        koef = 1;
    else
        units = '[cm]';
        koef = 1e2;
    end
    disp(sprintf('%s = %5.2f %s',ls{ii},val*koef, units))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% etc

disp(sprintf('\nPlasma dimensionless parameters\n'))

beta  = Vtp^2/Va^2;

if beta < 10/Mp_Me
    disp(sprintf('beta*Mp_Me = %1.5f',Vtp^2/Va^2*Mp_Me))
else
    disp(sprintf('beta  = %1.5f',Vtp^2/Va^2))
end
