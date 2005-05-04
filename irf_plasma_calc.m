function [Wpe,Wce,Wuh,Wpp,Wcp,WpO,WcO,Va,Vte,Le] = irf_plasma_calc(B,n,no,Te,Ti,noshow)
%IRF_PLASMA_CALC   Calculate basic plasma quantities
%
% irf_plasma_calc(B,n,no,Te,Ti)
%
% [Wpe,Wce,Wuh,Wpp,Wcp,WpO,WcO,Va,Vte,Le] = irf_plasma_calc(B,n,no,Te,Ti,flag);
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

% Copyright 1997-2005 Yuri Khotyaintsev
if nargin == 0, % take input from terminal
    B=irf_ask('Magnetic field in nT [%] >','B',10);
    n=irf_ask('H+ desity in cc [%] >','n',1);
    no=irf_ask('Oxygen density in percent [%] >','no',0);
    Te=irf_ask('Electron  temperature in eV [%] >','Te',100);
    Ti=irf_ask('Ion  temperature in eV [%] >','Ti',1000);
end
if nargin < 6
    noshow = 0;
end
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

Mp_Me = 1836.15;

Wpe = 8.973*sqrt(n);
Wce = 2.8e-2*B;
Wpp = 8.973*sqrt(n*(1-no/100)/Mp_Me);
WpO = 8.973*sqrt(n*(no/100)/Mp_Me/16);
Va = 21.8*B./sqrt(n.*(1+15*no/100));
Vte = 4.19*1e2*sqrt(Te);
Vtp = 9.79*sqrt(Ti);
Vts = 9.79*sqrt(Te);
VtO = Vtp/4;
Le = 5.3e3/sqrt(n);
Li = 3e8/(2*pi*Wpp*1e3);

Wpe = Wpe*1e3;
Wce = Wce*1e3;
Wuh = sqrt(Wce.^2+Wpe.^2);
Wpp = Wpp*1e3;
Wcp = Wce/Mp_Me;
WpO = WpO*1e3;
WcO = Wce/Mp_Me/16;
Wlh = sqrt((Wpp.^2).*Wce.^2./(Wce.^2+Wpe.^2)+Wcp.^2);

Roe = Vte./Wce*1e3; % in meters
Rop = Vtp./Wcp*1e3; % in meters
RoO = VtO./WcO*1e3;
Ros = Vts./Wcp*1e3;


% variables to return in case time series input
if strcmp(flag_time_series,'yes')
    var={...
        {'Wpe','Hz'},{'Wce','Hz'}...
        {'Wpp','Hz'},{'Wcp','Hz'}...
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequencies

freq = [Wpe Wce Wuh Wlh Wpp Wcp WpO WcO];
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

disp(sprintf('\nPlasma velcities\n'))
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
