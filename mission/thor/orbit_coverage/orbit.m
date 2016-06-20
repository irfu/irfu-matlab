function [tt xx yy] = orbit(r_per,r_ap,Ttot,varargin)
% ORBIT  Calculates Kepler orbit as a function of time.
%   Orbit starts at perihelion at t = 0.
%   
%   Example
%   [t,x,y] = orbit(r_per,r_ap,Ttot,method,'dt',dt);
%       Output:
%           t - time, s
%           x - x in GSE coordinates, RE
%           y - y in GSE coordinates, RE
%       Input:
%           r_per - perigee, RE
%           r_ap - apogee, RE
%           Ttot - total time to calculate the orbit for. If Ttot is larger
%                  than the period of the orbit, the orbit is only
%                  duplicated. It is assumed the orbit is fixed wrt the 
%                  star system, so a precession frequency (?) is added
%                  that displaces the orbit slightly.
%           method - 'E' - constant increase in eccentric anomaly, dt is
%                          not constant unless a 'dt'/dt is given, then the
%                          time series is interpolated to values dt
%                    't' - uses a constant timestep, given in 'dt'/dt, if
%                          not given, a default value of 2 min is used
%           dt - timestep, s, optional, if not given, the eccentric anomaly
%                is increased with constant step instead and dt is not
%                constant. If method 't' is given, but not 'dt'/dt, then a
%                default timestep of 120 s is used.

method = 'E'; % default
    
while ~isempty(varargin)
    if strcmp(varargin{1},'dt')        
        dt = varargin{2};
        varargin(1:2)=[];
    elseif strcmp(varargin{1},'t')
        method = 't';        
        varargin(1)=[];
    elseif strcmp(varargin{1},'E')
        method = 'E';        
        varargin(1)=[];
    end
end

% Physical parameters
G = 6.67384e-11; % m^3 kg^-1 s^-2 (N m^-2 kg^-2)
ME = 5.9722e24; % kg
RE = 6371*1e3; % m

r_ap = r_ap*RE; % m
r_per = r_per*RE; % m

% Orbital parameters
a = (r_ap + r_per)/2; % m, semi-major axis
mu = G*ME; % m^3 s^-2
e = (r_ap - r_per)/(r_ap + r_per); % eccentricity
n = sqrt(mu/a^3);
T = 2*pi/n; % orbital period
disp(['Apogee: ' num2str(r_ap/RE,'%.1f') ' RE,  Perigee: ' num2str(r_per/RE,'%.1f') ' RE'])
disp(['Orbit period: ' num2str(T/60/60,'%.2f') ' hours'])
disp(['Eccentricity: ' num2str(e,'%.2f')])

switch lower(method)
    case 'e' % constant eccentric anomaly step
        E=0:.01:2*pi; % eccentric anomaly
        %f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
        %f(f<0)=f(f<0)+2*pi;
        %f = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));
        f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); % true anomaly (angle from x,y=0)
        b = a*sqrt(1-e^2); % semi-minor axis?
        t = a*sqrt(a/mu)*(E-e*sin(E)); % time       
        x = a*(cos(E)-e);
        y = b*sin(E);
        r = sqrt(x.^2+y.^2);
    otherwise % constant timestep
        if ~exist('dt','var')
            dt = 120; % s, 120s=2min
        end
        if Ttot<T; t = 0:dt:Ttot;   
        else t = 0:dt:T;
        end
        nt = numel(t);
        tau = 0; % time of perihelion passage
        M = n*(t-tau); % mean anomaly;
        funEM = @(E,M) E-e*sin(E)-M;
        E = fsolve(@(E) funEM(E,M),M);
        r = a*(1-e*cos(E)); % radial distance from Earth
        f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); % true anomaly (angle from x,y=0)
        %f = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));
end        
                      
if Ttot<T; 
    rr = r;
    ff = f;
    tt = t;
else % copy orbit
    nT = ceil(Ttot/T);
    rr = repmat(r,1,nT);
    ff = repmat(f,1,nT);
    tt = [];
    for kk=1:nT;
        tt = [tt t+(kk-1)*T];        
    end
    %tt = linspace(0,nT*T,nT*nt);
end

% add precession of orbit as Earth goes around Sun
% 2*pi rad degrees in 60*60*24*365 seconds (one year)
ffplus = tt*2*pi/(60*60*24*365);
xx = rr.*cos(ff+ffplus);
yy = rr.*sin(ff+ffplus);

% remove again that extra odd part of an orbit
xx(tt>Ttot)=[];
yy(tt>Ttot)=[];
tt(tt>Ttot)=[];


% interpolate to constant dt
if exist('dt','var')
    tinterp = 1:dt:tt(end);
    xx = interp1(tt,xx,tinterp);
    yy = interp1(tt,yy,tinterp);
    tt = tinterp;
    
end