function [neSC,Iph0,Tph0,Iph1,Tph1] = scpot2ne(varargin)
% SCPOT2NE Compute number density from spacecraft potential. Function uses
% number density and electron temperature to compute electron thermal
% current. A best fit of the photoelectron current to the thermal current is
% used to determine the photoelectron currents and temperatures. These are
% then used to construct a new number density from the spacecraft potential.
% 
% [neSC,Iph0,Tph0,Iph1,Tph1] = mms.scpot2ne(SCpot,ne,Te[,Iasp]);
% 
% Inputs: 
%   SCpot - Spacecraft potential (TSeries format)
%   ne - electron number density (TSeries format)
%   Te - electron temperature (TSeries format). Function accepts scalar
%   temperature or tensor (scalar temperature is used in either case).
%   Iasp - ASPOC current in muA (TSeries format) [optional]
%
% Outputs:
%   neSC - number density estimated from SCpot, at the same resolution as
%   SCpot.
%   Iph0,Tph0,Iph1,Tph1 - Values of the photoelectron currents (muA) and
%   temperatures (eV). Two photoelectron populations are assumed. 
%
% Notes: 
%   * Usual assumptions are made for thermal and photoelectron current, vis.,
%   planar geometry for photoelectrons and spherical geometry for thermal
%   electrons. 
%   * Currently the calculation neglects the ion thermal current, secondary
%   electrons, and other sources of current.
%   * ASPOC on does not work very well.
% 
% Written by D. B. Graham

% Check input
SCpot = varargin{1};
ne = varargin{2};
Te = varargin{3};

nargin
ASPOCon = 0;
if nargin == 4
    Iasp = varargin{4};
    ASPOCon = 1;
    Iasp = Iasp.resample(ne);
end
    


% Check format of electron temperature
dimTe = size(squeeze(Te.data(1,:,:)));
dimTe = dimTe(1)*dimTe(2);
if dimTe == 9
    Te = irf.ts_tensor_xyz(Te.time,Te.data); 
    Te = Te.trace/3;
else
    if dimTe ~= 1
        irf.log('critical','Te format not recognized');
        return;
    end
end

% Define constants
Units = irf_units;
me = Units.me;
qe = Units.e;
Ssurf = 34;

veth = sqrt(2*qe*Te.data/me);

SCpotr = SCpot.resample(ne);

Ie = (1e12*qe*Ssurf/(2*sqrt(pi)))*ne.data.*veth.*(1+SCpotr.data./Te.data); % Thermal current in muA

% First a simple fit of Iph to Ie using 1 photoelectron population
if ASPOCon
    fsimp = @(x) sum(abs(Ie+Iasp.data-(x(1).*exp(-SCpotr.data./x(2)))),'omitnan');
else
    fsimp = @(x) sum(abs(Ie-(x(1).*exp(-SCpotr.data./x(2)))),'omitnan');
end

options = optimset('MaxFunEvals',5000);
[Xsimp,~] = fminsearch(@(x) fsimp(x),[500;3],options); % Nelder-Mead method
g1 = Xsimp(1);
g2 = Xsimp(2);

% Fit of Iph to Ie for two photoelectron populations
if ASPOCon
    f = @(x) sum(abs(Ie+Iasp.data-(x(1).*exp(-SCpotr.data./x(2)) + x(3).*exp(-SCpotr.data./x(4)))),'omitnan');
else
    f = @(x) sum(abs(Ie-(x(1).*exp(-SCpotr.data./x(2)) + x(3).*exp(-SCpotr.data./x(4)))),'omitnan');
end
[X,~] = fminsearch(@(x) f(x),[g1;g2;10;10],options);

Iph0 = X(1);
Tph0 = X(2);
Iph1 = X(3);
Tph1 = X(4);

Ter = Te.resample(SCpot);

veth = sqrt(2*qe*Ter.data/me);
neSC = 2*sqrt(pi)*(Iph0*exp(-SCpot.data/Tph0)+Iph1*exp(-SCpot.data/Tph1))*1e-12./(Ssurf*qe*veth.*(1 + SCpot.data./Ter.data));
neSC = irf.ts_scalar(SCpot.time,neSC);

end