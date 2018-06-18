function [J_probe, J_photo, J_plasma]=probe_current(probe,U_probe,R_sun,UV_factor,plasma)
% LP.PROBE_CURRENT calculate current to the probe
% J_probe=LP.PROBE_CURRENT(probe,U_probe,R_sun,UV_factor,plasma)
%
%   Calculates the total probe current to/from 
%   a cylindrical or spherical Langmuir probe.
%
% Input:
%  probe.type    - 'spherical','cylindrical','arbitrary'
%  probe.surface - 'themis','cassini' (one of flags in lp.photocurrent)
%  probe.cross_section_area - in m2
%  probe.total_area - in m2
%  U_probe    - probe potential (can be vector or matrix) 
%  R_sun      - distance from sun in AU
%  UV_factor  - default is 1
%  plasma     - describes plasma components (structure)
%    plasma.q - charge of species in e (the length of this vector corresponds to number of species)
%    plasma.m - mass of species in proton masses (0 corresponds to e- mass)
%    plasma.n - density of species [cc]
%    plasma.T - temperature [eV]
%    plasma.vsc - velocity of probe wrt. mmedia [m/s]
%
% [J_probe, J_photo, J_plasma]=LP.PROBE_CURRENT
%   Return current contributions from photoeletrons and all
%   the plasma components
%
% See also: LP.PHOTOCURRENT, LP.THERMAL_CURRENT

Units=irf_units;

if isempty(plasma) % calculate only photocurrent
    nSpecies=0;
else
    nSpecies=numel(plasma.q);
    J_plasma=cell(nSpecies,1);
    plasma.TK=plasma.T*Units.e/Units.kB;
end
if strcmpi(probe.type,'spherical'), probe_type=1;end
if strcmpi(probe.type,'cylindrical'), probe_type=2;end
if strcmpi(probe.type,'arbitrary'), probe_type=1;end

J_photo = -lp.photocurrent(probe.cross_section_area, U_probe, R_sun,probe.surface);
J_photo = J_photo .* UV_factor;
J_probe=J_photo; % initialize
for ii=1:nSpecies
    % density n
    q=plasma.q(ii);
    if numel(plasma.n)<nSpecies && ii > numel(plasma.n)
        n=plasma.n(end);
    else
        n=plasma.n(ii);
    end
    n=n*1e6; % to get density from cc to m3
    % temperature T
    if numel(plasma.TK)<nSpecies && ii > numel(plasma.TK)
        T=plasma.TK(end);
    else
        T=plasma.TK(ii);
    end
    % mass m
    if numel(plasma.m)<nSpecies && ii > numel(plasma.m)
        m=plasma.m(end);
    else
        m=plasma.m(ii);
    end
    if m==0 
        m=Units.me; 
    else
        m=Units.mp*m;
    end
    % velocity with respect to media
    if numel(plasma.vsc)<nSpecies && ii > numel(plasma.vsc)
        vsc=plasma.vsc(end);
    else
        vsc=plasma.vsc(ii);
    end
    J_thi = lp.thermal_current(probe_type,n,T,m,vsc,q,U_probe,probe.total_area);
    J_plasma{ii}=-sign(q)*J_thi; % positive current away from probe 
    J_probe=J_probe+J_plasma{ii};
end

