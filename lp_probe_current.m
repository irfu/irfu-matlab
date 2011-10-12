function [J_probe, J_photo, J_plasma]=lp_probe_current(probe_type, XA,Ap,U_probe,R_sun,UV_factor,plasma)
% LP_PROBE_CURRENT calculate current to the probe
% J_probe=LP_PROBE_CURRENT(probe_type, XA,Ap,U_probe,R_sun,UV_factor,plasma)
%
%   Calculates the total probe current to/from 
%   a cylindrical or spherical Langmuir probe.
%
% Input:
%  probe_type - spherical(1) or cylindrical (2).
%  XA         - cross section area
%  Ap         - probe area
%  U_probe    - probe potential (can be vector) 
%  R_sun      - distance from sun in AU
%  UV_factor  - default is 1
%  plasma     - describes plasma components (structure)
%    plasma.q - charge of species in e (the length of this vector corresponds to number of species)
%    plasma.m - mass of species in proton masses (0 corresponds to e- mass)
%    plasma.n - density of species [cc]
%    plasma.T - temperature [eV]
%    plasma.vsc - velocity of probe wrt. mmedia [m/s]
%
% [J_probe, J_photo, J_plasma]=LP_PROBE_CURRENT
%   Return current contributions from photoeletrons and all
%   the plasma components
%
% See also: LP_PHOTOCURRENT, LP_THERMAL_CURRENT

irf_units;

n_of_species=numel(plasma.q);
J_plasma=cell(n_of_species,1);
plasma.TK=plasma.T*Units.e/Units.kB;

J_photo = -lp_photocurrent( XA, U_probe, R_sun,'themis' );
J_photo = J_photo .* UV_factor;
J_probe=J_photo; % initialize
for ii=1:n_of_species,
    % density n
    q=plasma.q(ii);
    if numel(plasma.n)<n_of_species && ii > numel(plasma.n)
        n=plasma.n(end);
    else
        n=plasma.n(ii);
    end
    n=n*1e6; % to get density from cc to m3
    % temperature T
    if numel(plasma.TK)<n_of_species && ii > numel(plasma.TK)
        T=plasma.TK(end);
    else
        T=plasma.TK(ii);
    end
    % mass m
    if numel(plasma.m)<n_of_species && ii > numel(plasma.m)
        m=plasma.m(end);
    else
        m=plasma.m(ii);
    end
    if m==0, 
        m=Units.me; 
    else
        m=Units.mp*m;
    end
    % velocity with respect to media
    if numel(plasma.vsc)<n_of_species && ii > numel(plasma.vsc)
        vsc=plasma.vsc(end);
    else
        vsc=plasma.vsc(ii);
    end
    J_thi = lp_thermal_current( probe_type,n,T,m,vsc,q,U_probe,Ap);
    J_plasma{ii}=-sign(q)*J_thi; % positive current away from probe 
    J_probe=J_probe+J_plasma{ii};
end

