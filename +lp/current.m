function J=current(Probe,vectorU,rSunAU,factorUV,Plasma)
% LP.CURRENT calculate current to the probe
% J=LP.CURRENT(Probe,vectorU,rSunAU,factorUV,Plasma)
%
%   Calculates the total probe current to/from a Langmuir probe
%   consisting of cylinder/wire and sphere. Return current contributions:
%			J.total - total probe current
%     J.photo - photoeletron current
%    J.plasma - current of different plasma components
%
% Input:
%  Probe      - describes probe properties (structure or LP.LPROBE object)
%   Probe.Area.sunlit - total sunlit area [m^2]
%   Probe.Area.sphere - area of the sphere
%   Probe.Area.wire   - area of the wire
%   Probe.surface     - 'themis','cassini' (one of flags in lp.photocurrent)
%   Probe.areaSunlit  - use this if Area.sunlit not defined
%   Probe.areaTotal   - use this if Area.sphere and Area.wire are not
%                       defined, assume in this case sphere
%  vectorU    - probe potential (can be vector or matrix)
%  rSunAU     - distance from sun in AU
%  factorUV   - default is 1
%  Plasma     - describes plasma components (structure or LP.PLASMA object)
%    Plasma.qe - charge of species in e (the length of this vector corresponds to number of species)
%    Plasma.m  - mass of species in kg
%    Plasma.n  - density of species [m^-3]
%    Plasma.T  - temperature [eV]
%    Plasma.v  - velocity of probe wrt. media, zero if not given. [m/s]
%
% See also: LP.PHOTOCURRENT, LP.THERMAL_CURRENT

if isempty(Plasma) % calculate only photocurrent
  nPlasmaSpecies=0;
else
  nPlasmaSpecies=numel(Plasma.qe);
  J.plasma=cell(nPlasmaSpecies,1);
  %    plasma.TK=plasma.T*Units.e/Units.kB;
end

if isprop(Probe,'Area')
  areaSunlit = Probe.Area.sunlit;
else
  areaSunlit = Probe.areaSunlit;
end

if factorUV == 0
  J.photo = zeros(size(vectorU));
else
  J.photo = -lp.photocurrent(areaSunlit, vectorU, rSunAU,Probe.surfacePhotoemission);
  J.photo = J.photo .* factorUV;
end

J.total=J.photo; % initialize
for ii=1:nPlasmaSpecies
  % density n
  q=Plasma.qe(ii);
  if numel(Plasma.n)<nPlasmaSpecies && ii > numel(Plasma.n)
    n=Plasma.n(end);
  else
    n=Plasma.n(ii);
  end
  % temperature T
  if numel(Plasma.T)<nPlasmaSpecies && ii > numel(Plasma.T)
    T=Plasma.T(end);
  else
    T=Plasma.T(ii);
  end
  % mass m
  if numel(Plasma.m)<nPlasmaSpecies && ii > numel(Plasma.m)
    m=Plasma.m(end);
  else
    m=Plasma.m(ii);
  end
  % velocity with respect to media
  if numel(Plasma.v)<nPlasmaSpecies && ii > numel(Plasma.v)
    v=Plasma.v(end);
  else
    v=Plasma.v(ii);
  end
  J.plasma{ii}=thermal_current(Probe,n,T,m,v,q,vectorU);
  J.total=J.total+J.plasma{ii};
end
end

function jThermal = thermal_current(Lprobe,n,T,m,vsc,q,vectorU)
% calculates thermal current to Langmuir probe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function j_thermal = lp.thermal_current( p_type, N, T, m, V, Z, U, A )
%
%   Matlab function that calculates the thermal probe current to/from
%   a cylindrical or spherical body, e.g. a Langmuir probe or the a
%   spherical (cylindrical) S/C.
%
%   Input parameters:  N,T,m,Z    =  #density[m^-3], temperature[eV], mass[kg]
%                                    and charge [e] of current carrying
%                                    species.
%                      V          =  velocity of the body with respect
%                                    to the plasma [m/s].
%                      U          =  body potential [V]
%                      A          =  area of body [m^2]
%                      p_type       = spherical   (1) or
%                                   cylindrical (2).
%
%                      (all scalars, U may be a vector)

Units=irf_units;

% Initialize.
%%%%%%%%%%%%%
jThermal       = zeros(size(vectorU));
jThermalSphere = jThermal;
jThermalWire   = jThermal;

% If zero density return zero current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n==0,    return;end

% Is the body moving with a velocity, V, with
% respect to the plasma ?
% Criteria set such that it is considered
% important if V > 0.1 * V_th.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if vsc < 0.1 * sqrt(Units.e*T/m)

  % Ratio of potential to thermal energy.
  X = vectorU/T;

  % Total current to/from body.
  fluxIp = n*Units.e*sqrt( T*Units.e/(2.0*pi*m) );

else
  X = ( Units.e / (m*vsc^2/2 + Units.e*T) ) .* vectorU;
  fluxIp = n*Units.e*sqrt( vsc^2/16 + T*Units.e/(2.0*pi*m) );
end


indPositiveU = find( vectorU >= 0 );
indNegativeU = find( vectorU < 0 );

% Spherical body case.
%%%%%%%%%%%%%%%%%%%%%%
if isprop(Lprobe,'Area')
  A = Lprobe.Area.sphere;
else
  A = Lprobe.areaTotal;
end

Ip = A*fluxIp;

if q > 0
  jThermalSphere(indPositiveU) = Ip .* exp(-X(indPositiveU));
  jThermalSphere(indNegativeU) = Ip .* (1-X(indNegativeU));
elseif q < 0
  jThermalSphere(indPositiveU) = Ip .* (1+X(indPositiveU));
  jThermalSphere(indNegativeU) = Ip .* exp(X(indNegativeU));
end


% Cylindrical body case.
%%%%%%%%%%%%%%%%%%%%%%%%

if isprop(Lprobe,'Area')
  A = Lprobe.Area.wire;
  Ip = A*fluxIp;

  sq         = zeros(size(vectorU));
  %     erfv       = zeros( U_pts, 1 );

  sq(indNegativeU) = sqrt( abs(-X(indNegativeU)) );
  sq(indPositiveU) = sqrt( abs(+X(indPositiveU)) );
  erfv = erf( sq );

  if q > 0
    jThermalWire(indPositiveU) = Ip .* exp(-X(indPositiveU));
    jThermalWire(indNegativeU) = Ip .* ( (2/sqrt(pi)) .* sq(indNegativeU) ...
      + exp(-X(indNegativeU)) .* (1.0 - erfv(indNegativeU)) );
  elseif q < 0
    jThermalWire(indNegativeU) = Ip .* exp(X(indNegativeU));
    jThermalWire(indPositiveU) = Ip .* ( (2.0/sqrt(pi)) .* sq(indPositiveU) ...
      + exp(+X(indPositiveU)) .* (1.0 - erfv(indPositiveU)) );
  end
else
  jThermalWire=0;
end

jThermal = jThermalSphere + jThermalWire;
jThermal = -q*jThermal; % positive current is away from probe

end
