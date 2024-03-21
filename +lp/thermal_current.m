function [j_thermal] = thermal_current( p_type, N, T, m, V, Z, U, A )
% LP.THERMAL_CURRENT calculates thermal current to Langmuir probe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [j_thermal] = lp.thermal_current( p_type, N, T, m, V, Z, U, A )
%
%   Matlab function that calculates the thermal probe current to/from
%   a cylindrical or spherical body, e.g. a Langmuir probe or the a
%   spherical (cylindrical) S/C.
%
%   Input parameters:  N,T,m,Z    =  #density[m^-3], temperature[K], mass[kg]
%                                    and charge [+/-] of current carrying
%                                    species.
%                      V          =  velocity of the body with respect
%                                    to the plasma [m/s].
%                      U          =  body potential [V]
%                      A          =  area of body [m^2]
%                      p_type       = spherical   (1) or
%                                   cylindrical (2).
%
%                      (all scalars, U may be a vector)
%
%   Created by Jan-Erik Wahlund, Cornell University, October-1994.
%   Modified for use with isdat_2.6, J-E. Wahlund, IRF-Uppsala, 1999.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check # input/output parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
narginchk(8,8)
nargoutchk(1,1)


% Globals.
%%%%%%%%%%
Units=irf_units;

% Initialize.
%%%%%%%%%%%%%
j_thermal    = zeros(size(U));

% If zero density return zero current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if N==0 || T==0,    return;end

% Is the body moving with a velocity, V, with
% respect to the plasma ?
% Criteria set such that it is considered
% important if V > 0.1 * V_th.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if V < 0.1 * sqrt(Units.kB*T/m)

  % Ratio of potential to thermal energy.
  X = ( Units.e / (Units.kB*T) ) .* U;

  % Total current to/from body.
  Ip = A*N*Units.e*sqrt( T*Units.kB/(2.0*pi*m) );

else
  X = ( Units.e / (m*V^2/2 + Units.kB*T) ) .* U;
  Ip = A*N*Units.e*sqrt( V^2/16 + T*Units.kB/(2.0*pi*m) );
end


% Spherical body case.
%%%%%%%%%%%%%%%%%%%%%%
if p_type == 1

  pos_ind = find( U >= 0 );
  neg_ind = find( U < 0 );

  if Z > 0
    j_thermal(pos_ind) = Ip .* exp(-X(pos_ind));
    j_thermal(neg_ind) = Ip .* (1-X(neg_ind));
  elseif Z < 0
    j_thermal(pos_ind) = Ip .* (1+X(pos_ind));
    j_thermal(neg_ind) = Ip .* exp(X(neg_ind));
  end


  % Cylindrical body case.
  %%%%%%%%%%%%%%%%%%%%%%%%
elseif p_type == 2

  pos_ind = find( U >= 0 );
  neg_ind = find( U < 0 );

  sq         = zeros(size(U));
  %     erfv       = zeros( U_pts, 1 );

  sq(neg_ind) = sqrt( abs(-X(neg_ind)) );
  sq(pos_ind) = sqrt( abs(+X(pos_ind)) );
  erfv = erf( sq );

  if Z > 0
    j_thermal(pos_ind) = Ip .* exp(-X(pos_ind));
    j_thermal(neg_ind) = Ip .* ( (2/sqrt(pi)) .* sq(neg_ind) ...
      + exp(-X(neg_ind)) .* (1.0 - erfv(neg_ind)) );
  elseif Z < 0
    j_thermal(neg_ind) = Ip .* exp(X(neg_ind));
    j_thermal(pos_ind) = Ip .* ( (2.0/sqrt(pi)) .* sq(pos_ind) ...
      + exp(+X(pos_ind)) .* (1.0 - erfv(pos_ind)) );
  end

else
  disp('This probe type is not supported yet !');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



