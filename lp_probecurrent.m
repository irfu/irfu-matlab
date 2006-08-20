%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lp_probe_current
%
%   Matlab script that calculates the *total* probe current to/from 
%   a cylindrical or spherical Langmuir probe.
%   (all scalars, U may be a vector)
%  
%   For the time being, U_SC is not used.
%
%   Created by Jan-Erik Wahlund, Cornell University, October-1994.
%   Modified for use in isdat_2.6, J-E. Wahlund, IRF-Uppsala, 1999.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize.
%%%%%%%%%%%%%
  U_pts  = length( U_probe );

  J_the      = zeros(U_pts, 1);
  J_thi1     = zeros(U_pts, 1);
  J_thi2     = zeros(U_pts, 1);
  J_photo    = zeros(U_pts, 1);
  J_SC       = zeros(U_pts, 1);


% Photo-current.
%%%%%%%%%%%%%%%%

  % Probe photo-electrons.
  %%%%%%%%%%%%%%%%%%%%%%%%
  J_photo = lp_photocurrent( XA, U_probe, R_sun );

  % Magical photoelectron calibration constants, which 
  % probably are dependent on probe surface material
  % as well as the Solar UV radiation flux intensity.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  J_photo = J_photo .* UV_factor;


  % Rolf's S/C photo-electron cloud.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Nph     = 50e6;
  Wph     = 2.74;         % Typical photoelectron temperature.
  Tph     = Wph * 11605;
  Z       = -1;
  %U_cloud = 15;  % [V]
  U_cloud = 0;  % [V]
  U_eff   = U_probe - U_cloud;
  %J_SC = ThermalCurrent( probe_type, Nph, Tph, me, 0, Z, U_eff, Ap );
  J_SC = 0;


% Ion thermal current. (H+)
%%%%%%%%%%%%%%%%%%%%%%
  Z      = +1;
  mi1    = m_amu1 * mH;
  Ni1    = (1 - m2_fraction) * Ne;
  J_thi1 = lp_thermal_current( probe_type, Ni1, Ti, mi1, V_SC, Z, U_probe, Ap );


% Ion thermal current. (O+)
%%%%%%%%%%%%%%%%%%%%%%
  Z      = +1;
  mi2    = m_amu2 * mH;
  Ni2    = m2_fraction * Ne;
  J_thi2 = lp_thermal_current( probe_type, Ni2, Ti, mi2, V_SC, Z, U_probe, Ap );


% Electron thermal current.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Z = -1;
  J_the = lp_thermal_current( probe_type, Ne, Te, me, V_SC, Z, U_probe, Ap );


% Sum up !
%%%%%%%%%%
  J_probe = J_the - J_thi1 - J_thi2 - J_photo + J_SC; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

