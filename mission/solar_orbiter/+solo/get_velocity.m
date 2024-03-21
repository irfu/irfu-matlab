function soloVelocity = get_velocity(Tint, varargin)
%SOLO.GET_VELOCITY  Get the velocity of Solar Orbiter
%
% soloVelocity = solo.get_velocity(tint, PARAMS)
%
%
% A function to get the Velocity of Solar Orbiter for a given time interval
% Output is a TSeries object with velocity in the data variable "soloVelocity",
% if any positons found for the time interval "tint".
%
% Options:
%    'predicted' - get predicted position instead of flown
%    'frame' - frame of the data (see below, default: ECLIPJ2000)
%
% The following generic frames are defined:
%
%      SPICE Frame Name            Long-name
%      -------------------------   --------------------------------------------
%
%   SOLO mission specific generic frames:
%
%      SOLO_SUN_RTN                Sun Solar Orbiter Radial-Tangential-Normal
%      SOLO_SOLAR_MHP              S/C-centred mirror helioprojective
%      SOLO_IAU_SUN_2009           Sun Body-Fixed based on IAU 2009 report
%      SOLO_IAU_SUN_2003           Sun Body-Fixed based on IAU 2003 report
%      SOLO_GAE                    Geocentric Aries Ecliptic at J2000 (GAE)
%      SOLO_GSE                    Geocentric Solar Ecliptic at J2000 (GSE)
%      SOLO_HEE                    Heliocentric Earth Ecliptic at J2000 (HEE)
%
%   Heliospheric Coordinate Frames developed for the NASA STEREO mission:
%
%      SOLO_ECLIPDATE              Mean Ecliptic of Date Frame
%      SOLO_HCI                    Heliocentric Inertial Frame
%      SOLO_HEE_NASA               Heliocentric Earth Ecliptic Frame
%      SOLO_HEEQ                   Heliocentric Earth Equatorial Frame
%      SOLO_HEEQ                   Heliocentric Earth Equatorial Frame
%      SOLO_GEORTN                 Geocentric Radial Tangential Normal Frame
%
%   Heliocentric Generic Frames(*):
%
%      SUN_ARIES_ECL               Heliocentric Aries Ecliptic   (HAE)
%      SUN_EARTH_CEQU              Heliocentric Earth Equatorial (HEEQ)
%      SUN_EARTH_ECL               Heliocentric Earth Ecliptic   (HEE)
%      SUN_INERTIAL                Heliocentric Inertial         (HCI)
%
%   Geocentric Generic Frames:
%
%      EARTH_SUN_ECL   (*)         Geocentric Solar Ecliptic     (GSE)
%      EARTH_MECL_MEQX (*)         Earth Mean Ecliptic and Equinox of date
%                                  frame (Auxiliary frame for EARTH_SUN_ECL)
%      EARTH_MECL_MEQX_J2000       Earth Mean Ecliptic and Equinox at J2000
%                                  frame (Auxiliary frame for SOLO_GSE and
%                                  SOLO_HEE)
%
% Example:
% tint = irf.tint('2020-02-11T00:00:00Z/2020-08-01T15:00:00Z'); % time interval in TT2000 UTC
% %Get velocity in GSE
% soloVelocity = solo.get_velocity(Tint,'frame','SOLO_GSE');
% irf_plot(soloVelocity, '.');


frame = 'ECLIPJ2000';
flagFlownPredicted = 'flown';

if nargin > 1, have_options = 1; args = varargin;
else, have_options = 0;
end
while have_options
  l = 1;
  switch lower(args{1})
    case 'frame'
      frame = args{2}; l=2;
    case 'flown'
      flagFlownPredicted = 'flown';
    case 'predicted'
      flagFlownPredicted = 'predicted';
    otherwise
      irf.log('warning', ['Option ''' args{1} ''' not recognized']);
  end
  if length(args) > l, args = args(l+1:end);
  else, break
  end
end

t_step = 60*60; %3600 sec

solo.db_get_metakernel(flagFlownPredicted);

% Compute et (SPICE ephemeries time, make use of input in TT2000)
et = Tint.start.tts:t_step:Tint.stop.tts;

% The position of Solar Orbiter in units of km
posvel = cspice_spkezr('solo', et, frame, 'LT+s', 'Sun'); %This gives both position and velocity
% More frames (instead of 'ECLIPJ2000') described in/share/SPICE/Solar-Orbiter/kernels/fk/solo_ANC_soc-sci-fk_V06.tf

% Convert to utc and then to TT2000 in TSeries object
utc_tmp = cspice_et2utc(et, 'ISOC', 0);

% Note pos' since it is returned as 3xN but TSeries expects Nx3 (where N is number of records).
soloVelocity= irf.ts_vec_xyz(EpochTT(utc_tmp), posvel(4:6,:)');
soloVelocity.units = 'km/s'; % Add some metadata information (read when plotting in irf_plot)
soloVelocity.coordinateSystem=frame;
end
