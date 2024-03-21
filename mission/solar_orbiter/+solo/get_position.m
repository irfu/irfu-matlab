function soloPosition = get_position(Tint, varargin)
%SOLO.GET_POSITION  Get the position of Solar Orbiter
%
% soloPosition = solo.get_position(tint, PARAMS)
%
% A function to get the position of Solar Orbiter for a given time interval
% Output is a TSeries object with position in the data variable "pos",
% if any positons found for the time interval "tint".
%
% Note: Always returns coordinates with a sampling rate of one position per hour
%       (hardcoded).
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
%
% Tint = irf.tint('2020-07-31T00:00:00Z/2020-08-01T15:00:00Z'); % time interval in TT2000 UTC
% %Get predicted position in SOLO_GSE frame
% soloPosition = solo.get_position(Tint, 'predicted','frame','SOLO_GSE');
% irf_plot(soloPosition, '.');


frame = 'ECLIPJ2000';
flagFlownPredicted = 'flown';

if nargin > 1, has_options = 1; args = varargin;
else, has_options = 0;
end
while has_options
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

T_STEP_SEC = 60*60; % 3600 s

solo.db_get_metakernel(flagFlownPredicted);

% Compute et (SPICE ephemeries time, make use of input in TT2000)
et = Tint.start.tts:T_STEP_SEC:Tint.stop.tts;

% The position of Solar Orbiter in units of km
pos = cspice_spkpos('solo', et, frame, 'LT+s', 'Sun');
% More frames (instead of 'ECLIPJ2000') described in/share/SPICE/Solar-Orbiter/kernels/fk/solo_ANC_soc-sci-fk_V06.tf

% Convert to UTC and then to TT2000 in TSeries object
utc_tmp = cspice_et2utc(et, 'ISOC', 0);

% Note pos' since it is returned as 3xN but TSeries expects Nx3 (where N is number of records).
soloPosition = irf.ts_vec_xyz(EpochTT(utc_tmp), pos');
soloPosition.units = 'km'; % Add some metadata information (read when plotting in irf_plot)
soloPosition.coordinateSystem = frame;
end
