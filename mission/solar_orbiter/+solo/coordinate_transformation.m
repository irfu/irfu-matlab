function TSout = coordinate_transformation(TSin,from,to)
%SOLO.coordinate_transformtion  Transfrom between coordinate system 'from'
%to coordinate system 'to'
%supported coordinate systems: 'VSO' , 'SRF', 'SUN_RTN','SOLAR_MHP','GSE','HEE','HCI','GEORTN','HEEQ'

%
%
%      SOLO_VSO                    Venus solar orbiter
%      SOLO_SRF                    S/C-reference frame
%      SOLO_SUN_RTN                Sun Solar Orbiter Radial-Tangential-Normal
%      SOLO_SOLAR_MHP              S/C-centred mirror helioprojective
%      SOLO_GSE                    Geocentric Solar Ecliptic at J2000 (GSE)
%      SOLO_HEE                    Heliocentric Earth Ecliptic at J2000 (HEE)
%      SOLO_HCI                    Heliocentric Inertial Frame
%      SOLO_HEEQ                   Heliocentric Earth Equatorial Frame
%      SOLO_GEORTN                 Geocentric Radial Tangential Normal Frame


if nargin<3, error('This function takes input TSin, "from" and "to" coordinate systems'); end

CS={'VSO','SUN_RTN','SRF','SOLAR_MHP','GSE','HEE','HCI','GEORTN','HEEQ'};
if ~max(strcmpi(from,CS)), error([upper(from) ' coordinate system is not supported, see help for more info']);end
if ~max(strcmpi(to,CS)), error([upper(to) ' coordinate system is not supported, see help for more info']);end

solo.db_get_metakernel('flown');

% Compute et (SPICE ephemeries time, make use of input in TT2000)
et=irf_time(TSin.time,'EpochTT>tt')';

% The rotation matrix M (note that this is time dependent so it should be a
% 3x3xlength(time) matrix)

M=cspice_pxform(['SOLO_' upper(from)],['SOLO_' upper(to)],et);



out = zeros(size(TSin.data));
for idx=1:3
    out (:,idx) = ...
        squeeze(M(idx,1,:)).*TSin.data(:,1) + ...
        squeeze(M(idx,2,:)).*TSin.data(:,2) + ...
        squeeze(M(idx,3,:)).*TSin.data(:,3);
end

TSout=irf.ts_vec_xyz(TSin.time,out);
TSout.units = TSin.units;
TSout.siConversion = TSin.siConversion;

end


