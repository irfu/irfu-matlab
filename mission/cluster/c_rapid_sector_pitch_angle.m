function pitch_angles=c_rapid_sector_pitch_angle(diB,ic)
%pitch_angles=c_rapid_sector_pitch_angle(diB,ic);
%
% B - magnetic field sampled at times for which to calculate pitch angles
% in DSI ref frame
% ic - s/c number
% 
% pitch_anglse - matrix with 10 columns. 1 column time and others are pitch angles for sectors 1..9

narginchk(0,6)
if nargin==0, help c_rapid_sector_pitch_angle;return;end

t=diB(:,1);

a=c_load('A?',ic,'var'); 
phase=irf_resamp([a(:,1) unwrap(a(:,2)/180*pi)],t);
phase=mod(phase(:,2)*180/pi,360); % take away time column

phase_rapid=phase/180*pi + 60.167/180*pi; %#ok<NASGU> % rapid phase   

polar_angles=(170:-20:10)*pi/180; %#ok<NASGU>
sector_numbers=1:9;

c_eval('sector_?_dsi=[t sin(polar_angles(?))*cos(phase_rapid) -sin(polar_angles(?))*sin(phase_rapid) cos(polar_angles(?)*ones(size(t)))];',1:9);
c_eval('sector_?_dsi=irf_resamp(sector_?_dsi,diB);',1:9);
c_eval('sector_?_pitch_angle=irf_ang(sector_?_dsi,diB,1);',1:9)

pitch_angles=sector_1_pitch_angle;
c_eval('pitch_angles=[pitch_angles 180-sector_?_pitch_angle(:,2)];',2:9);

