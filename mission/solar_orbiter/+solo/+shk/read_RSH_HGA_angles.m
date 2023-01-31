%
% Read HGA (High Gain Antenna) angles from ROC's SolO HK XML files.
%
%
% DOCUMENTATION /COMMENTS FROM ROC
% ================================
% NCFT29T0 = "HGA Acquired Elevation in degrees"
% NCFT29S0 = "HGA Acquired Azimuth in degrees"
% """"""""
% HGA angles
%
% They are calibrated in degrees.
% The frame used is a bit weird.
%          -90 deg elevation is pointing in the spacecraft XY plane (i.e.
%          "spacecraft equatorial").
%
%
% Imagine starting from the stowed position (and recall the HGA is mounted on
% the minus -Z panel, on the "bottom" of the Spacecraft). There's a 90 deg
% rotation to deploy the antenna to a typical "equatorial" position. Then zero
% elevation corresponds to pointing the beam in the direction of SC -X, and the
% elevation motor is upside down, so motor elevation is back to front wrt
% right-handed SC +Z rotation.
%
% Azimuth is positive around -Z
% """""""" /https://confluence-lesia.obspm.fr/display/ROC/SOLO+HK+Parameter+data
%
%
% RETURN VALUES
% =============
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-10-12.
%
function [HgaElevation, HgaAzimuth] = read_RSH_HGA_angles(xmlFilePathsCa)
    % PROPOSAL: Automatic test code.
    %
    % PROPOSAL: Refactor into reading arbitrary element tags
    %           (e.g. EngineeringValue) for arbitrary list of mnemonics.
    % PROPOSAL: Save ~mnemonics in data structures.
    %   CON: Structs are no longer same size in all fields.


    % HK mnemonics
    % ------------
    % https://confluence-lesia.obspm.fr/display/ROC/SOLO+HK+Parameter+data
    % NCFT29T0 = HGA Acquired Elevation in degrees
	% NCFT29S0 = HGA Acquired Azimuth in degrees
    NAME_ELEVATION = 'NCFT29T0';
    NAME_AZIMUTH   = 'NCFT29S0';



    %t = tic();
    D = solo.shk.read_RSH_file_many(xmlFilePathsCa, {NAME_ELEVATION, NAME_AZIMUTH});
    %toc(t)    % TEMP

    Dt = datetime(...
        D.TimeStampAsciiA, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');

    % Sort arrays by time.
    % IMPLEMENTATION NOTE: Does this before separating arrays.
    [Dt, iSort] = sort(Dt);
    D.Name             = D.Name(iSort);
    D.EngineeringValue = D.EngineeringValue(iSort);

    % Sort data by HK Name.
    bElevation = strcmp(D.Name, NAME_ELEVATION);
    bAzimuth   = strcmp(D.Name, NAME_AZIMUTH);

    HgaElevation.Dt = Dt(bElevation);
    HgaAzimuth.Dt   = Dt(bAzimuth);
    HgaElevation.angleDeg = str2double(D.EngineeringValue(bElevation));
    HgaAzimuth.angleDeg   = str2double(D.EngineeringValue(bAzimuth));
end
