%
% Read solar panel/array angles from ROC's SolO HK XML files.
%
%
% RETURN VALUES
% =============
% Py, My
%       Structs with arrays for the respective time series.
%
%
% DOCUMENTATION /COMMENTS FROM ROC
% ================================
% """"""""
% Solar array angle
%
% # NCFT55P0 SADE-A calibrated PY array position (deg)
% # NCFT55U0 SADE-A calibrated MY array position (deg)
% # And equivalently NCFT56U0 and NCFT56Y0 if we are ever using the SADE-B
% # 4) APR temperature
%
% # This is a difficult one. Attachments below, including the list of PCDU
% thermistors from Airbus. You take your chances with these. The external URP is
% at least at the correct end of the PCDU, but it's external. The internal ones
% are (as far as I can see now) not close to the APR end.
% """""""" /https://confluence-lesia.obspm.fr/display/ROC/SOLO+HK+Parameter+data
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-09-07.
%
function [Py, My] = read_RSH_solar_array_angles(xmlFilePathsCa)
% PROPOSAL: Rename ~SADE-A.
%   PRO: Shorter.
%   PRO: More precise.
%   CON: Bad name if later have to use both SADE-A and SADE-B (e.g.
%        different for different time intervals, after one has failed).
%
% PROPOSAL: Refactor into reading arbitrary element tags
%           (e.g. EngineeringValue) for arbitrary list of mnemonics.


% HK mnemonics
% ------------
% https://confluence-lesia.obspm.fr/display/ROC/SOLO+HK+Parameter+data
% NCFT55P0 = SADE-A calibrated PY array position (deg)
% NCFT55U0 = SADE-A calibrated MY array position (deg)
% (SADE-A = Solar Array Drive Electronics, unit A.)
NAME_1 = 'NCFT55P0';
NAME_2 = 'NCFT55U0';



%t = tic();
D = solo.shk.read_RSH_file_many(xmlFilePathsCa, {NAME_1, NAME_2});
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
b1 = strcmp(D.Name, NAME_1);
b2 = strcmp(D.Name, NAME_2);

Py.Dt = Dt(b1);
My.Dt = Dt(b2);
Py.solarArrayAngleDeg = str2double(D.EngineeringValue(b1));
My.solarArrayAngleDeg = str2double(D.EngineeringValue(b2));
end
