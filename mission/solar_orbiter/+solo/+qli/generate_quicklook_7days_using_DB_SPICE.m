%
% Wrapper around solo.qli.generate_quicklook_7days(). The function obtains
% science data using solo.db_get_ts() ("DB"; indirectly) and geometry
% information using SPICE (indirectly).
%
%
% IMPLEMENTATION NOTE: This function is separate from
% solo.qli.generate_quicklook_7days() (1) partly for historical reasons, and (2)
% partly to facilitate (future) test code for
% solo.qli.generate_quicklook_7days() which then requires explicitly submitting
% all data as arguments to it (i.e. NOT use e.g. solo.db_get_ts(), or SPICE).
%
%
% ARGUMENTS
% =========
% Dt
%     UTC Datetime object to the midnight that begins a week.
% vhtFile6hPath
%     Path to file with 6h VHT data. Typically found at
%     brain:/data/solo/data_yuri/V_RPW.mat .
% outputDir1wPath
% irfLogoPath
%     Path to IRF logo. Empty ==> Do not plot any logo.
%     IMPLEMENTATION NOTE: The IRF logo should only be used for official
%     quicklooks.
%
%
% RETURN VALUES
% =============
% (None)
%
%
% Author: Based on code by Konrad Steinvall, IRF, Uppsala, Sweden. Modified by
% Erik P G Johansson, IRF, Uppsala, Sweden.
%
function generate_quicklook_7days_using_DB_SPICE(Dt, vhtFile6hPath, outputDir1wPath, irfLogoPath)
Tint = [
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt), ...
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt+caldays(7)), ...
  ];
solo.qli.utils.log_plot_function_time_interval(Tint)



Vht6h = load(vhtFile6hPath);

Data = [];

Data.Vrpw   = Vht6h.V_RPW.tlim(Tint);
% E-field:
Data.E      = solo.qli.utils.db_get_ts('solo_L3_rpw-bia-efield-10-seconds-cdag', 'EDC_SRF', Tint);
% RPW density:
Data.Ne     = solo.qli.utils.db_get_ts('solo_L3_rpw-bia-density-10-seconds-cdag', 'DENSITY', Tint);
% B-field:
Data.B      = solo.qli.utils.db_get_ts('solo_L2_mag-rtn-normal-1-minute', 'B_RTN', Tint);
% Proton & alpha temperature:
Data.Tpas   = solo.qli.utils.db_get_ts('solo_L2_swa-pas-grnd-mom', 'T', Tint);
% Proton & alpha velocity:
Data.Vpas   = solo.qli.utils.db_get_ts('solo_L2_swa-pas-grnd-mom', 'V_RTN', Tint);
% Proton & alpha density:
Data.Npas   = solo.qli.utils.db_get_ts('solo_L2_swa-pas-grnd-mom', 'N', Tint);
% Ion spectrum
Data.ieflux = solo.db_get_ts(          'solo_L2_swa-pas-eflux', 'eflux', Tint);
% TNR E-field
Data.Etnr   = solo.db_get_ts(          'solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);
% Solar Orbiter position
Data.soloPos = solo.qli.utils.get_SolO_position(Tint);

% Earth position (also uses SPICE)
EARTH_POS_DT_SEC = 60*60;
Data.earthPos    = solo.qli.utils.get_Earth_position(Tint, EARTH_POS_DT_SEC);

% ============================
% Plot data and save quicklook
% ============================
solo.qli.generate_quicklook_7days(Data, outputDir1wPath, Tint, irfLogoPath)
end



