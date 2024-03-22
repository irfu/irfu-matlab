%
% Wrapper around solo.qli.generate_quicklooks_24h_6h_2h(). The function obtains
% science data using solo.db_get_ts() ("DB"; indirectly) and geometry
% information using SPICE (indirectly).
%
%
% IMPLEMENTATION NOTE: This function is separate from
% solo.qli.generate_quicklooks_24h_6h_2h() (1) partly for historical reasons,
% and (2) partly to facilitate test code for
% solo.qli.generate_quicklooks_24h_6h_2h() which requires explicitly
% submitting all data as arguments to it (i.e. NOT use e.g. solo.db_get_ts() or
% SPICE which both read files on disk).
%
%
% ARGUMENTS
% =========
% Dt
%     UTC Datetime object to the midnight that begins the day for which
%     quicklooks should be made.
% vhtFile1hPath
%     Path to file with 1h VHT data. Typically found at
%     brain:/data/solo/data_yuri/V_RPW_1h.mat .
% OutputPaths
%     Struct with paths to output directories for different types of plots.
%     Has fields: .dir2h, .dir6h, .dir24h
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
function generate_quicklooks_24h_6h_2h_using_DB_SPICE(Dt, vhtFile1hPath, OutputPaths, irfLogoPath)
Tint = [
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt), ...
  solo.qli.utils.scalar_datetime_to_EpochTT(Dt+caldays(1))
  ];
solo.qli.utils.log_plot_function_time_interval(Tint)



Data = [];

Vht1h       = load(vhtFile1hPath);
Data.Vrpw   = Vht1h.V_RPW_1h.tlim(Tint);
% E-field
Data.E      = solo.qli.utils.db_get_ts('solo_L3_rpw-bia-efield-10-seconds-cdag', 'EDC_SRF', Tint);
% RPW density
Data.Ne     = solo.qli.utils.db_get_ts('solo_L3_rpw-bia-density-10-seconds-cdag', 'DENSITY', Tint);
% B-field
Data.B      = solo.qli.utils.db_get_ts('solo_L2_mag-rtn-normal', 'B_RTN', Tint);
% Proton & alpha temperature
Data.Tpas   = solo.qli.utils.db_get_ts('solo_L2_swa-pas-grnd-mom', 'T', Tint);
% Proton & alpha velocity
Data.Vpas   = solo.qli.utils.db_get_ts('solo_L2_swa-pas-grnd-mom', 'V_RTN', Tint);
% Proton & alpha density
Data.Npas   = solo.qli.utils.db_get_ts('solo_L2_swa-pas-grnd-mom', 'N', Tint);
% Ion spectrum
Data.ieflux = solo.qli.utils.db_get_ts('solo_L2_swa-pas-eflux', 'eflux', Tint);
Data.swaEnergyMetadata = solo.qli.utils.read_constant_metadata('solo_L2_swa-pas-eflux', 'Energy', Tint);

% TNR_BAND
% --------
% NOTE: Variable is not used very much. solo.qli.generate_quicklooks_24h_6h_2h()
% only checks if empty or not (sic!).
%
%      FIELDNAM        (CDF_CHAR/8): "TNR_BAND"
%      CATDESC         (CDF_CHAR/31): "TNR band of the current record "
%      VAR_NOTES       (CDF_CHAR/71): "TNR band of the current record. Possible values are: 1=A, 2=B, 3=C, 4=D"
% /solo_L2_rpw-tnr-surv-cdag_20240101_V02.cdf
Data.tnrBand = solo.qli.utils.db_get_ts('solo_L2_rpw-tnr-surv-cdag', 'TNR_BAND', Tint);

Data.Tnr = solo.read_TNR(Tint);



% Solar Orbiter position
Data.soloPos = solo.qli.utils.get_SolO_position(Tint);
% Earth position (also uses SPICE)
EARTH_POS_DT_SEC = 60*60;
Data.earthPos    = solo.qli.utils.get_Earth_position(Tint, EARTH_POS_DT_SEC);

if ~solo.qli.const.ENABLE_B
  Data.B = [];
end

% =============================
% Plot data and save quicklooks
% =============================
solo.qli.generate_quicklooks_24h_6h_2h(Data, OutputPaths, Tint, irfLogoPath)
end



% IMPLEMENTATION NOTE: Implementation is made to reproduce old behaviour
% (copy-paste). Can likely be cleaned up/simplified together with
% solo.read_TNR().
%
% function Tnr = read_TNR(Tint)
%
% try
%   Tnr = solo.read_TNR(Tint);
% catch Exc
%   if strcmp(Exc.identifier, 'read_TNR:FileNotFound')
%     Tnr = [];
%   end
% end
%
% end