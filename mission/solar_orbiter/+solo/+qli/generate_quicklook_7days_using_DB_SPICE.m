%
% Wrapper around solo.qli.generate_quicklook_7days(). The function obtains
% science data using solo.db_get_ts() ("DB"; indirectly) and geometry
% information using SPICE (indirectly).
%
%
% IMPLEMENTATION NOTE: This function is separate from
% solo.qli.generate_quicklook_7days() (1) partly for historical reasons, and (2)
% partly to facilitate test code for solo.qli.generate_quicklook_7days() which
% requires explicitly submitting all data as arguments to it (i.e. NOT use
% e.g. solo.db_get_ts() or SPICE which both read files from disk).
%
%
% ARGUMENTS
% =========
% Dt
%     UTC Datetime object to the midnight that begins the 7-day period for which
%     a quicklook should be made.
% vhtFile6hPath
%     Path to file with 6h VHT data. Typically found at
%     brain:/data/solo/data_yuri/V_RPW.mat .
% outputDir1wPath
%     Path to output directory.
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



Data = [];

Vht6h       = load(vhtFile6hPath);
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

Data.Tnr     = solo.read_TNR(Tint);



% Solar Orbiter position
Data.soloPos = solo.qli.utils.get_SolO_position(Tint);
% Earth position (also uses SPICE)
EARTH_POS_DT_SEC = 60*60;
Data.earthPos    = solo.qli.utils.get_Earth_position(Tint, EARTH_POS_DT_SEC);

if ~solo.qli.const.ENABLE_B
  Data.B = [];
end

% ============================
% Plot data and save quicklook
% ============================
solo.qli.generate_quicklook_7days(Data, outputDir1wPath, Tint, irfLogoPath)
end



% Calls solo.read_TNR() for multiple days and merges the result.
%
%
% RETURN VALUE
% ============
% TNR
%       Struct, if there are datasets.
%       [], if there are no datasets.
%
% IMPLEMENTATION NOTE: Implementation is made to reproduce old behaviour
% (copy-paste). Can likely be cleaned up/simplified together with
% solo.read_TNR().
%
% function Tnr = read_TNR(Tint)
%
% % UGLY! Obtains list of files to determine for which time intervals to call
% % solo.read_TNR()!!!
%
% TnrFileArray = solo.db_list_files('solo_L2_rpw-tnr-surv-cdag', Tint);
%
% tp  = [];
% pp  = [];
% Tnr = [];
% warning('off', 'fuzzy:general:warnDeprecation_Combine');
%
% % NOTE: Below loop takes most of the time. In each iteration,
% % solo.read_TNR() dominates the time consumption.
% for iFile = 1:length(TnrFileArray)
%   FileTint  = [TnrFileArray(iFile).start, TnrFileArray(iFile).stop];
%   [TnrFile] = solo.read_TNR(FileTint);    % Somewhat time-consuming.
%   if struct(TnrFile)
%     % NOTE: MATLAB documentation (R2019b):
%     % "combine will be removed in a future release"
%     Tnr.t = combine(tp, TnrFile.t);
%     tp    = Tnr.t;
%     Tnr.p = combine(pp, TnrFile.p);
%     pp    = Tnr.p;
%
%     % IMPLEMENTATION NOTE: Only read from TNRp from within this if
%     % clause, since it might not be a struct if read from elsewhere,
%     % even if it in principle means overwriting the value multiple
%     % times as for TNRp.f and TNRp.p_label.
%     Tnr.f       = TnrFile.f;
%     Tnr.p_label = TnrFile.p_label;
%   end
% end
%
% end