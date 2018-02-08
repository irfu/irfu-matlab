%
% Define constants used by and defined by TM_power_budget.engine.
%
% This function is meant to collect definitions of frequently used hardcoded values and do nothing else.
%
%
% Initially created 2018-01-31 by Erik Johansson, IRF Uppsala.
%
function Constants = default_constants()
% PROPOSAL: Merge default_instrument_modes into this one?

Constants = [];

[Constants.InsModeDescrList, Constants.RadModeDescrList] = TM_power_budget.default_instrument_modes();

Constants.SystemPrps = [];
Constants.SystemPrps.storageBytes = 0.25 * 1/8 * 2^40;   % 0.25 Tbit. Probably not certain. TBit or TByte?
Constants.SystemPrps.clusterSizeBytes    = 2 * 2^20;            % 2 MiB.
Constants.SystemPrps.survFileSizeBytes   = 1.0 * Constants.SystemPrps.clusterSizeBytes;   % Arbitrarily chosen
Constants.SystemPrps.richFileSizeBytes   = 0.5 * Constants.SystemPrps.clusterSizeBytes;   % Arbitrarily chosen

end
