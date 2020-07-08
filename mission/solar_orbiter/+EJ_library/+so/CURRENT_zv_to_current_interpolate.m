%
% Utility function for easily deriving (set) current at arbitrary timestamps almost directly from
% SOLO_L1_RPW-BIA-CURRENT dataset zVariables. Therefore has some extra features to handle the quirks of the dataset
% format.
%
% Could be useful for plotting.
% NOTE: "CURRENT" refers to DATASET_ID=SOLO_L1_RPW-BIA-CURRENT
% NOTE: Is used by BICAS, which want to detect anomaly (multiple instances of equivalent successive bias settings) and
% react to it itself.
%
%
% ARGUMENTS
% ==========
% t1        : 1D vector. Time. Sorted. Type unimportant (linear increasing with time). Can be TT2000.
%             NOTE: Must be floating-point since interp1 requires it.
% zvIBIASx1 : 1D vector. Same length as t1. New bias values.
%             NOTE: Must be floating-point since interp1 requires it.
% t2        : 1D vector. Time. Same type of time as t1.
%
% 
% RETURN VALUE
% ============
% zvIBIASx2         : Current values at t2. Double. NaN for timestamps before first input timestamp.
% duplicatesAnomaly : Whether detected known anomaly in SOLO_L1_RPW-BIA-CURRENT datasets.
%                     EJ_library.so.CURRENT_zv_to_current.
%                     RATIONALE: Caller (e.g. BICAS) can give error, warning.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
function [zvIBIASx2, duplicatesAnomaly] = CURRENT_zv_to_current_interpolate(t1, zvIBIASx1, t2)
    % PROPOSAL: Accept all three antennas at the same time.
    %   PRO: Can verify format better.
    %       CON: Only that every record contains data for exactly one antenna.
    % PROPOSAL: Configurable to accept exact CDF format.
    %   PRO: More useful outside BICAS.
    %   NOTE: Function would still need to know fill value.
    %
    % PROPOSAL: Function names?
    %   Reference to datasets
    
    % NOTE: interp1 requires double.
    t2 = double(t2);
    
    [t1b, zvIBIASx1b, duplicatesAnomaly] = EJ_library.so.CURRENT_zv_to_current(t1, zvIBIASx1);
    
    % IMPLEMENTATION NOTE: Bias currents are set VERY RARELY. Must therefore use interpolation method 'previous'.
    % NOTE: interp1 does NOT require t1 to be sorted.
    % NOTE: interp1 does not permit identical x1 values, not even if the coresponding y1 values are identical.
    zvIBIASx2 = interp1(t1b, zvIBIASx1b, t2, 'previous');
    
    % NOTE: interp1(... 'previous') returns NaN also for t2 > max(t1)!! Must therefore do this oneself.    
    zvIBIASx2(t2 > t1b(end)) = zvIBIASx1b(end);
end
