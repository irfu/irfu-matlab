%
% Convert LFR sampling frequency index value to LSF in Hz, for entire array.
% Useful for converting e.g. LFR zVariable FREQ to LSF.
%
% LSF = LFR Sampling Frequency (an explicit frequency value; not e.g. F0)
%
%
% ARGUMENTS
% =========
% iLsf
%       Numeric array. Arbitrary size. The LSF index, i.e. 1=LFR freq. F0, and
%       so on.
%       NOTE: Allows values 1,2,3,4 (not 0).
% freqHz
%       Sampling frequency in Hz.
%
%
% RETURN VALUES
% =============
% freqHz
%       Same size as iLsf.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-26, when broken out from other code.
%
function freqHz = get_LFR_frequency(iLsf)
    % PROPOSAL: Change name:
    %   get_LSF()
    %       CON: Too BICAS-specific abbreviation.
    %           CON: Already used in argument name, solo.hwzv.const.LSF_HZ,
    %                function documentation.
    %   get_LFR_sampling_frequency()

    LSF_HZ = solo.hwzv.const.LSF_HZ;

    % ASSERTION
    uniqueValues = unique(iLsf);
    if ~all(ismember(uniqueValues, [1:4]))
        % NOTE: Has to print without \n to keep all values on a single-line
        % string.
        %uniqueValuesStr = sprintf('%d ', uniqueValues);
        uniqueValuesStr = strjoin(...
            irf.str.sprintf_many('%d', uniqueValues), ', ');

        error(...
            ['Found unexpected values in LSF index (corresponding to', ...
            ' LFR FREQ+1). Unique values: %s.'], ...
            uniqueValuesStr)
    end

    % NOTE: Implementation that works for arrays of any size.
    freqHz = nan(size(iLsf));        % Allocate array and set default values.
    freqHz(iLsf==1) = LSF_HZ(1);
    freqHz(iLsf==2) = LSF_HZ(2);
    freqHz(iLsf==3) = LSF_HZ(3);
    freqHz(iLsf==4) = LSF_HZ(4);
end
