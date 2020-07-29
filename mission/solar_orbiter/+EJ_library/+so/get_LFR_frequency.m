%
% Convert LFR frequency index value to LSF in Hz, for entire array.
%
% Useful for converting e.g. LFR zVariable FREQ to LSF.
%
% LSF = LFR Sampling Frequency
%
%
% ARGUMENTS
% =========
% iLsf   : Numeric array. Arbitrary size. The LSF index, i.e. 1=LFR freq. F0, and so on.
%          NOTE: Allowes values 1,2,3,4 (not 0).
% freqHz : Frequency in Hz.
%
%
% RETURN VALUES
% =============
% freqHz : Same size as iLsf.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-06-26, when broken out from other code.
%
function freqHz = get_LFR_frequency(iLsf)
    
    LSF_HZ = EJ_library.so.constants.LSF_HZ;
    
    % PROPOSAL: Somehow avoid having an argument for lsfArrayHz. Have it as a constant somehow.
    
    % ASSERTION
    uniqueValues = unique(iLsf);
    if ~all(ismember(uniqueValues, [1,2,3,4]))
        uniqueValuesStr = sprintf('%d', uniqueValues);   % NOTE: Has to print without \n to keep all values on a single-line string.
        error('BICAS:proc_utils:Assertion:IllegalArgument:DatasetFormat', ...
            'Found unexpected values in LSF index (corresponding to LFR FREQ+1). Unique values: %s.', uniqueValuesStr)
    end
    
    % NOTE: Implementation that works for arrays of any size.
    freqHz = ones(size(iLsf)) * NaN;        % Allocate array and set default values.
    freqHz(iLsf==1) = LSF_HZ(1);
    freqHz(iLsf==2) = LSF_HZ(2);
    freqHz(iLsf==3) = LSF_HZ(3);
    freqHz(iLsf==4) = LSF_HZ(4);
end
