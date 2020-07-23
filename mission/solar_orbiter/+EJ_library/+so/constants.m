%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-06-26
%
classdef constants   % < handle
    
    properties(Constant)
        % LFR Sampling Frequencies (LSF): F0, F1, F2, F3
        % The variables names (F[0-3]) follow LFR's naming scheme.
        LSF_HZ         = [24576, 4096, 256, 16];
        LSF_NAME_ARRAY = {'F0', 'F1', 'F2', 'F3'};
        
        % Should at least refer to the "normal" LFR snapshot length that BICAS uses.
        % Notes imply that there may be other ones (calibration? LFR-HF? LFR-SCM?).
        LFR_SWF_SNAPSHOT_LENGTH = 2048;
        
        % Max absolute value of set current.
        % NOTE: Does not take into consideration that the actual min & max might be
        % slightly different due to that TM = -2^15 ... (2^15-1)., i.e.
        % max=(2^15-1)/2^15 * 60e-9 sampere.
        MAX_ABS_SAMPERE = 60e-6;
        
        TM_PER_SAMPERE = 32768 / 60e-6;
    end

end
