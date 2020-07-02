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
    end

end