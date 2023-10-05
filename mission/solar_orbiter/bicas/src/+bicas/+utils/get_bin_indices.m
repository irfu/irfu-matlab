%
% Given (1D) bin boundaries, and data (e.g. timestamps), find the indices to
% those data for every bin separately.
%
%
% RATIONALE FOR RETURN VALUES
% ===========================
% Can be used for finding e.g. the indices to CDF records for every bin
% separately, given specified bin boundaries in Epoch. Returned indices can then
% be used to collect zVar data for the respective bins.
%
%
% SPEED CONSIDERATIONS, IMPLEMENTATION NOTES
% ==========================================
% The original version of this code constituted the bulk of processing time when
% constructing downsampled datasets (before it was a separate function; roughly
% equivalent to internal implementation_RAW()). Optimizing this function for
% speed is therefore important.
% --
% NOTE: Since optimizing this function for speed is important, it is also
% important to be able conveniently investigate/experiment/debug it and design
% it, which means it is more important to keep it as a separate function and
% document it well.
%
%
% ARGUMENTS
% =========
% t
%       1D numeric array. Sorted increasing (non-strictly).
% bb
%       1D numeric array. Sorted increasing (non-strictly).
%       Bin boundaries (BB). There are no bins below or above the highest
%       boundaries (to infinities).
%       NOTE: Smaller boundary of a bin is inclusive, higher boundary is
%             exclusive.
%       NOTE: Can have zero-size bins, but no data will be assigned to those
%             bins (since there is no "t" which satisfies b =< t < b).
%       NOTE: Empty xBoundaries ==> iInBinCa == {}
% nBbThreshold
%       Scalar number, >=3.
%       Min number of bin boundaries for which to recursively split up the task
%       to speed up execution. Indirectly determines the number of recursive
%       calls.
%
%
% RETURN VALUES
% =============
% iInBinCa
%       1D cell array of 1D numerical arrays. {iBin}(iX) = Index into x, for x
%       index that has been assigned to bin iBin.
%       NOTE: Data outside the lowest and highest boundary are not represented.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-20.
%
function iInBinCa = get_bin_indices(t, bb, nBbThreshold)
    %
    % PROPOSAL: Implement using discretize().
    %   CON: Recursive implementation does not work with the special behaviour
    %        that the last bin includes the last edge/boundary (special
    %        behaviour)?
    %   CON: Has tried using discretize() and then either use find(), or not
    %        use find() but use index ZVs (quality variables only?) using logical
    %        indexing instead. Seems to take as much time again.
    %
    % PROPOSAL: Use parfor on "raw" implementation.
    %   CON: Does not seem to work for larger test dataset (brain)
    %        (L2-->L2+L3).
    %        solo_L2_rpw-lfr-surv-cwf-e_20200704_V01.cdf (99 MiB).
    %
    % PROPOSAL: Move to bicas.proc.dsr.
    %   CON: Potentially generic outside of BICAS.
    % --------------------------------------------------------------
    % Speed test for older implementation 2021-05-19, brain:
    % output1 : for
    % output2 : parfor
    % ==> output1/solo_L2_rpw-lfr-surv-cwf-e-1-second_20200704_V01.cdf.2021-05-19T19.55.50.log <==
    % 2021-05-19T19:58:18 -- DEBUG -- SPEED -- main_without_error_handling: 146.251 [s] (wall time)
    % 
    % ==> output1/solo_L3_rpw-bia-efield_20200704_V01.cdf.2021-05-19T19.58.18.log <==
    % 2021-05-19T19:58:38 -- DEBUG -- SPEED -- main_without_error_handling: 19.87 [s] (wall time)
    % 
    % ==> output2/solo_L2_rpw-lfr-surv-cwf-e-1-second_20200704_V01.cdf.2021-05-19T20.02.53.log <==
    % 2021-05-19T20:05:27 -- DEBUG -- SPEED -- main_without_error_handling: 152.514 [s] (wall time)
    % 
    % ==> output2/solo_L3_rpw-bia-efield_20200704_V01.cdf.2021-05-19T20.05.27.log <==
    % 2021-05-19T20:05:45 -- DEBUG -- SPEED -- main_without_error_handling: 18.5625 [s] (wall time)
    % --------------------------------------------------------------
    

    % BB = Bin Boundaries
    
    % ASSERTIONS
    assert(iscolumn(t),  'Argument t is not a column vector.')
    assert(iscolumn(bb), 'Argument binBoundaries is not a column vector.')
    
    % NOTE: REQUIRED by recursive implementation for adding back index offsets
    % when merging results from recursive calls.
    assert(issorted(t,  'ascend'), 'Argument t is not sorted and increasing.')
    
    % Does not need to check for STRICTLY ascending bin boundaries.
    % Algorithm works for non-strictly increasing values.
    assert(issorted(bb, 'ascend'), 'Argument bb is not sorted and increasing.')
    
    % Going below threshold leads to infinite recursion.
    assert(isscalar(nBbThreshold) && (nBbThreshold >= 3))
    
    % Slow for large vectors.
    % iInBinCa = implementation_RAW(t, binBoundaries);
    % Faster
    iInBinCa = implementation_RECURSIVE(t, bb, nBbThreshold);
end



% Faster recursive implementation.
%
% The basic implementation (implementation_RAW()) suffers (theoretically) from
% scaling time ~ nData*nBins. This implementation divides the task into two
% simpler tasks by dividing the bins and data into two groups, thereby
% disproportionally reducing the processing. Exact number of recursions is
% arbitrary.
%
function iInBinCa = implementation_RECURSIVE(t, bb, nBbThreshold)
    
    nBb = numel(bb);
    
    if nBb >= nBbThreshold
        i = round((1+nBb) / 2);
        
        % Group bins into two large bins.
        iInBinCa  = implementation_RAW(t, bb([1,i,end]));
        
        % NOTE: Separate variables which are useful for debugging (inspecting
        % values).
        t1  = t(iInBinCa{1});
        t2  = t(iInBinCa{2});
        bb1 = bb([1:i  ]);
        bb2 = bb([i:end]);
        
        % NOTE: RECURSIVE CALL
        iInBinCa1 = implementation_RECURSIVE(t1, bb1, nBbThreshold);
        iInBinCa2 = implementation_RECURSIVE(t2, bb2, nBbThreshold);
        
        % Add offset to the resulting indices from the SECOND recursive call
        % ------------------------------------------------------------------
        % ASSUMES: t is (non-strictly) increasing. Can otherwise not adjust
        % indices from second recursive call this easily since iInBinCa{1} would
        % not constitute a continuous range of indices starting at 1. Relaxing
        % the assumption would require more index magic.
        iInBinCa2 = cellfun(@(x) (x+numel(t1)), iInBinCa2, 'UniformOutput', false);
        
        iInBinCa = [iInBinCa1; iInBinCa2];
    else
        iInBinCa = implementation_RAW(t, bb);
    end

end



function iInBinCa = implementation_RAW(t, bb)
    % PROPOSAL: Special case for bb small: if-then
    
    nBins = numel(bb) - 1;

    % NOTE: Empty xBoundaries
    % ==> nBins    = -1
    % ==> iInBinCa = cell(0,1) (sic!),
    %     since cell(-1,1) == cell(0,1).
    iInBinCa = cell(nBins, 1);
    for iBin = 1:nBins
        bInBin = (bb(iBin) <= t) & (t < bb(iBin+1));

        % NOTE: find() does not always return a column vector for column vector
        % input.
        % find(zeros(0,1)) == <0x1>
        % find(zeros(1,1)) == <0x0>  (sic!)
        % find(zeros(2,1)) == <0x1>
        
        iInBin = find(bInBin);
        %iInBin = find([b; 0; 0]);   % Not a speed improvement.
        
        % Normalize result to always be a column vector.
        % One can also use
        %   b = [b; 0; 0];
        % before calling find() but that is much slower for non-recursive
        % implementation, and slightly slower for recursive implementation.
        if isempty(iInBin)
            iInBin = zeros(0,1);
        end
        
        iInBinCa{iBin} = iInBin;
    end
end
