function get_bin_indices___ATEST()
    % NOTE: Implementation might use internal hardcoded constant. Automatic
    % tests must include large enough input data to trigger alternate
    % implementation.
    
    % NOTE: Changes rows to columns to make typing arguments easier.
    function add_test(xRow, xBoundariesRow, iInBinCaRowRow)
        
        iInBinCa = cellfun(@(c) (c'), iInBinCaRowRow, 'UniformOutput', false)';
        
        % NOTE: Thresholds chosen to overlap with lengths of 
        for nBbThreshold = [3,4,5,6,10]
            tl{end+1} = EJ_library.atest.CompareFuncResult(...
                @bicas.utils.get_bin_indices, ...
                {xRow', xBoundariesRow', nBbThreshold}, {iInBinCa});
        end
    end
    
    tl = {};
    A1x0 = zeros(1,0);
    C1x0 = cell( 1,0);
    
    % Edge cases: No data, or no boundaries
    add_test(A1x0,    A1x0,   C1x0)
    add_test([3],     A1x0,   C1x0)
    add_test(A1x0,    [5],    C1x0)
    add_test(A1x0,    [5,10], {A1x0})
    add_test([3,13],  [5,10], {A1x0})
    
    add_test([6,7,8], [5,10], {[1,2,3]})
    
    add_test([3,6,7,8,11],  [5,10],      {[2,3,4]})
    add_test([3,6,7,8,11],  [0,5,10,15], {[1], [2,3,4], [5]})

    % Data on bin boundaries.
    add_test([5,10], [0,5,10,15], {A1x0, [1], [2]})
    
    % Data in zero-length bins (on boundaries).
    add_test([5,10], [0, 5,5, 10,10, 15], {A1x0, A1x0, [1], A1x0, [2]})
    
    % Zero-length bin boundaries.
    add_test([3,6,7,8,11],  [0,5,5,10,10,15], {[1], A1x0, [2,3,4], A1x0, [5]})
    
    % Data not sorted.
    %add_test([11,6,3,8,7],  [0,5,10,15], {[3], [2,4,5], [1]})

    
    % TODO: Data in zero-length bins.
    % TODO: Vary N_BB_MAX
    
    EJ_library.atest.run_tests(tl);
    
end

% tic; iRecordsInBinCa2 = bicas.utils.get_bin_indices(zvAllTt2000, boundariesTt2000); toc ; isequal(iRecordsInBinCa, iRecordsInBinCa2)