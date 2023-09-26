%
% Performance test for bicas.utils.SameRowsMap.setRows().
% The performance of this function is critical to the performance of LFR-SWF
% processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function SameRowsMap___setRows_SpeedTest
    clear classes
    
    N_DATA_POINTS    = 20;
    N_ARRAY_COLS     = 1024;
    N_ARRAY_ROWS_MIN = 1e1;
    N_ARRAY_ROWS_MAX = 3e3;
    FH = @speed_test_setRows;
    %FH = @speed_test_growing_array;
    %FH = @speed_test_preallocated_array;
    
    
    %close all
    nRowsArray = round( logspace(log10(N_ARRAY_ROWS_MIN), log10(N_ARRAY_ROWS_MAX), N_DATA_POINTS) );
    
    tSecArray = [];
    for nRows = nRowsArray
        fprintf('nRows = %i\n', nRows);
        tSec = FH(nRows, N_ARRAY_COLS);
        fprintf('tSec  = %f\n', tSec);
    
        tSecArray(end+1) = tSec;
    end
    
    yArray = tSecArray ./ nRowsArray;
    figure
    loglog(nRowsArray, yArray, '*-')
    title(irf.graph.escape_str(func2str(FH)))
    grid on
    xlabel('n = Array size')
    %ylabel('t [s] = Time per array')
    ylabel('t [s/row] = Time per array row')
end



function tSec = speed_test_setRows(nArrayRows, nArrayCols)
    
    bigArray   = repmat(linspace(1, nArrayRows, nArrayRows)', 1, nArrayCols);
    smallArray = NaN(1, nArrayCols);
    
    M1 = bicas.utils.SameRowsMap('char', nArrayRows, 'constant', bigArray,   {'K'});
    M2 = bicas.utils.SameRowsMap('char', 1,          'constant', smallArray, {'K'});
    
    t = tic();
    for i = 1:nArrayRows
        M1.setRows(M2, i)
    end
    tSec = toc(t);
    
    %assert(all(isnan(M1.get('K')), 'all'))
end



% 1e4, 1e4
function tSec = speed_test_growing_array(nArrayRows, nArrayCols)
    growingArray = zeros(0, nArrayCols);
    
    t = tic();
    for i = 1:nArrayRows
        growingArray(end+1, :) = NaN(1, nArrayCols);
    end
    tSec = toc(t);
    
    assert(isequal(size(growingArray), [nArrayRows, nArrayCols]))
end



% 1e6, 1e4
function tSec = speed_test_preallocated_array(nArrayRows, nArrayCols)
    preallocArray = zeros(nArrayRows, nArrayCols);
    
    t = tic();
    for i = 1:nArrayRows
        preallocArray(i, :) = NaN(1, nArrayCols);
    end
    tSec = toc(t);
    
    assert(isequal(size(preallocArray), [nArrayRows, nArrayCols]))
end
