%
% Performance test for bicas.utils.FillPositionsArray.subsasgn(), i.e. for
% assigning Fpa(...) = OtherFpa.
%
% The performance of that function could be critical to the performance of
% LFR-SWF processing if FPAs are ever used for updating the (preallocated)
% global array of samples. This code tests how modifying subsets of internal
% values scales with the size of the pre-allocated internal variable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function FillPositionsArray___subsasgn_SpeedTest
% PROPOSAL: Name that implies preallocation.

    clear classes
    
    N_DATA_POINTS    = 20;
    N_ARRAY_COLS     = 1024;
    N_ARRAY_ROWS_MIN = 1e1;
    N_ARRAY_ROWS_MAX = 5e3;
    FH = @speed_test_subsasgn;
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



function tSec = speed_test_subsasgn(nArrayRows, nArrayCols)
    
    bigArray   = repmat(linspace(1, nArrayRows, nArrayRows)', 1, nArrayCols);
    smallArray = NaN(1, nArrayCols);
    
    Fpa1 = bicas.utils.FillPositionsArray(bigArray,   'NO_FILL_POSITIONS');
    Fpa2 = bicas.utils.FillPositionsArray(smallArray, 'NO_FILL_POSITIONS');
    
    t = tic();
    for i = 1:nArrayRows
        Fpa1(i, :) = Fpa2;
    end
    tSec = toc(t);
    
    assert(all(isnan(Fpa1.get_data(-1)), 'all'))
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
