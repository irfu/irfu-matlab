%
% Performance test for bicas.utils.SameRowsMap.set_rows().
%
% The performance of that function is critical to the performance of LFR-SWF
% processing. This code tests how modifying subsets of internal values scales
% with the size of the pre-allocated internal variable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function SameRowsMap___set_rows_SpeedTest
clear classes

N_DATA_POINTS    = 20;
N_ARRAY_COLS     = 1024;
N_ARRAY_ROWS_MIN = 1e1;
N_ARRAY_ROWS_MAX = 1e4;
FH = @speed_test_set_rows;
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



function tSec = speed_test_set_rows(nArrayRows, nArrayCols)

bigArray   = repmat(linspace(1, nArrayRows, nArrayRows)', 1, nArrayCols);
smallArray = NaN(1, nArrayCols);

Srm1 = bicas.utils.SameRowsMap('string', nArrayRows, 'CONSTANT', bigArray,   {"K"});
Srm2 = bicas.utils.SameRowsMap('string', 1,          'CONSTANT', smallArray, {"K"});

t = tic();
for i = 1:nArrayRows
  Srm1.set_rows(Srm2, i)
end
tSec = toc(t);

assert(all(isnan(Srm1("K")), 'all'))
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
