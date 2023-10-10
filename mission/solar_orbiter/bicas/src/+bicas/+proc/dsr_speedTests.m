%
% Simple speed tests.
%
% OBSERVATIONS
% ============
% Eliminating inner function has effect, but only a small one (~10%).
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function dsr_speedTests()
    
    N_OSR_SAMPLES = 351103;
    % NOTE: N_BINS is not the same as DSR samples, if bin DSR output can be empty.
    N_BINS        = 3610;

    ZvOsrFpa        = bicas.utils.FillPositionsArray(ones(N_OSR_SAMPLES, 1), 'NO_FILL_POSITIONS');
    iRecordsInBinCa = cell(N_BINS, 1);
    
    for iBin = 1:N_BINS
        % Not entirely sure of calculation, but should be good enough.
        i1 = floor((iBin-1)/N_BINS * N_OSR_SAMPLES) + 1;
        i2 = ceil(  iBin   /N_BINS * N_OSR_SAMPLES);
        iRecordsInBinCa{iBin} = i1:i2;
    end
    
    test_speed(ZvOsrFpa, iRecordsInBinCa, @bicas.proc.dsr.downsample_ZV_minimum);
    test_speed(ZvOsrFpa, iRecordsInBinCa, @bicas.proc.dsr.downsample_ZV_minimum_W_INNER_ARRAYS)

    test_speed(ZvOsrFpa, iRecordsInBinCa, @bicas.proc.dsr.downsample_ZV_bitmask);
    test_speed(ZvOsrFpa, iRecordsInBinCa, @bicas.proc.dsr.downsample_ZV_bitmask_W_INNER_ARRAYS);
end



function test_speed(ZvOsrFpa, iRecordsInBinCa, fh)
    t = tic();
    [~] = fh(ZvOsrFpa, iRecordsInBinCa);
    tSec = toc(t);
    
    nOsrSamples = size(ZvOsrFpa, 1);
    nBins       = size(iRecordsInBinCa, 1);    
    funcName    = func2str(fh);
    
    fprintf('%-51s, %i [OSR samples], %i [bins]: %f [s], %d [s/OSR sample], %d [s/bin]\n', funcName, nOsrSamples, nBins, tSec, tSec/nOsrSamples, tSec/nBins);
end