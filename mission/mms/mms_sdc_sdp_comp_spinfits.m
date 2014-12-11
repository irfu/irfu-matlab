function fits = mms_sdc_sdp_comp_spinfits
% Function to calculate spinfits from data found in datamanager (DATAC).
%
% output: fits     - struct with the following fields
%           .epoch - Timestamps (in TT2000)
%           .sfit  - Spinfits (matrix with one fit per column or row) <--           CHECK which is which.
%           .sdev  - Standard deviation of each fit
%           .iter  - Numer of iterations required for each fit
%           .nBad  - Number of bad fits.
%

% No arguments in, as of yet.
narginchk(0,0);

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
fits = MMS_CONST.Error;

% Some default settings
maxIt = 3;      % Maximum of iterations to run fit
nTerms = 3;     % Number of terms to fit, Y = A + B*sin(wt) + C*cos(wt) +..., must be odd.
minPts = 4;     % Minumum of points for one fit
fitEvery = 5*10^9;   % Fit every X nanoseconds.
fitInterv = 20*10^9; % Fit over X nanoseconds interval.

procId = mms_sdc_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.l2pre}
    hk_101 = mms_sdc_sdp_datamanager('hk_101');
    if isnumeric(hk_101) && numel(hk_101)==1 && hk_101==MMS_CONST.Error,
      irf.log('warning','Bad hk_101 input'); return
    end
    dce = mms_sdc_sdp_datamanager('dce');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad dce input'); return
    end
    phase = mms_sdc_sdp_datamanager('phase');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad phase input'); return
    end
    
    % Call mms_spinfit_mx mex file with loaded data
    [time, sfit, sdev, iter, nBad] = mms_spinfit_mx(maxIt, ...
      minPts, nTerms, double(dce.time)', double(dce.e12.data)', phase.data', fitEvery, fitInterv);
    
    % Store output. BUGFIX REMEMBER TO CHANGE ROW / COLUMNS in
    % mms_spinfit_mx
    fits = struct('time', int64(time'), 'sfit', single(sfit'), 'sdev', single(sdev'), 'iter', single(iter'),...
      'nBad', single(nBad'));
    
  otherwise
    irf.log('warning','Only L2Pre supported for spinfits calculation.'); return
end

end
