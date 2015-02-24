function fits = mms_sdc_sdp_comp_spinfits
% Function to calculate spinfits from data found in datamanager (DATAC).
%
% output: fits     - struct with the following fields
%           .epoch - Timestamps (in TT2000)
%              and for each sdp pair (e12 or e34)
%           .sfit.e12  - Spinfits (matrix with one fit per row) for e12
%           .sfit.e34  - Same for pair e34
%           .sdev.e12  - Standard deviation of each fit for e12
%           .sdev.e34  - Same for pair e34
%           .iter.e12  - Numer of iterations required for each fit for e12
%           .iter.e34  - Same for pair e34
%           .nBad.e12  - Number of bad fits for e12.
%           .nBad.e34  - Same for pair e34
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

sdpPair = {'e12', 'e34'};
time=[];
sfit = struct(sdpPair{1}, [],sdpPair{2}, []);
sdev = struct(sdpPair{1}, [],sdpPair{2}, []);
iter = struct(sdpPair{1}, [],sdpPair{2}, []);
nBad = struct(sdpPair{1}, [],sdpPair{2}, []);

procId = mms_sdc_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.l2pre}
    % FIXME: REPLACE HK 101 with DEFATT for L2pre/L2a etc.
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
    
    % Calculate first timestamp of spinfits to be after start of dce time 
    % and evenly divisable with fitEvery. 
    % I.e. if fitEvery = 5 s, then spinfit timestamps would be like 
    % [00.00.00; 00.00.05; 00.00.10; 00.00.15;] etc.
    % For this one must rely on breakdowntt2000 as the TT2000 (int64)
    % includes things like leap seconds.
    t1 = spdfbreakdowntt2000(dce.time(1)); % Start time in format [YYYY MM DD HH MM ss mm uu nn]
    % Evenly divisable timestamp with fitEvery after t1, in ns.
    t2 = ceil((t1(6)*10^9+t1(7)*10^6+t1(8)*10^3+t1(9))/fitEvery)*fitEvery;
    % Note; spdfcomputett2000 can handle any column greater than expected,
    % ie "62 seconds" are re-calculated to "1 minute and 2 sec".
    t3.sec = floor(t2/10^9);
    t3.ms  = floor((t2-t3.sec*10^9)/10^6);
    t3.us  = floor((t2-t3.sec*10^9-t3.ms*10^6)/10^3);
    t3.ns  = floor(t2-t3.sec*10^9-t3.ms*10^6-t3.us*10^3);
    % Compute what TT2000 time that corresponds to, using computeTT2000.
    t0 = spdfcomputett2000([t1(1) t1(2) t1(3) t1(4) t1(5) t3.sec t3.ms t3.us t3.ns]);

    if( (dce.time(1)<=t0) && (t0<=dce.time(end)))
      for iPair=1:numel(sdpPair)
        % Call mms_spinfit_mx mex file with loaded data
        [time, sfit.(sdpPair{iPair}), sdev.(sdpPair{iPair}), ...
          iter.(sdpPair{iPair}), nBad.(sdpPair{iPair})] = mms_spinfit_mx(maxIt, ...
          minPts, nTerms, double(dce.time)', double(dce.(sdpPair{iPair}).data)', ...
          phase.data' + rad2deg(MMS_CONST.Phaseshift.(sdpPair{iPair})) ,...
          fitEvery, fitInterv, double(t0));

        % Replace non valid values -159e7 (hardcoded in sfit.h used by C/mex)
        % with NaN in Matlab.
        sfit.(sdpPair{iPair})(sfit.(sdpPair{iPair})==-159e7)=NaN;
        sdev.(sdpPair{iPair})(sdev.(sdpPair{iPair})==-159e7)=NaN;
        iter.(sdpPair{iPair})(iter.(sdpPair{iPair})==-159e7)=NaN;
        nBad.(sdpPair{iPair})(nBad.(sdpPair{iPair})==-159e7)=NaN;

        % Change row/column written in mms_spinfit_mx and convert to single
        sfit.(sdpPair{iPair})=single(sfit.(sdpPair{iPair}))';
        sdev.(sdpPair{iPair})=single(sdev.(sdpPair{iPair}))';
        iter.(sdpPair{iPair})=single(iter.(sdpPair{iPair}))';
        nBad.(sdpPair{iPair})=single(nBad.(sdpPair{iPair}))';
      end
    else
      warnStr = sprintf('Too short time series, no data cover first spinfit timestamp (t0=%i)',t0);
      irf.log('warning', warnStr);
    end
    % Store output.
    fits = struct('time', int64(time'), 'sfit', sfit,...
      'sdev', sdev, 'iter', iter, 'nBad', nBad);
    
  otherwise
    irf.log('warning','Only L2Pre supported for spinfits calculation.'); return
end

end
