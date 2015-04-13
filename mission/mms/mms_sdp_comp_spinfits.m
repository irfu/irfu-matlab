function fits = mms_sdp_comp_spinfits
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
%minPts = 4;     % Minumum of points for one fit
minFraction = 0.20; % Minumum fraction of points required for one fit (minPts = minFraction * fitInterv [s] * samplerate [smpl/s] )
fitEvery = 5*10^9;   % Fit every X nanoseconds.
fitInterv = 20*10^9; % Fit over X nanoseconds interval.

sdpPair = {'e12', 'e34'};
time=[];
sfit = struct(sdpPair{1}, [],sdpPair{2}, []);
sdev = struct(sdpPair{1}, [],sdpPair{2}, []);
iter = struct(sdpPair{1}, [],sdpPair{2}, []);
nBad = struct(sdpPair{1}, [],sdpPair{2}, []);

dce = mms_sdp_datamanager('dce');
if mms_is_error(dce)
  errStr='Bad DCE input, cannot proceed.';
  irf.log('critical',errStr); error(errStr);
end
phase = mms_sdp_datamanager('phase');
if mms_is_error(phase)
  errStr='Bad PHASE input, cannot proceed.';
  irf.log('critical',errStr); error(errStr);
end
samplerate = mms_sdp_datamanager('samplerate');
if mms_is_error(samplerate)
  errStr='Bad SAMPLERATE input, cannot proceed.';
  irf.log('critical',errStr); error(errStr);
end
% Calculate minumum number of points req. for one fit covering fitInterv
minPts = minFraction * samplerate * fitInterv/10^9; % "/10^9" as fitInterv is in [ns].

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
    sigE = sdpPair{iPair};
    probePhaseRad = unwrap(phase.data*pi/180) - MMS_CONST.Phaseshift.(sigE);
    dataIn = dce.(sigE).data;
    bits = bitor(MMS_CONST.Bitmask.SIGNAL_OFF,MMS_CONST.Bitmask.SWEEP_DATA);
    dataIn = mask_bits(dataIn, dce.(sigE).bitmask, bits);
    idxBad = isnan(dataIn); dataIn(idxBad) = [];
    timeIn = dce.time; timeIn(idxBad) = [];
    probePhaseRad(idxBad) = [];
    
    % Call mms_spinfit_m, .m interface file for the mex compiled file
    % XXX FIXME: converting time here to double reduces the precision.
    % It would be best if the function accepted time as seconds from 
    % the start of the day or t0
    [time, sfit.(sigE), sdev.(sdpPair{iPair}), iter.(sigE), nBad.(sigE)] = ...
      mms_spinfit_m(maxIt, minPts, nTerms, double(timeIn), double(dataIn), ...
      probePhaseRad, fitEvery, fitInterv, t0);
    
    % Change to single
    sfit.(sigE) = single(sfit.(sigE));
    sdev.(sigE) = single(sdev.(sigE));
    iter.(sigE) = single(iter.(sigE));
    nBad.(sigE) = single(nBad.(sigE));
  end
else
  warnStr = sprintf(['Too short time series:'...
    ' no data cover first spinfit timestamp (t0=%i)'],t0);
  irf.log('warning', warnStr);
end
% Store output.
fits = struct('time', int64(time), 'sfit', sfit,...
  'sdev', sdev, 'iter', iter, 'nBad', nBad);

end
