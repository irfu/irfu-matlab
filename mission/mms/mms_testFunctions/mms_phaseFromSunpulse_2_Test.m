% MMS_PHASEFROMSUNPULSE_2_TEST is a unit testing framework for testing
% mms_phaseFromSunpulse calculation of spin phase from sunpulse data
% (found in HK_101).
%   results = MMS_PHASEFROMSUNPULSE_2_TEST creates a unit testing
%   framework for testing phase calculation, and make sure that they are as
%   expected. Including testing of flags identifying extrapolations, as
%   well as data gaps, etc.
%   These functions require Matlab R2013b or later.
%
%       Example:
%               results = mms_phaseFromSunpulse_2_Test
%               results.run
%
%       See also MATLAB.UNITTEST.

function tests = mms_phaseFromSunpulse_2_Test
    % Verify version of Matlab (R2013b or later).
    if( verLessThan('matlab', '8.3') )
        error('Require R2013b or later to run this test. Please upgrade.');
    end
    % Setup basics, such as logging.
    global ENVIR MMS_CONST;
    if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
    ENVIR = mms_sdc_sdp_init;
    % Run the tests
    tests = functiontests(localfunctions);
end

function testCompleteSunpulses(testCase)
  % This test will try to calculate phase for timestamps which are exactly
  % halfway between the sunpulses.

  % Get inputs (reqTimes start halfway between two sunpulses and continues
  % with 1 second resolution). Sunpulses are computed with perfect spin
  % period of 20 sec.
  [reqTimes, hk_101] = getInputDataForTests;

  % Compute phase and corresponding flags.
  [phase, flags] = mms_sdp_phase_2(hk_101, reqTimes);

  % Expected phase, first reqTimes was exactly halfway between the two
  % first distinct sunpulses => phase = 180 deg from sunpulse sensor
  % (76 deg off from BCS-X).
  expPhase = mod((180-76):18:(180-76+18*(length(reqTimes)-1)),360)';

  % Phase computed after the last pulse (ie not in between two pulses)
  % should be marked by an extra "10" added to its flag. Otherwise it
  % should be "0" indicating "perfect" phase.
  expFlags = zeros(size(expPhase), 'int16');
  tmp = unique(sort(hk_101.sunpulse)); tmp = tmp(end)-diff(tmp(end-1:end))/2;
  expFlags(reqTimes>=tmp) = int16(10);

  % Verify phase, within rounding error of the static values listed above.
  verifyEqual(testCase, phase, expPhase, 'AbsTol', 10^-12);
  % Verify flags
  verifyEqual(testCase, flags, expFlags);
end


function [reqTimes, hk_101] = getInputDataForTests
  % Return "prestine" sunpulse data and time stamps to compute phase for.
  hk_101 = struct('time', [], 'sunpulse', [], 'sunssps', [], 'iifsunper', []);
  % Time interval: '2016-01-01T00:00:00.000000000Z' to '2016-01-02T00:00:00.000000000Z'
  timeInterval = int64([504878468184000000; 504964868184000000]);
  % Every 10 seconds (adjust by 1 sec => hk timestamps should be after
  % corresponding pulse was measured).
  hk_101.time = (timeInterval(1):10^10:timeInterval(end) + int64(10^9))';
  % Compute sunpulses with assumed perfect 20 second spin period
  tmp1 = timeInterval(1):2.0*10^10:timeInterval(end);
  tmp1 = tmp1 - int64(1.0*10^10); % Half a spin before
  hk_101.sunpulse = sort([tmp1; tmp1]);
  hk_101.sunpulse = hk_101.sunpulse(1:length(hk_101.time))';
  hk_101.sunssps = zeros(size(hk_101.time), 'uint8'); % All zeros => perfect sunpulse data
  hk_101.iifsunper = 20000000*ones(size(hk_101.time), 'uint32');

  % Request time evenly from start to end every 1 second, this will cover
  % some times after the last pulse (for which extrapolation will be needed).
  reqTimes = (timeInterval(1):10^9:timeInterval(end))';

end