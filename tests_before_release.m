function tests_before_release
% This is a combined function for running multiple tests (matlab.unittest)
% before creating a new stable release of irfu-matlab.

% It should be possible to construct some form of continous-integration
% system running this file to verify new commits does not break existing
% code. And thereby running this combine tests on multiple release versions
% of Matlab. As of 2016/10/ this would include at least R2013a and upwards.

% http://blogs.mathworks.com/developer/2015/01/20/the-other-kind-of-continuous-integration/

failed = false;
DEBUG = false; % Display some extra information or not.

% Setup paths etc.
irf;

% Tests to run (ie file name of file containing a "matlab.unittest.TestCase")
testsToRun = {...
  'TestTimeArray', ...              % IRF generic
  'test_irf_time', ...
  'irf.test_geocentric_coordinate_transformation', ...
  'testC4', ...                     % Cluster specific
  'test_mms_defatt_phase', ...      % MMS specific
  'test_mms_spinfit', ...
  'test_mms_dsl2gse', ...
  'mms_phaseFromSunpulse_2_Test'};

for ii = 1:length(testsToRun)
  try
    test = eval([testsToRun{ii},';']);
    results = test.run;
    if DEBUG, display(results); end
    if any([results.Failed]), failed = true; end
  catch ME
    display(getReport(ME, 'extended'));
    failed = true;
    continue
  end
end

if(failed)
  disp('================================================================');
  warning('Test failed! Do not release irfu-matlab until all issues solved.');
  disp('================================================================');
  %exit(1);
else
  disp('================================================================');
  disp('     All automatic tests completed without issues               ');
  disp('================================================================');
  %exit(0);
end

end