function tests_before_release
% This is a combined function for running multiple tests (matlab.unittest)
% before creating a new stable release of irfu-matlab.

% It should be possible to construct some form of continous-integration
% system running this file to verify new commits does not break existing
% code. And thereby running this combine tests on multiple release versions
% of Matlab. As of 2016/10/ this would include at least R2013a and upwards.

% https://blogs.mathworks.com/developer/2015/01/20/the-other-kind-of-continuous-integration/


% Setup paths etc.
irf;

% Tests to run (ie file name of file containing a "matlab.unittest.TestCase")
testsToRun = {...
  'TestTimeArray', ...              % IRF generic
  'test_irf_time', ...
  'TestTSeries', ...
  'irf.test_geocentric_coordinate_transformation', ...
  'test_irf_plot', ...
  'testC4', ...                     % Cluster specific
  'test_mms_defatt_phase', ...      % MMS specific
  'test_mms_spinfit', ...
  'test_mms_dsl2gse', ...
  'mms_phaseFromSunpulse_2_Test'};

import matlab.unittest.plugins.TestReportPlugin;
runner = matlab.unittest.TestRunner.withTextOutput;
runner.addPlugin(TestReportPlugin.producingHTML('Verbosity',3));
runner.run(testsuite(testsToRun));

% CHECK output, if any problems do not release new version of irfu-matlab!

end
