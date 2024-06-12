function tests_before_release
% This is a combined function for running multiple tests (matlab.unittest)
% before creating a new stable release of irfu-matlab.

% It should be possible to construct some form of continous-integration
% system running this file to verify new commits does not break existing
% code. And thereby running this combine tests on multiple release versions
% of Matlab. As of 2016/10/ this would include at least R2013a and upwards.

% https://blogs.mathworks.com/developer/2015/01/20/the-other-kind-of-continuous-integration/
% https://blogs.mathworks.com/developer/2020/12/15/cloud-ci-services/
%
% 2021/04: Possibly move to Matlab with GitHub Actions directly instead of going via Travis-ci.
% (Travis-CI is free for opensource but initially limited to a total of 10'000 minutes CPU time,
% Github actions is free for opensource but use a 2'000 minutes/month limit and is preferred by
% many communities.)
% Details: https://github.com/matlab-actions/overview (note: not yet available on github marketplace)

% Setup paths etc.
irf;

% Tests to run (i.e. file name of file containing a "matlab.unittest.TestCase")
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
  'mms_phaseFromSunpulse_2_Test', ...
  'whamp.test_whamp_module'};             % +whamp specific
suite = testsuite(testsToRun);

% Add tests for MATLAB packages in which automated test files can be
% automatically found via MATLAB's filenaming convention.
%
% NOTE: Excludes BICAS which requires MATLAB R2019b.
% NOTE: Searches for tests recursively.
for pkgPathCa = {'irf.fs', 'irf.str', 'irf.utils', 'solo'}
  suite = [ ...
    suite, matlab.unittest.TestSuite.fromPackage(...
    pkgPathCa{1}, 'IncludingSubpackages', true) ...
    ];
end

% Generate test report (PDF), if not pre-existing.
import matlab.unittest.plugins.TestReportPlugin;
runner = matlab.unittest.TestRunner.withTextOutput;
ciPath = fullfile(irf('path'), 'ciPath');
if ~exist(ciPath, 'dir')
  mkdir(ciPath);
end
runner.addPlugin(TestReportPlugin.producingPDF(...
  fullfile(ciPath, 'report.pdf'), 'Verbosity', 3 ...
));

% RUN TESTS
% ---------
% CHECK output. If there are any problems do not release new version of irfu-matlab!
% assertSuccess should allow CI/CD Matlab Actions to report error if any
% test failed.
% NOTE: assertSuccess() only exists in MATLAB R2020a and later.
assertSuccess(runner.run(suite));

end
