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
explicitTestFiles = {...
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
  'whamp.test_whamp_module', ...    % +whamp specific
  'solo.db_list_files___UTEST'};
TestArray = testsuite(explicitTestFiles);

% Add tests for MATLAB packages for which automated test files can be
% found AUTOMATICALLY and RECURSIVELY via MATLAB's filenaming convention.
%
% NOTE: Excluding solo.qli since (1) solo.qli.batch tests fail due to not being
% able to launch/access processing pool for unknown reasons (tested code only
% uses processing pool pool for logging some if its parameters), and (2) some
% of those particular tests are slow. See GitHub CI/CD log for commit where CI
% failed: ecbb4b13a Erik P G Johansson (2024-09-16 16:37:07 +0200) Merge branch
% 'devel' into SOdevel.
% --
% 2024-09-23: Testing to include solo.qli after possibly relevant modification:
% commit 669019b4c Erik P G Johansson (2024-09-23 13:22:23 +0200).
% solo.qli.batch: delete(gcp('nocreate')) before creating parpool()
for pkgPathCa = {...
    'irf', ...
    'solo.BSACT_utils', 'solo.adm', 'solo.hwzv', 'solo.shk', 'solo.qli'}
  TestArray = [ ...
    TestArray, matlab.unittest.TestSuite.fromPackage(...
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
assertSuccess(runner.run(TestArray));

end
