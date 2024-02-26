function tests_before_release_R2021a
% Function that runs multiple tests (matlab.unittest) under MATLAB R2021a
% specifically.
%
% RATIONALE
% ---------
% The primary use case for this file is BICAS (SolO/RPW/BIAS calibration code)
% which must be able to run on MATLAB R2019b (sic!) per agreement with LESIA/ROC
% which runs it. However, MATLAB's "Action for Setting Up MATLAB on
% GitHub-Hosted Runner" (https://github.com/matlab-actions/setup-matlab/) only
% supports MATLAB R2021a and later. Therefore, it is not possible to run CI test
% under R2019b. Therefore running BICAS tests on GitHub under MATLAB R2021a,
% despite it being suboptimal.

% Setup paths etc.
irf;

% Create empty list which is added to by later code.
suite = testsuite();

% Add tests for MATLAB packages in which automated test files can be
% automatically found via MATLAB's filenaming convention.
for pkgPathCa = {'bicas', 'solo.hwzv'}
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
  fullfile(ciPath, 'report_R2021a.pdf'), 'Verbosity', 3 ...
  ));

% RUN TESTS
assertSuccess(runner.run(suite));

end
