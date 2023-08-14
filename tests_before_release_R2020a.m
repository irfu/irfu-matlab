function tests_before_release_R2020a
% Function that runs multiple tests (matlab.unittest) under MATLAB R2020a
% specifically.
%
% RATIONALE
% ---------
% The primary use case for this file is BICAS (SolO/RPW/BIAS calibration code)
% which requires MATLAB R2019b (sic!) per agreement with LESIA/ROC which runs
% it. However, it is not possible to run tests under R2019b since MATLAB's
% "Action for Setting Up MATLAB on GitHub-Hosted Runner"
% (https://github.com/matlab-actions/setup-matlab/) only supports MATLAB R2020a
% and later. Therefore running BICAS tests on GitHub under MATLAB R2020a,
% despite it being suboptimal.

% Setup paths etc.
irf;

% Create empty list which is added to by later code.
suite = testsuite();

% Add tests for MATLAB packages in which automated test files can be
% automatically found via MATLAB's filenaming convention.
for pkgPathCa = {'bicas'}
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
  fullfile(ciPath, 'report_R2019b.pdf'), 'Verbosity', 3 ...
));

% RUN TESTS
assertSuccess(runner.run(suite));

end