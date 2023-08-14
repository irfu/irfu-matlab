function tests_before_release_R2019b
% Function that runs multiple tests (matlab.unittest) under MATLAB R2019b
% specifically. The primary example of this is BICAS (SolO/RPW/BIAS calibration
% code) which requires MATLAB R2019b per agreement with LESIA/ROC, but there
% could be other cases too.

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
% ---------
% NOTE: Can not run assertSuccess() since it is not supported on MATLAB R2019b.
% Below assertion should implement the same function.
TestResultArray = runner.run(suite);
assert(all([TestResultArray(:).Passed]))

end