function tests_before_release_BICAS
% Function which runs multiple tests (matlab.unittest) for BICAS specifically.
%
%
% RATIONALE FOR HAVING SEPARATE TESTS
% -----------------------------------
% (1) BICAS (SolO/RPW/BIAS Calibration Software) is required to support a
% specific MATLAB version per agreement with LESIA/ROC which runs it. As of
% 2024-05-28, this MATLAB version is planned to be MATLAB R2024a, and therefore
% tests should run on that version.
% (2) This code can also be run more seldomly (e.g. only for branch SOdevel) to
% save on limited server resources for CI.

% Setup paths etc.
irf;

% Create empty list which is added to by later code.
suite = testsuite();

% Add tests for MATLAB packages in which automated test files can be
% automatically found via MATLAB's filenaming convention.
%
% NOTE: Searches for tests recursively.
% NOTE: Includes tests for non-BICAS packages which BICAS uses.
for pkgPathCa = {'bicas', 'solo.adm', 'solo.hwzv'}
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
  fullfile(ciPath, 'report_BICAS.pdf'), 'Verbosity', 3 ...
  ));

% RUN TESTS
assertSuccess(runner.run(suite));

end
