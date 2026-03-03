%
% Run all of BICAS's automatic tests. This function is not called by BICAS
% proper. It is only intended to be called manually during development to
% trigger all automated tests related to BICAS.
%
% NOTE: Will ALWAYS only run in the CURRENT irfu-matlab git repo, not
% necessarily the irfu-matlab git repo used for running BICAS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-19.
%
function run_UTESTs()
tic

% IMPLEMENTATION NOTE: runtests('bicas', 'IncludeSubpackages', true); used to
% work, but stopped working after creating a file src/bicas.m, presumably
% because it clashes with src/+bicas/. Must therefore use another command.
% testResults = runtests(bicas.utils.get_BICAS_source_dir(), 'IncludeSubfolders', true);
testAr0 = matlab.unittest.TestSuite.fromPackage('bicas', 'IncludingSubpackages', true);



% IMPLEMENTATION NOTE: matlab.unittest.TestSuite.fromPackage('bicas', ...)
% empirically excludes the bicas___UTEST.m files which is outside of
% the +bicas/ directory. Must therefore specify it manually.
bicasTestPath = fullfile(bicas.utils.get_BICAS_source_dir(), "bicas___UTEST.m");
testAr1 = matlab.unittest.TestSuite.fromFile(bicasTestPath);



% NOTE: This path will not work for the directory configuration used when
% delivering BICAS to ROC.
psp2neTestPath = fullfile( ...
  bicas.utils.get_BICAS_root_dir(), "..", "+solo", "psp2ne___UTEST.m");
testAr2 = matlab.unittest.TestSuite.fromFile(psp2neTestPath);



testAr = [testAr0, testAr1, testAr2];

assertSuccess(matlab.unittest.TestRunner.withTextOutput.run(testAr));

toc
end
