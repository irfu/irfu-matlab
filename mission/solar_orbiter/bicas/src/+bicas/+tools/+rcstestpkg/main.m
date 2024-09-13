%
% Create an (unzipped) RCS test package for BICAS.
%
% Every official BICAS delivery to ROC must, by agreement, be accompanied with
% a "test package" containing CDF files corresponding to the BICAS input &
% output for all official SWMs, i.e. L1R->L2.
%
% RCS test data packages are described in "ROC Engineering Guidelines For
% External Users", ROC-GEN-SYS-NTT-00019-LES, Iss.02, Rev.01 draft, Section 2.3
% (which is quite different from iss2rev0).
%
% See readme.txt for example description of the directory structure in an RCS
% test data package.
%
%
% ARGUMENTS
% =========
% outputParentDir
% letterVersion
%       One capital letter to be used in the test package version.
%       NOTE: The BICAS version is automatically used.
% configFile
%       JSON file with paths to input CDFs. See
%       bicas.tools.rcstestpkg.main___UTEST for an "example".
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
function main(outputParentDir, letterVersion, configFile)
% PROPOSAL: Better name for function.
%     create
%     create_package
% PROPOSAL: Better name for package/tool?
%     BICAS
%     ROC
%     "RCS test data package" (~formal term)
%     create
%     --
%     rcstestpkg
%     rcstdpkg  (TD = Test Data)
%     rtdp      (= RCS Test Data Package)
%
bicas.tools.rcstestpkg.misc.create_RCS_test_pkg(...
  outputParentDir, letterVersion, configFile, false)
end
