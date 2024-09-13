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
%       JSON file with (1) paths to input CDFs, and (2) expected BICAS source
%       directory. See bicas.tools.rcstestpkg.misc___UTEST.create_config_file()
%       for an "example".
%
%
% RETURN VALUES
% =============
% (none)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
function main(outputParentDir, letterVersion, configFile)
% bicas.tools.rcstestpkg.main('/nonhome_data/SOLAR_ORBITER/bicas_test_packages/temp', 'A', '/nonhome_data/work_files/SOLAR_ORBITER/rcstestpkg_config.json')
%
% PROPOSAL: Better name for function.
%     create
%     create_package
%     create_RCS_test_package
%     create_RCS_test_data_package
%     create_RTDP
% PROPOSAL: Better name for package/tool?
%     BICAS
%     ROC
%     "RCS test data package" (~formal term)
%     create
%     --
%     rcstestpkg
%     testpkg
%     rcstdpkg  (TD = Test Data)
%     rtdp      (= RCS Test Data Package)

bicas.tools.rcstestpkg.misc.create_RCS_test_pkg(...
  outputParentDir, letterVersion, configFile, false)
end
