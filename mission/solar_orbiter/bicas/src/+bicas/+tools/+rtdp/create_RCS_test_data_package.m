% Function for creating
% (1) an (unzipped) RTDP, and
% (2) a zipped RDTP.
%
% Every official BICAS delivery to ROC must, by agreement, be accompanied with
% an "RCS test data package" (RDTP) containing CDF files corresponding to the
% BICAS input & output for all official SWMs, i.e. L1R->L2.
%
% RTDPs are described in "ROC Engineering Guidelines For
% External Users", ROC-GEN-SYS-NTT-00019-LES, Iss.02, Rev.01 draft, Section 2.3
% (which is quite different from iss2rev0).
%
% See readme.txt for example description of the directory structure in an RCS
% test data package.
%
% ~BUG: For letter versions other than "A", release_notes.txt will only contain
%       a reference to the specified version and not previous ones. For this
%       case, the file must be added manually and the zipped RDTP must be
%       created manually.
%
% NOTE: Never calls bash wrapper script "bicas", which that it is also never
%       tested (should also not be necessary).
%
%
% ARGUMENTS
% =========
% outputParentDir
% letterVersion
%       One capital letter to be used in the test package version. The first
%       RTDP for any given BICAS version must be "A", the second "B" and so on.
%       NOTE: The BICAS version does not need to specified since it is
%       automatically obtained from BICAS (bicas.const).
% configFile
%       JSON file with (1) paths to input CDFs, and (2) expected BICAS source
%       directory. See bicas.tools.rtdp.misc___UTEST.create_config_file()
%       for an "example".
%
%
% RETURN VALUES
% =============
% rtdpDir
% rtdpZippedFile
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
function [rtdpDir, rtdpZipFile] = create_RCS_test_data_package(...
  outputParentDir, letterVersion, configFile)
% bicas.tools.rtdp.create_RCS_test_data_package('/nonhome_data/SOLAR_ORBITER/bicas_test_packages/temp', 'A', '/nonhome_data/work_files/SOLAR_ORBITER/rtdp_config.json')
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
%     rtdp
%     testpkg
%     rcstdpkg  (TD = Test Data)
%     rtdp      (= RCS Test Data Package)

tic
[rtdpDir, rtdpZipFile] = bicas.tools.rtdp.misc.create_RCS_test_data_package(...
  outputParentDir, letterVersion, configFile, false);
toc
end
