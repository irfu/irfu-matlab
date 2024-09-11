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
%
function main(outputParentDir, letterVersion, configFile)
% TODO-DEC: Name for package/tool?
%     BICAS
%     ROC
%     RCS test data package
%     create
%     --
%     rcstestpkg
%     rcstdpkg  (TD = Test Data)
%     rtdp      (= RCS Test Data Package)
%
% PROPOSAL: Zip package.
%   CON: Can not manually update readme.txt, release_notes.txt
%   NOTE: There is zip support in MATLAB.
%         https://se.mathworks.com/help/matlab/ref/zip.html
%   PROPOSAL: Separate command.
%     CON: ~Superfluous?
%
% PROPOSAL: Use bicas.tools.batch functionality.
%   Ex: bicas.tools.batch.autocreate_input_BPCIs()
%
% PROPOSAL: Check the MATLAB version when calling BICAS.
%   TODO-NI/TODO-DEC: Where is this authoritatively specified where?
% PROPOSAL: Check the git repo version when calling BICAS.
%   PROPOSAL: Specify in config file.
%
% TODO-DEC: How handle BICAS call in tests?
%   PROPOSAL: Mock object.
%     CON: Overkill for such a simple application. Needs abstract object+2
%          subclasses.
%   PROPOSAL: Switch/flag for whether to call BICAS or not.
%
% PROPOSAL: Call test for SWD file.
%   CON: Related to BICAS deliveries, but unrelated to RCS test data packages.
% PROPOSAL: Verify existence of all .txt files.
%
% PROPOSAL: Check that using the correct directory with source code (bicas_ROC
%           git repo). Specify in config file.
%   NOTE: Must be able to run the TEST code in arbitrary irfu-matlab directory.
%         Tests can specify the current irfu-matlab directory when generating
%         the config file.
%
% PROBLEM: Not checking BICAS bash file.
%
% =============================================================================
% EXAMPLE DIRECTORY TREE
% =============================================================================
% TESTDATA_RODP_BICAS_V8.2.1A
% ├── LFR-SBM1-CWF-E
% │   ├── expected_outputs
% │   │   ├── manifest.txt
% │   │   └── solo_L2_rpw-lfr-sbm1-cwf-e-cdag_20210715T234148-20210715T235548_V02.cdf
% │   └── inputs
% │       ├── manifest.txt
% │       ├── solo_HK_rpw-bia_20210715_V05.cdf
% │       ├── solo_L1R_rpw-lfr-sbm1-cwf-e-cdag_20210715T234148-20210715T235548_V02.cdf
% │       └── solo_L1_rpw-bia-current-cdag_20210701-20210731_V34.cdf
% ├── LFR-SBM2-CWF-E
% │   ├── expected_outputs
% │   │   ├── manifest.txt
% │   │   └── solo_L2_rpw-lfr-sbm2-cwf-e-cdag_20220922T134335-20220922T154536_V01.cdf
% │   └── inputs
% │       ├── manifest.txt
% │       ├── solo_HK_rpw-bia_20220922_V06.cdf
% │       ├── solo_L1R_rpw-lfr-sbm2-cwf-e-cdag_20220922T134335-20220922T154536_V01.cdf
% │       └── solo_L1_rpw-bia-current-cdag_20220901-20220930_V37.cdf
% ├── LFR-SURV-CWF-E
% │   ├── expected_outputs
% │   │   ├── manifest.txt
% │   │   └── solo_L2_rpw-lfr-surv-cwf-e-cdag_20200213_V10.cdf
% │   └── inputs
% │       ├── manifest.txt
% │       ├── solo_HK_rpw-bia_20200213_V07.cdf
% │       ├── solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200213_V10.cdf
% │       └── solo_L1_rpw-bia-current-cdag_20200211-20200229_V03.cdf
% ├── LFR-SURV-SWF-E
% │   ├── expected_outputs
% │   │   ├── manifest.txt
% │   │   └── solo_L2_rpw-lfr-surv-swf-e-cdag_20200213_V10.cdf
% │   └── inputs
% │       ├── manifest.txt
% │       ├── solo_HK_rpw-bia_20200213_V07.cdf
% │       ├── solo_L1R_rpw-lfr-surv-swf-e-cdag_20200213_V10.cdf
% │       └── solo_L1_rpw-bia-current-cdag_20200211-20200229_V03.cdf
% ├── TDS-LFM-CWF-E
% │   ├── expected_outputs
% │   │   ├── manifest.txt
% │   │   └── solo_L2_rpw-tds-lfm-cwf-e-cdag_20200225_V07.cdf
% │   └── inputs
% │       ├── manifest.txt
% │       ├── solo_HK_rpw-bia_20200225_V06.cdf
% │       ├── solo_L1R_rpw-tds-lfm-cwf-e-cdag_20200225_V07.cdf
% │       └── solo_L1_rpw-bia-current-cdag_20200211-20200229_V02.cdf
% ├── TDS-LFM-RSWF-E
% │   ├── expected_outputs
% │   │   ├── manifest.txt
% │   │   └── solo_L2_rpw-tds-lfm-rswf-e-cdag_20200409_V08.trimmed.cdf
% │   └── inputs
% │       ├── manifest.txt
% │       ├── solo_HK_rpw-bia_20200409_V08.cdf
% │       ├── solo_L1R_rpw-tds-lfm-rswf-e-cdag_20200409_V08.trimmed.cdf
% │       └── solo_L1_rpw-bia-current-cdag_20200401-20200430_V03.cdf
% ├── readme.txt
% └── release_notes.txt



Bso = bicas.create_default_BSO();
Bso.make_read_only();

Swml = bicas.swm.get_SWML(Bso);

Config = bicas.tools.rcstestpkg.Config(configFile);

% Create root directory.
pkgDirName = bicas.tools.rcstestpkg.misc.create_test_package_directory_name(letterVersion);
pkgDir = bicas.tools.rcstestpkg.misc.mkdir(outputParentDir, pkgDirName);

bicas.tools.rcstestpkg.misc.create_readme_file(pkgDir)
bicas.tools.rcstestpkg.misc.create_release_notes_file(pkgDir, letterVersion)

for iSwm = 1:numel(Swml.List)
  Swm = Swml.List(iSwm);

  bicas.tools.rcstestpkg.misc.create_SWM_directory(pkgDir, Swm, Config)
end

end
