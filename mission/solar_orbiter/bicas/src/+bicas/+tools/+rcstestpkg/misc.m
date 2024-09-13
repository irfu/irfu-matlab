%
% Collection of various functions. Name is likely temporary. Exact set of
% functions is likely temporary until a better organization of code is found.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef misc



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    % PROPOSAL: Copy from bicas.const.SWD_METADATA('SWD.release.author').
    CONTACT_PERSON = bicas.const.SWD_METADATA('SWD.release.author');
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function create_RCS_test_pkg(outputParentDir, letterVersion, configFile, automatedTestRun)
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

      assert(islogical(automatedTestRun))



      Bso = bicas.create_default_BSO();
      Bso.make_read_only();

      Swml = bicas.swm.get_SWML(Bso);

      Config = bicas.tools.rcstestpkg.Config(configFile);

      % Create root directory.
      pkgDirName = bicas.tools.rcstestpkg.misc.create_test_package_directory_name(letterVersion);
      pkgDir     = bicas.tools.rcstestpkg.misc.mkdir(outputParentDir, pkgDirName);

      bicas.tools.rcstestpkg.misc.create_readme_file(pkgDir)
      bicas.tools.rcstestpkg.misc.create_release_notes_file(pkgDir, letterVersion)

      for iSwm = 1:numel(Swml.List)
        Swm = Swml.List(iSwm);

        bicas.tools.rcstestpkg.misc.create_SWM_directory(pkgDir, Swm, Config)
      end

    end



    function create_SWM_directory(parentDir, Swm, Config)
      swmDir = bicas.tools.rcstestpkg.misc.mkdir(parentDir, Swm.cliOption);

      inputsDir  = bicas.tools.rcstestpkg.misc.mkdir(swmDir, 'inputs');
      outputsDir = bicas.tools.rcstestpkg.misc.mkdir(swmDir, 'expected_outputs');

      bicasArgsCa        = {Swm.cliOption};

      BpciInputDsmdArray = solo.adm.DSMD.empty(0, 1);
      inputFilenamesCa   = cell(0, 1);
      inputCohbCa        = cell(0, 1);
      for i = 1:numel(Swm.inputsList)
        InputDataset = Swm.inputsList(i);
        cohb         = InputDataset.cliOptionHeaderBody;

        inputSrcFile = Config.get_input_file(Swm.cliOption, cohb);
        irf.assert.file_exists(inputSrcFile)

        bicasArgsCa = [bicasArgsCa; {['--', cohb]}; {inputSrcFile}];

        inputFilenamesCa{end+1, 1} = irf.fs.get_name(inputSrcFile);
        inputCohbCa     {end+1, 1} = cohb;

        % Store information for later being able to automatically construct the
        % output filename.
        InputDsmd = solo.adm.paths_to_DSMD_array({inputSrcFile});
        % Check that filename can be parsed.
        assert(isscalar(InputDsmd), ...
          'Can not parse the filename of file "%s".', inputSrcFile)
        BpciInputDsmdArray(end+1, 1) = InputDsmd;

        copyfile(inputSrcFile, inputsDir)
      end
      bicas.tools.rcstestpkg.misc.create_manifest_file(...
        inputsDir, inputFilenamesCa, inputCohbCa)



      outputFilenamesCa = cell(0, 1);
      outputCohbCa      = cell(0, 1);
      for i = 1:numel(Swm.outputsList)
        OutputDataset = Swm.outputsList(i);
        cohb          = OutputDataset.cliOptionHeaderBody;

        % Automatically construct output dataset filename.
        % NOTE: Uses bicas.tools.batch function.
        outputFile    = bicas.tools.batch.get_BPCI_output_path2(...
          BpciInputDsmdArray, solo.adm.DSMD.empty(0, 1), ...
          OutputDataset.dsi, 'HIGHEST_USED', outputsDir, false);

        irf.assert.path_is_available(outputFile)

        bicasArgsCa = [bicasArgsCa; {['--', cohb]}; {outputFile}];

        outputFilenamesCa{end+1, 1} = irf.fs.get_name(inputSrcFile);
        outputCohbCa     {end+1, 1} = cohb;
      end
      bicas.tools.rcstestpkg.misc.create_manifest_file(...
        outputsDir, outputFilenamesCa, outputCohbCa)



      %errorCode = bicas.main(bicasArgsCa{:})
      %assert(errorCode == 0)
    end



    function pkgDirName = create_test_package_directory_name(letterVersion)
      irf.assert.castring_regexp(letterVersion, '[A-Z]')

      % Ex: TESTDATA_RODP_BICAS_V8.2.1A
      bicasVerStr = bicas.const.SWD_METADATA('SWD.release.version');

      pkgDirName = sprintf('TESTDATA_RODP_BICAS_V%s%s', bicasVerStr, letterVersion);
    end



    function newDirPath = mkdir(parentDir, dirName)
      newDirPath = fullfile(parentDir, dirName);

      % IMPLEMENTATION NOTE: Checking for pre-existing directory for preventing
      % the user from overwriting(?) a pre-existing test package.
      irf.assert.path_is_available(newDirPath)

      [success, message, errorMessageId] = mkdir(parentDir, dirName);
      assert(success, errorMessageId, message)
    end



    function create_manifest_file(datasetDir, inputFilenamesCa, cohbCa)
      assert(iscell(inputFilenamesCa))
      assert(iscell(cohbCa))
      nDatasets = irf.assert.sizes(...
        inputFilenamesCa, [-1, 1], ...
        cohbCa,           [-1, 1]);

      manifestFile = fullfile(datasetDir, "manifest.txt");

      fid = fopen(manifestFile, "w");
      for i = 1:nDatasets
        inputFilename = inputFilenamesCa{i};
        cohb          = cohbCa{i};
        fprintf(fid, "%s\t%s\n", inputFilename, cohb);
      end
      fclose(fid);
    end



    function create_readme_file(parentDir)
      s = sprintf('Contact person; %s', bicas.tools.rcstestpkg.misc.CONTACT_PERSON);

      readmePath = fullfile(parentDir, "readme.txt");
      irf.fs.write_file(readmePath, uint8(s)')
    end



    function create_release_notes_file(parentDir, letterVersion)
      dateStr = string(datetime('now','TimeZone','local','Format','yyyy-mm-d'));
      s = sprintf('Version %s; %s; %s:', letterVersion, dateStr, bicas.tools.rcstestpkg.misc.CONTACT_PERSON);

      readmePath = fullfile(parentDir, "release_notes.txt");
      irf.fs.write_file(readmePath, uint8(s)')
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)
  end    % methods(Static, Access=private)



end
