%
% matlab.unittest automatic test code for bicas.ga.get_output_dataset_GAs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_output_dataset_GAs___UTEST < matlab.unittest.TestCase



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
  end



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and teardown
    % methods which store/read their own data from the testCase object.
    % dir
    L
    Bso
  end



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    % Technically, additional properties of testCase objects with cell array
    % default values. Test methods with arguments with the same name will be
    % called once for every element in the cell arrays.
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      % Fixture = testCase.applyFixture(...
      %   matlab.unittest.fixtures.TemporaryFolderFixture);
      % % NOTE: The same fixture should always return the same directory.
      % testCase.dir = Fixture.Folder;
      testCase.L   = bicas.Logger('NO_STDOUT', false);
      testCase.Bso = bicas.create_default_BSO();
      testCase.Bso.make_read_only()
    end



  end



  %##########
  %##########
  % TEARDOWN
  %##########
  %##########
  methods(TestMethodTeardown)



    function teardown(testCase)
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      %=================
      % 3x InputDataset
      %=================

      InputDatasetsMap = containers.Map();

      InputDatasetsMap('HK_cdf') = bicas.InputDataset( ...
        struct(), struct(), struct(), struct(), ...
        '/subdir_HK/solo_L1_rpw-bia-current-cdag_20200211-20200229_V02.cdf');

      InputDatasetsMap('CUR_cdf') = bicas.InputDataset( ...
        struct(), struct(), struct(), struct(), ...
        '/subdir_CUR/solo_L1_rpw-bia-current-cdag_20200211-20200229_V02.cdf');

      InputDatasetsMap('SCI_cdf') = bicas.InputDataset( ...
        struct(), struct(), struct(), struct(), ...
        '/subdir_DCI/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200213_V07.cdf');

      %===============
      % OutputDataset
      %===============
      Epoch = int64([0; 1]);
      Zv    = struct('Epoch', Epoch);
      Ga    = struct('OBS_ID', {' '}, 'SOOP_TYPE', {' '});
      % NOTE: LFR RCT does not have values for CAL_ENTITY_NAME,
      % CAL_ENTITY_AFFILIATION, CAL_EQUIPMENT. Code should handle this
      % gracefully.
      RctdBias = bicas.proc.L1L2.cal.rct.RctDataTest( ...
        'solo_CAL_rpw-bias_20200210-20991231_V01.cdf', ...
        '01', 'BIAS team', 'Swedish Institute of Space Physics (IRF)', 'BIAS');
      RctdLfr = bicas.proc.L1L2.cal.rct.RctDataTest( ...
        'SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf', ...
        ' ', [], [], []);
      RctdCa = {RctdBias; RctdLfr};
      OutputDataset = bicas.OutputDataset(Zv, Ga, RctdCa);

      outputFilename = 'solo_L2_rpw-lfr-surv-cwf-e_20200213_V01.cdf';
      outputDsi      = 'SOLO_L2_RPW-LFR-SURV-CWF-E';

      %##################
      % CALL TESTED CODE
      %##################
      ActOutGaSubset = bicas.ga.get_output_dataset_GAs(...
        InputDatasetsMap, OutputDataset, outputFilename, outputDsi, ...
        testCase.Bso, testCase.L);

      %######################################################################
      % Examining return value
      % NOTE: Not checking most of the return value. Can add more as needed.
      %######################################################################

      function test_cell_column(value)
        testCase.assertTrue(iscolumn(value))
        testCase.assertTrue(iscell(  value))
      end
      function test_RctdCa_size(gaValue)
        test_cell_column(gaValue)
        testCase.assertEqual(size(gaValue), size(RctdCa))
      end

      test_RctdCa_size(ActOutGaSubset.CAL_ENTITY_AFFILIATION)
      test_RctdCa_size(ActOutGaSubset.CAL_ENTITY_NAME)
      test_RctdCa_size(ActOutGaSubset.CAL_EQUIPMENT)
      test_RctdCa_size(ActOutGaSubset.CALIBRATION_TABLE)
      test_RctdCa_size(ActOutGaSubset.CALIBRATION_VERSION)

      test_cell_column(ActOutGaSubset.Parents)
      testCase.assertEqual(size(ActOutGaSubset.Parents), [3,1])

      testCase.assertEqual(ActOutGaSubset.CALIBRATION_VERSION{1}, '01')
      testCase.assertEqual(ActOutGaSubset.CALIBRATION_VERSION{2}, ' ')
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)
  end    % methods(Access=private)



end
