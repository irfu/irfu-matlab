%
% Implementation for automatic tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BicasProcessingAccessTest < bicas.tools.batch.BicasProcessingAccessAbstract
  % PROPOSAL: Support raising exception.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=private)
    SwmArray
    errorInputFilesCaCa
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % errorInputFilesCaCa
    %       Data structure for specifying which calls should end with an
    %       emulated error.
    %       {i}{j} = Char array path to input file j. Order within {i} does not
    %       matter.
    %       When the set of input paths is equal to {i}, then the class will
    %       emulate an BICAS error.
    %       IMPLEMENTATION NOTE: Must base the criterion for error on the input
    %       since the class must work also for parallelized code. Using only the
    %       call order number ("the n'th call") would not work.
    %
    function obj = BicasProcessingAccessTest(SwmArray, errorInputFilesCaCa)
      % ASSERTIONS
      assert(isa(SwmArray, 'bicas.swm.SoftwareMode') & iscolumn(SwmArray))
      assert(iscell(errorInputFilesCaCa) & iscolumn(errorInputFilesCaCa))
      for i = 1:numel(errorInputFilesCaCa)
        errorInputFilesCa = errorInputFilesCaCa{i};
        assert(iscell(errorInputFilesCa) & iscolumn(errorInputFilesCa))
      end

      obj.SwmArray            = SwmArray;

      obj.errorInputFilesCaCa = errorInputFilesCaCa;
    end



    % OVERRIDE
    function [varargout] = bicas_main(obj, varargin)
      % ASSERTIONS
      for i = 1:numel(varargin)
        assert(ischar(varargin{i}))
      end

      [inputPathsCa, outputPathsCa] = obj.parse_CLI_arguments(varargin);

      % ASSERTION: Input dataset exists.
      for i = 1:numel(inputPathsCa)
        irf.assert.file_exists(inputPathsCa{i})
      end


      if obj.is_non_zero_error_call(inputPathsCa)
        %===========================
        % CASE: Non-zero error code
        %===========================
        % Do no processing (generate no output files)
        % -------------------------------------------
        % NOTE: One could imagine simulating an error after or between output
        % datasets are generated but that should be overkill.

        [varargout{1}] = 1;
      else
        %=====================================
        % CASE: No error. Emulate processing.
        %=====================================

        % Output datasets: Create empty output datasets
        for i = 1:numel(outputPathsCa)
          irf.fs.write_empty_file({outputPathsCa{i}});
        end

        [varargout{1}] = 0;
      end

    end    % bicas_main



  end    % methods(Access=public)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function [inputPathsCa, outputPaths] = parse_CLI_arguments(obj, cliArgCa)
      %==============
      % Identify SWM
      %==============
      swmCliOption = cliArgCa{1};
      iSwm = find(strcmp(swmCliOption, {obj.SwmArray(:).cliOption}));
      assert(isscalar(iSwm), ...
        'Did not find exactly one SWM which matches CLI arg. "%s".', ...
        swmCliOption)
      Swm  = obj.SwmArray(iSwm);

      %================
      % Input datasets
      %================
      inputPathsCa = cell(0, 1);
      for iInput = 1:numel(Swm.inputsList)
        cohb      = Swm.inputsList(iInput).cliOptionHeaderBody;
        iCoh      = find(strcmp(['--', cohb], cliArgCa));
        assert(isscalar(iCoh), ...
          'Can not identify exactly one BICAS argument with cohb="%s"', cohb)
        inputPathsCa{end+1, 1} = cliArgCa{iCoh+1};
      end

      %=================
      % Output datasets
      %=================
      outputPaths = cell(0, 1);
      for iOutput = 1:numel(Swm.outputsList)
        cohb                  = Swm.outputsList(iOutput).cliOptionHeaderBody;
        iCoh                  = find(strcmp(['--', cohb], cliArgCa));
        outputPaths{end+1, 1} = cliArgCa{iCoh+1};

        % NOTE: Can not assert that output datasets do not pre-exist,
        % since bicas.batch.main permits overwriting files.
      end
    end



    function isNonZeroErrorCall = is_non_zero_error_call(obj, inputPathsCa)
      for i = 1:numel(obj.errorInputFilesCaCa)
        errorInputFilesCa = obj.errorInputFilesCaCa{i};
        if isempty(setdiff(inputPathsCa, errorInputFilesCa))
          isNonZeroErrorCall = true;
          return
        end
      end

      isNonZeroErrorCall = false;
    end



  end    % methods(Access=private)



end
