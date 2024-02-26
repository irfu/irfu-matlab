%
% Implementation for automatic tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BicasProcessingAccessTest < bicas.tools.batch.BicasProcessingAccessAbstract
  % PROPOSAL: Support returning non-zero error code.
  % PROPOSAL: Support raising exception.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=private)
    SwmArray
    callNonZeroErrorArray
  end    % properties(SetAccess=immutable)
  properties(Access=private)
    % The number of times method "bicas_main" has been called.
    nCalls
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % callNonZeroErrorArray
    %       Array of numbers. Call numbers for when method bicas_main()
    %       should return non-zero error code and simulate failure.
    %       1=First call.
    %
    function obj = BicasProcessingAccessTest(SwmArray, callNonZeroErrorArray)
      assert(isa(SwmArray, 'bicas.swm.SoftwareMode') & iscolumn(SwmArray))
      assert(isnumeric(callNonZeroErrorArray))
      assert(iscolumn(callNonZeroErrorArray) || isempty (callNonZeroErrorArray))

      obj.SwmArray              = SwmArray;
      obj.callNonZeroErrorArray = callNonZeroErrorArray;
      obj.nCalls                = 0;
    end



    % OVERRIDE
    function [varargout] = bicas_main(obj, varargin)

      obj.nCalls = obj.nCalls + 1;
      iCall = obj.nCalls;
      if ismember(iCall, obj.callNonZeroErrorArray)
        %=============================
        % CASE: Return non-zero error
        %=============================

        % Do no processing (generate no output files)
        % -------------------------------------------
        % NOTE: One could imagine simulating an error after or between
        % output datasets are generated but that should be overkill.

        [varargout{1}] = 1;
      else
        %=======================================
        % CASE: Error code zero. Do processing.
        %=======================================

        for i = 1:numel(varargin)
          assert(ischar(varargin{i}))
        end

        swmCliOption = varargin{1};

        iSwm = find(strcmp(swmCliOption, {obj.SwmArray(:).cliOption}));
        assert(isscalar(iSwm), 'Did not find exactly one SWM which matches CLI arg. "%s".', swmCliOption)
        Swm  = obj.SwmArray(iSwm);



        %=========================================
        % Input datasets: Assert that files exist
        %=========================================
        for iInput = 1:numel(Swm.inputsList)
          cohb      = Swm.inputsList(iInput).cliOptionHeaderBody;
          iCoh      = find(strcmp(['--', cohb], varargin));
          assert(isscalar(iCoh), ...
            'Can not identify exactly one BICAS argument with cohb="%s"', cohb)
          inputPath = varargin{iCoh+1};

          % ASSERTION: Input dataset exists.
          irf.assert.file_exists(inputPath)
        end

        %===============================================
        % Output datasets: Create empty output datasets
        %===============================================
        for iOutput = 1:numel(Swm.outputsList)
          cohb           = Swm.outputsList(iOutput).cliOptionHeaderBody;
          iCoh           = find(strcmp(['--', cohb], varargin));
          outputFilename = varargin{iCoh+1};

          % NOTE: Can not assert that output datasets do not pre-exist,
          % since bicas.batch.main permits overwriting files.

          % CREATE EMPTY OUTPUT DATASET.
          irf.fs.create_empty_file(outputFilename);
        end

        [varargout{1}] = 0;
      end



    end



  end    % methods(Access=public)



end
