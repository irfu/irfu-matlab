%
% Info required for processing datasets using one single BICAS call, but not the
% argument list as such.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-05-27.
%
classdef BicasProcessingCallInfo
    % PROPOSAL: Better name. "info" is too generic.
    %   NOTE: Compare bicas.tools.batch.BicasProcessingCallSummary.
    %   NOTE: Only includes information on SWM + input + output, but excludes
    %         custom settings. Excludes --version, --swdescriptor, --config etc.
    %   PROPOSAL: Only imply specifying datasets (paths) in & out, not
    %             specifying the entire call.
    %   --
    %   BICAS
    %   Datasets
    %   Paths
    %   Input/output (datasets)
    %   I/O
    %   Processing (as opposed to non-processing)
    %   Data, Info, Record
    %   --
    %   BicasCallPathsIO, BicasCallDatasetsIO
    %
    % PROPOSAL: Be able to generate BPCI from SWM for testing purposes.
    %   See bicas.tools.batch.autocreate_one_SWM_BPCI().
    %   CON: Still needs to be able to associate COBHs with paths and filenames.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        swmCliOption
        inputsArray
        outputsArray
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = BicasProcessingCallInfo(swmCliOption, inputsArray, outputsArray)
            assert(ischar(swmCliOption))
            assert(iscolumn(inputsArray)  & isa(inputsArray,  'bicas.tools.batch.BpciInput'))
            assert(iscolumn(outputsArray) & isa(outputsArray, 'bicas.tools.batch.BpciOutput'))

            obj.swmCliOption = swmCliOption;
            obj.inputsArray  = inputsArray;
            obj.outputsArray = outputsArray;
        end



        function filenameCa = get_output_filenames(obj)
            filenameCa = cellfun(...
                @irf.fs.get_name, {obj.outputsArray.path}', ...
                'UniformOutput', false);
            filenameCa = filenameCa(:);
        end



        % Create modified copy of object. Prepend all output paths with a
        % specified path.
        function obj = prepend_output_directory(obj, outputDir)
            outputsArray = bicas.tools.batch.BpciOutput.empty(0, 1);
            for i = 1:numel(obj.outputsArray)
                o = obj.outputsArray(i);
                outputsArray(i, 1) = bicas.tools.batch.BpciOutput(...
                    o.cohb, o.dsi, fullfile(outputDir, o.path));
            end

            obj = bicas.tools.batch.BicasProcessingCallInfo(...
                obj.swmCliOption, obj.inputsArray, outputsArray);
        end



    end    % methods(Access=public)



end
