%
% Store for the content of one config file.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Config



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=private)
    JsonStruct
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = Config(configFile)
      uint8Array     = irf.fs.read_file(configFile);
      obj.JsonStruct = jsondecode(char(uint8Array)');
    end



    function inputFile = get_input_dataset(obj, swmCliOption, inputCohb)
      jsonSwmCliOption = obj.JSON_key_str(swmCliOption);
      if ~isfield(obj.JsonStruct.inputDatasets, jsonSwmCliOption)
        error('Can not find SWM "%s" in configuration.', swmCliOption)
      end
      SwmJson = obj.JsonStruct.inputDatasets.(jsonSwmCliOption);

      jsonInputCohb = obj.JSON_key_str(inputCohb);
      if ~isfield(SwmJson, jsonInputCohb)
        error('Can not find input dataset "%s" for SWM "%s" in configuration.', ...
          inputCohb, swmCliOption)
      end
      inputFile = SwmJson.(jsonInputCohb);
    end



    function bicasRootDir = get_BICAS_root_dir(obj)
      bicasRootDir = obj.JsonStruct.bicasRootDir;
    end



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)
  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Convert arbitrary string to string that will be used by MATLAB as a
    % corresponding field name in the structs which represent the content of a
    % JSON file. Since MATLAB represents the content of JSON files as structs,
    % the JSON keys can not be arbitrary strings, but are converted to strings
    % which are legal struct field names.
    function jsonKeyStr = JSON_key_str(s)
      s(s=='-')  = '_';
      jsonKeyStr = s;
    end



  end    % methods(Static, Access=private)



end
