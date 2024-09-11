%
% Store for the content of the config file.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Config
  % PROPOSAL: Automatic test code.



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



    function inputFile = get_input_file(obj, swmCliOption, inputCohb)
      SwmJson = obj.JsonStruct.(obj.JSON_key_str(swmCliOption));

      inputFile = SwmJson.(obj.JSON_key_str(inputCohb));
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



    function jsonKeyStr = JSON_key_str(s)
      s(s=='-') = '_';
      jsonKeyStr = s;
    end



  end    % methods(Static, Access=private)



end
