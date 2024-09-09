%
% RCTD class for test purposes.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctDataTest < bicas.proc.L1L2.cal.rct.RctDataAbstract



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = RctDataTest(varargin)
      obj@bicas.proc.L1L2.cal.rct.RctDataAbstract(varargin{:});
    end



  end    % methods(Access=public)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Implementation of abstract method.
    function log_RCT(obj, L)
    end



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Implementation of abstract method.
    function [RctRawData] = read_RCT(filePath)
      RctRawData = [];
    end



  end    % methods(Static)



end
