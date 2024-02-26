%
% Abstract class for handling the communication with BICAS for the purpose of
% processing data.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) BicasProcessingAccessAbstract < handle



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract)



    % Call bicas.main() with the exact same arguments and return value(s).
    [varargout] = bicas_main(obj, varargin);



  end    % methods(Access=public)



end
