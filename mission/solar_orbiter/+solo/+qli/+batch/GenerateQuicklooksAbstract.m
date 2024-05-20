%
% Abstract code which, thorugh subclasses, provides access to code which
% generates actual quicklooks (files). This class exists so that the dependence
% on this code can be exchanged for custom implementations for tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef GenerateQuicklooksAbstract < handle
  % IMPLEMENTATION NOTE: "Must" be a handle class to store updateable
  % information on calls to method.



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract)

    generate_quicklooks_24h_6h_2h_using_DB_SPICE(obj, varargin)

    generate_quicklook_7days_using_DB_SPICE(obj, varargin)

  end



end
