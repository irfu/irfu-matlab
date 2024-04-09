%
% Implementation of abstract superclass that passes on method calls to the
% nominal code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef GenerateQuicklooksImplementation < solo.qli.batch.GenerateQuicklooksAbstract



  %##########################################################
  %##########################################################
  % PUBLIC INSTANCE METHODS THAT OVERRIDE SUPERCLASS METHODS
  %##########################################################
  %##########################################################
  methods(Access=public)



    function generate_quicklooks_24h_6h_2h_using_DB_SPICE(obj, varargin)
      solo.qli.generate_quicklooks_24h_6h_2h_using_DB_SPICE(varargin{:})
    end



    function generate_quicklook_7days_using_DB_SPICE(obj, varargin)
      solo.qli.generate_quicklook_7days_using_DB_SPICE(varargin{:})
    end



  end    % methods(Access=public)



end
