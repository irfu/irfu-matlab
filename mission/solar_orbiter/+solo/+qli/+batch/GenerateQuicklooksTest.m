%
% Implementation of abstract superclass for the purpose of tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef GenerateQuicklooksTest < solo.qli.batch.GenerateQuicklooksAbstract



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private)
    % 1D column arrays of the dates for which
    % generate_quicklooks_24h_6h_2h_using_DB_SPICE()
    % /generate_quicklook_7days_using_DB_SPICE()
    % was called, in the order it was called.
    UmdDt24h6h2hCallsArray
    UmdDt7daysCallsArray

    % 1D column arrays of the dates for which
    % generate_quicklooks_24h_6h_2h_using_DB_SPICE()
    % /generate_quicklook_7days_using_DB_SPICE()
    % shall raise exception.
    UmdDt24h6h2ExceptionArray
    UmdDt7daysExceptionArray
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % Syntax 1: (no arguments)
    %       Object will raise no exceptions.
    % Syntax 2:
    %       varargin{1}
    %           UMD datetime array for which method
    %           generate_quicklooks_24h_6h_2h_using_DB_SPICE() should raise
    %           exception.
    %       varargin{2}
    %           UMD datetime array for which method
    %           generate_quicklook_7days_using_DB_SPICE() should raise
    %           exception.
    function obj = GenerateQuicklooksTest(varargin)
      obj.UmdDt24h6h2hCallsArray = solo.qli.const.EMPTY_DT_ARRAY;
      obj.UmdDt7daysCallsArray   = solo.qli.const.EMPTY_DT_ARRAY;

      switch(nargin)
        case 0
          obj.UmdDt24h6h2ExceptionArray = solo.qli.const.EMPTY_DT_ARRAY;
          obj.UmdDt7daysExceptionArray  = solo.qli.const.EMPTY_DT_ARRAY;
        case 2
          obj.UmdDt24h6h2ExceptionArray = varargin{1};
          obj.UmdDt7daysExceptionArray  = varargin{2};
        otherwise
          error('Illegal number of arguments.')
      end

      irf.dt.assert_UTC_midnight(obj.UmdDt24h6h2ExceptionArray)
      irf.dt.assert_UTC_midnight(obj.UmdDt7daysExceptionArray)
      assert(iscolumn(obj.UmdDt24h6h2ExceptionArray))
      assert(iscolumn(obj.UmdDt7daysExceptionArray))
    end



    function generate_quicklooks_24h_6h_2h_using_DB_SPICE(obj, Dt, varargin)
      obj.UmdDt24h6h2hCallsArray = [obj.UmdDt24h6h2hCallsArray; Dt];
      if ismember(Dt, obj.UmdDt24h6h2ExceptionArray)
        error('Test error')
      end
    end



    function generate_quicklook_7days_using_DB_SPICE(obj, Dt, varargin)
      obj.UmdDt7daysCallsArray = [obj.UmdDt7daysCallsArray; Dt];
      if ismember(Dt, obj.UmdDt7daysExceptionArray)
        error('Test error')
      end
    end



  end    % methods(Access=public)



end
