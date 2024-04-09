%
% Implementation of abstract superclass for tests.
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
    % generate_quicklooks_24h_6h_2h_using_DB_SPICE()/generate_quicklook_7days_using_DB_SPICE()
    % was called, in the order it was called.
    Dt24h6h2hArray
    Dt7daysArray

    % 1D column arrays of the dates for which
    % generate_quicklooks_24h_6h_2h_using_DB_SPICE()/generate_quicklook_7days_using_DB_SPICE()
    % shall raise exception.
    Dt24h6h2ExceptionArray
    Dt7daysExceptionArray
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = GenerateQuicklooksTest(varargin)
      obj.Dt24h6h2hArray = solo.qli.const.EMPTY_DT_ARRAY;
      obj.Dt7daysArray   = solo.qli.const.EMPTY_DT_ARRAY;

      switch(nargin)
        case 0
          obj.Dt24h6h2ExceptionArray = solo.qli.const.EMPTY_DT_ARRAY;
          obj.Dt7daysExceptionArray  = solo.qli.const.EMPTY_DT_ARRAY;
        case 2
          obj.Dt24h6h2ExceptionArray = varargin{1};
          obj.Dt7daysExceptionArray  = varargin{2};
        otherwise
          error('Illegal number of arguments.')
      end
    end



    function generate_quicklooks_24h_6h_2h_using_DB_SPICE(obj, Dt, varargin)
      obj.Dt24h6h2hArray = [obj.Dt24h6h2hArray; Dt];
      if ismember(Dt, obj.Dt24h6h2ExceptionArray)
        error('Test error')
      end
    end



    function generate_quicklook_7days_using_DB_SPICE(obj, Dt, varargin)
      obj.Dt7daysArray = [obj.Dt7daysArray; Dt];
      if ismember(Dt, obj.Dt7daysExceptionArray)
        error('Test error')
      end
    end



  end    % methods(Access=public)



end
