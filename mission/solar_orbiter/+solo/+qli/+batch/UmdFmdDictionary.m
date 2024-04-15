%
% Dictionary for day (midnight, UTC) --> timestamp (no timezone time zone).
% Intended for storing one FMD per day of data, without assuming continuous time
% ranges of data.
%
% Keys are referred to as "Day", and values as FMD.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef UmdFmdDictionary
  % PROPOSAL: Better name,
  %   day, datasets, FMD, midnight, UTC
  %   DayFmdDictionary = DFD
  %   DataDayFmdDictionary = DDFD
  %   UmdFmdDictionary  = UFDictionary = UFD
  %
  % PROBLEM: Property "Dict" is not private despite that it should be.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Access=private)
    Dict
  end
  properties(Access=public)
    DaysDtArray
    FmdDtArray
    n
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Asserts that there are no duplicate constructor input data days.
    function obj = UmdFmdDictionary(varargin)

      if numel(varargin) == 0
        DatasetDaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;
        FmdDtArray         = datetime(cell(0, 1));
      elseif numel(varargin) == 2
        DatasetDaysDtArray = varargin{1};
        FmdDtArray         = varargin{2};
      else
        error('Illegal number of arguments.')
      end

      % ASSERTIONS
      solo.qli.utils.assert_UTC_midnight_datetime(DatasetDaysDtArray)
      assert(isa(FmdDtArray, 'datetime'))
      assert(strcmp(FmdDtArray.TimeZone, ''))
      irf.assert.sizes(...
        DatasetDaysDtArray, [-1], ...
        FmdDtArray,         [-1])
      % Unique days.
      assert(numel(unique(DatasetDaysDtArray)) == numel(DatasetDaysDtArray))

      obj.Dict = dictionary(DatasetDaysDtArray, FmdDtArray);
      % NOTE: dictionary does not store TimeZone for empty arrays, i.e.
      % obj.Dict.keys.TimeZone == ''.
    end



    function isKey = is_key(obj, DayDt)
      isKey = obj.Dict.isKey(DayDt);
    end



    function varargout = subsref(obj, S)

      switch S(1).type
        case '()'
          assert(isscalar(S))
          assert(isscalar(S.subs))

          key = S.subs{1};
          solo.qli.utils.assert_UTC_midnight_datetime(key);

          varargout = {obj.Dict(key)};

        case '.'
          % Call method (sic!)
          [varargout{1:nargout}] = builtin('subsref', obj, S);

        otherwise
          error('BICAS:Assertion', 'Does not support operation.')
      end

    end



    function obj = subsasgn(obj, S, FmdDt)
      assert(isscalar(S))
      assert(strcmp(S.type, '()'))
      assert(isscalar(S.subs))

      key = S.subs{1};

      solo.qli.utils.assert_UTC_midnight_datetime(key);
      assert(isa(FmdDt, 'datetime'))
      assert(strcmp(FmdDt.TimeZone, ''))

      obj.Dict(key) = FmdDt;
    end



    function isEqual = isequal(obj, otherUfd)
      if ~isa(obj, 'solo.qli.batch.UmdFmdDictionary')
        isEqual = false;
      elseif ~isa(otherUfd, 'solo.qli.batch.UmdFmdDictionary')
        isEqual = false;
      else
        isEqual = isequal(obj.Dict, otherUfd.Dict);
      end
    end



    % ~Generic dictionary utility function.
    %
    % Set one dictionary key value. In case of key collision, use the smaller of
    % the old and new value.
    %
    function obj = set_if_smaller(obj, DayDt, FmdDt)
      S = struct('type', '()', 'subs', {{DayDt}});

      if obj.is_key(DayDt)
        oldValue = obj.subsref(S);
        if oldValue > FmdDt
          obj = obj.subsasgn(S, FmdDt);
        end
      else
        obj = obj.subsasgn(S, FmdDt);
      end
    end



    % ~Generic dictionary utility function.
    %
    % Set one dictionary key value. In case of key collision, use the larger of
    % the old and new value.
    %
    function obj = set_if_greater(obj, DayDt, FmdDt)
      S = struct('type', '()', 'subs', {{DayDt}});

      if obj.is_key(DayDt)
        oldValue = obj.subsref(S);
        if oldValue < FmdDt
          obj = obj.subsasgn(S, FmdDt);
        end
      else
        obj = obj.subsasgn(S, FmdDt);
      end
    end



  end    % methods(Access=public)
  methods



    % NOTE: Does not guarantee order, except that it is consistent with
    % get.FmdDtArray().
    function DaysDtArray = get.DaysDtArray(obj)
      DaysDtArray = obj.Dict.keys;

      % Normalize TimeZone (empty for empty array).

      % IMPLEMENTATION NOTE: dictionary does not store TimeZone for empty
      % arrays, i.e. obj.Dict.keys.TimeZone == ''.
      if isempty(DaysDtArray)
        DaysDtArray.TimeZone = 'UTCLeapSeconds';
      end
    end



    % NOTE: Does not guarantee order, except that it is consistent with
    % get.DaysDtArray().
    function FmdDtArray = get.FmdDtArray(obj)
      FmdDtArray = obj.Dict.values;
    end



    function n = get.n(obj)
      n = obj.Dict.numEntries;
    end



  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Merge arbitrary number of UFDs into one. Only keep the greatest value
    % when there are key collisions.
    %
    function Ufd = merge_max(UfdCa)
      assert(iscell(UfdCa))

      OutputDfmmd = solo.qli.batch.UmdFmdDictionary();

      for iUfd = 1:numel(UfdCa)
        InputUfd  = UfdCa{iUfd};

        DayDtArray = InputUfd.DaysDtArray();
        FmdDtArray = InputUfd.FmdDtArray();

        for iKey = 1:numel(DayDtArray)
          OutputDfmmd = OutputDfmmd.set_if_greater(DayDtArray(iKey), FmdDtArray(iKey));
        end
      end

      Ufd = OutputDfmmd;
    end



  end    % methods(Static)



end
