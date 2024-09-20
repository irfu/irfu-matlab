%
% Collection of shorter ~generic utility functions that could conceivably be
% used all over BICAS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-27, with moved from bicas.proc.utils.
%
classdef utils
  % PROPOSAL: More automatic test code.
  % PROPOSAL: Rename.
  %   PRO: "utils" implies that code is generic, while the code seems to *not* be.
  %
  % PROPOSAL: Replace functions returning paths with constants as far as is
  %           possible (all except bicas.utils.get_BICAS_root_dir()?).
  %   CON: Abandons any idea of temporarily changing the location of the default
  %        config file etc. for tests.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Get path to the root of the BICAS directory structure.
    function bicasRootDir = get_BICAS_root_dir()
      % ASSUMES: The current file is in the <BICAS>/src/+bicas/ directory.
      % Use path of the current MATLAB file.
      [matlabSrcDir, ~, ~] = fileparts(mfilename('fullpath'));
      bicasRootDir         = fullfile(matlabSrcDir, '..', '..');
      bicasRootDir         = irf.fs.get_abs_path(bicasRootDir);
    end



    function swdFile = get_SWD_file()
      swdFile = fullfile(...
        bicas.utils.get_BICAS_root_dir(), ...
        bicas.const.SWD_FILENAME);
    end



    function bicasConfigDir = get_BICAS_config_dir()
      bicasConfigDir = fullfile(...
        bicas.utils.get_BICAS_root_dir(), ...
        bicas.const.DEFAULT_CONFIG_DIR_RPATH);
    end



    function bicasDefaultConfigFile = get_BICAS_default_config_file()
      bicasDefaultConfigFile = fullfile(...
        bicas.utils.get_BICAS_root_dir(), ...
        bicas.const.DEFAULT_CONFIG_DIR_RPATH, ...
        bicas.const.DEFAULT_CONFIG_FILENAME);
    end



    % Whether two sets (arrays) of arbitrary objects are set equal.
    %
    % NOTE: Compare objects, not handles. NaN == NaN.
    % NOTE: Ignores duplicated objects within arrays.
    function equal = object_sets_isequaln(Ar1, Ar2)
      % PROPOSAL: Better name

      % IMPLEMENTATION NOTE: Should not be most optimal implementation,
      % but good enough.
      equal = bicas.utils.is_subset_isequaln(Ar1, Ar2) ...
        && bicas.utils.is_subset_isequaln(Ar2, Ar1);
    end



    function utcStr = TT2000_to_UTC_str(zvTt2000, nSecondDecimals)
      bicas.utils.assert_ZV_Epoch(zvTt2000)
      utcStr = irf.cdf.TT2000_to_UTC_str(zvTt2000, nSecondDecimals);
    end



    % Log human readable summary of a set of zVar-like variables.
    %
    % NOTE: Ignores string ZVs.
    % NOTE: Can not handle FPAs.
    %
    %
    % ARGUMENTS
    % =========
    % Zvs
    %       Struct with ~zVariables.
    %       NOTE: Uses field name to determine whether field is Epoch-like
    %             or not.
    %
    function log_ZVs(Zvs, Bso, L)
      % PROBLEM: Can not manually specify which variables are Epoch-like.
      % PROBLEM: Can not manually specify variable name strings.
      %   Ex: process_HK_CDF_to_HK_on_SCI_TIME: Print different versions
      %       of time for comparison. Want whitespace
      %
      % PROPOSAL: For min-max values, also print difference.
      %   Ex: Time difference for Epoch.
      %       TODO-DEC: How print time difference?
      %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13:02:34
      %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13h02m34s

      LL = 'debug';

      fnList     = fieldnames(Zvs);
      ColumnStrs = irf.ds.empty_struct([0,1], ...
        'name', 'size', 'nNan', 'percentageNan', ...
        'nUniqueValues', 'values');

      for iFn = 1:numel(fnList)
        zvName  = fnList{iFn};
        zvValue = Zvs.(zvName);

        if iscolumn(zvValue) && isa(zvValue, 'int64') ...
            && any(irf.str.regexpf(...
            zvName, {'Epoch.*', '.*Epoch', '.*tt2000.*'}))
          % CASE: Epoch-like variable.

          ColumnStrs(end+1) = bicas.utils.get_array_statistics_strings(...
            zvName, zvValue, 'Epoch_LIKE_ZV', Bso);

        elseif isnumeric(zvValue)
          % CASE: Non-Epoch-like numeric variable.

          ColumnStrs(end+1) = bicas.utils.get_array_statistics_strings(...
            zvName, zvValue, 'NUMERIC_ZV', Bso);

        elseif ischar(zvValue)

          % Example of string valued (but irrelevant) CDF zVariables:
          % ACQUISITION_TIME_LABEL
          % Ignore. Do nothing.

        else
          error('BICAS:Assertion', ...
            'Can not handle zVar "%s" ("%s").', zvName, class(zvValue))
        end
      end

      HEADER_STRS = {'Name', 'Size', '#NaN', '%NaN', '#Uniq', 'Values'};
      dataStrs = {};
      dataStrs(:,1) = {ColumnStrs(:).name}';
      dataStrs(:,2) = {ColumnStrs(:).size}';
      dataStrs(:,3) = {ColumnStrs(:).nNan}';
      dataStrs(:,4) = {ColumnStrs(:).percentageNan}';
      dataStrs(:,5) = {ColumnStrs(:).nUniqueValues}';
      dataStrs(:,6) = {ColumnStrs(:).values}';
      columnAdjustments = [{'left', 'left'}, repmat({'right'}, 1,3), {'left'}];
      [HEADER_STRS, dataStrs, columnWidths] = ...
        irf.str.assist_print_table(...
        HEADER_STRS, dataStrs,  columnAdjustments);

      L.log(LL, strjoin(HEADER_STRS, ' '))
      L.log(LL, repmat('=', 1, sum(columnWidths) + numel(HEADER_STRS) - 1))
      for iRow = 1:numel(ColumnStrs)
        L.log(LL, strjoin(dataStrs(iRow, :), ' '))
      end
      L.logf(LL, [...
        '    #NaN = Number of NaN\n', ...
        '    #Uniq = Number of unique values incl.', ...
        ' NaN which counts as equal to itself.\n', ...
        '    Mm = min-max\n', ...
        '    Us = Unique values (explicitly listed)\n'])
    end



    %############
    % ASSERTIONS
    %############



    % Assert that variable is a "zVar Epoch-like" variable.
    function assert_ZV_Epoch(zvEpoch)
      % NOTE: No check for monotonically increasing timestamps. Done in
      % other locations. Universally? Slow?

      % PROPOSAL: Move to irf.
      %   PRO: Used by solo.sp.summary_plot.

      assert(iscolumn(zvEpoch), ...
        'BICAS:Assertion:IllegalArgument', ...
        'Argument is not a column vector')
      assert(isa(zvEpoch, 'int64'), ...
        'BICAS:Assertion:IllegalArgument', ...
        'Argument has the wrong class.')

      % Use?!!! Too processing heavy?!
      %validateattributes(Epoch, {'numeric'}, {'increasing'})
    end



    % Assert that variable is a "zVar ACQUISITION_TIME-like" variable.
    function assert_ZV_ACQUISITION_TIME(ACQUISITION_TIME)

      EMID = 'BICAS:Assertion:IllegalArgument';

      assert(isa(  ACQUISITION_TIME, 'uint32'), ...
        EMID, 'ACQUISITION_TIME is not uint32.')
      irf.assert.sizes(ACQUISITION_TIME, [NaN, 2])
      assert(all(  ACQUISITION_TIME(:, 1) >= 0), ...
        EMID, 'ACQUISITION_TIME has negative number of integer seconds.')
      % IMPLEMENTATION NOTE: Does not need to check for negative values
      % due to uint32.
      assert(all(  ACQUISITION_TIME(:, 2) < 65536), ...
        EMID, 'ACQUISITION_TIME subseconds out of range.')
    end



    % NOTE: Not implemented as assertion function in order to make it
    % possible to return proper error message on fail.
    function success = validate_ZV_QUALITY_FLAG(QUALITY_FLAG)
      if ~isa(QUALITY_FLAG, 'uint8')
        success = false;
      elseif ~all(...
          bicas.const.QUALITY_FLAG_MIN <= QUALITY_FLAG & ...
          bicas.const.QUALITY_FLAG_MAX >= QUALITY_FLAG ...
          )
        success = false;
      else
        success = true;
      end
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Whether Ar1 is a subset of Ar2.
    %
    % NOTE: Compare objects, not handles. NaN == NaN.
    % NOTE: Ignores duplicated objects within arrays.
    function isSubset = is_subset_isequaln(Ar1, Ar2)
      % PROPOSAL: Convert into generic function.

      for i = 1:numel(Ar1)
        if isempty(bicas.utils.find_first_isequaln(Ar1(i), Ar2))
          % CASE: Ar1(i) not found in Ar2
          isSubset = false;
          return
        end
      end

      isSubset = true;
    end



    % First index into Ar for which isequaln(x, Ar(i)).
    function i = find_first_isequaln(x, Ar)
      for i = 1:numel(Ar)
        if isequaln(x, Ar(i))
          return
        end
      end
      i = [];
    end



    function ColumnStrs = get_array_statistics_strings(...
        varName, varValue, varType, Bso)
      %
      % Derive statistics on the contents of a numeric variable (any
      % dimensionality) and return it so that it can easily be logged, e.g. in
      % a table:
      %   ** Array size
      %   ** Number of and percentage NaN,
      %   ** unique values, min-max.
      % Primarily intended for zVariables and derivatives thereof. Can be
      % useful for knowing which settings are used (e.g. DIFF_GAIN),
      % constant/varying bias current, suspect input datasets.
      %
      % Effectively a utility function for bicas.utils.log_zVar().
      %
      %
      % ARGUMENTS
      % =========
      % varName
      % varValue
      % varType
      %       String constant. 'NUMERIC_ZV' or 'Epoch_LIKE_ZV'.
      %       Determines how varValue is interpreted.
      %
      %
      % RETURN VALUE
      % ============
      % ColumnStrs
      %       Struct with fields corresponding to different column values for
      %       one row in a table.
      %

      % PROPOSAL: Test code.
      % PROPOSAL: Move to +utils.
      % PROPOSAL: Special log function for ZVs. Can print CDF type (implicitly range).
      % PROPOSAL: Print MATLAB class (implicitly range).
      % PROPOSAL: Better function name. Should imply that it generates strings for logging, not the logging
      %           itself.
      % PROPOSAL: Handle fill/pad value?
      % PROPOSAL: Only work on FPAs.

      % ASSERTION
      assert(isnumeric(varValue))

      uniqueValues  = bicas.utils.unique_values_NaN(varValue);
      nUniqueValues = numel(uniqueValues);
      nValues       = numel(varValue);

      %=================================
      % Construct string: variable size
      %=================================
      % Create comma-separated list of numbers.
      sizeStr = strjoin(arrayfun(...
        @(n) num2str(n), size(varValue), 'UniformOutput', 0),',');
      sizeStr = sprintf('(%s)', sizeStr);

      switch(varType)

        case 'NUMERIC_ZV'
          % ASSERTION
          assert(ndims(varValue) <= 3, ...
            'BICAS:Assertion:IllegalArgument', ...
            'v is not numerical with max 3 dimensions.')

          nNan             = sum(isnan(varValue(:)));
          nNanStr          = num2str(nNan);
          if nValues == 0
            percentageNanStr = '-';
          else
            percentageNan    = round((nNan/nValues)*100);
            percentageNanStr = sprintf('%i%%', percentageNan);
          end

          %===================================
          % Construct string: range of values
          %===================================
          if nUniqueValues > Bso.get_fv('LOGGING.MAX_NUMERIC_UNIQUES_PRINTED')
            vMin = min(min(min(varValue)));
            vMax = max(max(max(varValue)));

            % IMPLEMENTATION NOTE: Space around "--" to make it
            % easier to spot minus sign in a negative max number.
            valuesStr = sprintf('Mm: %d -- %d', vMin, vMax);
          else
            if nUniqueValues == 0
              valuesStr = '';
            else
              valuesStr = ['Us: ', sprintf('%d ', uniqueValues)];
            end
          end

        case 'Epoch_LIKE_ZV'
          % ASSERTIONS
          bicas.utils.assert_ZV_Epoch(varValue)

          nNanStr          = '-';
          percentageNanStr = '- ';   % NOTE: Extra whitespace.

          if nUniqueValues > Bso.get_fv('LOGGING.MAX_TT2000_UNIQUES_PRINTED')
            epochMinStr = bicas.utils.TT2000_to_UTC_str(min(varValue), 9);
            epochMaxStr = bicas.utils.TT2000_to_UTC_str(max(varValue), 9);
            valuesStr   = sprintf('Mm: %s -- %s', epochMinStr, epochMaxStr);
          elseif nValues >= 1
            bicas.utils.assert_ZV_Epoch(uniqueValues)
            valueStrs = irf.cdf.TT2000_to_UTC_str_many(uniqueValues, 9);
            valuesStr = ['Us: ', strjoin(valueStrs, ', ')];
          else
            valuesStr = '-';
          end

        otherwise
          error('BICAS:Assertion', ...
            'Illegal argument varType="%s"', varType)
      end

      % Assemble the final strings.
      ColumnStrs.name             = varName;
      ColumnStrs.size             = sizeStr;
      ColumnStrs.nNan             = nNanStr;
      ColumnStrs.percentageNan    = percentageNanStr;
      ColumnStrs.nUniqueValues    = num2str(nUniqueValues);
      ColumnStrs.values           = valuesStr;
    end



  end    % methods(Static, Access=private)



end
