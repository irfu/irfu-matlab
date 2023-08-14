%
% Collection of shorter ~generic utility functions that could conceivably be
% used all over BICAS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-27, with moved from bicas.proc.utils.
%
classdef utils
    % PROPOSAL: Automatic test code.



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)


        % Get path to the root of the BICAS directory structure.
        function bicasRootPath = get_BICAS_root_path()
            % ASSUMES: The current file is in the <BICAS>/src/+bicas/ directory.
            % Use path of the current MATLAB file.
            [matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));
            bicasRootPath         = irf.fs.get_abs_path(...
                fullfile(matlabSrcPath, '..', '..'));
        end



        function utcStr = TT2000_to_UTC_str(zvTt2000)
        % Convert tt2000 value to UTC string with nanoseconds.

            bicas.utils.assert_ZV_Epoch(zvTt2000)

            utcStr = irf.cdf.TT2000_to_UTC_str(zvTt2000);
        end



        % Log human readable summary of a set of zVar-like variables.
        % NOTE: Ignores string ZVs.
        %
        %
        % ARGUMENTS
        % =========
        % Zvs : Struct with ~zVariables.
        %       NOTE: Uses field name to determine whether field is Epoch-like
        %             or not.
        %
        function log_ZVs(Zvs, SETTINGS, L)
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
                        zvName, zvValue, 'Epoch', SETTINGS);

                elseif isnumeric(zvValue)
                    % CASE: Non-Epoch-like numeric variable.

                    ColumnStrs(end+1) = bicas.utils.get_array_statistics_strings(...
                        zvName, zvValue, 'numeric', SETTINGS);

                elseif ischar(zvValue)

                    % Example of string valued (but irrelevant) CDF zVariables:
                    % ACQUISITION_TIME_LABEL
                    % Ignore. Do nothing.

                else
                    error('BICAS:Assertion', ...
                        'Can not handle zVar "%s".', zvName)
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



    end    % methods(Static)

    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        function ColumnStrs = get_array_statistics_strings(...
                varName, varValue, varType, SETTINGS)
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
        %       String constant. 'numeric' or 'Epoch'.
        %       Determines how varValue is interpreted.
        %
        %
        % RETURN VALUE
        % ============
        % ColumnStrs
        %       Struct with fields corresponding to different column values for
        %       one row in a table.
        %

            % PROPOSAL: Handle fill/pad value?
            % PROPOSAL: Move to +utils.
            % PROPOSAL: Special log function for ZVs. Can print CDF type (implicitly range).
            % PROPOSAL: Print MATLAB class (implicitly range).
            % PROPOSAL: Better function name. Should imply that it generates strings for logging, not the logging
            %           itself.

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

                case 'numeric'
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
                    if nUniqueValues > SETTINGS.get_fv('LOGGING.MAX_NUMERIC_UNIQUES_PRINTED')
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

                case 'Epoch'
                    % ASSERTIONS
                    bicas.utils.assert_ZV_Epoch(varValue)

                    nNanStr          = '-';
                    percentageNanStr = '- ';   % NOTE: Extra whitespace.

                    if nUniqueValues > SETTINGS.get_fv('LOGGING.MAX_TT2000_UNIQUES_PRINTED')
                        epochMinStr = bicas.utils.TT2000_to_UTC_str(min(varValue));
                        epochMaxStr = bicas.utils.TT2000_to_UTC_str(max(varValue));
                        valuesStr   = sprintf('Mm: %s -- %s', epochMinStr, epochMaxStr);
                    elseif nValues >= 1
                        bicas.utils.assert_ZV_Epoch(uniqueValues)
                        valueStrs = irf.cdf.TT2000_to_UTC_str_many(uniqueValues);
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
