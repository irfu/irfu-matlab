%
% Collection of minor utility functions (in the form of static methods) used for
% data processing.
%
% proc_utils = processing utilities
%
%
% TERMINOLOGY
% ===========
% SPR = Samples Per (CDF-like) Record
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-10-10
%
classdef proc_utils
%============================================================================================================
% PROPOSAL: Split up in separate files?!
% PROPOSAL: Move some functions to "utils".
%   Ex: log_array, log_struct_array
%
% PROPOSAL: Write test code for ACQUISITION_TIME_to_TT2000 and its inversion.
%
% N-->1 sample/record
%    NOTE: Time conversion may require moving the zero-point within the snapshot/record.
%    PROPOSAL: : All have nSamplesPerOldRecord as column vector.
%       PRO: LFR.
%    PROPOSAL: First convert column data to 2D data (with separate functions), then reshape to 1D with one common function.
%       CON: Does not work for ACQUISITION_TIME since two columns.
%
% PROPOSAL: Split up SPR functions.
%   convert_N_to_1_SPR_redistribute     -- Keep
%   convert_N_to_1_SPR_Epoch            -- increment_in_record + convert_N_to_1_SPR_redistribute
%   convert_N_to_1_SPR_ACQUISITION_TIME -- Keep
%   convert_1_to_1_SPR_by_repeating     -- convert_1_to_N_SPR_by_repeating + convert_N_to_1_SPR_redistribute



    methods(Static, Access=public)

        
        
        function c2 = select_row_range_from_cell_comps(c1, iFirst, iLast)
        % For every cell in a cell array, select an index range in the first
        % dimension for every cell array component.
            
            % ASSERTIONS
            bicas.proc_utils.assert_cell_array_comps_have_same_N_rows(c1)
            
            for i = 1:numel(c1)
                c2{i} = c1{i}(iFirst:iLast, :, :,:,:,:);
            end
        end



        function S = set_struct_field_rows(S, SNew, iRowsArray)
        % Generic utility function.
        % Overwrite struct fields at specific field rows using other struct
        % fields.
        %
        % ARGUMENTS
        % =========
        % S          : Struct. Only numeric fields.
        %              All fields have same number of rows.
        % SNew       : Struct. Only numeric fields. 
        %              All fields have same number of rows. Same fields as S.
        % iRowsArray : 1D array. Same length as number of rows in SNew fields.

            % ASSERTIONS
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S);
            nRowsSa = bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(SNew);
            assert(numel(iRowsArray) == nRowsSa)
            EJ_library.assert.castring_sets_equal(fieldnames(S), fieldnames(SNew))
            
            
            
            fieldNamesList = fieldnames(SNew);
            for i=1:length(fieldNamesList)
                fn = fieldNamesList{i};
                assert(isnumeric(S.(fn)))
                assert(isnumeric(SNew.(fn)))
                
                S.(fn)(iRowsArray, :) = SNew.(fn)(:, :);
            end
        end
        


        function tt2000 = ACQUISITION_TIME_to_TT2000(ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC)
        % Convert time in from ACQUISITION_TIME to tt2000 which is used for
        % Epoch in CDF files.
        % 
        % 
        % ARGUMENTS
        % =========
        % ACQUSITION_TIME            : NOTE: Can be negative since it is uint32.
        % ACQUISITION_TIME_EPOCH_UTC :
        %               Numeric row vector. The time in UTC at
        %               which ACQUISITION_TIME == [0,0], expressed as
        %               [year, month, day, hour, minute, second, 
        %                millisecond, microsecond(0-999), nanoseconds(0-999)].
        %
        %
        % RETURN VALUE
        % ============
        % tt2000 : NOTE: int64
        %
        
            bicas.proc_utils.assert_ACQUISITION_TIME(ACQUISITION_TIME)
            
            % at = ACQUISITION_TIME
            ACQUISITION_TIME = double(ACQUISITION_TIME);
            atSeconds = ACQUISITION_TIME(:, 1) + ACQUISITION_TIME(:, 2) / 65536;
            % NOTE: spdfcomputett2000 returns int64 (as it should).
            tt2000 = spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) + int64(atSeconds * 1e9);
        end
        

        
        function ACQUISITION_TIME = TT2000_to_ACQUISITION_TIME(tt2000, ACQUISITION_TIME_EPOCH_UTC)
        % Convert from tt2000 to ACQUISITION_TIME.
        %
        % t_tt2000 : Nx1 vector. Tequired to be int64 like the real zVar Epoch.
        % ACQUISITION_TIME : Nx2 vector. uint32.
        %       NOTE: ACQUSITION_TIME can not be negative since it is uint32.
        %
        
            % ASSERTIONS
            bicas.proc_utils.assert_zv_Epoch(tt2000)

            % NOTE: Important to type cast to double because of multiplication
            atSeconds = double(int64(tt2000) - spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC)) * 1e-9;    % at = ACQUISITION_TIME
            
            % ASSERTION: ACQUISITION_TIME must not be negative.
            if any(atSeconds < 0)
                error('BICAS:proc_utils:Assertion:IllegalArgument:DatasetFormat', 'Can not produce ACQUISITION_TIME (uint32) with negative number of integer seconds.')
            end
            
            atSeconds = round(atSeconds*65536) / 65536;
            atSecondsFloor = floor(atSeconds);
            
            ACQUISITION_TIME = uint32([]);
            ACQUISITION_TIME(:, 1) = uint32(atSecondsFloor);
            ACQUISITION_TIME(:, 2) = uint32((atSeconds - atSecondsFloor) * 65536);
            % NOTE: Should not be able to produce ACQUISITION_TIME(:, 2)==65536
            % (2^16) since atSeconds already rounded (to parts of 2^-16).
        end



        function filteredData = filter_rows(data, bRowFilter)
        % Function intended for filtering out data from a zVariable by setting
        % parts of it to NaN. Also useful for constructing aonymous functions.
        %
        %
        % ARGUMENTS
        % =========
        % data         : Numeric array with N rows.                 
        % bRowFilter   : Numeric/logical column vector with N rows.
        %
        %
        % RETURN VALUE
        % ============
        % filteredData :
        %         Array of the same size as "data", such that
        %         filteredData(i,:,:) == NaN,         for rowFilter(i)==0.
        %         filteredData(i,:,:) == data(i,:,:), for rowFilter(i)~=0.

            % ASSERTIONS
            assert(islogical(bRowFilter))    % Mostly to make sure the caller knows that it represents true/false.
            assert(isfloat(data), ...
                'BICAS:proc_utils:Assertion:IllegalArgument', ...
                'Argument "data" is not a floating-point class (can not represent NaN).')
            % Not really necessary to require row vector, only 1D vector.
            EJ_library.assert.sizes(...
                data,       [-1, NaN, NaN], ...
                bRowFilter, [-1])



            % Copy all data
            filteredData = data;
            
            % Overwrite data that should not have been copied with NaN
            % --------------------------------------------------------
            % IMPLEMENTATION NOTE: Command works empirically for filteredData
            % having any number of dimensions. However, if rowFilter and
            % filteredData have different numbers of rows, then the final array
            % may get the wrong dimensions (without triggering error!) since new
            % array components (indices) are assigned. ==> Having a
            % corresponding ASSERTION is important!
            filteredData(bRowFilter, :) = NaN;
        end

        
        
        function ACQUISITION_TIME_2 = convert_N_to_1_SPR_ACQUISITION_TIME(...
            ACQUISITION_TIME_1, nSpr, freqWithinRecords, ACQUISITION_TIME_EPOCH_UTC)
        % Function intended for converting ACQUISITION_TIME (always one time per
        % record) from many samples/record to one sample/record. See
        % convert_N_to_1_SPR_Epoch which is analogous.
        % 
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % ACQUISITION_TIME_1         : Nx2 vector.
        % freqWithinRecords          : Nx2 vector.
        % ACQUISITION_TIME_2         : Nx2 vector.
        % ACQUISITION_TIME_EPOCH_UTC : UTC as 1x9 row vector.
        %
        % NOTE: Theoretically, the function should be independent of the exact
        % value of ACQUISITION_TIME_EPOCH_UTC.

        % Command-line algorithm "test code":
        % clear; t_rec = [1;2;3;4]; f = [5;1;5;20]; N=length(t_rec); M=5; I_sample=repmat(0:(M-1), [N, 1]); F=repmat(f, [1,M]); T_rec = repmat(t_rec, [1,M]); T = T_rec + I_sample./F; reshape(T', [numel(T), 1])
            
            % ASSERTIONS
            bicas.proc_utils.assert_ACQUISITION_TIME(ACQUISITION_TIME_1)

            tt2000_1           = bicas.proc_utils.ACQUISITION_TIME_to_TT2000(ACQUISITION_TIME_1, ACQUISITION_TIME_EPOCH_UTC);
            tt2000_2           = EJ_library.so.convert_N_to_1_SPR_Epoch(     tt2000_1,           nSpr, freqWithinRecords);
            ACQUISITION_TIME_2 = bicas.proc_utils.TT2000_to_ACQUISITION_TIME(tt2000_2,           ACQUISITION_TIME_EPOCH_UTC);
        end
        
        
        
        function zv = set_NaN_after_snapshots_end(zv, snapshotLengths)
            % ASSERTIONS
            [nRecords, snapshotMaxLength] = EJ_library.assert.sizes(...
                zv,              [-1, -2], ...
                snapshotLengths, [-1]);
            assert(snapshotMaxLength >= max([snapshotLengths; 0]))
            % Add zero to vector so that max gives sensible value for empty snapshotLengths.
                        
            % IMPLEMENTATION
            for iRecord = 1:nRecords
                zv(iRecord, (snapshotLengths(iRecord)+1):end) = NaN;
            end
        end


        
        % Convert 2D array --> 1D cell array of 1D arrays, one per source row.
        %
        % ARGUMENTS
        % =========
        % M                     : 2D matrix
        % nCopyColsPerRowVec    : 1D column vector.
        %                         {i}=Number of elements to copy from M{i,:}.
        %
        % RETURN VALUE
        % ============
        % ca                    : Column cell array of 1D vectors.
        function ca = convert_matrix_to_cell_array_of_vectors(M, nCopyColsPerRowArray)
            EJ_library.assert.vector(nCopyColsPerRowArray)
            nRows = EJ_library.assert.sizes(M, [-1, NaN], nCopyColsPerRowArray, [-1, 1]);
            
            ca = cell(size(M, 1), 1);
            for iRow = 1:nRows
                ca{iRow} = M(iRow, 1:nCopyColsPerRowArray(iRow));
            end
        end
        

        
        % ARGUMENTS
        % =========
        % ca                    : Column cell array of 1D vectors.
        % nMatrixColumns        : Scalar. Number of columns in M.
        % M                     : Numeric 2D matrix.
        %                         NOTE: Sets unset elements to NaN.
        % nCopyColsPerRowVec    : 1D vector. {i}=Length of ca{i}=Number of elements copyied to M{i,:}.
        function [M, nCopyColsPerRowVec] = convert_cell_array_of_vectors_to_matrix(ca, nMatrixColumns)
            assert(iscell(ca))
            EJ_library.assert.vector(ca)
            assert(isscalar(nMatrixColumns))
            EJ_library.assert.vector(nMatrixColumns)
            
            nCopyColsPerRowVec = zeros(numel(ca), 1);   % Always column vector.
            M                  = nan(  numel(ca), nMatrixColumns);
            for iRow = 1:numel(nCopyColsPerRowVec)
                nCopyColsPerRowVec(iRow)            = numel(ca{iRow});
                M(iRow, 1:nCopyColsPerRowVec(iRow)) = ca{iRow};
            end
            
        end

        
        
        function zv_DELTA_PLUS_MINUS = derive_DELTA_PLUS_MINUS(freqHz, nSpr)
        %
        % Derive value for zVar DELTA_PLUS_MINUS.
        %
        % NOTE: All values on any given row (CDF record) of DELTA_PLUS_MINUS are
        % identical. Not sure why multiple values per row are needed but it is
        % probably intentional, as per YK's instruction.
        %
        %
        % ARGUMENTS
        % =========
        % freqHz : Frequency column vector in s^-1.
        %          Can not handle freqHz=NaN since the output is an
        %          integer (assertion).
        % nSpr   : Number of samples/record.
        %
        % 
        % RETURN VALUE
        % ============
        % DELTA_PLUS_MINUS : Analogous to BIAS zVariable. CDF_INT8=int64.
        %                    NOTE: Unit ns.
        %
        
            ZV_DELTA_PLUS_MINUS_DATA_TYPE = 'CDF_INT8';
            
            % ASSERTIONS
            nRecords = EJ_library.assert.sizes(freqHz, [-1]);
            assert(isfloat(freqHz) && all(isfinite(freqHz)), ...
                'BICAS:proc_utils:Assertion:IllegalArgument', ...
                'Argument "freqHz" does not consist of non-NaN floats.')
            assert(isscalar(nSpr), ...
                'BICAS:proc_utils:Assertion:IllegalArgument', ...
                'Argument "nSpr" is not a scalar.')
            
            
            
            zv_DELTA_PLUS_MINUS = zeros([nRecords, nSpr]);
            %DELTA_PLUS_MINUS = zeros([nRecords, 1]);    % Always 1 sample/record.
            for i = 1:nRecords
                % NOTE: Converts [s] (1/freqHz) --> [ns] (DELTA_PLUS_MINUS) so
                % that the unit is the same as for Epoch.
                % NOTE: Seems to work for more than 2D.
                zv_DELTA_PLUS_MINUS(i, :) = 1./freqHz(i) * 1e9 * 0.5;    % Unit: nanoseconds
            end
            zv_DELTA_PLUS_MINUS = cast(zv_DELTA_PLUS_MINUS, ...
                EJ_library.cdf.convert_CDF_type_to_MATLAB_class(...
                    ZV_DELTA_PLUS_MINUS_DATA_TYPE, 'Only CDF data types'));
        end



        function utcStr = TT2000_to_UTC_str(zvTt2000)
        % Convert tt2000 value to UTC string with nanoseconds.
            
            bicas.proc_utils.assert_zv_Epoch(zvTt2000)
            
            utcStr = EJ_library.cdf.TT2000_to_UTC_str(zvTt2000);
        end
        
        
        
        function ColumnStrs = log_array(varName, varValue, varType, SETTINGS)
            % Logs statistics on the contents of a numeric variable (any
            % dimensionality):
            %   ** Array size
            %   ** Number of and percentage NaN,
            %   ** unique values, min-max.
            % Primarily intended for zVariables and derivatives thereof. Can be
            % useful for knowing which settings are used (e.g. DIFF_GAIN),
            % constant/varying bias current, suspect input datasets.
            %
            % IMPLEMENTATION NOTE: Deliberately short function name to not
            % clutter the log.
            %
            %
            % ARGUMENTS
            % =========
            % varName  :
            % varValue :
            % varType  : String constant. 'numeric' or 'Epoch'.
            %            Determines how varValue is interpreted.
            %
            %
            % RETURN VALUE
            % ============
            % ColumnStrs : Struct with fields corresponding to different column
            %              values for one row in a table.
            %
            
            % PROPOSAL: Handle fill/pad value?
            % PROPOSAL: Move to +utils.
            % PROPOSAL: Special log function for zVars. Can print CDF type (implicitly range).
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
            sizeStr = strjoin(arrayfun(@(n) num2str(n), size(varValue), 'UniformOutput', 0),',');
            sizeStr = sprintf('(%s)', sizeStr);
            
            switch(varType)
                
                case 'numeric'
                    % ASSERTION
                    assert(ndims(varValue) <= 3, ...
                        'BICAS:proc_utils:Assertion:IllegalArgument', ...
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
                    bicas.proc_utils.assert_zv_Epoch(varValue)
                    
                    nNanStr          = '-';
                    percentageNanStr = '- ';   % NOTE: Extra whitespace.
                    
                    if nUniqueValues > SETTINGS.get_fv('LOGGING.MAX_TT2000_UNIQUES_PRINTED')
                        epochMinStr = bicas.proc_utils.TT2000_to_UTC_str(min(varValue));
                        epochMaxStr = bicas.proc_utils.TT2000_to_UTC_str(max(varValue));
                        valuesStr   = sprintf('Mm: %s -- %s', epochMinStr, epochMaxStr);
                    elseif nValues >= 1
                        bicas.proc_utils.assert_zv_Epoch(uniqueValues)
                        valueStrs = EJ_library.cdf.TT2000_to_UTC_str_many(uniqueValues);
                        valuesStr = ['Us: ', strjoin(valueStrs, ', ')];
                    else
                        valuesStr = '-';
                    end
                    
                otherwise
                    error('BICAS:proc_utils', 'Illegal argument varType="%s"', varType)
            end
            
            % Assemble the final strings.
            ColumnStrs.name             = varName;
            ColumnStrs.size             = sizeStr;
            ColumnStrs.nNan             = nNanStr;
            ColumnStrs.percentageNan    = percentageNanStr;
            ColumnStrs.nUniqueValues    = num2str(nUniqueValues);
            ColumnStrs.values           = valuesStr;
        end



        % Log human readable summary of a set of zVar-like variables.
        % NOTE: Ignores string zVars.
        % 
        % ARGUMENTS
        % =========
        % Zvs : Struct with ~zVariables.
        %       NOTE: Uses field name to determine whether field is Epoch-like
        %             or not.
        %
        function log_zVars(Zvs, SETTINGS, L)
            % PROBLEM: Can not manually specify which variables are Epoch-like.
            % PROBLEM: Can not manually specify variable name strings.
            %   Ex: process_HK_CDF_to_HK_on_SCI_TIME: Print different versions of time for comparison. Want whitespace
            %
            % PROPOSAL: For min-max values, also print difference.
            %   Ex: Time difference for Epoch.
            %       TODO-DECISION: How print time difference?
            %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13:02:34
            %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13h02m34s
            
            LL = 'debug';

            fnList     = fieldnames(Zvs);
            ColumnStrs = EJ_library.utils.empty_struct([0,1], ...
                'name', 'size', 'nNan', 'percentageNan', 'nUniqueValues', 'values');
            
            for iFn = 1:numel(fnList)
                zvName  = fnList{iFn};
                zvValue = Zvs.(zvName);
                
                if iscolumn(zvValue) && isa(zvValue, 'int64') ...
                        && any(EJ_library.str.regexpf(zvName, {'Epoch.*', '.*Epoch', '.*tt2000.*'}))
                    % CASE: Epoch-like variable.
                    
                    ColumnStrs(end+1) = bicas.proc_utils.log_array(zvName, zvValue, 'Epoch', SETTINGS);
                    
                elseif isnumeric(zvValue)
                    % CASE: Non-Epoch-like numeric variable.
                    
                    ColumnStrs(end+1) = bicas.proc_utils.log_array(zvName, zvValue, 'numeric', SETTINGS);
                    
                elseif ischar(zvValue)
                    
                    % Example of string valued (but irrelevant) CDF zVariables: ACQUISITION_TIME_LABEL
                    % Ignore. Do nothing.
                    
                else
                    error('BICAS:proc_utils:Assertion', 'Can not handle zVar "%s".', zvName)
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
            [HEADER_STRS, dataStrs, columnWidths] = EJ_library.str.assist_print_table(...
                HEADER_STRS, dataStrs,  columnAdjustments);

            L.log(LL, strjoin(HEADER_STRS, ' '))
            L.log(LL, repmat('=', 1, sum(columnWidths) + numel(HEADER_STRS) - 1))
            for iRow = 1:numel(ColumnStrs)
                L.log(LL, strjoin(dataStrs(iRow, :), ' '))
            end
            L.logf(LL, [...
                '    #NaN = Number of NaN\n', ...
                '    #Uniq = Number of unique values incl. NaN which counts as equal to itself.\n', ...
                '    Mm = min-max\n', ...
                '    Us = Unique values (explicitly listed)\n'])
        end



        % Assert that variable is an "zVar Epoch-like" variable.
        function assert_zv_Epoch(zvEpoch)

            assert(iscolumn(zvEpoch),     'BICAS:proc_utils:Assertion:IllegalArgument', 'Argument is not a column vector')
            assert(isa(zvEpoch, 'int64'), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Argument has the wrong class.')

            % Use?!!! Too processing heavy?!
            %validateattributes(Epoch, {'numeric'}, {'increasing'})
        end

        
        
        function assert_ACQUISITION_TIME(ACQUISITION_TIME)
        % Assert that variable is an "zVar ACQUISITION_TIME-like" variable.
        
            EMID = 'BICAS:proc_utils:Assertion:IllegalArgument';
        
            assert(isa(  ACQUISITION_TIME, 'uint32'),     EMID, 'ACQUISITION_TIME is not uint32.')
            EJ_library.assert.sizes(ACQUISITION_TIME, [NaN, 2])
            assert(all(  ACQUISITION_TIME(:, 1) >= 0),    EMID, 'ACQUISITION_TIME has negative number of integer seconds.')
            % IMPLEMENTATION NOTE: Does not need to check for negative values due to uint32.
            assert(all(  ACQUISITION_TIME(:, 2) < 65536), EMID, 'ACQUISITION_TIME subseconds out of range.')
        end
        
        
        
        function nRows = assert_struct_num_fields_have_same_N_rows(S)
        % Assert that data structure have the same number of rows in its
        % constituent parts.
        %
        % Useful for structs where all fields represent CDF zVariables and/or
        % derivatives thereof, the size in the first index (number of CDF
        % record) should be equal.
        %
        %
        % ARGUMENTS
        % =========
        % S : Struct
        %       Fields may (ony) be of the following types. Number of rows must
        %       be identical for the specified data structure components
        %       (right-hand side).
        %       (1) numeric/logical fields  : The field itself.
        %       (2) cell fields             : Cell array components (not
        %                                     the cell array itself!).
        %       (3) struct field (in struct): The inner struct's fields (not
        %                                     recursive).

        % NOTE: Function name somewhat bad.
        % PROPOSAL: Make recursive?!
        % PROPOSAL: Implement using new features in EJ_library.assert.sizes.
        
            fieldNamesList1 = fieldnames(S);
            nRowsArray = [];
            for iFn1 = 1:length(fieldNamesList1)
                fieldValue = S.(fieldNamesList1{iFn1});
                
                if isnumeric(fieldValue) || islogical(fieldValue)
                    
                    nRowsArray(end+1) = size(fieldValue, 1);
                    
                elseif iscell(fieldValue)
                    
                    for iCc = 1:numel(fieldValue)
                        nRowsArray(end+1) = size(fieldValue{iCc}, 1);
                    end
                    
                elseif isstruct(fieldValue)
                    
                    fieldNamesList2 = fieldnames(fieldValue);
                    for iFn2 = 1:length(fieldNamesList2)
                        nRowsArray(end+1) = size(fieldValue.(fieldNamesList2{iFn2}), 1);
                    end
                    
                else
                    
                    error('BICAS:proc_utils:Assertion', ...
                        'Can not handle this type of struct field.')
                    
                end
            end
            
            % NOTE: Empty vector if nRowsArray is empty.
            nRows = unique(nRowsArray);   
            
            % NOTE: length==0 valid for struct containing zero numeric fields.
            if length(unique(nRowsArray)) > 1    
                error('BICAS:proc_utils:Assertion', ...
                    ['Numeric fields and cell array components in struct do not have the same number', ...
                    ' of rows (likely corresponding to CDF zVar records).'])
            end
        end
        
        
        
        % Assert that cell array components all have the same number of rows.
        %
        % c : Cell array.
        function assert_cell_array_comps_have_same_N_rows(c)
            nRowsArray = cellfun(@(v) (size(v,1)), c, 'UniformOutput', true);
            EJ_library.assert.all_equal( nRowsArray )
        end



    end   % Static
    
    
    
end
