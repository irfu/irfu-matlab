% Collection of minor utility functions (in the form of static methods) used for data processing.
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
%   Ex: log_array, log_struct_array, log_tt2000_array (uses bicas.proc_utils_assert_zv_Epoch)
%
% PROPOSAL: Write test code for ACQUISITION_TIME_to_tt2000 and its inversion.
%
% PROPOSAL: Replace find_last_same_subsequence with function that returns list of sequences.
%   PRO: Can naturally handle zero records.
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
        % For every cell in a cell array, select an index range in the first dimension for every cell array component.
            
            % ASSERTIONS
            bicas.proc_utils.assert_cell_array_comps_have_same_N_rows(c1)
            
            for i = 1:numel(c1)
                c2{i} = c1{i}(iFirst:iLast, :, :,:,:,:);
            end
        end



%         function S = add_rows_to_struct_fields(S, SAmendment)
%         % Generic utility function.
%         % Add values to every struct field by adding components after their highest row index (let them grow in the row
%         % index).
%         %
%         % NOTE 2020-07-29: Strong indications that using this function is inefficient if called for one added record at
%         % a time.
%         %
%         % NOTE: Keep function for a while for potential speed comparisons.
% 
%             bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S);
%             bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(SAmendment);
%             
%             fieldNamesList = fieldnames(SAmendment);
%             for i=1:length(fieldNamesList)
%                 fn = fieldNamesList{i};
%                 
%                 S.(fn) = [S.(fn) ; SAmendment.(fn)];
%             end
%         end

        
        
        function S = set_struct_field_rows(S, SAmendment, iRowsArray)
        % Generic utility function.
        % Set values in every struct field.

            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S);
            nRowsSa = bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(SAmendment);
            assert(numel(iRowsArray) == nRowsSa)
            EJ_library.assert.castring_sets_equal(fieldnames(S), fieldnames(SAmendment))
            
            fieldNamesList = fieldnames(SAmendment);
            for i=1:length(fieldNamesList)
                fn = fieldNamesList{i};
                assert(isnumeric(S.(fn)))
                
                S.(fn)(iRowsArray, :) = SAmendment.(fn)(:, :);
            end
        end
        


        function tt2000 = ACQUISITION_TIME_to_tt2000(ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC)
        % Convert time in from ACQUISITION_TIME to tt2000 which is used for Epoch in CDF files.
        % 
        % 
        % ARGUMENTS
        % =========
        % ACQUSITION_TIME            : NOTE: Can be negative since it is uint32.
        % ACQUISITION_TIME_EPOCH_UTC : Numeric row vector. The time in UTC at which ACQUISITION_TIME is [0,0] as
        %                              Year-month-day-hour-minute-second-millisecond-mikrosecond(0-999)-nanoseconds(0-999)
        %
        % RETURN VALUE
        % ============
        % tt2000 : NOTE: int64
        
            bicas.proc_utils.assert_ACQUISITION_TIME(ACQUISITION_TIME)
            
            ACQUISITION_TIME = double(ACQUISITION_TIME);
            atSeconds = ACQUISITION_TIME(:, 1) + ACQUISITION_TIME(:, 2) / 65536;    % at = ACQUISITION_TIME
            tt2000 = spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) + int64(atSeconds * 1e9);   % NOTE: spdfcomputett2000 returns int64 (as it should).
        end
        

        
        function ACQUISITION_TIME = tt2000_to_ACQUISITION_TIME(tt2000, ACQUISITION_TIME_EPOCH_UTC)
        % Convert from tt2000 to ACQUISITION_TIME.
        %
        % t_tt2000 : Nx1 vector. Tequired to be int64 like the real zVar Epoch.
        % ACQUISITION_TIME : Nx2 vector. uint32.
        %       NOTE: ACQUSITION_TIME can not be negative since it is uint32.
        
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
            % NOTE: Should not be able to produce ACQUISITION_TIME(:, 2)==65536 (2^16) since atSeconds already rounded (to parts of 2^-16).
        end
        
        
        
%         % Find sequences of constant value for a set of non-empty N-D vectors of identical length. Return sequences in
%         % the format of indices to "edges", here defined as the set union of
%         % (1) The first index
%         % (2) The last index+1
%         % (3) Every index which is the first index in a sequence of unchanging values for all vectors
%         % NOTE: NaN counts as equal to itself.
%         % NOTE: Needs to work for NaN in order to handle demultipexer mode and diff gain being NaN (unknown).
%         %
%         %
%         % ARGUMENTS AND RETURN VALUE
%         % ==========================
%         % varargin  : Matrices with same size in the first index. Must be at least one argument. Max 2-D.
%         % iEdgeList : 1D vector. Minimum-length 2.
%         % --
%         % RATIONALE: The function uses varargin (and the possibility to submit many separate 1D vectors) instead of one
%         % matrix argument (the caller merges) to make it possible to have vectors of multiple variable types (MATLAB
%         % classes) and different dimensions (column, row etc).
%         % --
%         % RATIONALE: The return format is chosen such that it is easy to merge it with other lists of edges from other
%         % sources.
%         %
%         %
%         % EXECUTION SPEED
%         % ===============
%         % Empirically, it can be useful to have a fast implementation of this function.
%         % IMPLEMENTATION NOTE: One can re-implement using subsref to handle higher-dimensional matrices, but this slows
%         % down the function, by factor of ~20. Presently hard-coded to permit up to 2-D matrices by the number indexing
%         % (hard-coded number of colons), but this can be increased to an arbitrary finite limit.
%         % Has kept multiple implementations to be able to compare speeds.
%         %
%         function iEdgeList = find_constant_sequences(varargin)
%             % PROPOSAL: Replace using EJ_library.utils.split_by_change.
%             
%             nArgs = numel(varargin);
%             
%             % ASSERTION
%             assert(nArgs >= 1, 'BICAS:proc_utils:Assertion:IllegalArgument', 'Must have at least one argument.')
% 
%             nRows  = size(varargin{1},1);
%             % Pre-allocate. Should be same size for all arguments and therefore does not need to be re-initialized/cleared.
%             diff_v = zeros(nRows-1, 1);
% 
% 
% 
%             if 1
%                 %=========================================================================
%                 % IMPLEMENTATION 1a
%                 % * One call to "isequalnan" for every row and argument.
%                 % NOTE: Slightly faster to index only once, and save the result (v1, v2).
%                 %=========================================================================
%                 
%                 iEdgeListArray = cell(nArgs, 1);    % Initialize empty variable.
%                 for iArg = 1:nArgs
%                     v = varargin{iArg};   % Argument before assertions.
%                     
%                     % ASSERTIONS
%                     assert(~isempty(v))
%                     assert(nRows == size(v, 1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Arguments have different number of rows.')
%                     assert(ndims(v) <= 2)
%                     
%                     v1 = v(1, :);
%                     for iRow = 1:nRows-1
%                         
%                         v2 = v(iRow+1, :);
%                         
%                         %==================================================================================================
%                         % Compare two slices of v, and subsequent in the first index of v
%                         % ----------------------------------------------------------------
%                         % IMPLEMENTATION NOTE: Uses "isequaln" to treat NaN as equal to itself. A side effect is that also
%                         % Inf equals itself, and that one can (untested) have arrays of non-numeric data (structs, chars,
%                         % objects etc).
%                         %==================================================================================================
%                         %diff_v(iRow) = ~isequaln(v(iRow,:), v(iRow+1,:));
%                         diff_v(iRow) = ~isequaln(v1, v2);
%                         
%                         v1 = v2;
%                     end
%                     
%                     iEdgeListArray{iArg} = [1; 1+find(diff_v); nRows+1];
%                 end
%                 iEdgeList = bicas.proc_utils.merge_index_edge_lists(iEdgeListArray{:});
%                 
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if 0
%                 %======================================================
%                 % IMPLEMENTATION 1b
%                 % * One call to "isequalnan" for every row and argument.
%                 % * Uses subsref ==> Makes it really slow.
%                 %======================================================
%                 
%                 iEdgeListArray = cell(nArgs, 1);    % Initialize empty variable.
%                 S = struct('type', '()', 'subs', {repmat({':'}, ndims(varargin{1}), 1)});
%                 
%                 for iArg = 1:nArgs
%                     v = varargin{iArg};   % Argument before assertions.
%                     
%                     % ASSERTIONS
%                     assert(~isempty(v))
%                     assert(nRows == size(v, 1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Arguments have different number of rows.')
%                     assert(ndims(v) <= 2)
%                     
%                     S.subs{1} = 1;
%                     v1        = subsref(v, S);
%                     for iRow = 1:nRows-1
%                         
%                         S.subs{1} = iRow + 1;
%                         v2        = subsref(v, S);
%                         
%                         %==================================================================================================
%                         % Compare two slices of v, and subsequent in the first index of v
%                         % ----------------------------------------------------------------
%                         % IMPLEMENTATION NOTE: Uses "isequaln" to treat NaN as equal to itself. A side effect is that also
%                         % Inf equals itself, and that one can (untested) have arrays of non-numeric data (structs, chars,
%                         % objects etc).
%                         %==================================================================================================
%                         diff_v(iRow) = ~isequaln(v1, v2);
%                         
%                         v1 = v2;
%                     end
%                     
%                     iEdgeListArray{iArg} = [1; 1+find(diff_v); nRows+1];
%                 end
%                 iEdgeList = bicas.proc_utils.merge_index_edge_lists(iEdgeListArray{:});
%                 
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if 0
%                 %========================================================
%                 % IMPLEMENTATION 2
%                 % * One call to "isequalnan" for every row (not argument).
%                 % Seems slower than not.
%                 %========================================================
%                 v1 = cell(nArgs, 1);
%                 v2 = cell(nArgs, 1);
%                 for iArg = 1:nArgs
%                     v = varargin{iArg};   % Argument before assertions.
%                     
%                     % ASSERTIONS
%                     assert(~isempty(v))
%                     assert(nRows == size(v, 1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Arguments have different number of rows.')
%                     assert(ndims(v) <= 2)
%                     
%                     v1{iArg} = v(1, :);
%                 end
%                 for iRow = 1:nRows-1
%                     for iArg = 1:nArgs
%                         v = varargin{iArg};
%                         v2{iArg} = v(iRow+1, :);
%                     end
%                     diff_v(iRow) = ~isequaln(v1, v2);
%                     v1 = v2;
%                 end
%                 iEdgeList = [1; 1+find(diff_v); nRows+1];
%                 
%             end
%         end
%         
%         
%         
%         function iEdgeList = merge_index_edge_lists(varargin)
%             iEdgeList = [];
%             for i = 1:numel(varargin)
%                 v = varargin{i};
%                 
%                 % ASSERTIONS
%                 EJ_library.assert.vector(v)
%                 % Verifies that it is an edge list, and that the "convention" for what is an edge list (include
%                 % beginning and end) has not changed.
%                 assert(v(1) == 1)    
%                 
%                 % NOTE: Works with (forces) column vectors to make concatenations reliable
%                 iEdgeList = [iEdgeList; varargin{i}(:)];
%             end
%             iEdgeList = sort(unique(iEdgeList));
%         end
%         
%         
%         
%         % EXPERIMENTAL
%         %
%         % Convert a list of edges (indexes) into adjacent sequences of indices, represented by lists of the first and
%         % last index for each sequence. Each sequence begins and ends with an edge.
%         %
%         % iEdgeList  : Sorted numeric 1D vector. If empty or scalar, then empty vectors are returned.
%         % iFirstList, iLastList : Vectors with first and last index for each sequence.
%         %
%         function [iFirstList, iLastList] = index_edges_2_first_last(iEdgeList)
%             assert(issorted(iEdgeList))
%             EJ_library.assert.vector(iEdgeList)
%             
%             iFirstList = iEdgeList(1:end-1);
%             iLastList  = iEdgeList(2:end) - 1;
%         end
        
        
        
        function filteredData = filter_rows(data, bRowFilter)
        % Function intended for filtering out data from a zVariable by setting parts of it to NaN. Also useful for
        % constructing aonymous functions.
        %
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % data         : Numeric array with N rows.                 (Intended to represent a zVariable with N records.)
        % bRowFilter   : Numeric/logical column vector with N rows. (Intended to represent a zVariable with N records.)
        % filteredData : Array of the same size as "data", such that
        %                filteredData(i,:,:, ...) == NaN,              for rowFilter(i)==0.
        %                filteredData(i,:,:, ...) == data(i,:,:, ...), for rowFilter(i)~=0.

            % ASSERTIONS
            assert(islogical(bRowFilter))    % Mostly to make sure the caller knows that it represents true/false.
            assert(isfloat(data), ...
                'BICAS:proc_utils:Assertion:IllegalArgument', 'Argument "data" is not a floating-point class (can not represent NaN).')
            % Not really necessary to require row vector, only 1D vector.
            assert(iscolumn(bRowFilter), ...
                'BICAS:proc_utils:Assertion:IllegalArgument', 'Argument "rowFilter" is not a column vector.')
            assert(size(bRowFilter, 1) == size(data, 1), ...
                'BICAS:proc_utils:Assertion:IllegalArgument', 'Numbers of records do not match.')



            % Copy all data
            filteredData = data;
            
            % Overwrite data that should not have been copied with NaN
            % --------------------------------------------------------
            % IMPLEMENTATION NOTE: Command works empirically for filteredData having any number of dimensions. However,
            % if rowFilter and filteredData have different numbers of rows, then the final array may get the wrong
            % dimensions (without triggering error!) since new array components (indices) are assigned. ==> Having a
            % corresponding ASSERTION is important!
            filteredData(bRowFilter, :) = NaN;
        end

        
        
        function ACQUISITION_TIME_2 = convert_N_to_1_SPR_ACQUISITION_TIME(...
            ACQUISITION_TIME_1, nSpr, freqWithinRecords, ACQUISITION_TIME_EPOCH_UTC)
        % Function intended for converting ACQUISITION_TIME (always one time per record) from many samples/record to one
        % sample/record. See convert_N_to_1_SPR_Epoch which is analogous.
        % 
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % ACQUISITION_TIME_1         : Nx2 vector.
        % freqWithinRecords          : Nx2 vector.
        % ACQUISITION_TIME_2         : Nx2 vector.
        % ACQUISITION_TIME_EPOCH_UTC : UTC as 1x9 row vector.
        %
        % NOTE: Theoretically, the function should be independent of the exact value of ACQUISITION_TIME_EPOCH_UTC.

        % Command-line algorithm "test code":
        % clear; t_rec = [1;2;3;4]; f = [5;1;5;20]; N=length(t_rec); M=5; I_sample=repmat(0:(M-1), [N, 1]); F=repmat(f, [1,M]); T_rec = repmat(t_rec, [1,M]); T = T_rec + I_sample./F; reshape(T', [numel(T), 1])
            
            % ASSERTIONS
            bicas.proc_utils.assert_ACQUISITION_TIME(ACQUISITION_TIME_1)

            tt2000_1           = bicas.proc_utils.ACQUISITION_TIME_to_tt2000(ACQUISITION_TIME_1, ACQUISITION_TIME_EPOCH_UTC);
            tt2000_2           = EJ_library.so.convert_N_to_1_SPR_Epoch(     tt2000_1,           nSpr, freqWithinRecords);
            ACQUISITION_TIME_2 = bicas.proc_utils.tt2000_to_ACQUISITION_TIME(tt2000_2,           ACQUISITION_TIME_EPOCH_UTC);
        end
        
        
        
        function zv = set_NaN_after_snapshots_end(zv, snapshotLengths)
            assert(iscolumn(snapshotLengths))
            
            assert(ndims(zv) == 2)   % NOTE: ndims always returns at least 2.
            assert(size(zv,1) == numel(snapshotLengths))
            assert(size(zv,2) >= max([snapshotLengths; 0]))   % Add zero so that max gives sensible value for empty snapshotLengths.
            
            for iRecord = 1:numel(snapshotLengths)
                zv(iRecord, (snapshotLengths(iRecord)+1):end) = NaN;
            end
        end


        
        % Convert 2D array --> 1D cell array of 1D arrays, one per source row.
        %
        % ARGUMENTS
        % =========
        % M                     : 2D matrix
        % nCopyColsPerRowVec    : 1D column vector. {i}=Number of elements to copy from M{i,:}.
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
            M                  = zeros(numel(ca), nMatrixColumns) * NaN;
            for iRow = 1:numel(nCopyColsPerRowVec)
                nCopyColsPerRowVec(iRow)            = numel(ca{iRow});
                M(iRow, 1:nCopyColsPerRowVec(iRow)) = ca{iRow};
            end
            
        end

        
        
        function DELTA_PLUS_MINUS = derive_DELTA_PLUS_MINUS(freqHz, nSpr)
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % freqHz           : Frequency column vector in s^-1. Can not handle freqHz=NaN since the output is an integer.
        % nSpr             : Number of samples/record.
        % DELTA_PLUS_MINUS : Analogous to BIAS zVariable. CDF_INT8=int64. NOTE: Unit ns.
        
            ZV_DELTA_PLUS_MINUS_DATA_TYPE = 'CDF_INT8';
            
            if ~iscolumn(freqHz) || ~isfloat(freqHz) || any(isnan(freqHz))
                error('BICAS:proc_utils:Assertion:IllegalArgument', '"freqHz" is not a column vector of non-NaN floats.')
            elseif ~isscalar(nSpr)
                error('BICAS:proc_utils:Assertion:IllegalArgument', '"nSpr" is not a scalar.')
            end
            
            nRecords = size(freqHz, 1);
            DELTA_PLUS_MINUS = zeros([nRecords, nSpr]);
            for i = 1:length(freqHz)
                DELTA_PLUS_MINUS(i, :) = 1./freqHz(i) * 1e9 * 0.5;      % Seems to work for more than 2D.
            end
            DELTA_PLUS_MINUS = cast(DELTA_PLUS_MINUS, EJ_library.cdf.convert_CDF_type_to_MATLAB_class(...
                ZV_DELTA_PLUS_MINUS_DATA_TYPE, 'Only CDF data types'));
        end
        
        
        
        % mSize : [nRows, nColumns, ...] so that the return value from the size() function can be used.
        function M = create_NaN_array(mSize)
            assert(numel(mSize) >= 2)
            
            M = NaN * zeros(mSize);
        end


        
        function utcStr = tt2000_to_UTC_str(zvTt2000)
        % Convert tt2000 value to UTC string with nanoseconds.
            
            bicas.proc_utils.assert_zv_Epoch(zvTt2000)
            
            utcStr = EJ_library.cdf.tt2000_to_UTC_str(zvTt2000);
        end
        
        
        
        function ColumnStrs = log_array(varName, varValue, varType, SETTINGS)
            % Logs statistics on the contents of a numeric variable (any dimensionality):
            %   ** Array size
            %   ** Number of and percentage NaN,
            %   ** unique values, min-max.
            % Primarily intended for zVariables and derivatives thereof. Can be useful for knowing which settings are
            % used (e.g. DIFF_GAIN), constant/varying bias current, suspect input datasets.
            %
            % IMPLEMENTATION NOTE: Deliberately short function name to not clutter the log.
            %
            %
            % ARGUMENTS
            % =========
            % varName  :
            % varValue :
            % varType  : String constant. 'numeric' or 'Epoch'. Determines how varValue is interpreted.
            %
            %
            % RETURN VALUE
            % ============
            % ColumnStrs : Struct with fields corresponding to different column values for one row in a table.
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
                    assert(ndims(varValue) <= 3, 'BICAS:proc_utils:Assertion:IllegalArgument', 'v is not numerical with max 3 dimensions.')

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
                        
                        % IMPLEMENTATION NOTE: Space around "--" to make it easier to spot minus sign in a negative max number.
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
                        epochMinStr = bicas.proc_utils.tt2000_to_UTC_str(min(varValue));
                        epochMaxStr = bicas.proc_utils.tt2000_to_UTC_str(max(varValue));
                        valuesStr   = sprintf('Mm: %s -- %s', epochMinStr, epochMaxStr);
                    elseif nValues >= 1
                        bicas.proc_utils.assert_zv_Epoch(uniqueValues)
                        valueStrs = EJ_library.cdf.tt2000_to_UTC_str_many(uniqueValues);
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
        %       NOTE: Uses field name to determine whether field is Epoch-like or not.
        %
        function log_zVars(Zvs, SETTINGS, L)
            % PROBLEM: Can not manually specify which variables are Epoch-like.
            % PROBLEM: Can not manually specify variable name strings.
            %   Ex: process_HK_to_HK_on_SCI_TIME: Print different versions of time for comparison. Want whitespace
            %
            % PROPOSAL: For min-max values, also print difference.
            %   Ex: Time difference for Epoch.
            %       TODO-DECISION: How print time difference?
            %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13:02:34
            %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13h02m34s
            
            LL = 'debug';

            fnList     = fieldnames(Zvs);
            ColumnStrs = EJ_library.utils.empty_struct([0,1], 'name', 'size', 'nNan', 'percentageNan', 'nUniqueValues', 'values');
            
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
            
            headerStrs = {'Name', 'Size', '#NaN', '%NaN', '#Uniq', 'Values'};
            dataStrs = {};
            dataStrs(:,1) = {ColumnStrs(:).name}';
            dataStrs(:,2) = {ColumnStrs(:).size}';
            dataStrs(:,3) = {ColumnStrs(:).nNan}';
            dataStrs(:,4) = {ColumnStrs(:).percentageNan}';
            dataStrs(:,5) = {ColumnStrs(:).nUniqueValues}';
            dataStrs(:,6) = {ColumnStrs(:).values}';
            columnAdjustments = [{'left', 'left'}, repmat({'right'}, 1,3), {'left'}];
            [headerStrs, dataStrs, columnWidths] = EJ_library.str.assist_print_table(...
                headerStrs, dataStrs,  columnAdjustments);

            L.log(LL, strjoin(headerStrs, ' '))
            L.log(LL, repmat('=', 1, sum(columnWidths) + numel(headerStrs) - 1))
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
        % Assert that data structure have the same number of rows in its constituent parts.
        %
        % Useful for structs where all fields represent CDF zVariables and/or derivatives thereof, the size in the first
        % index (number of CDF record) should be equal.
        %
        % ARGUMENTS
        % =========
        % S : Struct
        %       Fields may (ony) be of the following types. Number of rows must be identical for the specified
        %       data structure components (right-hand side).
        %       (1) numeric/logical fields          : The field itself
        %       (2) cell fields                     : Cell array components (not the cell array itself!)
        %       (3) struct field (inside the struct): The inner struct's fields (not recursive)

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
                    
                    error('BICAS:proc_utils:Assertion', 'Can not handle this type of struct field.')
                    
                end
            end
            
            nRows = unique(nRowsArray);   % NOTE: Empty vector if nRowsArray is empty.
            
            if length(unique(nRowsArray)) > 1    % NOTE: length==0 valid for struct containing zero numeric fields.
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
