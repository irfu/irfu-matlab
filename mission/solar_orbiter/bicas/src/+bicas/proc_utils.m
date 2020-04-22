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
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
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



        function S = add_rows_to_struct_fields(S, SAmendment)
        % Generic utility function.
        % Add values to every struct field by adding components after their highest row index (let them grow in the row
        % index).
 

            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S);
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(SAmendment);
            
            fieldNamesList = fieldnames(SAmendment);
            for i=1:length(fieldNamesList)
                fn = fieldNamesList{i};
                
                S.(fn) = [S.(fn) ; SAmendment.(fn)];
            end
        end

        

        function freqHz = get_LFR_frequency(iLsf, lsfArrayHz)
        % Convert LFR zVariable FREQ values to Hz. The usefulness of this function stems from how the LFR
        % datasets are defined.
        %
        % ARGUMENTS
        % =========
        % iLsf   : The LSF index, i.e. 1=LFR freq. F0, and so on.
        % freqHz : Frequency in Hz.
        
        % PROPOSAL: Somehow avoid having an argument for lsfArrayHz. Have it as a constant somehow.

            % ASSERTION
            uniqueValues = unique(iLsf);
            if ~all(ismember(uniqueValues, [1,2,3,4]))
                uniqueValuesStr = sprintf('%d', uniqueValues);   % NOTE: Has to print without \n to keep all values on a single-line string.
                error('BICAS:proc_utils:Assertion:IllegalArgument:DatasetFormat', ...
                    'Found unexpected values in LSF index (corresponding to LFR FREQ+1). Unique values: %s.', uniqueValuesStr)
            end
            
            % NOTE: Implementation that works for arrays of any size.
            freqHz = ones(size(iLsf)) * NaN;        % Allocate array and set default values.
            freqHz(iLsf==1) = lsfArrayHz(1);
            freqHz(iLsf==2) = lsfArrayHz(2);
            freqHz(iLsf==3) = lsfArrayHz(3);
            freqHz(iLsf==4) = lsfArrayHz(4);
        end
        
        
        
        function Rx = get_LFR_Rx(R0, R1, R2, iLsf)
        % Return the relevant value of LFR CDF zVariables R0, R1, or R2, or a hypothetical but analogous "R3" which is always 1.
        %
        % ARGUMENTS
        % =========
        % R0, R1, R2, FREQ : LFR CDF zVariables. All must have identical array sizes.
        %                    FREQ(i) == 0 (not 1) ==> F0 and so on.
        % Rx               : Same size array as R0, R1, R2, FREQ. The relevant values are copied, respectively, from
        %                    R0, R1, R2, or an analogous hypothetical "R3" that is a constant (=1) depending on
        %                    the value of FREQ in the corresponding component.
        %                    NOTE: Not MATLAB class "logical".
        %
        % NOTE: Works for all array sizes.
            
            Rx = NaN * ones(size(iLsf));        % Set to NaN (should always be overwritten if code works).
            
            I = (iLsf==1);   Rx(I) = R0(I);
            I = (iLsf==2);   Rx(I) = R1(I);
            I = (iLsf==3);   Rx(I) = R2(I);
            I = (iLsf==4);   Rx(I) = 1;      % The value of a hypothetical (non-existant, constant) analogous zVariable "R3".
            
            assert(all(~isnan(Rx)))
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
        
        
        
        % Find sequences of constant value for a set of non-empty N-D vectors of identical length. Return sequences in
        % the format of indices to "edges", here defined as the set union of
        % (1) The first index
        % (2) The last index+1
        % (3) Every index which is the first index in a sequence of unchanging values for all vectors
        % NOTE: NaN counts as equal to itself.
        % NOTE: Needs to work for NaN in order to handle demultipexer mode and diff gain being NaN (unknown).
        %
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % varargin  : Matrices with same size in the first index. Must be at least one argument. Max 2-D.
        % iEdgeList : 1D vector. Minimum-length 2.
        % --
        % RATIONALE: The function uses varargin (and the possibility to submit many separate 1D vectors) instead of one
        % matrix argument (the caller merges) to make it possible to have vectors of multiple variable types (MATLAB
        % classes) and different dimensions (column, row etc).
        % --
        % RATIONALE: The return format is chosen such that it is easy to merge it with other lists of edges from other
        % sources.
        %
        %
        % EXECUTION SPEED
        % ===============
        % Empirically, it can be useful to have a fast implementation of this function.
        % IMPLEMENTATION NOTE: One can re-implement using subsref to handle higher-dimensional matrices, but this slows
        % down the function, by factor of ~20. Presently hardcoded to permit up to 2-D matrices by the number indexing
        % (hardcoded number of colons), but this can be increased to an arbitrary finite limit.
        % Has kept multiple implementations to be able to compare speeds.
        %
        function iEdgeList = find_constant_sequences(varargin)
            % PROPOSAL: Rename to imply that the function finds edges (separating sequences), not sequences.
            %tic
            nArgs = numel(varargin);
            
            % ASSERTION
            assert(nArgs >= 1, 'BICAS:proc_utils:Assertion:IllegalArgument', 'Must have at least one argument.')
                        
            nRows  = size(varargin{1},1);
            % Pre-allocate. Should be same size for all arguments and therefore does not need to be re-initialized/cleared.
            diff_v = zeros(nRows-1, 1);



            if 1
                %=========================================================================
                % IMPLEMENTATION 1a
                % * One call to "isequalnan" for every row and argument.
                % NOTE: Slightly faster to index only once, and save the result (v1, v2).
                %=========================================================================
                
                iEdgeListArray = cell(nArgs, 1);    % Initialize empty variable.
                for iArg = 1:nArgs
                    v = varargin{iArg};   % Argument before assertions.
                    
                    % ASSERTIONS
                    assert(~isempty(v))
                    assert(nRows == size(v, 1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Arguments have different number of rows.')
                    assert(ndims(v) <= 2)
                    
                    v1 = v(1, :);
                    for iRow = 1:nRows-1
                        
                        v2 = v(iRow+1, :);
                        
                        %==================================================================================================
                        % Compare two slices of v, and subsequent in the first index of v
                        % ----------------------------------------------------------------
                        % IMPLEMENTATION NOTE: Uses "isequaln" to treat NaN as equal to itself. A side effect is that also
                        % Inf equals itself, and that one can (untested) have arrays of non-numeric data (structs, chars,
                        % objects etc).
                        %==================================================================================================
                        %diff_v(iRow) = ~isequaln(v(iRow,:), v(iRow+1,:));
                        diff_v(iRow) = ~isequaln(v1, v2);
                        
                        v1 = v2;
                    end
                    
                    iEdgeListArray{iArg} = [1; 1+find(diff_v); nRows+1];
                end
                iEdgeList = bicas.proc_utils.merge_index_edge_lists(iEdgeListArray{:});
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if 0
                %======================================================
                % IMPLEMENTATION 1b
                % * One call to "isequalnan" for every row and argument.
                % * Uses subsref ==> Makes it really slow.
                %======================================================
                
                iEdgeListArray = cell(nArgs, 1);    % Initialize empty variable.
                S = struct('type', '()', 'subs', {repmat({':'}, ndims(varargin{1}), 1)});
                
                for iArg = 1:nArgs
                    v = varargin{iArg};   % Argument before assertions.
                    
                    % ASSERTIONS
                    assert(~isempty(v))
                    assert(nRows == size(v, 1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Arguments have different number of rows.')
                    assert(ndims(v) <= 2)
                    
                    S.subs{1} = 1;
                    v1        = subsref(v, S);
                    for iRow = 1:nRows-1
                        
                        S.subs{1} = iRow + 1;
                        v2        = subsref(v, S);
                        
                        %==================================================================================================
                        % Compare two slices of v, and subsequent in the first index of v
                        % ----------------------------------------------------------------
                        % IMPLEMENTATION NOTE: Uses "isequaln" to treat NaN as equal to itself. A side effect is that also
                        % Inf equals itself, and that one can (untested) have arrays of non-numeric data (structs, chars,
                        % objects etc).
                        %==================================================================================================
                        diff_v(iRow) = ~isequaln(v1, v2);
                        
                        v1 = v2;
                    end
                    
                    iEdgeListArray{iArg} = [1; 1+find(diff_v); nRows+1];
                end
                iEdgeList = bicas.proc_utils.merge_index_edge_lists(iEdgeListArray{:});
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if 0
                %========================================================
                % IMPLEMENTATION 2
                % * One call to "isequalnan" for every row (not argument).
                % Seems slower than not.
                %========================================================
                v1 = cell(nArgs, 1);
                v2 = cell(nArgs, 1);
                for iArg = 1:nArgs
                    v = varargin{iArg};   % Argument before assertions.
                    
                    % ASSERTIONS
                    assert(~isempty(v))
                    assert(nRows == size(v, 1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'Arguments have different number of rows.')
                    assert(ndims(v) <= 2)
                    
                    v1{iArg} = v(1, :);
                end
                for iRow = 1:nRows-1
                    for iArg = 1:nArgs
                        v = varargin{iArg};
                        v2{iArg} = v(iRow+1, :);
                    end
                    diff_v(iRow) = ~isequaln(v1, v2);
                    v1 = v2;
                end
                iEdgeList = [1; 1+find(diff_v); nRows+1];
                
            end
        end
        
        
        
        % EXPERIMENTAL
        function iEdgeList = merge_index_edge_lists(varargin)
            iEdgeList = [];
            for i = 1:numel(varargin)
                v = varargin{i};
                
                % ASSERTIONS
                EJ_library.assert.vector(v)
                assert(v(1) == 1)    % Verifies that it is an edge list, and that the "convention" for what is an edge list (include beginning and end) has not changed.
                
                % NOTE: Works with (forces) column vectors to make concatenations reliable
                iEdgeList = [iEdgeList; varargin{i}(:)];
            end
            iEdgeList = sort(unique(iEdgeList));
        end
        
        
        
        % EXPERIMENTAL
        %
        % Convert a list of edges (indexes) into adjacent sequences of indices, represented by lists of the first and
        % last index for each sequence. Each sequence begins and ends with an edge.
        %
        % iEdgeList  : Sorted numeric 1D vector. If empty or scalar, then empty vectors are returned.
        % iFirstList, iLastList : Vectors with first and last index for each sequence.
        %
        function [iFirstList, iLastList] = index_edges_2_first_last(iEdgeList)
            assert(issorted(iEdgeList))
            EJ_library.assert.vector(iEdgeList)
            
            
            iFirstList = iEdgeList(1:end-1);
            iLastList  = iEdgeList(2:end) - 1;
        end
        
        
        
        function filteredData = filter_rows(data, rowFilter)
        % Function intended for filtering out data from a zVariable by setting parts of it to NaN.
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % data         : Numeric array with N rows.                 (Intended to represent a zVariable with N records.)
        % rowFilter    : Numeric/logical column vector with N rows. (Intended to represent a zVariable with N records.)
        % filteredData : Array of the same size as "data", such that
        %                filteredData(i,:,:, ...) == NaN,              for rowFilter(i)==0.
        %                filteredData(i,:,:, ...) == data(i,:,:, ...), for rowFilter(i)~=0.

            % ASSERTIONS
            if ~iscolumn(rowFilter)     % Not really necessary to require row vector, only 1D vector.
                error('BICAS:proc_utils:Assertion:IllegalArgument', '"rowFilter" is not a column vector.')
            elseif size(rowFilter, 1) ~= size(data, 1)
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Numbers of records do not match.')
            elseif ~isfloat(data)
                error('BICAS:proc_utils:Assertion:IllegalArgument', '"data" is not a floating-point class (can not represent NaN).')
            end
            
            
            
            % Copy all data
            filteredData = data;
            
            % Overwrite data that should not have been copied with NaN
            % --------------------------------------------------------
            % IMPLEMENTATION NOTE: Command works empirically for filteredData having any number of dimensions. However,
            % if rowFilter and filteredData have different numbers of rows, then the final array may get the wrong
            % dimensions (without triggering error!) since new array components (indices) are assigned. ==> Having a
            % corresponding ASSERTION is important!
            filteredData(rowFilter==0, :) = NaN;
        end



        function zv2 = convert_N_to_1_SPR_redistribute(zv1)
        % Convert zVariable-like variable from N samples/record to 1 sample/record (from a matrix to a column vector).
        % Increases number of records.
        %
        % ARGUMENT
        % ========
        % zv1     : (iRecord, iSnapshotSample, iChannel)
        %
        % RETURN VALUE
        % ============
        % zv2     : (iRecord, iChannel). Same number of components, nRecords2 = nRecords1*nSpr1

            EJ_library.assert.size(zv1, [NaN, NaN, NaN])
            
            zv  = permute(zv1, [2,1,3]);
            zv2 = reshape(zv, size(zv,1) * size(zv,2), size(zv,3));
            
            EJ_library.assert.size(zv2, [NaN, NaN])
        end



        function zv2 = convert_1_to_1_SPR_by_repeating(zv1, nRepeatsPerRecord1)
        % Two steps:
        % (1) Convert zVariable-like variable from 1 value/record to N values/record (same number of records) by
        %     repeating within record.
        % (2) Convert zVariable-like variable from N values/record to 1 value/record by redistributing values.
            
            % ASSERTIONS
            assert(iscolumn(zv1),                'BICAS:proc_utils:Assertion:IllegalArgument', 'zv1 is not a column vector')
            assert(isscalar(nRepeatsPerRecord1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'nRepeatsPerRecord1 is not a scalar')
            
            zv2 = bicas.proc_utils.convert_1_to_N_SPR_by_repeating(zv1, nRepeatsPerRecord1);
            zv2 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(zv2);
        end
        
        
        
        function zv2 = convert_1_to_N_SPR_by_repeating(zv1, nRepeatsPerRecord1)
            % NOTE: Maybe somewhat unnecessary function.
            
            assert(iscolumn(zv1),                'BICAS:proc_utils:Assertion:IllegalArgument', 'zv1 is not a column vector')
            assert(isscalar(nRepeatsPerRecord1), 'BICAS:proc_utils:Assertion:IllegalArgument', 'nRepeatsPerRecord1 is not a scalar')
            
            zv2 = repmat(zv1, [1,nRepeatsPerRecord1]);
        end

        
        
        function newTt2000 = convert_N_to_1_SPR_Epoch( oldTt2000, nSpr, freqHzWithinRecords )
        % Convert time series zVariable (column) equivalent to converting N-->1 samples/record, assuming time increments
        % with frequency within each snapshot.
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % oldTt2000         : Nx1 vector.
        % nSpr              : Positive integer. Scalar. Number of values/samples per record (SPR).
        % freqWithinRecords : Nx1 vector. Frequency of samples within a subsequence (CDF record). Unit: Hz.
        % newTt2000         : Nx1 vector. Like oldTt2000 but each single time (row) has been replaced by a constantly
        %                     incrementing sequence of times (rows). Every such sequence begins with the original value,
        %                     has length nSpr with frequency freqWithinRecords(i).
        %                     NOTE: There is no check that the entire sequence is monotonic. LFR data can have snapshots
        %                           (i.e. snapshot records) that overlap in time!
            
        % PROPOSAL: Turn into more generic function, working on number sequences in general.
        % PROPOSAL: N_sequence should be a column vector.
        %    NOTE: TDS-LFM-RSWF, LFR-SURV-CWF have varying snapshot lengths.
        %    PRO: Could be useful for converting N->1 samples/record for calibration with transfer functions.
        %       NOTE: Then also needs function that does the reverse.
        % PROPOSAL: Replace by some simpler(?) algorithm that uses column/matrix multiplication.
            
            % ASSERTIONS
            bicas.proc_utils.assert_zv_Epoch(oldTt2000)
            if numel(nSpr) ~= 1
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'nSpr not scalar.')
            elseif size(freqHzWithinRecords, 1) ~= size(oldTt2000, 1)
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'freqWithinRecords and oldTt2000 do not have the same number of rows.')
            end
            
            nRecords = numel(oldTt2000);
            
            % Express frequency as period length in ns (since tt2000 uses ns as a unit).
            % Use the same MATLAB class as tt.
            % Unique frequency per record.
            periodNsColVec = int64(1e9 ./ freqHzWithinRecords);   
            periodNsMatrix = repmat(periodNsColVec, [1, nSpr]);
                        
            % Conventions:
            % ------------
            % Time unit: ns (as for tt2000)            
            % Algorithm should require integers to have a very predictable behaviour (useful when testing).
            
            % Times for the beginning of every record.
            tt2000RecordBeginColVec = oldTt2000;
            tt2000RecordBeginMatrix = repmat(tt2000RecordBeginColVec, [1, nSpr]);
            
            % Indices for within every record (start at zero for every record).
            iSampleRowVec = int64(0:(nSpr-1));
            iSampleMatrix = repmat(iSampleRowVec, [nRecords, 1]);
            
            % Unique time for every sample in every record.
            tt2000Matrix = tt2000RecordBeginMatrix + iSampleMatrix .* periodNsMatrix;
            
            % Convert to 2D matrix --> 1D column vector.
            %newTt2000 = reshape(tt2000Matrix', [nRecords*nSpr, 1]);
            newTt2000 = bicas.proc_utils.convert_N_to_1_SPR_redistribute(tt2000Matrix);
        end
        
        
        
        function ACQUISITION_TIME_2 = convert_N_to_1_SPR_ACQUISITION_TIME(  ACQUISITION_TIME_1, nSpr, freqWithinRecords, ACQUISITION_TIME_EPOCH_UTC )
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
            tt2000_2           = bicas.proc_utils.convert_N_to_1_SPR_Epoch(  tt2000_1,           nSpr, freqWithinRecords);
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


        
        % M                     : 2D matrix
        % nCopyColsPerRowVec    : 1D vector. {i}=Number of elements to copy from M{i,:}.
        % ca                    : Column cell array of 1D vectors.
        function ca = convert_matrix_to_cell_array_of_vectors(M, nCopyColsPerRowVec)
            EJ_library.assert.vector(nCopyColsPerRowVec)
            assert(ismatrix(M))
            assert(size(M, 1) == length(nCopyColsPerRowVec))
            
            ca = cell(size(M, 1), 1);
            for iRow = 1:numel(nCopyColsPerRowVec)
                ca{iRow} = M(iRow, 1:nCopyColsPerRowVec(iRow));
            end
        end
        

        
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
            DELTA_PLUS_MINUS = cast(DELTA_PLUS_MINUS, EJ_library.utils.convert_CDF_type_to_MATLAB_class(...
                ZV_DELTA_PLUS_MINUS_DATA_TYPE, 'Only CDF data types'));
        end
        
        
        
        % 2020-03-10: Function seems to not be used any more. DSAMP_TIME abolished from BIAS datasets?
%         function SAMP_DTIME = derive_SAMP_DTIME(freqHz, nSpr)
%         %
%         % ARGUMENTS
%         % =========
%         % freqHz     : Frequency column vector in s^-1. Can not handle freq=NaN since the output is an integer.
%         % nSpr       : Number of samples per record (SPR).
%         %
%         % RETURN VALUE
%         % ============
%         % SAMP_DTIME : Analogous to BIAS zVariable with CDF_UINT4=uint32. NOTE: Unit ns.
%         %
%         % ~BUG: The LFR/TDS/BIAS dataset skeletons specify that zVariable SAMP_DTIME is CDF_UINT4 in unit ns which should
%         % have a too small a range for some snapshots. Therefore, this conversion will eliminate most Example: LFR
%         % 2048/256 Hz = 8e9 ns > 2^31 ns ~ 2e9 ns.
%         % 2017-03-08: Xavier Bonnin (LESIA) and Bruno Katra (RPW/LFR) are aware of this and seem to have implemented it
%         % that way intentionally!!!
%         %
%         % IMPLEMENTATION NOTE: Algorithm should require integers in order to have a very predictable behaviour (useful
%         % when testing).
%         
%             % ASSERTIONS
%             if ~iscolumn(freqHz) || ~isfloat(freqHz) || any(isnan(freqHz))
%                 error('BICAS:proc_utils:Assertion:IllegalArgument', '"freqHz" is not a column vector of non-NaN floats.')
%             elseif ~isscalar(nSpr)
%                 error('BICAS:proc_utils:Assertion:IllegalArgument', '"nSpr" is not a scalar.')
%             end
%             
%             nRecords = size(freqHz, 1);
%             
%             % Express frequency as period length in ns (since tt2000 uses ns as a unit).
%             % Use the same MATLAB class as tt
%             % Unique frequency per record.
%             periodNsColVec = int64(1e9 ./ freqHz);   % Ns = ns = nanoseconds
%             periodNsMatrix = repmat(periodNsColVec, [1, nSpr]);
%                         
%             % Conventions:
%             % ------------
%             % Time unit: ns (as for tt2000)            
%             % Algorithm should require integers to have a very predictable behaviour (useful when testing).
%             
%             % Indices for within every record (start at zero for every record).
%             iSampleRowVec = int64(0:(nSpr-1));
%             iSampleMatrix = repmat(iSampleRowVec, [nRecords, 1]);
%             
%             % Unique time for every sample in every record.
%             SAMP_DTIME = iSampleMatrix .* periodNsMatrix;
%             
%             SAMP_DTIME = cast(SAMP_DTIME, EJ_library.utils.convert_CDF_type_to_MATLAB_class('CDF_UINT4',  'Only CDF data types'));
%         end
        
        
        
        % mSize : [nRows, nColumns, ...] so that the return value from the size() function can be used.
        function M = create_NaN_array(mSize)
            assert(numel(mSize) >= 2)
            
            M = NaN * zeros(mSize);
        end


        
        function y2 = nearest_interpolate_float_records(zvTt2000_1, y1, zvTt2000_2, method)
            % Interpolate ~zVariable to other points in time using nearest neighbour interpolation.
            % Values outside the interval covered by the old time series will be set to NaN.
            %
            % This is intended for interpolating HK and current values to SCI record times.
            %
            % IMPLEMENTATION NOTE: interp1 does seem to require y1 to be float. Using NaN as a "fill value" for the
            % return value imples that it too has to be a float.
            
            % PROPOSAL: Assertion for whether interpolating is possible.
            %   NOTE: Depends on method of interpolation / purpose.
            %       Ex: Currents. Currents can always be extrapolated forward in time, but not backward.
            %   PROPOSAL: method='previous' ==> check min
            
            switch(method)
                case 'nearest'
                    assert(...
                        EJ_library.utils.is_range_subset(zvTt2000_2, zvTt2000_1), ...
                        'BICAS:proc_utils:interpolate_float_records:Assertion:IllegalArgument', ...
                        'Can not interpolate data since the time range of zvTt2000_2 is not a subset of zvTt2000_1.')

                case 'previous'
                    % IMPLEMENTATION NOTE: Used for currents which can be extrapolate forward, but not backward.
                    assert(min(zvTt2000_1) <= min(zvTt2000_2), ...
                        'BICAS:proc_utils:nearest_interpolate_float_records:Assertion:IllegalArgument', ...
                        'Can not interpolate data since the time range of zvTt2000_2 does not begin after zvTt2000_1 begins.')

                otherwise
                    error(...
                        'BICAS:proc_utils:nearest_interpolate_float_records:Assertion:IllegalArgument', ...
                        'Illegal argument method="%s".', ...
                        method)
            end
            
            bicas.proc_utils.assert_zv_Epoch(zvTt2000_1)
            bicas.proc_utils.assert_zv_Epoch(zvTt2000_2)
            y2 = interp1(double(zvTt2000_1), y1, double(zvTt2000_2), method, NaN);
        end
        


        function utcStr = tt2000_to_UTC_str(zvTt2000)
        % Convert tt2000 value to UTC string with nanoseconds.
            
            bicas.proc_utils.assert_zv_Epoch(zvTt2000)
            
            utcStr = EJ_library.utils.CDF_tt2000_to_UTC_str(zvTt2000);
        end
        
        
        
        function ColumnStrs = log_array(varName, varValue, varType)
            % Logs statistics on the contents of a numeric variable (any dimensionality): Number of & percentage NaN, unique
            % values, min-max. Primarily intended for zVariables and derivatives thereof. Can be useful for knowing which
            % settings are used (e.g. DIFF_GAIN), constant/varying bias current, suspect input datasets.
            %
            % IMPLEMENTATION NOTE: Deliberately short function name to not clutter the log.
            %
            %
            % ARGUMENTS
            % =========
            % Alternative 1: 'explanation' (string literal) : Prints explanation of the (very condensed) log messages.
            % Alternative 2: variableName, variableValue    : Print statistics on line
            %
            % ASSUMPTIONS
            % ===========
            % variableValue is numeric
            
            % PROPOSAL: Handle fill/pad value?
            % PROPOSAL: Move to +utils.
            % PROPOSAL: Special log function for zVars. Can print CDF type (implicitly range).
            % PROPOSAL: Print MATLAB class (implicitly range).
            
            MAX_EPOCH_UNIQUES_PRINTED = 2;
            
            global SETTINGS
            
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
                    if nUniqueValues > SETTINGS.get_fv('LOGGING.MAX_UNIQUES_PRINTED')
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
                    
                    if nUniqueValues > MAX_EPOCH_UNIQUES_PRINTED
                        epochMinStr = bicas.proc_utils.tt2000_to_UTC_str(min(varValue));
                        epochMaxStr = bicas.proc_utils.tt2000_to_UTC_str(max(varValue));
                        valuesStr   = sprintf('Mm: %s -- %s', epochMinStr, epochMaxStr);
                    elseif nValues >= 1
                        for i = 1:numel(uniqueValues)
                            valueStrs{end+1} = bicas.proc_utils.tt2000_to_UTC_str(uniqueValues{i});
                        end
                        valuesStr = ['Us: ', strjoin(valueStrs, ', ')];
                    else
                        valuesStr = '-';
                    end
                    
                otherwise
                    error('BICAS:proc_utils', 'Illegal argument varType="%s"', varType)
            end
            
            % Assemble the final string
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
        function log_zVars(Zvs, L)
            % PROBLEM: Can not manually specify which variables are Epoch-like.
            % PROBLEM: Can not manually specify variable name strings.
            %   Ex: process_HK_to_HK_on_SCI_TIME: Print different versions of time for comparison. Want whitespace
            %
            % PROPOSAL: For min-max values, also print difference.
            %   Ex: Time difference for Epoch.
            %       TODO-DECISION: How print time difference?
            %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13:02:34
            %           PROPOSAL: Days-hours-minutes-seconds, e.g. 56 days, 13h02m34s
            
            LOG_LEVEL = 'debug';

            fnList     = fieldnames(Zvs);
            ColumnStrs = EJ_library.utils.empty_struct([0,1], 'name', 'size', 'nNan', 'percentageNan', 'nUniqueValues', 'values');
            
            for iFn = 1:numel(fnList)
                zvName  = fnList{iFn};
                zvValue = Zvs.(zvName);
                
                if iscolumn(zvValue) && isa(zvValue, 'int64') ...
                        && any(EJ_library.utils.regexpf(zvName, {'Epoch.*', '.*Epoch', '.*tt2000.*'}))
                    % CASE: Epoch-like variable.
                    
                    ColumnStrs(end+1) = bicas.proc_utils.log_array(zvName, zvValue, 'Epoch');
                    
                elseif isnumeric(zvValue)
                    % CASE: Non-Epoch-like numeric variable.
                    
                    ColumnStrs(end+1) = bicas.proc_utils.log_array(zvName, zvValue, 'numeric');
                    
                elseif ischar(zvValue)
                    
                    % Example of string valued (but irrelevant) CDF zVariables: ACQUISITION_TIME_LABEL
                    ;   % Ignore
                    
                else
                    error('BICAS:proc_utils:Assertion', 'Can not handle zVar "%s".', zvName)
                end
            end
            
            headerStrs = {'Name', 'Size', '#NaN', '%NaN', '#Uniq', 'Values'};
            tableStrs = {};
            tableStrs(:,1) = {ColumnStrs(:).name}';
            tableStrs(:,2) = {ColumnStrs(:).size}';
            tableStrs(:,3) = {ColumnStrs(:).nNan}';
            tableStrs(:,4) = {ColumnStrs(:).percentageNan}';
            tableStrs(:,5) = {ColumnStrs(:).nUniqueValues}';
            tableStrs(:,6) = {ColumnStrs(:).values}';
            tableColumnAdjustments = [{'left', 'left'}, repmat({'right'}, 1,3), {'left'}];
            [headerStrs, tableStrs, columnWidths] = EJ_library.utils.assist_print_table(headerStrs, tableStrs,  tableColumnAdjustments);

            L.log(LOG_LEVEL, strjoin(headerStrs, ' '))
            L.log(LOG_LEVEL, repmat('=', 1, sum(columnWidths) + numel(headerStrs) - 1))
            for iRow = 1:numel(ColumnStrs)
                L.log(LOG_LEVEL, strjoin(tableStrs(iRow, :), ' '))
            end
            L.logf(LOG_LEVEL, [...
                '    #NaN = Number of NaN\n', ...
                '    #Uniq = Number of unique values incl. NaN which counts as equal to itself.\n', ...
                '    Mm = min-max\n', ...
                '    Us = Unique values (explicitly listed)\n'])
        end



        % Assert that variable is an "zVar Epoch-like" variable.
        function assert_zv_Epoch(Epoch)
        % PROPOSAL: Change name: assert_zv_Epoch_zvar
        % PROPOSAL: Separate functions: assert_zv_Epoch_zvar, assert_zv_Epoch.
        
            if ~iscolumn(Epoch)
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Argument is not a column vector')   % Right ID?                
            elseif ~isa(Epoch, 'int64')
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Argument has the wrong class.')   % Right ID?
            end
            
            % Use?!!! Too processing heavy?!
            %validateattributes(Epoch, {'numeric'}, {'increasing'})
        end

        
        
        function assert_ACQUISITION_TIME(ACQUISITION_TIME)
        % Assert that variable is an "zVar ACQUISITION_TIME-like" variable.
        
            if ~isa(ACQUISITION_TIME, 'uint32')
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME is not uint32.')
            elseif ndims(ACQUISITION_TIME) ~= 2
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME is not 2D.')
            elseif size(ACQUISITION_TIME, 2) ~= 2
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME does not have two columns.')
            elseif any(ACQUISITION_TIME(:, 1) < 0)
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME has negative number of integer seconds.')
            elseif any(65536 <= ACQUISITION_TIME(:, 2))    % Does not need to check for negative values due to uint32.
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME subseconds out of range.')
            end
        end
        
        
        
        function assert_struct_num_fields_have_same_N_rows(S)
        % Assert that the below variables inside the struct have the same size in the first index.
        % (1) numeric fields                  : The field itself
        % (2) cell fields                     : Cell array components
        % (3) struct field (inside the struct): The inner struct's fields (not recursive)
        % Asserts that there are no other types of variables.
        %
        % Useful for structs where all fields represent CDF zVariables and/or derivatives thereof, the size in the first
        % index (number of CDF record) should be equal.
        
        % NOTE: Function name somewhat bad.
        % PROPOSAL: Make recursive?!
        % PROPOSAL: Implement using new features in EJ_library.assert.size.
        
            fieldNamesList1 = fieldnames(S);
            nRows = [];
            for iFn1 = 1:length(fieldNamesList1)
                fieldValue = S.(fieldNamesList1{iFn1});
                
                if isnumeric(fieldValue)
                    nRows(end+1) = size(fieldValue, 1);
                elseif iscell(fieldValue)
                    for iCc = 1:numel(fieldValue)
                        nRows(end+1) = size(fieldValue{iCc}, 1);
                    end
                elseif isstruct(fieldValue)
                    fieldNamesList2 = fieldnames(fieldValue);
                    for iFn2 = 1:length(fieldNamesList2)
                        nRows(end+1) = size(fieldValue.(fieldNamesList2{iFn2}), 1);
                    end
                else
                    error('BICAS:proc_utils:Assertion', 'Can not handle this type of struct field.')
                end
            end
            if length(unique(nRows)) > 1    % NOTE: length==0 valid for struct containing zero numeric fields.
                error('BICAS:proc_utils:Assertion', ...
                    'Numeric fields and cell array components in struct do not have the same number of rows (likely corresponding to CDF zVar records).')
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

