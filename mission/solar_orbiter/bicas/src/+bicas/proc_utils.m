% Collections of minor utility functions (in the form of static methods) used for data processing.
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
classdef proc_utils

%============================================================================================================
% PROPOSAL: Split up in separate files?!
% PROPOSAL: Move some functions to "utils".
%   Ex: add_rows_to_struct_fields, select_row_range_from_struct_fields
%   Ex: log_array, log_struct_array, log_tt2000_array (uses bicas.proc_utils_assert_Epoch)
%
% PROPOSAL: Write test code for ACQUISITION_TIME_to_tt2000 and its inversion.
% PROPOSAL: Reorg select_row_range_from_struct_fields into returning a list of intervals instead.
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
%   convert_N_to_1_SPR_repeat           -- repeat_in_record    + convert_N_to_1_SPR_redistribute
%   convert_N_to_1_SPR_Epoch            -- increment_in_record + convert_N_to_1_SPR_redistribute
%   convert_N_to_1_SPR_ACQUISITION_TIME -- Keep



    methods(Static, Access=public)

%         function S = select_row_range_from_struct_fields(S, iFirst, iLast)
%         % Given a struct, select a subset of that struct defined by a range of ROW indices for every field.
%         % Generic utility function.
%         %
%         % Compare add_rows_to_struct_fields. Name chosen in analogy.
%         
%             bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S);
%         
%             fieldNameList = fieldnames(S);
%             nRows = NaN;                   % Initial non-sensical value which is later replaced.
%             for i=1:length(fieldNameList)
%                 fn = fieldNameList{i};
%                 
%                 % ASSERTIONS
%                 if isnan(nRows)
%                     nRows = size(S.(fn), 1);
%                     if (nRows < iFirst) || (nRows < iLast)
%                         error('BICAS:proc_utils:Assertion', 'iFirst or iLast outside of interval of indices (rows).')
%                     end
%                 elseif nRows ~= size(S.(fn), 1)
%                    error('BICAS:proc_utils:Assertion', 'Not all struct fields have the same number of rows.')
%                 end
%                 
%                 S.(fn) = S.(fn)(iFirst:iLast, :, :);
%             end
%         end



        function c2 = select_row_range_from_cell_comps(c1, iFirst, iLast)
            % ASSERTIONS
            bicas.proc_utils.assert_cell_array_comps_have_same_N_rows(c1)
            
            for i = 1:numel(c1)
                c2{i} = c1{i}(iFirst:iLast, :, :);
            end
        end

        

        function S = add_rows_to_struct_fields(S, SAmendment)
        % Generic utility function.
        % Add values to every struct field by adding components after their highest row index (let them grow in
        % the row index).
        %
        % Compare select_row_range_from_struct_fields. Name chosen in analogy.
 

            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S);
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(SAmendment);
            
            fieldNamesList = fieldnames(SAmendment);
            for i=1:length(fieldNamesList)
                fn = fieldNamesList{i};
                
                S.(fn) = [S.(fn) ; SAmendment.(fn)];
            end
        end

        

        function freqHz = get_LFR_frequency(iLsf)
        % Convert LFR zVariable FREQ values to Hz. The usefulness of this function stems from how the LFR
        % datasets are defined.
        %
        % ARGUMENTS
        % =========
        % iLsf   : The LSF index, i.e. 1=LFR freq. F0, and so on.
        % freqHz : Frequency in Hz.

            global SETTINGS

            % ASSERTION
            uniqueValues = unique(iLsf);
            if ~all(ismember(uniqueValues, [1,2,3,4]))
                uniqueValuesStr = sprintf('%d', uniqueValues);   % NOTE: Has to print without \n to keep all values on a single-line string.
                error('BICAS:proc_utils:Assertion:IllegalArgument:DatasetFormat', ...
                    'Found unexpected values in LSF index (corresponding to LFR FREQ+1). Unique values: %s.', uniqueValuesStr)
            end
            
            % NOTE: Implementation that works for arrays of any size.
            freqHz = ones(size(iLsf)) * NaN;        % Allocate array and set default values.
            freqHz(iLsf==1) = SETTINGS.get_fv('PROCESSING.LFR.F0_HZ');
            freqHz(iLsf==2) = SETTINGS.get_fv('PROCESSING.LFR.F1_HZ');
            freqHz(iLsf==3) = SETTINGS.get_fv('PROCESSING.LFR.F2_HZ');
            freqHz(iLsf==4) = SETTINGS.get_fv('PROCESSING.LFR.F3_HZ');
        end
        
        
        
%         function GAMMA = get_simple_demuxer_gamma(DIFF_GAIN)
%         % Translate a scalar zVariable value DIFF_GAIN to an actual scalar gamma used in simplified calibration.
%         % NaN translates to NaN.
%         
%             global SETTINGS
%             
%             % ASSERTION
%             if numel(DIFF_GAIN) ~= 1
%                 error('BICAS:proc_utils:Assertion:IllegalArgument', 'Illegal argument value "DIFF_GAIN". Must be scalar (not array).')
%             end
% 
%             switch(DIFF_GAIN)
%                 case 0    ; GAMMA = SETTINGS.get_fv('PROCESSING.CALIBRATION.SCALAR.GAMMA_LOW_GAIN');
%                 case 1    ; GAMMA = SETTINGS.get_fv('PROCESSING.CALIBRATION.SCALAR.GAMMA_HIGH_GAIN');
%                 otherwise
%                     if isnan(DIFF_GAIN)
%                         GAMMA = NaN;
%                     else
%                         error('BICAS:proc_utils:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value "DIFF_GAIN"=%d.', DIFF_GAIN)                    
%                     end
%             end
%         end
        
        
        
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
%             tt2000 = spdfcomputett2000(SETTINGS.get_fv('PROCESSING.ACQUISITION_TIME_EPOCH_UTC')) + int64(atSeconds * 1e9);   % NOTE: spdfcomputett2000 returns int64 (as it should).
            tt2000 = spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) + int64(atSeconds * 1e9);   % NOTE: spdfcomputett2000 returns int64 (as it should).
        end
        

        
        function ACQUISITION_TIME = tt2000_to_ACQUISITION_TIME(tt2000, ACQUISITION_TIME_EPOCH_UTC)
        % Convert from tt2000 to ACQUISITION_TIME.
        %
        % t_tt2000 : Nx1 vector. Tequired to be int64 like the real zVar Epoch.
        % ACQUISITION_TIME : Nx2 vector. uint32.
        %       NOTE: ACQUSITION_TIME can not be negative since it is uint32.
        
            % ASSERTIONS
            bicas.proc_utils.assert_Epoch(tt2000)

            % NOTE: Important to type cast to double because of multiplication
%             atSeconds = double(int64(tt2000) - spdfcomputett2000(SETTINGS.get_fv('PROCESSING.ACQUISITION_TIME_EPOCH_UTC'))) * 1e-9;    % at = ACQUISITION_TIME
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
                EJ_library.utils.assert.vector(v)
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
            EJ_library.utils.assert.vector(iEdgeList)
            
            
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
        %
        % ARGUMENT
        % ========
        % v     : (iRecord, iSnapshotSample, iChannel)
        %
        % RETURN VALUE
        % ============
        % v     : (iRecord, iChannel)

            EJ_library.utils.assert.size(zv1, [NaN, NaN, NaN])
            
            zv  = permute(zv1, [2,1,3]);
            zv2 = reshape(zv, size(zv,1) * size(zv,2), size(zv,3));
            
            EJ_library.utils.assert.size(zv2, [NaN, NaN])
        end



        function newData = convert_N_to_1_SPR_repeat(oldData, nRepeatsPerOldRecord)
        % (1) Convert zVariable-like variable from 1 value/record to N values/record (same number of records) by repeating within record.
        % (2) Convert zVariable-like variable from N values/record to 1 value/record by redistributing values.
            
            % ASSERTIONS
            assert(iscolumn(oldData),              'BICAS:proc_utils:Assertion:IllegalArgument', 'oldData is not a column vector')
            assert(isscalar(nRepeatsPerOldRecord), 'BICAS:proc_utils:Assertion:IllegalArgument', 'nSamplesPerOldRecord is not a scalar')
            
            newData = repmat(oldData, [1,nRepeatsPerOldRecord]);
            newData = bicas.proc_utils.convert_N_to_1_SPR_redistribute(newData);
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
            bicas.proc_utils.assert_Epoch(oldTt2000)
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
            EJ_library.utils.assert.vector(nCopyColsPerRowVec)
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
            EJ_library.utils.assert.vector(ca)
            assert(isscalar(nMatrixColumns))
            EJ_library.utils.assert.vector(nMatrixColumns)
            
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
            DELTA_PLUS_MINUS = cast(DELTA_PLUS_MINUS, bicas.utils.convert_CDF_type_to_MATLAB_class('CDF_INT8',  'Only CDF data types'));
        end
        
        
        
        function SAMP_DTIME = derive_SAMP_DTIME(freqHz, nSpr)
        % freqHz     : Frequency column vector in s^-1. Can not handle freq=NaN since the output is an integer.
        % nSpr       : Number of samples per record (SPR).
        % SAMP_DTIME : Analogous to BIAS zVariable with CDF_UINT4=uint32. NOTE: Unit ns.
        %
        % ~BUG: The LFR/TDS/BIAS dataset skeletons specify that zVariable SAMP_DTIME is CDF_UINT4 in unit ns which should
        % have a too small a range for some snapshots. Therefore, this conversion will eliminate most Example: LFR
        % 2048/256 Hz = 8e9 ns > 2^31 ns ~ 2e9 ns.
        % 2017-03-08: Xavier Bonnin (LESIA) and Bruno Katra (RPW/LFR) are aware of this and seem to have implemented it
        % that way intentionally!!!
        %
        % IMPLEMENTATION NOTE: Algorithm should require integers in order to have a very predictable behaviour (useful
        % when testing).
        
            % ASSERTIONS
            if ~iscolumn(freqHz) || ~isfloat(freqHz) || any(isnan(freqHz))
                error('BICAS:proc_utils:Assertion:IllegalArgument', '"freqHz" is not a column vector of non-NaN floats.')
            elseif ~isscalar(nSpr)
                error('BICAS:proc_utils:Assertion:IllegalArgument', '"nSpr" is not a scalar.')
            end
            
            nRecords = size(freqHz, 1);
            
            % Express frequency as period length in ns (since tt2000 uses ns as a unit).
            % Use the same MATLAB class as tt
            % Unique frequency per record.
            periodNsColVec = int64(1e9 ./ freqHz);   % Ns = ns = nanoseconds
            periodNsMatrix = repmat(periodNsColVec, [1, nSpr]);
                        
            % Conventions:
            % ------------
            % Time unit: ns (as for tt2000)            
            % Algorithm should require integers to have a very predictable behaviour (useful when testing).
            
            % Indices for within every record (start at zero for every record).
            iSampleRowVec = int64(0:(nSpr-1));
            iSampleMatrix = repmat(iSampleRowVec, [nRecords, 1]);
            
            % Unique time for every sample in every record.
            SAMP_DTIME = iSampleMatrix .* periodNsMatrix;
            
            SAMP_DTIME = cast(SAMP_DTIME, bicas.utils.convert_CDF_type_to_MATLAB_class('CDF_UINT4',  'Only CDF data types'));
        end
        
        
        
        % mSize : [nRows, nColumns, ...] so that the return value from the size() function can be used.
        function M = create_NaN_array(mSize)
            assert(numel(mSize) >= 2)
            
            M = NaN * zeros(mSize);
        end



        function newData = nearest_interpolate_float_records(oldData, oldTt2000, newTt2000)
        % Interpolate ~zVariable to other points in time using nearest neighbour interpolation.
        % Values outside the interval covered by the old time series will be set to NaN.
        %
        % This is intended for interpolating HK values to SCI record times.
        %
        % IMPLEMENTATION NOTE: interp1 does seem to require oldData to be float. Using NaN as a "fill value" for the
        % return value imples that it too has to be a float.
            
            bicas.proc_utils.assert_Epoch(oldTt2000)
            bicas.proc_utils.assert_Epoch(newTt2000)
            newData = interp1(double(oldTt2000), oldData, double(newTt2000), 'nearest', NaN);
        end



        function utcStr = tt2000_to_UTC_str(tt2000)
        % Convert tt2000 value to UTC string with nanoseconds.
        %
        % Example: 2016-04-16T02:26:14.196334848
        % NOTE: This is the inverse to spdfparsett2000.
            
            bicas.proc_utils.assert_Epoch(tt2000)
            
            %  spdfbreakdowntt2000 converts the CDF TT2000 time, nanoseconds since
            %               2000-01-01 12:00:00 to UTC date/time.
            %
            %     OUT = spdfbreakdowntt2000(tt2000) returns the UTC date/time from CDF TT2000
            %     time. OUT is an array with each row having nine (9) numerical values
            %     for year, month, day, hour, minute, second, millisecond, microsecond
            %     and nanosecond.
            v = spdfbreakdowntt2000(tt2000);
            utcStr = sprintf('%04i-%02i-%02iT%02i:%02i:%2i.%03i%03i%03i', v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9));
        end
        
        
        
        % NOTE: Does not recognize HK datasets.
        % NOTE: Only classifies input datasets. (Is there a good reason for this?)
        % NOTE: Function deliberately ignores Skeleton_version.
        % IMPLEMENTATION NOTE: Still recognizes old ROC-SGSE datasets since they may be found in global attribute
        % DATASET_ID in old test files.
        %
        function C = classify_DATASET_ID(datasetId)
            % PROPOSAL: Use regexp instead.
            %   PRO: Can more easily handle old ROC-SGSE datasets.
            %
            % PROPOSAL: Implement assertsion on DATASET_ID via this function.
            %   Ex: bicas.assert_DATASET_ID
            %   Ex: bicas.swmde_defs.assert_DATASET_ID
            %   CON: Requires strict matching.
            %   PRO: Does not spread out the knowledge of DATASET_IDs.
            %   PROPOSAL: Flag for obsoleted DATASET_IDs that may be found in input datasets. Caller decides how to
            %       respond.
            % NEED?!: Some way of determining whether an obsoleted and current DATASET_ID are equivalent.
            
            EJ_library.utils.assert.castring(datasetId)
            
            C.isLfrSbm1 = 0;
            C.isLfrSbm2 = 0;
            C.isLfrSwf  = 0;
            C.isTdsCwf  = 0;
            C.isTdsRswf = 0;
            C.isL1      = 0;
            C.isL1R     = 0;
            
            % Shorten function name. MF = Matching Function
            mf = @(regexpPatternList) (any(EJ_library.utils.regexpf(datasetId, regexpPatternList)));
            
            if     mf({'(ROC-SGSE|SOLO)_L1R_RPW-LFR-SBM1-CWF-E', ...
                                   'SOLO_L1_RPW-LFR-SBM1-CWF'})
                C.isLfrSbm1  = 1;
            elseif mf({'(ROC-SGSE|SOLO)_L1R_RPW-LFR-SBM2-CWF-E', ...
                                   'SOLO_L1_RPW-LFR-SBM2-CWF'})
                C.isLfrSbm2  = 1;
            elseif mf({'(ROC-SGSE|SOLO)_L1R_RPW-LFR-SURV-CWF-E', ...
                                   'SOLO_L1_RPW-LFR-SURV-CWF'})
                % Do nothing
            elseif mf({'(ROC-SGSE|SOLO)_L1R_RPW-LFR-SURV-SWF-E', ...
                        '(ROC-SGSE|SOLO)_L1_RPW-LFR-SURV-SWF'})
                C.isLfrSwf  = 1;
            elseif mf({'(ROC-SGSE|SOLO)_L1R_RPW-TDS-LFM-CWF-E', ...
                                   'SOLO_L1_RPW-TDS-LFM-CWF'})
                C.isTdsCwf  = 1;
            elseif mf({'(ROC-SGSE|SOLO)_L1R_RPW-TDS-LFM-RSWF-E', ...
                                   'SOLO_L1_RPW-TDS-LFM-RSWF'})
                C.isTdsRswf = 1;
            else
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Illegal DATASET_ID. datasetId="%s"', datasetId)
            end
            
            if     mf({'(ROC-SGSE|SOLO)_L1_RPW-.*'})
                C.isL1      = 1;
            elseif mf({'(ROC-SGSE|SOLO)_L1R_RPW-.*'})
                C.isL1R     = 1;
            else
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Illegal DATASET_ID. datasetId="%s"', datasetId)
            end
            
            EJ_library.utils.assert.struct2(C, {...
                'isLfrSbm1', ...
                'isLfrSbm2', ...
                'isLfrSwf', ...
                'isTdsCwf', ...
                'isTdsRswf', ...
                'isL1', ...
                'isL1R'}, {})
        end
        
        
        
        function log_array(varargin)
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
        %
        
            global SETTINGS                
        
            if nargin == 1
                % ASSERTION
                if ~strcmp(varargin{1}, 'explanation')
                    error('BICAS:proc_utils:Assertion:IllegalArgument', 'Wrong number of arguments')
                end
                    
                EXPLANATION_STRING = 'Explanation for variable log messages: (x,y, ...)=size of variable; #=number of ...; Us=unique values (incl. NaN which counts as equal to itself); Mm=min-max';
                bicas.log('info', EXPLANATION_STRING)
                return
                
            elseif nargin == 2
                
                variableName  = varargin{1};
                variableValue = varargin{2};
                
                % ASSERTIONS
                if ~isnumeric(variableValue) || ndims(variableValue) > 3   % NOTE: min-max limit number of dimensions.
                    error('BICAS:proc_utils:Assertion:IllegalArgument', 'v is not numerical with max 3 dimensions.')
                end
                
                nValues       = numel(variableValue);
                nUniqueValues = length(bicas.utils.unique_values_NaN(variableValue));
                nNan          = sum(isnan(variableValue(:)));
                
                %============================
                % Construct string: NaN info
                %============================
                if nValues == 0
                    nanStr = '';
                else
                    nanStr = sprintf('#NaN=%3d%%=%d', round((nNan/numel(variableValue))*100), nNan);
                end
                
                %=================================
                % Construct string: variable size
                %=================================
                % Create comma-separated list of numbers.
                sizeStr = strjoin(arrayfun(@(n) num2str(n), size(variableValue), 'UniformOutput', 0),',');
                
                %===================================
                % Construct string: range of values
                %===================================
                if nUniqueValues > SETTINGS.get_fv('LOGGING.MAX_UNIQUES_PRINTED')
                    vMin = min(min(min(variableValue)));
                    vMax = max(max(max(variableValue)));
                    valuesStr = sprintf('Mm: %d--%d', vMin, vMax);
                else
                    if nUniqueValues == 0
                        valuesStr = '';
                    else
                        valuesStr = ['Us: ', sprintf('%d ', bicas.utils.unique_values_NaN(variableValue))];
                    end
                end
                
                %========================================================================================
                % Assemble the final string
                % -------------------------
                % Examples for choosing column sizes:
                % Long variable names:       HK_BIA_MODE_BIAS3_ENABLED
                %                            <V03_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E>.QUALITY_BITMASK
                %                            <V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF>.BIAS_MODE_BIAS1_ENABLED
                % Long variable size string: (90,90,2048)
                %========================================================================================
                outputStr = sprintf('%-61s (%-10s): #Us=%5d (%-16s) %s', variableName, sizeStr, nUniqueValues, nanStr, valuesStr);
                
                bicas.log('info', outputStr)
            else
                
                % ASSERTION
                error('BICAS:proc_utils:Assertion:IllegalArgument', 'Wrong number of arguments')
                
            end
        end
        
        
        
        function log_struct_arrays(variableName, variableValue)
            
            bicas.proc_utils.log_array('explanation')
            log_struct_arrays_INNER(variableName, variableValue)
        
            function log_struct_arrays_INNER(variableName, variableValue)
                % Call log_array recursively for struct.
                %
                % NOTE: Special case for variables/fields named "Epoch" of type int64.
                
                if iscolumn(variableValue) && isa(variableValue, 'int64') && ~isempty(regexp(variableName, 'Epoch$', 'once'))
                    
                    bicas.proc_utils.log_tt2000_array(variableName, variableValue);
                    
                elseif isnumeric(variableValue)
                    
                    bicas.proc_utils.log_array(variableName, variableValue)
                    
                elseif isstruct(variableValue)
                    
                    fieldNamesList = fieldnames(variableValue);
                    for i = 1:length(fieldNamesList)
                        fieldName = fieldNamesList{i};
                        
                        % NOTE: RECURSIVE CALL
                        log_struct_arrays_INNER(...
                            [variableName, '.', fieldName], ...
                            variableValue.(fieldName))
                    end
                    
                elseif ischar(variableValue)
                    
                    % Example of string valued (but irrelevant) CDF zVariables: ACQUISITION_TIME_LABEL
                    ;   % Ignore
                    
                else
                    
                    error('BICAS:proc_utils:Assertion', 'variableValue is neither numeric nor struct.')
                    
                end
            end
            
        end
        
        
        
        function log_tt2000_array(variableName, tt2000)
        % Log summary of series of times.
        %
        % tt2000 : A vector of tt2000 values.
        %
        % NOTE: Assumes that t is sorted in time, increasing.
        % NOTE: Can handle zero values.
        
        % PROPOSAL: Move to +utils.
        
            bicas.proc_utils.assert_Epoch(tt2000)
            
            if ~isempty(tt2000)
                strFirst = bicas.proc_utils.tt2000_to_UTC_str(tt2000(1));
                strLast  = bicas.proc_utils.tt2000_to_UTC_str(tt2000(end));
                bicas.logf('info', '%s: %s -- %s', variableName, strFirst, strLast)
            else
                bicas.logf('info', '%s: <empty>', variableName)
            end
        end

        
        
        % Assert that variable is an "zVar Epoch-like" variable.
        function assert_Epoch(Epoch)
        % PROPOSAL: Change name: assert_Epoch_zvar
        % PROPOSAL: Separate functions: assert_Epoch_zvar, assert_Epoch.
        
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
            EJ_library.utils.assert.all_equal( nRowsArray )
        end
        
        
        
        function doOverlap = ranges_overlap(v1, v2)
            EJ_library.utils.assert.vector(v1)
            EJ_library.utils.assert.vector(v2)
            
            doOverlap = (min(v2) <= max(v1)) && (min(v1) <= max(v2));
        end
        
        
        
    end   % Static
    
    
    
end

