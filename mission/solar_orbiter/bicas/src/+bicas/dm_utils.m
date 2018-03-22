classdef dm_utils
% Collections of minor utility functions (in the form of static methods) used by data_manager.
% The functions are collected here to reduce the size of data_manager.
%
% dm_utils = data_manager utilities
%
% SPR = Samples per record
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-10

%============================================================================================================
% PROPOSAL: Split up in separate files?!
% PROPOSAL: Move some functions to "utils".
%   Ex: add_components_to_struct, select_subset_from_struct
%   Ex: log_array, log_struct_array, log_tt2000_array (uses bicas.dm_utils_assert_Epoch)
%
% PROPOSAL: Write test code for ACQUISITION_TIME_to_tt2000 and its inversion.
% PROPOSAL: Reorg select_subset_from_struct into returning a list of intervals instead.
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

        function s = select_subset_from_struct(s, iFirst, iLast)
        % Given a struct, select a subset of that struct defined by a range of ROW indices for every field.
        % Generic utility function.
        
        % PROPOSAL: Use ~assert_unvaried_N_rows.
        % PROPOSAL: Better name. select_struct_columns_subset/range/interval.
        
            fieldNameList = fieldnames(s);
            nRows = NaN;                   % Initial non-sensical value which is later replaced.
            for i=1:length(fieldNameList)
                fn = fieldNameList{i};
                
                % ASSERTIONS
                if isnan(nRows)
                    nRows = size(s.(fn), 1);
                    if (nRows < iFirst) || (nRows < iLast)
                        error('BICAS:dm_utils:Assertion', 'iFirst or iLast outside of interval of indices (rows).')
                    end
                elseif nRows ~= size(s.(fn), 1)
                   error('BICAS:dm_utils:Assertion', 'Not all struct fields have the same number of rows.')
                end
                
                s.(fn) = s.(fn)(iFirst:iLast, :, :);
            end
        end
        
        

        function s = add_components_to_struct(s, structAmendment)
        % Add values to every struct field by adding components after their highest row index (let them grow in
        % the row index).
        
        % PROPOSAL: Better name. ~rows, ~fields, ~records
        %   Ex: add_row_components_to_struct_fields
            
            % Generic utility function.
            fieldNamesList = fieldnames(structAmendment);
            for i=1:length(fieldNamesList)
                fn = fieldNamesList{i};
                
                s.(fn) = [s.(fn) ; structAmendment.(fn)];
            end
        end



        function freq = get_LFR_frequency(FREQ)
        % Convert LFR zVariable FREQ values to Hz. The usefulness of this function stems from how the LFR
        % datasets are defined.
        %
        % FREQ : The FREQ zVariable in LFR CDFs (contains constants representing frequencies, themselves NOT being frequencies).
        % freq : Frequency in Hz.
            
            global SETTINGS
            
            % ASSERTION
            unique_values = unique(FREQ);
            if ~all(ismember(unique_values, [0,1,2,3]))
                unique_values_str = sprintf('%d', unique_values);   % NOTE: Has to print without \n to keep all values on a single-line string.
                error('BICAS:dm_utils:Assertion:IllegalArgument:DatasetFormat', 'Found unexpected values in (LFR) FREQ (unique values: %s).', unique_values_str)
            end
            
            % NOTE: Implementation that works for arrays of any size.
            freq = ones(size(FREQ)) * -1;        % Allocate array and set default values.
            freq(FREQ==0) = SETTINGS.get_fv('LFR.F0');
            freq(FREQ==1) = SETTINGS.get_fv('LFR.F1');
            freq(FREQ==2) = SETTINGS.get_fv('LFR.F2');
            freq(FREQ==3) = SETTINGS.get_fv('LFR.F3');
        end
        
        
        
        function GAMMA = get_simple_demuxer_gamma(DIFF_GAIN)
        % Translate a scalar zVariable value DIFF_GAIN to an actual scalar gamma used in simplified calibration.
        % NaN translates to NaN.
        
            global SETTINGS
            
            % ASSERTION
            if numel(DIFF_GAIN) ~= 1
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Illegal argument value "DIFF_GAIN". Must be scalar (not array).')
            end
            
            switch(DIFF_GAIN)
                case 0    ; GAMMA = SETTINGS.get_fv('SIMPLE_DEMUXER.GAMMA_LOW_GAIN');
                case 1    ; GAMMA = SETTINGS.get_fv('SIMPLE_DEMUXER.GAMMA_HIGH_GAIN');
                otherwise
                    if isnan(DIFF_GAIN)
                        GAMMA = NaN;
                    else
                        error('BICAS:dm_utils:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value "DIFF_GAIN"=%d.', DIFF_GAIN)                    
                    end
            end
        end
        
        
        
        function Rx = get_LFR_Rx(R0, R1, R2, FREQ)
        % Return the relevant value of LFR CDF zVariables R0, R1, or R2, or a hypothetical but analogous "R3" which is always 1.
        %
        % R0, R1, R2, FREQ : LFR CDF zVariables. All must have identical array sizes.
        % Rx               : Same size array as R0, R1, R2, FREQ. The relevant values are copied, respectively, from
        %                    R0, R1, R2, or an analogous hypothetical "R3" that is a constant (=1) depending on
        %                    the value of FREQ in the corresponding component.
        %
        % NOTE: Works for all array sizes.
            
            Rx = -ones(size(FREQ));        % Set to -1 (should always be overwritten).
            
            I = (FREQ==0); Rx(I) = R0(I);
            I = (FREQ==1); Rx(I) = R1(I);
            I = (FREQ==2); Rx(I) = R2(I);
            I = (FREQ==3); Rx(I) = 1;      % The value of a hypothetical (non-existant, constant) analogous zVariable "R3".
        end



        function tt2000 = ACQUISITION_TIME_to_tt2000(ACQUISITION_TIME)
        % Convert time in from ACQUISITION_TIME to tt2000 which is used for Epoch in CDF files.
        % 
        % NOTE: t_tt2000 is in int64.
        % NOTE: ACQUSITION_TIME can not be negative since it is uint32.
        
            global SETTINGS
            
            bicas.dm_utils.assert_ACQUISITION_TIME(ACQUISITION_TIME)
            
            ACQUISITION_TIME = double(ACQUISITION_TIME);
            atSeconds = ACQUISITION_TIME(:, 1) + ACQUISITION_TIME(:, 2) / 65536;   % at = ACQUISITION_TIME
            tt2000 = spdfcomputett2000(SETTINGS.get_fv('ACQUISITION_TIME_EPOCH_UTC')) + int64(atSeconds * 1e9);   % NOTE: spdfcomputett2000 returns int64 (as it should).
        end
        
        
        
        function ACQUISITION_TIME = tt2000_to_ACQUISITION_TIME(tt2000)
        % Convert from tt2000 to ACQUISITION_TIME.
        %
        % t_tt2000 : Nx1 vector. Tequired to be int64 like the real zVar Epoch.
        % ACQUISITION_TIME : Nx2 vector. uint32.
        %       NOTE: ACQUSITION_TIME can not be negative since it is uint32.
        
            global SETTINGS
            
            % ASSERTIONS
            bicas.dm_utils.assert_Epoch(tt2000)

            % NOTE: Important to type cast to double because of multiplication
            atSeconds = double(int64(tt2000) - spdfcomputett2000(SETTINGS.get_fv('ACQUISITION_TIME_EPOCH_UTC'))) * 1e-9;    % at = ACQUISITION_TIME
            
            % ASSERTION: ACQUISITION_TIME must not be negative.
            if any(atSeconds < 0)
                error('BICAS:dm_manager:Assertion:IllegalArgument:DatasetFormat', 'Can not produce ACQUISITION_TIME (uint32) with negative number of integer seconds.')
            end
            
            atSeconds = round(atSeconds*65536) / 65536;
            atSecondsFloor = floor(atSeconds);
            
            ACQUISITION_TIME = uint32([]);
            ACQUISITION_TIME(:, 1) = uint32(atSecondsFloor);
            ACQUISITION_TIME(:, 2) = uint32((atSeconds - atSecondsFloor) * 65536);
            % NOTE: Should not be able to produce ACQUISITION_TIME(:, 2)==65536 (2^16) since atSeconds already rounded (to parts of 2^-16).
        end



        function [iFirstList, iLastList] = find_sequences(varargin)
        % For a non-empty set of column vectors, find all subsequences of continuously constant values in all the vectors.
        % Useful for finding continuous sequences of (CDF) records with identical settings.
        % NOTE: NaN counts as equal to itself.
        %
        % ARGUMENTS
        % ==================
        % Arguments (varargin{i}) : Column arrays of the same size. Can have zero rows (but must still have one column;
        %                           size 0xN). Does not have to be numeric as long as isequaln can handle it.
        % iFirst, iLast           : Vectors with indices to the first and last index of each sequence.
        
            % ASSERTIONS & variable extraction
            if isempty(varargin)
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'There are no vectors to look for sequences in.')
            end
            nRecords = size(varargin{1}, 1);
            for kArg = 1:length(varargin)
                if ~iscolumn(varargin{kArg}) || nRecords ~= size(varargin{kArg}, 1)
                    error('BICAS:dm_utils:Assertion:IllegalArgument', 'varargins are not all same-size column vectors.')
                end
            end                
            
            iFirstList = [];
            iLastList  = [];
            iFirst = 1;
            iLast = iFirst;
            while iFirst <= nRecords
                
                while iLast+1 <= nRecords       % For as long as there is another row...
                    foundLast = false;
                    for kArg = 1:length(varargin)
                        if ~isequaln(varargin{kArg}(iFirst), varargin{kArg}(iLast+1))    % NOTE: Use "isequaln" that treats NaN as any other value.
                            % CASE: This row is different from the previous one.
                            foundLast = true;
                            break
                        end
                    end
                    if foundLast
                        break
                    end
                    
                    iLast = iLast + 1;                    
                end

                iFirstList(end+1) = iFirst;
                iLastList(end+1)  = iLast;

                iFirst = iLast + 1;
            end
        end
        
        
        
        function filteredData = filter_rows(data, rowFilter)
        % Function intended for filtering out (copying selectively) data from a zVariable. Not copied values are set to
        % NaN.
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % data         : Numeric array with N rows.                 (Intended to represent zVariables with N records.)
        % rowFilter    : Numeric/logical column vector with N rows. (Intended to represent zVariables with N records.)
        % filteredData : Array of the same size as "records", such that
        %                filteredData(i,:,:, ...) == NaN,                 for record_filter(i)==0.
        %                filteredData(i,:,:, ...) == records(i,:,:, ...), for record_filter(i)~=0.

            % ASSERTIONS
            if ~iscolumn(rowFilter)     % Not really necessary to require row vector, only 1D vector.
                error('BICAS:dm_utils:Assertion:IllegalArgument', '"rowFilter" is not a column vector.')
            elseif size(rowFilter, 1) ~= size(data, 1)
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Numbers of records do not match.')
            elseif ~isfloat(data)
                error('BICAS:dm_utils:Assertion:IllegalArgument', '"data" is not a floating-point class (can not represent NaN).')
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



        function newData = convert_N_to_1_SPR_redistribute(oldData)
        % Convert zVariable-like variable from N samples/record to 1 sample/record (from a matrix to a column vector).
        
            % ASSERTIONS
            % NOTE: ndims always returns at least two, which is exactly what we want, also for empty and scalars, and row vectors.
            if ndims(oldData) > 2
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'oldData has more than two dimensions.')
            end
            
            newData = reshape(oldData', numel(oldData), 1);
        end
        
        
        
        function newData = convert_N_to_1_SPR_repeat(oldData, nRepeatsPerOldRecord)
        % (1) Convert zVariable-like variable from 1 value/record to N values/record (same number of records) by repeating within record.
        % (2) Convert zVariable-like variable from N values/record to 1 value/record by redistributing values.
            
            % ASSERTIONS
            if ~(iscolumn(oldData))
                error('BICAS:dm_utils:Assertion', 'oldData is not a column vector')
            elseif ~isscalar(nRepeatsPerOldRecord)
                error('BICAS:dm_utils:Assertion', 'nSamplesPerOldRecord is not a scalar')
            end
            
            newData = repmat(oldData, [1,nRepeatsPerOldRecord]);
            %newData = reshape(newData', [numel(newData), 1]);     % NOTE: Must transpose first.
            newData = bicas.dm_utils.convert_N_to_1_SPR_redistribute(newData);
        end
        
        
        
        function newTt2000 = convert_N_to_1_SPR_Epoch( oldTt2000, nSpr, freqHzWithinRecords )
        % Convert time series zVariable (column) equivalent to converting N-->1 samples/record, assuming time increments
        % with frequency in each snapshot.
        %
        % oldTt2000  : Nx1 vector.
        % newTt2000  : Nx1 vector. Like oldTt2000 but each single time (row) has been replaced by a constantly
        %              incrementing sequence of times (rows). Every such sequence begins with the original value, has
        %              length nSpr with frequency freqWithinRecords(i).
        %              NOTE: There is no check that the entire sequence is monotonic. LFR data can have snapshots (i.e.
        %              snapshot records) that overlap in time!
        % nSpr                    : Positive integer. Scalar. Number of values/samples per record (SPR).
        % freqWithinRecords  : Nx1 vector. Frequency of samples within a subsequence (CDF record). Unit: Hz.
            
        % PROPOSAL: Turn into more generic function, working on number sequences in general.
        % PROPOSAL: N_sequence should be a column vector.
        %    NOTE: TDS-LFM-RSWF, LFR-SURV-CWF have varying snapshot lengths.
        %    PRO: Could be useful for converting N->1 samples/record for calibration with transfer functions.
        %       NOTE: Then also needs function that does the reverse.
        % PROPOSAL: Replace by some simpler(?) algorithm that uses column/matrix multiplication.
            
            % ASSERTIONS
            bicas.dm_utils.assert_Epoch(oldTt2000)
            if numel(nSpr) ~= 1
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'nSpr not scalar.')
            elseif size(freqHzWithinRecords, 1) ~= size(oldTt2000, 1)
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'freqWithinRecords and oldTt2000 do not have the same number of rows.')
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
            newTt2000 = bicas.dm_utils.convert_N_to_1_SPR_redistribute(tt2000Matrix);
        end
        
        
        
        function ACQUISITION_TIME_2 = convert_N_to_1_SPR_ACQUISITION_TIME(  ACQUISITION_TIME_1, nSpr, freqWithinRecords  )
        % Function intended for converting ACQUISITION_TIME (always one time per record) from many samples/record to one
        % sample/record. See convert_N_to_1_SPR_Epoch which is analogous.
        % 
        % ACQUISITION_TIME_1 : Nx2 vector.
        % ACQUISITION_TIME_2 : Nx2 vector.

        % Command-line algorithm "test code":
        % clear; t_rec = [1;2;3;4]; f = [5;1;5;20]; N=length(t_rec); M=5; I_sample=repmat(0:(M-1), [N, 1]); F=repmat(f, [1,M]); T_rec = repmat(t_rec, [1,M]); T = T_rec + I_sample./F; reshape(T', [numel(T), 1])
            
            % ASSERTIONS
            bicas.dm_utils.assert_ACQUISITION_TIME(ACQUISITION_TIME_1)

            tt2000_1 = bicas.dm_utils.ACQUISITION_TIME_to_tt2000(ACQUISITION_TIME_1);
            tt2000_2 = bicas.dm_utils.convert_N_to_1_SPR_Epoch(tt2000_1, nSpr, freqWithinRecords);
            ACQUISITION_TIME_2 = bicas.dm_utils.tt2000_to_ACQUISITION_TIME(tt2000_2);
        end
        
        
        
        function DELTA_PLUS_MINUS = derive_DELTA_PLUS_MINUS(freqHz, nSpr)
        % freqHz  : Frequency column vector in s^-1. Can not handle freqHz=NaN since the output is an integer.
        % nSpr    : Number of samples/record.
        % DELTA_PLUS_MINUS : Analogous to BIAS zVariable. CDF_INT8=int64. NOTE: Unit ns.
            
            if ~iscolumn(freqHz) || ~isfloat(freqHz) || any(isnan(freqHz))
                error('BICAS:dm_utils:Assertion:IllegalArgument', '"freqHz" is not a column vector of non-NaN floats.')
            elseif ~isscalar(nSpr)
                error('BICAS:dm_utils:Assertion:IllegalArgument', '"nSpr" is not a scalar.')
            end
            
            nRecords = size(freqHz, 1);
            DELTA_PLUS_MINUS = zeros([nRecords, nSpr]);
            for i = 1:length(freqHz)
                DELTA_PLUS_MINUS(i, :) = 1/freqHz(i) * 1e9 * 0.5;      % Seems to work for more than 2D.
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
                error('BICAS:dm_utils:Assertion:IllegalArgument', '"freqHz" is not a column vector of non-NaN floats.')
            elseif ~isscalar(nSpr)
                error('BICAS:dm_utils:Assertion:IllegalArgument', '"nSpr" is not a scalar.')
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
        
        
        
        function M = create_NaN_array(nRowsColumns)
            % Input on the format [nRows, nColumns] so that the return value from the size() function can be used.
            
            M = NaN * zeros(nRowsColumns(1), nRowsColumns(2));
        end



        function newData = nearest_interpolate_float_records(oldData, oldTt2000, newTt2000)
        % Interpolate ~zVariable to other points in time using nearest neighbour interpolation.
        % Values outside the interval covered by the old time series will be set to NaN.
        %
        % This is intended for interpolating HK values to SCI record times.
        %
        % IMPLEMENTATION NOTE: interp1 does seem to require oldData to be float. Using NaN as a "fill value" for the
        % return value imples that it too has to be a float.
            
            bicas.dm_utils.assert_Epoch(oldTt2000)
            bicas.dm_utils.assert_Epoch(newTt2000)
            newData = interp1(double(oldTt2000), oldData, double(newTt2000), 'nearest', NaN);
        end



        function utcStr = tt2000_to_UTC_str(tt2000)
        % Convert tt2000 value to UTC string with nanoseconds.
        %
        % Example: 2016-04-16T02:26:14.196334848
        % NOTE: This is the inverse to spdfparsett2000.
            
            bicas.dm_utils.assert_Epoch(tt2000)
            
            v = spdfbreakdowntt2000(tt2000);
            utcStr = sprintf('%04i-%02i-%02iT%02i:%02i:%2i.%03i%03i%03i', v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9));
        end
        
        
        
        function uniqueValues = unique_values_NaN(A)
        % Return number of unique values in array, treating +Inf, -Inf, and NaN as equal to themselves (separately).
        % (MATLAB's "unique" function does not do this for NaN.)
        %
        % NOTE: Should work for all dimensionalities.
           
        % PROPOSAL: Move to +utils.
            
            % NOTE: "unique" has special behaviour which must be taken into account:
            % 1) Inf and -Inf are treated as equal to themselves.
            % 2) NaN is treated as if it is NOT equal itself. ==> Can thus return multiple instances of NaN.
            % 3) "unique" always puts NaN at the then of the vector of unique values (one or many NaN).
            uniqueValues = unique(A);
            
            % Remove all NaN unless it is found in the last component (save one legitimate occurrence of NaN, if there is any).
            % NOTE: Does work for empty matrices.
            uniqueValues(isnan(uniqueValues(1:end-1))) = [];
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
                    error('BICAS:dm_utils:Assertion:IllegalArgument', 'Wrong number of arguments')
                end
                    
                EXPLANATION_STRING = 'Explanation for variable log messages: (x,y, ...)=size of variable; #=Number of ...; Us=Unique values (incl. NaN which counts as equal to itself); Mm=Min-max';
                irf.log('n', EXPLANATION_STRING)
                return
                
            elseif nargin == 2
                
                variableName  = varargin{1};
                variableValue = varargin{2};
                
                % ASSERTIONS
                if ~isnumeric(variableValue) || ndims(variableValue) > 3   % NOTE: min-max limit number of dimensions.
                    error('BICAS:dm_utils:Assertion:IllegalArgument', 'v is not numerical with max 3 dimensions.')
                end
                
                nValues       = numel(variableValue);
                nUniqueValues = length(bicas.dm_utils.unique_values_NaN(variableValue));
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
                        valuesStr = ['Us: ', sprintf('%d ', bicas.dm_utils.unique_values_NaN(variableValue))];
                    end
                end
                
                %======================================================
                % Assemble the final string
                % -------------------------
                % Examples for choosing column sizes:
                % Long variable names:       HK_BIA_MODE_BIAS3_ENABLED
                %                            <L2S_LFR-SURV-CWF-E_V02>.QUALITY_BITMASK
                %                            <L2R_LFR-SURV-CWF_V01>.BIAS_MODE_BIAS1_ENABLED
                % Long variable size string: (90,90,2048)
                %======================================================
                outputStr = sprintf('%-46s (%-10s): #Us=%5d (%-16s) %s', variableName, sizeStr, nUniqueValues, nanStr, valuesStr);
                
                irf.log('n', outputStr)
            else
                
                % ASSERTION
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Wrong number of arguments')
                
            end
        end
        
        
        
        function log_struct_arrays(variableName, variableValue)
            
            bicas.dm_utils.log_array('explanation')
            log_struct_arrays_INNER(variableName, variableValue)
        
            function log_struct_arrays_INNER(variableName, variableValue)
                % Call log_array recursively for struct.
                %
                % NOTE: Special case for variables/fields named "Epoch" of type int64.
                
                if iscolumn(variableValue) && isa(variableValue, 'int64') && ~isempty(regexp(variableName, 'Epoch$'))
                    
                    bicas.dm_utils.log_tt2000_array(variableName, variableValue);
                    
                elseif isnumeric(variableValue)
                    
                    bicas.dm_utils.log_array(variableName, variableValue)
                    
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
                    
                    error('BICAS:dm_utils:Assertion', 'variableValue is neither numeric nor struct.')
                    
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
        
            bicas.dm_utils.assert_Epoch(tt2000)
            
            if ~isempty(tt2000)
                strFirst = bicas.dm_utils.tt2000_to_UTC_str(tt2000(1));
                strLast  = bicas.dm_utils.tt2000_to_UTC_str(tt2000(end));
                irf.log('n', sprintf('%s: %s -- %s', variableName, strFirst, strLast))
            else
                irf.log('n', sprintf('%s: <empty>', variableName))
            end
        end

        
        
        function assert_Epoch(Epoch)
        % Assert that variable is an "zVar Epoch-like" variable.
        
            if ~iscolumn(Epoch)
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Argument is not a column vector')   % Right ID?                
            elseif ~isa(Epoch, 'int64')
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Argument has the wrong class.')   % Right ID?
            end
        end

        
        
        function assert_ACQUISITION_TIME(ACQUISITION_TIME)
        % Assert that variable is an "zVar ACQUISITION_TIME-like" variable.
        
            if ~isa(ACQUISITION_TIME, 'uint32')
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME is not uint32.')
            elseif ndims(ACQUISITION_TIME) ~= 2
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME is not 2D.')
            elseif size(ACQUISITION_TIME, 2) ~= 2
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME does not have two columns.')
            elseif any(ACQUISITION_TIME(:, 1) < 0)
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME has negative number of integer seconds.')
            elseif any(65536 <= ACQUISITION_TIME(:, 2))    % Does not need to check for negative values due to uint32.
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'ACQUISITION_TIME subseconds out of range.')
            end
        end
        
        
        
        function assert_unvaried_N_rows(s)
        % Assert that all NUMERIC fields in a structure have the same number of rows.
        %
        % Useful since in data_manager, much code assumes that struct fields represent CDF zVar records which should
        % have the same number of rows.
        %
        % s : A struct to be tested.
            
            % PROPOSAL: Better name.
            %   Ex: _equal_rows, _equal_N_rows, _same_N_rows, _equal_nbr_of_rows
            
            fieldNamesList = fieldnames(s);
            nRows = [];
            for i = 1:length(fieldNamesList)
                fn = fieldNamesList{i};
                
                if isnumeric(s.(fn))
                    nRows(end+1) = size(s.(fn), 1);
                end
            end
            if length(unique(nRows)) > 1    % NOTE: length==0 valid for struct containing zero numeric fields.
                error('BICAS:dm_utils:Assertion', 'Numeric fields in struct do not have the same number of rows (likely corresponding to CDF zVar records).')
            end
        end
        
    end   % Static
    
    
    
end

