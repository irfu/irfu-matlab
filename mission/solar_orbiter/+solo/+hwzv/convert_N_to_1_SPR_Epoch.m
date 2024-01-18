%
% Convert time series zVariable (column) equivalent to converting
% N-->1_samples/record, assuming constant time increments with sampling
% frequency within each snapshot.
%
% NOTE: Function is in particular meant for converting a zVar Epoch with one
% snapshot per record, to a zVar-like Epoch with one (snapshot) sample per
% record.
%
% NOTE: Function can be meaningful also for scalar TT2000 values when applied to
% one single snapshot.
%
%
% ARGUMENTS
% =========
% oldTt2000         : Nx1 vector.
% nSpr              : Positive integer. Scalar. Number of values/samples per
%                     record (SPR).
% freqWithinRecords : Nx1 vector. Frequency of samples within a subsequence
%                     (CDF record). Unit: Hz.
%
%
% RETURN VALUE
% ============
% newTt2000 : Nx1 vector. Like oldTt2000 but each single time (row) has been
%             replaced by a constantly incrementing sequence of times (rows).
%             Every such sequence begins with the original value, has length
%             nSpr with frequency freqWithinRecords(i).
%             NOTE: There is no check that the entire sequence is monotonic. LFR
%             data can have snapshots (i.e. snapshot records) that overlap in
%             time!
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-28, based on older code.
%
function newTt2000 = convert_N_to_1_SPR_Epoch( oldTt2000, nSpr, freqHzWithinRecords )

% PROPOSAL: Turn into more generic function, working on number sequences in general.
%   CON: Much TT2000-special code.
% PROPOSAL: N_sequence should be a column vector.
%    NOTE: TDS-LFM-RSWF, LFR-SURV-CWF have varying snapshot lengths.
%    PRO: Could be useful for converting N->1 samples/record for calibration with transfer functions.
%       NOTE: Then also needs function that does the reverse.
% PROPOSAL: Replace by some simpler(?) algorithm that uses column/matrix multiplication.
%
% PROPOSAL: Rename Epoch-->TT2000.
% PROPOSAL: Rename Epoch

% ASSERTIONS
if ~iscolumn(oldTt2000)
  error(...
    'convert_N_to_1_SPR_Epoch:Assertion:IllegalArgument', ...
    'Argument is not a column vector')   % Right ID?
elseif ~isa(oldTt2000, 'int64')
  error(...
    'convert_N_to_1_SPR_Epoch:Assertion:IllegalArgument', ...
    'Argument has the wrong class.')   % Right ID?
end
if numel(nSpr) ~= 1
  error(...
    'convert_N_to_1_SPR_Epoch:Assertion:IllegalArgument', ...
    'nSpr not scalar.')
elseif size(freqHzWithinRecords, 1) ~= size(oldTt2000, 1)
  error(...
    'convert_N_to_1_SPR_Epoch:Assertion:IllegalArgument', ...
    'freqWithinRecords and oldTt2000 do not have the same number of rows.')
end
assert(iscolumn(freqHzWithinRecords))

nRecords = numel(oldTt2000);

% Express frequency as period length in ns (since tt2000 uses ns as a unit).
% Use the same MATLAB class as tt.
% Unique frequency per record.
periodNsColVec = int64(1e9 ./ freqHzWithinRecords);
periodNsMatrix = repmat(periodNsColVec, [1, nSpr]);

% Conventions:
% ------------
% Time unit: ns (as for tt2000)
% Algorithm should require integers to have a very predictable behaviour
% (useful when testing).

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
newTt2000 = solo.hwzv.convert_N_to_1_SPR_redistribute(tt2000Matrix);
end
