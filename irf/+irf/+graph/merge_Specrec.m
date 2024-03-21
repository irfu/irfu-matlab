%
% Take multiple Specrec structs as returned from irf_powerfft() (irfu-matlab git
% repo) and merge them into one Specrec struct. The input structs may have
% different sets of frequencies, but should not overlap in time.
%
% NOTE: THIS FUNCTION IS A BIT OF A HACK to make it possible to apply
% irf_powerfft() to data with (wildly) changing sampling frequencies, yet call
% irf_spectrogram() only once.
% The intended use is to
% (1) Split the underlying data (samples) into time intervals with one
%     approximately constant sampling frequency per time interval.
% (2) For each time interval, call irf_powerfft().
% (3) Merge the resulting Specrec structs using this function.
% (4) Call irf_spectrogram() ONCE using the result from step (3).
%
% The resulting Specrec struct describes a larger spectrum, with potentially
% large parts legitimately set to NaN. It contains the set union of the source
% structs' frequencies and timestamps.
%
% One advantage with this is to reduce processing and speed up(?) scripts.
%
% NOTE: Memory use could be a potential problem since internal data size should
% be on the same order as data set zVars. Has not yet observed that to be a
% problem though. /2020-08-16
%
% NOTE: Currently only supports SpecrecCa{i}.p with one cell array component.
%
%
% ARGUMENTS
% =========
% SpecrecCa
%       1D cell array of "Specrec" structs as returned by irf_powerfft(), except
%       that it may also have had .dt added to it.
%       NOTE: Must only contain one sample per timestamp.
%
%
% RETURN VALUES
% =============
% Specrec
%       "Specrec" struct. Consists of the merger of structs in SpecrecCa.
%       Specrec.p{1}(i,j)==NaN for values/indices not assigned by any argument
%       to this function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-08-14.
%
function Specrec = merge_Specrec(SpecrecCa)
% PROPOSAL: Automatic test code.
%
% INCOMPLETE: Does not complement NaN between non-NaN values in the frequency coordinate.
%   PROPOSAL: Generic function x,y-->x,y such that NaN values are replaced by nearest value, if not farther away
%   than threeshold.
%       CON: When function is applied, coordinates have already been merged, and it is hard to calculate the max
%            nearest distance for replacing NaN to use.
%
% PROPOSAL: Move to irfu-matlab's solo.sp?!!
%
% PROPOSAL: Have time_pad_NaN() handle having only one timestamp? Zero (timestamps?)
%
% PROPOSAL: Not add Nan values in time (use time_pad_NaN()). Set .dt instead. -- DONE
%   PRO: Intended use of .dt.
%   CON: Some Specrec already have .dt values.
%       PROPOSAL: Keep without modifying.
%
% PROPOSAL: assert_nonempty_Specrec() as separate function.
%
% PRESUMED BUG:
%   F0-F2, V1_DC, 00:00-05:30. Drawn out, "homogenous" spectrum.
%   solo_L3_rpw-lfr-surv-swf-e_20210102_V04.png  /modif 2021-03-06
%   One snapshot/spectrum at beginning of day which is drawn out? -- FIXED
%   PROPOSAL: Better way of deducing time between snapshots.
%
% PROPOSAL: Permit merging to one empty Specrec.
%   TODO-DEC: Should that be represented as []?
%   TODO-DEC: If represented as a struct, should it have .dt?



assert(iscell(SpecrecCa), 'SpecrecCa is not a cell array.')

% Remove empty Specrecs, which indeed are legal return results from
% irf_powerfft().
bEmpty    = cellfun(@isempty, SpecrecCa);
SpecrecCa = SpecrecCa(~bEmpty);

N = numel(SpecrecCa);

assert(numel(SpecrecCa) >= 1, ...
  'SpecrecCa is empty (after removing individually empty Specrecs).')

for iS = 1:N
  S = SpecrecCa{iS};

  if iS == 1
    hasDt = isfield(S, 'dt');
  end

  % ASSERTION
  assert_nonempty_Specrec(S, hasDt)

  %         if N >= 2
  %             % Pad with extra NaN samples in time.
  %             %S = time_pad_NaN(S);
  %         end

  if iS == 1
    % Create Specrec "STot" which will grow as content from other
    % Specrecs (S) are added.
    % IMPLEMENTATION NOTE: Must create this AFTER the first Specrec has
    % been padded.
    STot = S;
  else
    [tf, pArray] = irf.ds.merge_coordinated_arrays(...
      NaN, ...
      {STot.t,  STot.f},  STot.p{1}, ...
      {S.t, S.f}, S.p{1});

    STot.t = tf{1};
    STot.f = tf{2};
    STot.p = {pArray};
    if hasDt
      STot.dt = [STot.dt; S.dt(:)];
    end
  end

end

% Fill in samples=NaN (use nearest value) placed between the original
% non-NaN samples over frequencies (same timestamp/spectrum)).
%
% TODO-NI: Why is this needed?!
for iTime = 1:numel(STot.t)

  STot.p{1}(iTime, :) = use_nearest_nonNaN(STot.f, STot.p{1}(iTime, :));
end



% ASSERTION
assert_nonempty_Specrec(STot, hasDt)

Specrec = STot;
end



% NOTE: Does not permit S == [].
% NOTE: Only permits scalar S.p ("nonempty").
function assert_nonempty_Specrec(S, hasDt)
fnCa = {'t', 'p', 'f'};
if hasDt
  fnCa{end+1} = 'dt';
end

irf.assert.struct(S, fnCa, {})
assert(isscalar(S.p))
assert(iscell(S.p))
irf.assert.sizes(...
  S.t,    [-1], ...
  S.p{1}, [-1, -2], ...
  S.f,    [-2]);

if hasDt
  irf.assert.sizes(...
    S.t,    [-1], ...
    S.dt,   [-1])
end

assert(issorted(S.t, 'strictascend'), 'S.t is not monotonically increasing.')
assert(issorted(S.f, 'strictascend'), 'S.f is not monotonically increasing.')
end



% Replace NaN values with nearest non-NaN value, unless outside range of non-NaN
% values.
%
% ARGUMENTS
% =========
% x, y : Same-sized 1D arrays.
%
function y = use_nearest_nonNaN(x,y)
bFinite = ~isnan(y);
x2 = x(bFinite);
y2 = y(bFinite);
if any(bFinite)
  y = interp1(x2, y2, x, 'nearest');   % No extrapolation.
end
end



% Add spectrum samples NaN at timestamps before first time stamp, and after last
% timestamp.
% function Specrec = time_pad_NaN(Specrec)
%
%     t = Specrec.t;
%
%     % NOTE: Requires numel(t) >= 2.
%     assert(numel(t) >= 2, ...
%         'Specrec.t contains less than two timestamps. Can not handle this case.')
%
%     C = 0.01;   % Good choice? Use 1?
%     t1 = t(1)   - (t(  2) - t(    1))*C;
%     t2 = t(end) + (t(end) - t(end-1))*C;
%     tPadded = [t1; t; t2];
%
%     Specrec.t = tPadded;
%
%     Specrec.p{1} = padarray(Specrec.p{1}, [1, 0], NaN, 'both');
%     % NOTE: .f unchanged
% end
