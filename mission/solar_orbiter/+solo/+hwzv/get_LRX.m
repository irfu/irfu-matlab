%
% Return the relevant value of LFR CDF zVariables R0, R1, or R2, or a
% hypothetical but analogous "R3" which is always 1.
%
% The abbreviation LRX is defined in BICAS
% (irfu-matlab/mission/solar_orbiter/bicas/readme.txt).
%
%
% ARGUMENTS
% =========
% zvR0, zvR1, zvR2, iLsf
%       LFR CDF zVariable-like variables. All must have identical array sizes.
%
%
% RETURN VALUE
% ============
% zvLrx
%       Same size array as arguments. The relevant values are copied,
%       respectively, from R0, R1, R2, or an analogous hypothetical "R3" that is
%       a constant (=1), depending on the value of iLsf in the corresponding
%       component.
%       iLsf(i) == 1 ==> Rx(i) == zvR0(i),
%       iLsf(i) == 2 ==> Rx(i) == zvR1(i), and so on.
%       NOTE: Numeric (like R0, R1, R2). Not MATLAB class "logical".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created earliest 2016-10-10 and latest 2019.
%
function zvLrx = get_LRX(zvR0, zvR1, zvR2, iLsf)
% TODO-DEC: Should convey iLsf=NaN as NaN, or assert that iLsf ~= NaN?

irf.assert.sizes(...
  zvR0, [-1, 1], ...
  zvR1, [-1, 1], ...
  zvR2, [-1, 1], ...
  iLsf, [-1, 1]);

% Set to NaN (should always be overwritten if code works) and iLsf has
% correct values.
zvLrx = nan(size(iLsf));

b = (iLsf==1);   zvLrx(b) = zvR0(b);
b = (iLsf==2);   zvLrx(b) = zvR1(b);
b = (iLsf==3);   zvLrx(b) = zvR2(b);
b = (iLsf==4);   zvLrx(b) = 1;
% Last value for (iLsf==4) is the value of a hypothetical (non-existent,
% constant) analogous zVariable "R3".

% NOTE: Prevents that iLsf=NaN ==> NaN. Desirable?
assert(all(~isnan(zvLrx)), ...
  ['Likely that argument iLsf has illegal values (not an integer', ...
  ' 1-4). Are you using zVar FREQ (integers 0-3)?'])
end
