%
% Return the relevant value of LFR CDF zVariables R0, R1, or R2, or a
% hypothetical but analogous "R3" which is always 1.
%
%
% ARGUMENTS
% =========
% R0, R1, R2, iLsf : LFR CDF zVariable-like. All must have identical array sizes.
%                    iLsf(i) == 1 (not 0) ==> F0 and so on.
%
%
% RETURN VALUE
% ============
% Rx : Same size array as arguments. The relevant values are copied,
%      respectively, from R0, R1, R2, or an analogous hypothetical "R3" that is
%      a constant (=1) depending on the value of iLsf in the corresponding
%      component.
%      NOTE: Numeric (like R0, R1, R2). Not MATLAB class "logical".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created earliest 2016-10-10 and latest 2019.
%
function Rx = get_LFR_Rx(R0, R1, R2, iLsf)
    % TODO-DEC: Should convey iLsf=NaN as NaN, or assert that iLsf ~= NaN?
    
    EJ_library.assert.sizes(...
        R0,   [-1, 1], ...
        R1,   [-1, 1], ...
        R2,   [-1, 1], ...
        iLsf, [-1, 1]);
    
    % Set to NaN (should always be overwritten if code works) and iLsf has
    % correct values.
    Rx = nan(size(iLsf));
    
    b = (iLsf==1);   Rx(b) = R0(b);
    b = (iLsf==2);   Rx(b) = R1(b);
    b = (iLsf==3);   Rx(b) = R2(b);
    b = (iLsf==4);   Rx(b) = 1;
    % Last one is the value of a hypothetical (non-existant, constant) analogous
    % zVariable "R3".
    
    % NOTE: Prevents that iLsf=NaN ==> NaN. Desirable?
    assert(all(~isnan(Rx)), 'Likely that argument iLsf has illegal values (not an integer 1-4). Are you using zVar FREQ (integers 0-3)?')
end
