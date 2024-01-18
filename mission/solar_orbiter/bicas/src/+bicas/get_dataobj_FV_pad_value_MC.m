%
%
% ARGUMENTS
% =========
% Do
%       dataobj object.
%
%
% RETURN VALUES
% =============
% fv
%       Empty if there is no fill value.
% padValue
% mc
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [fv, padValue, mc] = get_dataobj_FV_pad_value_MC(Do, zvName)
    % NOTE: Uncertain how it handles the absence of a fill value. (Or is fill value mandatory?)
    % PROPOSAL: Remake into general-purpose function.
    % PROPOSAL: Remake into just using the do.Variables array?
    %    NOTE: Has to derive CDF variable type from do.Variables too.
    % PROPOSAL: Incorporate into dataobj?! Isn't there a buggy function/method there already?
    %
    % PROPOSAL: Use irfu-matlab's getfillval() instead.
    %   CON: This function also returns pad value. There is no known getpadval().
    %   CON: Can not handle TT2000 (?)

    % ===========================
    % Fill value and MATLAB class
    % ===========================

    % Obtain tentative values
    % NOTE: Special function for dataobj.
    fv = getfillval(Do, zvName);
    mc = Do.data.(zvName).type;

    % NOTE: For unknown reasons, the fill value for tt2000 zVariables (or at
    %       least "Epoch") is stored as a UTC(?) string.
    % NOTE: dataobj (probably) has a special case for "type" for "Epoch" or
    %       TT2000.
    if strcmp(mc, 'tt2000')
        % NOTE: Uncertain if this is the correct conversion function.
        fv = spdfparsett2000(fv);
        mc = 'int64';
    end

    % =========
    % Pad value
    % =========
    iZv      = strcmp(Do.Variables(:,1), zvName);
    padValue = Do.Variables{iZv, 9};
    % Comments in "spdfcdfinfo.m" should indirectly imply that column 9 is pad
    % values since the structure/array commented on should be identical.
end
