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
% fillValue
%       Empty if there is no fill value.
% padValue
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function [fillValue, padValue] = get_fill_pad_values(Do, zvName)
    % NOTE: Uncertain how it handles the absence of a fill value. (Or is fill value mandatory?)
    % PROPOSAL: Remake into general-purpose function.
    % PROPOSAL: Remake into just using the do.Variables array?
    %    NOTE: Has to derive CDF variable type from do.Variables too.
    % PROPOSAL: Incorporate into dataobj?! Isn't there a buggy function/method there already?
    %
    % PROPOSAL: Use irfu-matlab's getfillval() instead.
    %   CON: This function also returns pad value. There is no known getpadval().
    %   CON: Can not handle TT2000 (?)
    
    % NOTE: Special function for dataobj.
    fillValue = getfillval(Do, zvName);
    % NOTE: For unknown reasons, the fill value for tt2000 zVariables (or at
    % least "Epoch") is stored as a UTC(?) string.
    if strcmp(Do.data.(zvName).type, 'tt2000')
        % NOTE: Uncertain if this is the correct conversion function.
        fillValue = spdfparsett2000(fillValue);
    end
    
    iZVariable = strcmp(Do.Variables(:,1), zvName);
    padValue   = Do.Variables{iZVariable, 9};
    % Comments in "spdfcdfinfo.m" should indirectly imply that column 9 is pad
    % values since the structure/array commented on should be identical.
end
