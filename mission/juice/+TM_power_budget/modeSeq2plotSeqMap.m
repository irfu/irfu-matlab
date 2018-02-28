%
% Function for converting a mode sequence to a list (map) of arrays which are suitable for plotting (in-situ/radio) mode
% sequences.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% ModeSeq         : 1D struct array with fields "id" (string) and "beginSec" (number). Intended to be either an in-situ
%                   or radio mode sequence as used in a call to "engine". "beginSeq" is assumed to be monotonically
%                   increasing.
% globalBeginSec, 
% globalEndSec    : Scalars. Specify interval which the returned arrays cover exactly.
% ModeMap         : containers.Map. Keys=ID strings, values = 1D arrays on the form the form
%                   [tBeginSec0, tEndSec0, NaN, tBeginSec1, tEndSec1, NaN, ...] and describing the intervals when
%                   the specified mode runs. Convenient for plotting when a mode runs (can plot as one line with holes
%                   in it).
% 
%
% Created 2018-02-28 by Erik Johansson, IRF Uppsala.
%
function ModeMap = modeSeq2plotSeqMap(ModeSeq, globalBeginSec, globalEndSec)
% NOTE: Questionable if efficient.
% IMPLEMENTATION NOTE: Can not specify Inf as lastEndSec since plotting does not handle Inf correctly (point is ignored,
% like for NaN). Can therefore not set it to that as a hardcoded constant.
% PROPOSAL: Change function name.


ModeMap = containers.Map;

for i = 1:numel(ModeSeq)
    id = ModeSeq(i).id;
    
    beginSec = ModeSeq(i).beginSec;
    if i < numel(ModeSeq)
        endSec = ModeSeq(i+1).beginSec;
    else
        endSec = max(beginSec, globalEndSec);
    end
    
    beginSec = max(beginSec, globalBeginSec);
    endSec   = min(endSec,   globalEndSec);
    if (endSec - beginSec <= 0)
        continue
    end
    
    if ~ModeMap.isKey(id)
        plotArray = [];
    else
        plotArray = ModeMap(id);
    end
    plotArray(end+1:end+3) = [beginSec, endSec, NaN];
    
    ModeMap(id) = plotArray;
end

end
