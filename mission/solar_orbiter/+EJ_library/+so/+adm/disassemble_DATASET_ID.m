%
% Given a DATASET_ID, return the constituent parts.
%
%
% NOTES
% =====
% """"""""
% 7.1.2 ROC dataset identifier naming convention
% Each RPW data set must be identified with a unique uppercase string of the following form:
% [Source_name]_[Level]_[Descriptor]
% """"""""
% /"ROC Data Products", ROC-PRO-DAT-NTT-00006-LES, 01/02 draft
%
%
% ARGUMENTS
% =========
% datasetId : Any RPW DATASET_ID that could possibly be legally defined, i.e.
%             including all RPW DATASET_IDs in use, including ones unrelated to
%             BIAS.
%
%
% RETURN VALUES
% =============
% sourceName : String.
% level      : String. NOTE: Includes "L" in e.g. "L2".
% descriptor : String.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-09-29.
%
function [sourceName, level, descriptor] = disassemble_DATASET_ID(datasetId)
    
    subStrList = EJ_library.str.regexp_str_parts(...
        datasetId, { ...
        '(ROC-SGSE|SOLO)', '_', ...
        '[^_]*',   '_', ...
        'RPW-[A-Z0-2-]*'}, ...
        'assert match');
    
    sourceName = subStrList{1};
    level      = subStrList{3};
    descriptor = subStrList{5};
    
    EJ_library.so.adm.assert_dataset_level(level)
end
