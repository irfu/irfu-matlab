%
% Given a DATASET_ID, return the constituent parts.
%
%
% NOTES
% =====
% """"""""
% 7.1.2 ROC dataset identifier naming convention
% Each RPW data set must be identified with a unique uppercase string of the
% following form:
% [Source_name]_[Level]_[Descriptor]
% """"""""
% /"ROC Data Products", ROC-PRO-DAT-NTT-00006-LES, 01/02 draft
%
%
% ARGUMENTS
% =========
% datasetId
%       Any DATASET_ID that could possibly be legally defined, including ones
%       unrelated to BIAS and for instruments other than RPW.
%
%
% RETURN VALUES
% =============
% sourceName
%       String. "ROC-SGSE" or "SOLO".
% level
%       String. NOTE: Includes "L" in e.g. "L2".
% descriptor
%       String.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-29.
%
function [sourceName, level, descriptor] = disassemble_DATASET_ID(datasetId)
% PROPOSAL: Automatic test code.
% PROPOSAL: Rename to use other word than "disassemble".
%   PRO: Unconventional.
%   PROPOSAL:
%     ~parse
%     ~interpret

subStrCa = irf.str.regexp_str_parts(...
  datasetId, { ...
  '(ROC-SGSE|SOLO)', '_', ...
  '[^_]*',   '_', ...
  '[A-Z0-2-]*'}, ...
  'assert match');

sourceName = subStrCa{1};
level      = subStrCa{3};
descriptor = subStrCa{5};

solo.adm.assert_dataset_level(level)
end
