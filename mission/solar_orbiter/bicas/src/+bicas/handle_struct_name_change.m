%
% To be used to analyze the result of use of
% irf.ds.normalize_struct_fieldnames(). React depending on BSO.
%
%
% ARGUMENTS
% =========
% fnChangeList
%       Returned from irf.ds.normalize_struct_fieldnames
% msgFunc
%       Function handle: msgStr = func(oldFieldname, newFieldname)
%       NOTE: Return value is passed to bicas.Logger.log (not logf), i.e.
%       multi-row messages must end with line feed.
% varargin
%       List of pairs of arguments.
%       varargin{2*m + 1} : Fieldname (new/after change) for which to react.
%       varargin{2*m + 2} : BSO key which determines the policy. Must have
%                           value WARNING or ERROR.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-09
%
function handle_struct_name_change(...
  fnChangeList, Bso, L, anomalyDescrMsgFunc, varargin)
%
% PROPOSAL: Somehow generalize to something that can handle "all" fnChangeList returned from
% normalize_struct_fieldnames.
%   PROPOSAL: Submit function that returns error/warning/log message. Accepts arguments for old+new
%             fieldname. Can thus handle e.g. both zVar names and global attributes.
%
% PROPOSAL: Argument for error ID to pass to bicas.default_anomaly_handling
%
% PROPOSAL: New name more analogous to irf.ds.normalize_struct_fieldnames (which function is to be used
%           with).
%   PROPOSAL: handle_struct_normalization_anomaly
%
% PROBLEM: Concept of changing fieldnames of actual struct is bad when struct is so large that it may cause memory
%          problems. May want to avoid modifying struct in order to help MATLAB's code optimization (prevent
%          temporary copies).

while numel(varargin) >= 2    % Iterate over pairs of varargin components.
  newFn      = varargin{1};
  settingKey = varargin{2};
  varargin   = varargin(3:end);

  % NOTE: i==0 <==> no match.
  [~, i] = ismember(newFn, {fnChangeList(:).newFieldname});
  if i > 0
    % CASE: Found a fieldname change to react to.
    [settingValue, settingKey] = Bso.get_fv(settingKey);
    anomalyDescrMsg = anomalyDescrMsgFunc(...
      fnChangeList(i).oldFieldname, ...
      fnChangeList(i).newFieldname);

    assert(...
      isempty(fnChangeList(i).ignoredCandidateFieldnames), ...
      ['Function not designed for handling non-empty', ...
      ' .ignoredCandidateFieldnames.'])

    bicas.default_anomaly_handling(L, ...
      settingValue, settingKey, ...
      'ERROR_WARNING_ILLEGAL_SETTING', anomalyDescrMsg, 'BICAS:Assertion')
  end
end

assert(numel(varargin) == 0)
end
