
%
%
% ARGUMENTS
% =========
% fnChangeList : Returned from EJ_library.utils.normalize_struct_fieldnames
% msgFunc  : Function handle: msgStr = func(oldFieldname, newFieldname)
%            NOTE: Return value is passed to bicas.logger.log (not logf), i.e. multi-row messages must end with line feed.
% varargin : List of pairs of arguments.
%            varargin{2*m + 1} : Fieldname (new/after change) for which to react.
%            varargin{2*m + 2} : SETTINGS key which determines the policy. Must have value PERMIT, WARNING, or
%                                ERROR.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
function handle_struct_name_change(fnChangeList, SETTINGS, L, msgFunc, varargin)
    % PROPOSAL: Somehow generalize to something that can handle "all" fnChangeList returned from
    % normalize_struct_fieldnames.
    %   PROPOSAL: Submit function that returns error/warning/log message. Accepts arguments for old+new
    %             fieldname. Can thus handle e.g. both zVar names and global attributes.
    
    while numel(varargin) >= 2    % Iterate over pairs of varargin components.
        newFn      = varargin{1};
        settingKey = varargin{2};
        varargin   = varargin(3:end);
        
        [~, i] = ismember(newFn, {fnChangeList(:).newFieldname});   % NOTE: i==0 <==> no match.
        if i > 0
            [settingValue, settingKey] = SETTINGS.get_fv(settingKey);
            anomalyDescrMsg = anomalyDescrMsgFunc(fnChangeList(i).oldFieldname, fnChangeList(i).newFieldname);
            
            assert(isempty(fnChangeList(i).ignoredCandidateFieldnames), ...
                'Function not designed for handling non-empty .ignoredCandidateFieldnames.')
            
            bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', anomalyDescrMsg, 'BICAS:Assertion')
        end
    end
    
    assert(numel(varargin) == 0)
end



