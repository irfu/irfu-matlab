% Standard handling of fnChangeList returned from bicas.utils.normalize_struct_fieldnames based on standardized
% SETTINGS values.
%
%
% ARGUMENTS
% =========
% varargin : List of pairs of arguments.
%            varargin{2*m + 1} : Fieldname (new/after change) for which to react.
%            varargin{2*m + 2} : SETTINGS key which determines the policy. Must have value PERMIT, WARNING, or
%                                ERROR.
% msgFunc  : Function handle: msgStr = func(oldFieldname, newFieldname)
%            NOTE: Return value is passed to bicas.log (not logf), i.e. multi-row messages must end with line feed.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
function handle_struct_name_change(fnChangeList, SETTINGS, msgFunc, varargin)
    % PROPOSAL: Somehow generalize to something that can handle "all" fnChangeList returned from
    % normalize_struct_fieldnames.
    %   PROPOSAL: Submit function that returns error/warning/log message. Accepts arguments for old+new
    %             fieldname. Can thus handle e.g. both zVar names and global attributes.
    
    while numel(varargin) >= 2
        newFn       = varargin{1};
        settingsKey = varargin{2};
        varargin    = varargin(3:end);
        
        [~, i] = ismember(newFn, {fnChangeList(:).newFieldname});   % NOTE: i==0 <==> no match.
        if i > 0
            settingZvNameChangePolicy = SETTINGS.get_fv(settingsKey);
            msg = msgFunc(fnChangeList(i).oldFieldname, fnChangeList(i).newFieldname);
            
            switch(settingZvNameChangePolicy)
                case 'PERMIT'
                    % Do nothing
                case 'WARNING'
                    bicas.log('warning', msg)
                    warning(msg)
                case 'ERROR'
                    error('BICAS:Assertion', msg)
                otherwise
                    error('BICAS:proc_sub:Assertion', 'Illegal setting for %s="%s"', settingsKey, settingZvNameChangePolicy)
            end
        end
    end
    
    assert(numel(varargin) == 0)
end



