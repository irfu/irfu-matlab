%
% Convert a settings value (from bicas.Settings) to  one-row string that can be
% displayed, e.g. in a log message.
%
% NOTE: This function should ideally be able to handle all settings values, but
% not necessarily constrain them, i.e. it may handle a larger range of values.
%   Ex: Recursive array or cell arrays, 2D arrays (not impl.), varying types in
%       same cell array.
%
%
% ARGUMENTS
% =========
% value
%       A settings value from an instance of bicas.Settings, i.e. the value in a
%       key-value pair (one of the versions of the value).
%
%
% RETURN VALUES
% =============
% displayStr
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-12, by breaking out code from bicas.sprint_SETTINGS().
%
function displayStr = settings_value_to_display_str(value)
    % PROPOSAL: Shorter name.
    
    if ischar(value)

        displayStr = ['"', value, '"'];

    elseif islogical(value)

        assert(isscalar(value))
        if value
            displayStr = 'true';
        else
            displayStr = 'false';
        end

    elseif isnumeric(value)

        assert(isvector(value))
        
        if isscalar(value)
            displayStr = sprintf('%g', value);
        else
            % RECURSIVE CALL
            displayStr = sprintf('[%s]', many_display_str(num2cell(value)));
        end

    elseif iscell(value)

        assert(isvector(value))
        
        % RECURSIVE CALL
        displayStr = sprintf('{%s}', many_display_str(value));
    else

        error(...
            'BICAS:Assertion', ...
            ['can not convert SETTINGS value (overriden or not)'])
    end
end



function displayStr = many_display_str(ca)
    % PROPOSAL: Better name.
    
    assert(isvector(ca))    
    
    displayStrCa = cell(numel(ca), 1);
    for i = 1:numel(ca)

        % RECURSIVE CALL
        displayStrCa{i} = bicas.settings_value_to_display_str(ca{i});
    end
    
    displayStr = strjoin(displayStrCa, ', ');
end
