%
% Automated test code for interpret_config_file.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-25
%
function regexp_str_parts_TEST()

    KEY_VALUE_REGEXP_LIST         = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *'};
    KEY_VALUE_COMMENT_REGEXP_LIST = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *', '(#.*)?'};

    args = {};
    exp = {};
     args{end+1} = {'key',  {'[a-z]+'}};
     exp{end+1}  = {'key'}';
     args{end+1} = {'key=value',  {'[a-z]+', '=', '.*'}};
     exp{end+1}  = {'key', '=', 'value'}';
     
     args{end+1} = {'key = "value"', KEY_VALUE_REGEXP_LIST};
     exp{end+1}  = {'key', ' = ', '"', 'value', '"', ''}';
     args{end+1} = {'key=   "value"   ', KEY_VALUE_REGEXP_LIST};
     exp{end+1}  = {'key', '=   ', '"', 'value', '"', '   '}';
     args{end+1} = {'key_name.subset   =" .-/ "   ', KEY_VALUE_REGEXP_LIST};
     exp{end+1}  = {'key_name.subset', '   =', '"', ' .-/ ', '"', '   '}';
     
    args{end+1} = {'key = "value"', KEY_VALUE_COMMENT_REGEXP_LIST};
    exp{end+1}  = {'key', ' = ', '"', 'value', '"', '', ''}';
    args{end+1} = {'key = "value"   # Comment', KEY_VALUE_COMMENT_REGEXP_LIST};
    exp{end+1}  = {'key', ' = ', '"', 'value', '"', '   ', '# Comment'}';
    
    %args{end+1} = {'ab', {'a', 'b', 'c'}};    % Generates function error, as expected.
    %exp{end+1}  = {'a', 'b', ''}';
    %args{end+1} = {'ac', {'a', 'b', 'c'}};    % Generates function error, as expected.
    %exp{end+1}  = {'a', '', 'c'}';
    %args{end+1} = {'a', {'a', 'b', 'c'}};    % Generates function error, as expected.
    %exp{end+1}  = {'a', '', ''}';
     
    %args{end+1} = {'ac', {'a', 'b'}};    % Generates function error, as expected.
    %exp{end+1}  = {'a'}';
    
    %args{end+1} = {'key = value', KEY_VALUE_REGEXP_LIST};  % Generates function error, as expected.
    %exp{end+1}  = {'key', ' = ', '', 'value', '', ''}';
    
    
    for k = 1:length(args)
        res = bicas.utils.regexp_str_parts(args{k}{:});
        if ~isequaln(res, exp{k})
            k
            args{k}{:}
            exp{k}
            res
            error('FAIL')
        end
        disp('TEST OK')
    end
end
