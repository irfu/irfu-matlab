%
% Automated test code for interpret_config_file.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-25
%
function regexp_str_parts___ATEST()

    KEY_VALUE_REGEXP_LIST         = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *'};
    KEY_VALUE_COMMENT_REGEXP_LIST = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *', '(#.*)?'};
    
    new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.utils.regexp_str_parts, inputs, outputs));
    tl = {};
    
    tl{end+1} = new_test({'key',  {'[a-z]+'}}, {{'key'}'});
    tl{end+1} = new_test({'key=value',         {'[a-z]+', '=', '.*'}}, {{'key', '=', 'value'}'});
    tl{end+1} = new_test({'key = "value"',     KEY_VALUE_REGEXP_LIST}, {{'key', ' = ', '"', 'value', '"', ''}'});
    tl{end+1} = new_test({'key=   "value"   ', KEY_VALUE_REGEXP_LIST}, {{'key', '=   ', '"', 'value', '"', '   '}'});
    tl{end+1} = new_test({'key_name.subset   =" .-/ "   ', KEY_VALUE_REGEXP_LIST}, {{'key_name.subset', '   =', '"', ' .-/ ', '"', '   '}'});
    tl{end+1} = new_test({'key = "value"',                 KEY_VALUE_COMMENT_REGEXP_LIST}, {{'key', ' = ', '"', 'value', '"', '', ''}'});
    tl{end+1} = new_test({'key = "value"   # Comment',     KEY_VALUE_COMMENT_REGEXP_LIST}, {{'key', ' = ', '"', 'value', '"', '   ', '# Comment'}'});
    
    
    tl{end+1} = new_test({'ab', {'a', 'b', 'c'}}, 'MException');
    tl{end+1} = new_test({'ac', {'a', 'b', 'c'}}, 'MException');
    tl{end+1} = new_test({'a',  {'a', 'b', 'c'}}, 'MException');
    tl{end+1} = new_test({'ac', {'a', 'b'     }}, 'MException');
    tl{end+1} = new_test({'key = value', KEY_VALUE_REGEXP_LIST}, 'MException');
    
    
    EJ_library.atest.run_tests(tl)    
end
