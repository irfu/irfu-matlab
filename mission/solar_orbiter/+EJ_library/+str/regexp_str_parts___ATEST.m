%
% Automated test code for interpret_config_file.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-25
%
function regexp_str_parts___ATEST()

    KEY_ANY_VALUE_REGEXP_LIST            = {'[a-z]+', '=', '.*'};   % Not whitespace before or after =.
    KEY_QUOTED_VALUE_REGEXP_LIST         = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *'};    % Require quoted value.
    KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST = {'[a-zA-Z0-9._]+', ' *= *', '"', '[^"]*', '"', ' *', '(#.*)?'};

    ES = char(zeros(1,0));    % Emty string, size 1x0.

    new_test   = @(inputs, outputs)        (EJ_library.atest.CompareFuncResult(@EJ_library.str.regexp_str_parts, inputs, outputs));
    new_test_A = @(arg1, arg2, outputs) (new_test({arg1, arg2, 'assert match'    }, outputs));   % A = Assert match
    new_test_P = @(arg1, arg2, outputs) (new_test({arg1, arg2, 'permit non-match'}, outputs));   % P = Permit non-match
    
    
    
    tl = {};

    tl{end+1} = new_test_A('',     {''},       {{''}',    ES, true});

    % Exact match
    % -----------
    % Non-empty regexp
    tl{end+1} = new_test_A('abc', {'a', 'b', 'c'}, {{'a' 'b' 'c'}', ES, true});
    tl{end+1} = new_test_P('abc', {'a', 'b', 'c'}, {{'a' 'b' 'c'}', ES, true});
    % Include empty regexp.
    tl{end+1} = new_test_A('abc', {'a', 'b', '', 'c'}, {{'a' 'b' '' 'c'}', ES, true});
    tl{end+1} = new_test_P('abc', {'a', 'b', '', 'c'}, {{'a' 'b' '' 'c'}', ES, true});
    
    % Matching, until running out of string
    % -------------------------------------
    tl{end+1} = new_test_A('ab', {'a', 'b', 'c'}, 'MException');
    tl{end+1} = new_test_P('ab', {'a', 'b', 'c'}, {{'a', 'b'}', ES, false});
    
    % First regexp matches, but second does not
    % -----------------------------------------
    % Non-empty remaining string.
    tl{end+1} = new_test_A('ac', {'a', 'b', 'c'}, 'MException');          
    tl{end+1} = new_test_P('ac', {'a', 'b', 'c'}, {{'a'}, 'c', false});
    tl{end+1} = new_test_A('ac', {'a', 'b'     }, 'MException');       
    tl{end+1} = new_test_P('ac', {'a', 'b'     }, {{'a'}, 'c', false});
    % Running out of string
    tl{end+1} = new_test_A('a',  {'a', 'b', 'c'}, 'MException');   
    tl{end+1} = new_test_P('a',  {'a', 'b', 'c'}, {{'a'}, ES,  false});



    tl{end+1} = new_test_A('word',  {'[a-z]+'}, {{'word'}', ES, true});
    tl{end+1} = new_test_A('key = value', KEY_QUOTED_VALUE_REGEXP_LIST, 'MException');

    tl{end+1} = new_test_A('key=value',                     KEY_ANY_VALUE_REGEXP_LIST,            {{'key',                '=',         'value'            }', ES, true});
    tl{end+1} = new_test_A('key = "value"',                 KEY_QUOTED_VALUE_REGEXP_LIST,         {{'key',               ' = ',   '"', 'value', '"', ''   }', ES, true});
    tl{end+1} = new_test_A('key=   "value"   ',             KEY_QUOTED_VALUE_REGEXP_LIST,         {{'key',                '=   ', '"', 'value', '"', '   '}', ES, true});
    tl{end+1} = new_test_A('key_name.subset   =" .-/ "   ', KEY_QUOTED_VALUE_REGEXP_LIST,         {{'key_name.subset', '   =',    '"', ' .-/ ', '"', '   '}', ES, true});
    tl{end+1} = new_test_A('key = "value"',                 KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {{'key',               ' = ',   '"', 'value', '"', '', ''            }', ES, true});
    
    % Match
    tl{end+1} = new_test_A('key = "value"   # Comment',     KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {{'key',               ' = ',   '"', 'value', '"', '   ', '# Comment'}', ES, true});
    tl{end+1} = new_test_P('key = "value"   # Comment',     KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {{'key',               ' = ',   '"', 'value', '"', '   ', '# Comment'}', ES, true});
    % Non-match
    tl{end+1} = new_test_A('key = "value"   % Comment',     KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, 'MException');
    tl{end+1} = new_test_P('key = "value"   % Comment',     KEY_QUOTED_VALUE_COMMENT_REGEXP_LIST, {{'key',               ' = ',   '"', 'value', '"', '   ', ''}', '% Comment', false});



    EJ_library.atest.run_tests(tl)    
end
