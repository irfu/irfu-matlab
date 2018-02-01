%
% Automated test code fore interpret_config_file.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-25
%
function interpret_config_file_TEST()

    args = {};
    exp = {};
    %args{end+1} = {{'# Comment'}};
    %exp{end+1}  = containers.Map('KeyType', 'char', 'ValueType', 'char');    
    %args{end+1} = {{'key="value"'}};    
    %exp{end+1}  = containers.Map({'key'}, {'value'});
    %args{end+1} = {{'key="value"   # Comment'}};    
    %exp{end+1}  = containers.Map({'key'}, {'value'});

    args{end+1} = {{...
        '# Comment', ...
        '', ...
        '   ', ...
        'key.1="value1"', ...
        'key_2="value2"   # Comment', ...
        'key-3  =   ""' ...
        }};
    exp{end+1}  = containers.Map({'key.1', 'key_2', 'key-3'}, {'value1', 'value2', ''});
    
    for k = 1:length(args)
        res = bicas.interpret_config_file(args{k}{:});
        if ~isequaln(res, exp{k})
            args{k}{:}
            exp{k}
            res
            error('FAIL')
        end
        disp('TEST OK')
    end
end
