% Test code for bicas.utils.parse_CLI_flags
function parse_CLI_flags_TEST

FlagsConfigMap = containers.Map(...
        {'a', 'b', 'c'}, ...
        {...
            struct('cliFlagString', '-a', 'occurranceRequirement', '0-1', 'nValues', 0), ...
            struct('cliFlagString', '-b', 'occurranceRequirement', '1',   'nValues', 1), ...
            struct('cliFlagString', '-c', 'occurranceRequirement', '0-N', 'nValues', 2)...
        });

    
    
exp =  {};
args = {};

args{end+1} = {strsplit('-b 123'), FlagsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{}, {{'123'}}, {}});

args{end+1} = {strsplit('-a -b 123'), FlagsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{cell(1,0)}, {{'123'}}, {}});

args{end+1} = {strsplit('-a -b 123 -c 8 9'), FlagsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{cell(1,0)}, {{'123'}}, {{'8', '9'}}});

args{end+1} = {strsplit('-c 6 7 -a -b 123 -c 8 9'), FlagsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{cell(1,0)}, {{'123'}}, {{'6', '7'}, {'8', '9'}}});

args{end+1} = {strsplit('-c 6 7 -b 123 -c 8 9'), FlagsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{}, {{'123'}}, {{'6', '7'}, {'8', '9'}}});


    
for k = 1:length(args)
    result = bicas.utils.parse_CLI_flags(args{k}{:});
    if ~isequaln( exp{k}, result )
        exp{k}
        error('FAIL')
    end
end

end