%
% Automatic test code for bicas.utils.parse_CLI_options.
%
function parse_CLI_options_TEST

OptionsConfigMap = containers.Map(...
        {'a', 'b', 'c'}, ...
        {...
            struct('optionHeader', '-a', 'occurrenceRequirement', '0-1',   'nValues', 0), ...
            struct('optionHeader', '-b', 'occurrenceRequirement', '1',     'nValues', 1), ...
            struct('optionHeader', '-c', 'occurrenceRequirement', '0-inf', 'nValues', 2)...
        });


    
exp =  {};
args = {};

args{end+1} = {strsplit('-b 123'), OptionsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{}, {{'123'}}, {}});

args{end+1} = {strsplit('-a -b 123'), OptionsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{cell(1,0)}, {{'123'}}, {}});

args{end+1} = {strsplit('-a -b 123 -c 8 9'), OptionsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{cell(1,0)}, {{'123'}}, {{'8', '9'}}});

args{end+1} = {strsplit('-c 6 7 -a -b 123 -c 8 9'), OptionsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{cell(1,0)}, {{'123'}}, {{'6', '7'}, {'8', '9'}}});

args{end+1} = {strsplit('-c 6 7 -b 123 -c 8 9'), OptionsConfigMap};
exp{end+1}  = containers.Map({'a', 'b', 'c'}, {{}, {{'123'}}, {{'6', '7'}, {'8', '9'}}});


    
for k = 1:length(args)
    result = bicas.utils.parse_CLI_options(args{k}{:});
    if ~isequaln( exp{k}, result )
        exp{k}
        error('FAIL')
    end
    disp('TEST OK')
end

end
