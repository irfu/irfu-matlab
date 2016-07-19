% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
% Parses list of command-line arguments assuming it is a list of flags and pairs of flag+value.
% Function tries to give accurate user-friendly errors for non-compliant arguments and absence of
% required arguments.
%
% Flag = A predefined (hardcoded, more or less) string meant to match a single argument, e.g. "--version"
% Value = An argument that specifies a more or less arbitrary value and comes after a flag, e.g. "--file <value>".
%
% flags : containers.Map, number/string-->struct
%    <keys>   : Arbitrary unique strings to identify the flags in the return value.
%    <values> : Information about each flag (syntax).
%    .CLI_str        : The command-line flag string (e.g. "--version"), including any prefix (e.g. dash).
%    .is_required    : Whether the flag (and any specified value) is required as opposed to optional.
%    .expects_value  : Whether the flag expects the following argument to be a value connected to the flag.
%
% flag_results : containers.Map, number/string-->string/number
%    <keys>
%    <values>  : For flag without value: 0=not set, 1=set
%                For flag with value:
%                   ~ischar(..) ==> the flag was not set (no such argument). (In reality numeric zero.)
%                    ischar(..) ==> the flag was set (incl. empty string). The string is the value.
%                NOTE: Can not use/implement isempty(..) as criterion since empty strings can be
%                legitimate argument values. Note that isempty('')==0.
%
function flag_results = parse_CLI_flags(arguments, flags)
%
% IMPLEMENTATION NOTE: Reasons for using containers.Map (instead of arrays, array of structs, cell
% array, structure of structures).
% 1) Caller can easily build up list of flags by amending list.
% 2) Can use key strings to identify flags (rather than the CLI flag strings themselves)
% 3) Can use more key strings (more characters) than for structure field names.
% 4) Easy for caller to group subsets of flags by keeping track of sets of keys. The caller can
% merge groups of flags (before submitting as one parameter), and can split the returned result into
% the groups of flags.
% Ex: Input files as opposed to output directory, log directory, config file.
% NOTE: The caller can easily(?) convert result into struct (one field per key-value pair).
%
% PROPOSAL: Change name: parse_arguments? parse_CLI_arguments?
%
% QUESTION: The function permits flags without value. Is this functionality really needed?
%

global ERROR_CODES



% ASSERTION CHECKS
if ~iscell(arguments)
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Parameter is not a cell array.');
elseif ~isa(flags, 'containers.Map')
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Parameter is not a containers.Map.');
end



CLI_strs = {};
flag_results = containers.Map;
for skey = flags.keys
    s = flags(skey{1});
    CLI_strs{end+1} = s.CLI_str;
    
    % Create return structure. Default: No flags found
    % NOTE: Applies to both flags with and without values!
    flag_results(skey{1}) = 0;         
end



% ASSERTION CHECKS: Check that there are no flag duplicates.
if length(CLI_strs) ~= length(unique(CLI_strs))
    errorp(ERROR_CODES.ASSERTION_ERROR, ...
        'The code is configured to accept multiple IDENTICAL command-line flags. This indicates a bug.')
end



ia = 1;
while ia <= length(arguments)    % ia = i_argument
    arg = arguments{ia};
    
    % Find matching flag (and key)
    flag = [];
    for skey = flags.keys
        key = skey{1};
        current_flag = flags(key);
        if strcmp(arg, current_flag.CLI_str)
            flag = current_flag;
            break
        end
    end    
    if isempty(flag)
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret command-line argument "%s". There is no such flag.', arg)
    end
    
    if ~(isnumeric(flag_results(key)) && (flag_results(key) == 0))
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'The command-line flag "%s" was specified (at least) twice.', arg)
    end
    
    if flag.expects_value
        if ia >= length(arguments)
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, ...
                'Can not find the argument that is expected to follow command-line flag "%s".', arg)
        end
        ia = ia + 1;
        flag_results(key) = arguments{ia};
    else
        flag_results(key) = 1;
    end
    
    ia = ia + 1;
end   % while



% Check that all required flags were set.
for skey = flags.keys
    flag = flags(skey{1});
    flag_result = flag_results(skey{1});
    if flag.is_required && (isnumeric(flag_result) && (flag_result == 0))
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Could not find required command-line flag "%s".', flag.CLI_str)
    end
end



% Convert flags (array) to structure.
% -----------------------------------
% IMPLEMENTATION NOTE: Up until here the algorithm works with the function argument describing the
% flags, and amends the results to it. First here does it copy the data to a variable that will be
% returned. This function has been modified several times to modify the output data structure and it
% is best to keep this conversion separate so that it is easy to modify the output data structure of
% the function again. The algorithm is also simpler(?) when working with arrays.
% flag_results = [];
% for jf = 1:length(flags)    
%     %flag_results.(flag.return_field_name).is_set = flags(jf).is_set;
%     %flag_results.(flag.return_field_name).value  = flags(jf).value;
%     
%     flag = flags(jf);
%     flag_results.(flag.return_field_name) = rmfield(flag, 'return_field_name');
%     
% end

end
