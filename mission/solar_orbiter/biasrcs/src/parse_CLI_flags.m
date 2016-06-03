% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
% Parses list of command-line arguments assuming it is a list of flags and pairs of flag+value.
% Function tries to give accurate user-friendly errors for non-compliant arguments and absence of
% required arguments.
%
% Flag = A predefined (hardcoded, more or less) string meant to match a single argument.
% Value = An argument that specifies a more or less arbitrary value and comes after a flag.
%
% flags : array of structs.
% .return_field_name : Name of corresponding field in the return struct.
% .CLI_name          : The command-line flag string.
% .is_required       : ("Required" as opposed to "optional".)
% .expects_value
%
% parsed_args : struct with fields corresponding to flags(i).return_field_name for every flag. Every
% such field contains a copy of an element of flags with some fields added.
% .(<return_field_name>).is_set
% .(<return_field_name>).value
% (+field inherited from flags(i))
%
function parsed_args = parse_CLI_flags(arguments, flags)
%
% TODO/PROPOSAL: Change name: parse_arguments? parse_CLI_arguments?
%
% IMPLEMENTATION NOTE:
% Reasons for implementing function parameter "flags" as array with one struct per flag (instead of
% struct with structs):
% (1) Make it possible for the caller to successively add flags with the same ID and have the
% function check for name collisions.
%
% QUESTION: The function permits flags without value. Is this functionality really needed?
%
global ERROR_CODES



% ASSERTION CHECKS
if ~isstruct(flags)
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Parameter is not struct.');
end



% Add necessary fields to "flags".
for i = 1:length(flags)
    flags(i).is_set = 0;
    flags(i).value = [];
end



% ASSERTION CHECKS: Check that are no configuration doubles.
CLI_names          = {flags.CLI_name};
return_field_names = {flags.return_field_name};
if length(CLI_names) ~= length(unique(CLI_names))
    errorp(ERROR_CODES.ASSERTION_ERROR, ...
        'The code is configured to accept multiple IDENTICAL command-line flags. This indicates a bug.')
elseif length(return_field_names) ~= length(unique(return_field_names))
    errorp(ERROR_CODES.ASSERTION_ERROR, ...
        'The code is configured to return multiple IDENTICAL field names. This indicates a bug.')
end



ia = 1;
while ia <= length(arguments)    % ia = i_argument
    arg = arguments{ia};
    
    flag = [];
    for jf = 1:length(flags)
        if strcmp(arg, flags(jf).CLI_name)
            flag = flags(jf);
            break
        end
    end
    if isempty(flag)
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret argument "%s".', arg)
    end
    
    if flag.is_set ~= 0
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Specified the command-line flag "%s" (at least) twice.', arg)
    end
    flag.is_set = 1;
    
    if flag.expects_value
        if ia >= length(arguments)
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, ...
                'Can not find the argument that is expected to follow command-line flag "%s".', arg)
        end
        ia = ia + 1;
        flag.value = arguments{ia};
    end
    
    flags(jf) = flag;

    ia = ia + 1;
end   % while


% Check that all required flags have been set.
for jf = 1:length(flags)
    flag = flags(jf);
    if flag.is_required && ~flag.is_set
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Missing required command-line flag "%s".', flag.CLI_name)
    end
end


% Convert flags (array) to structure.
% -----------------------------------
% IMPLEMENTATION NOTE: Up until here the algorithm works with the function argument, and amends the
% results to it. First here does it copies the data to a variable that will be returned. This
% function has been modified several times and it is best to keep this conversion separate so that
% it is easy modify the output of the function again. The algorithm is also simpler(?) working with
% arrays.
parsed_args = [];
for jf = 1:length(flags)
    flag = flags(jf);
    parsed_args.(flag.return_field_name) = rmfield(flag, 'return_field_name');
    
end



%------------------------------

end
