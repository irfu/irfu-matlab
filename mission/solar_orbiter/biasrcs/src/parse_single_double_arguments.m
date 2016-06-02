% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
% Parse list of command line arguments as if it was a list of "single arguments" and
% "double arguments". Flags are only permitted if they are contained in the submitted lists.
%
% Made-up terminology used in this file:
% Flag        = A command-line argument that (intentionally) matches a predefined (hardcoded)
%               string value (e.g. "--version", "--log").
% Flag value  = A (relatively) arbitrary command-line argument following a flag that it is associated with.
% Single flag = A standalone flag, (sort of) unrelated to other command-line arguments.
% Double flag = The combination of a flag and the following flag value (e.g. "--log ~/temp/").
%               Effectively a way to set a string variable via command-line arguments.
%
% NOTE: Checks preemptively for permitting the same flag twice.
% NOTE: Checks for arguments setting the same flag twice (counting both single & double flags together).
%
% opt = optional
% sgl = single
% dbl = double
function set_flags = parse_single_double_arguments(arguments, opt_sgl_flags, opt_dbl_flags)
%
% QUESTION: What terminology for flags/options/command-line arguments should one use?
% PROPOSAL: Add MANDATORY double flags (not just optional)?

global ERROR_CODES

set_flags.singles = {};
set_flags.doubles = {};
set_flags.double_values = {};



% ASSERTION CHECK
% Note: {C1{:} C2{:}} is a trick to merge cell matrices regardless of their sizes (e.g. row/column vector).
opt_flags = {opt_sgl_flags{:}, opt_dbl_flags{:}};
if length(opt_flags) ~= length(unique(opt_flags))
    errorp(ERROR_CODES.ASSERTION_ERROR, ...
        'The code is configured to accept multiple IDENTICAL command-line flags. This indicates a bug.')
end



i = 1;
while i <= length(arguments)
    arg = arguments{i};
    found_sf = length(find(strcmp(arg, opt_sgl_flags)));   % sf = single flag
    found_df = length(find(strcmp(arg, opt_dbl_flags)));   % df = double flag
    
    if (found_sf == 1) && (found_df == 0)
        
        set_flags.singles{end+1} = arg;
        
    elseif (found_sf == 0) && (found_df == 1)
        
        set_flags.doubles{end+1} = arg;
        if i < length(arguments)
            i = i + 1;
            set_flags.double_values{end+1} = arguments{i};
        else
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not find value expected to follow command-line flag "%s".', arg)
        end
        
    else
        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret argument "%s".', arg)
    end
    
    i = i + 1;
end   % while



% Note: {C1{:} C2{:}} is a trick to merge cell matrices regardless of their sizes (e.g. row/column vector).
all_set_flags = {set_flags.singles{:}, set_flags.doubles{:}};
if length(all_set_flags) ~= length(unique(all_set_flags))
    errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Specified multiple identical CLI argument flags.')
end


end
