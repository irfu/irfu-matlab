% Informal, non-automatic test code

clear
CONSTANTS = bicas.constants([]);


% flags.a.cliString    = '-a'
% flags.a.isRequired   = 0;
% flags.a.expectsValue = 0;
% 
% flags.b.cliString    = '-b';
% flags.b.isRequired   = 0;
% flags.b.expectsValue = 1;
% 
% flags.c.cliString    = '-c';
% flags.c.isRequired   = 1;
% flags.c.expectsValue = 1;
% 
% %pa1 = bicas.utils.parse_CLI_flags(strsplit('-a'), flags)
% %pa1 = bicas.utils.parse_CLI_flags(strsplit('-a -b'), flags)
% pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123'), flags)
% %pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c'), flags)
% pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999'), flags)
% 
% pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999 -b'), flags)


flags = containers.Map;
flags('a') = struct('cliString', '-a', 'isRequired', 0, 'expectsValue', 0);
flags('b') = struct('cliString', '-b', 'isRequired', 0, 'expectsValue', 1);
flags('c') = struct('cliString', '-c', 'isRequired', 1, 'expectsValue', 1);


%pa1 = bicas.utils.parse_CLI_flags(strsplit('-a'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-a -b'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c'), flags);
pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999 -b'), flags);
