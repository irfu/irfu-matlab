% Informal, non-automatic test code

clear
init_global_constants

% flags.a.CLI_name = '-a'
% flags.a.is_required = 0;
% flags.a.expects_value = 0;
% 
% flags.b.CLI_name = '-b';
% flags.b.is_required = 0;
% flags.b.expects_value = 1;
% 
% flags.c.CLI_name = '-c';
% flags.c.is_required = 1;
% flags.c.expects_value = 1;
% 
% %pa1 = bicas.utils.parse_CLI_flags(strsplit('-a'), flags)
% %pa1 = bicas.utils.parse_CLI_flags(strsplit('-a -b'), flags)
% pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123'), flags)
% %pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c'), flags)
% pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999'), flags)
% 
% pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999 -b'), flags)


flags = containers.Map;
flags('a') = struct('CLI_name', '-a', 'is_required', 0, 'expects_value', 0);
flags('b') = struct('CLI_name', '-b', 'is_required', 0, 'expects_value', 1);
flags('c') = struct('CLI_name', '-c', 'is_required', 1, 'expects_value', 1);


%pa1 = bicas.utils.parse_CLI_flags(strsplit('-a'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-a -b'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c'), flags);
pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999'), flags);
%pa1 = bicas.utils.parse_CLI_flags(strsplit('-b 123 -c 999 -b'), flags);


