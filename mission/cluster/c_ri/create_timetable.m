function [timetable_b, timetable_n] = create_timetable(from,to,cl_id)
%CREATE_TIMETABLE  create a list of FGM data available
%
% [timetable_b, timetable_n] = create_timetable(from,to,cl_id)
%
%     Create a chronoligically ordered list of times
%     for which FGM data is available
%
% Input:
%     from, to - time in epoch
%     cl_id -which cluster satellite
%         NOTE: from and to must be within the same day
%
% Output:
%     timetable - [from | to | mode] , where mode (b/n): burst or normal
%
% See also CREATE_FILE, GET_TIMETABLE
%

%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
narginchk(3,3)
if ~(isnumeric(from) && from >0), error('FROM must be epoch'), end
if ~(isnumeric(to) && to>from), error('TO must be larger then FROM'), end
if ~(isnumeric(cl_id) && any(cl_id==(1:4))), error('CL_ID must be 1..4'), end

file_p_n_b = create_file(from,cl_id,'b');
file_p_n_n = create_file(from,cl_id,'n');

timetable_b = get_timetable(from,to,file_p_n_b);
timetable_n = get_timetable(from,to,file_p_n_n);

unix_command = sprintf('rm %s %s', file_p_n_b ,file_p_n_n);
[s,h] = unix(unix_command);
if s~=0, warning('IRFU:unix','output from %s:\n%s', unix_command, h), end
