function file_path_and_name = create_file(time,mode,cl_nr)
%
%Input:
% time -time in epoch
% mode -b/n, means burst or normal
%
%Output:
% file_path_and_name -where the file was saved
%
%Descrition of the function:
% Creates a header file fot the day of "time" and saves it the disk
%
%Using:
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

time_u = fromepoch(time);

%should be in the form ex:020302*fb.01
y = datestr(time_u,11);
m = datestr(time_u,5);
d = datestr(time_u,7);
path = sprintf('/data/cluster/DDS/');

%the filename 
filename = sprintf('%s%s%s*f%s.0%d',y,m,d,mode,cl_nr);
p_f = sprintf('%s%s',path,filename);

%the output file 
output_file = sprintf('t_%s%s%sf%s.0%d',y,m,d,mode,cl_nr);
output_path = sprintf('');
output_p_f = sprintf('%s%s',output_path,output_file);

%the call to the unix function ddsls for the burstmode
%unix_command = sprintf('/home/scb/fgm/bin86/ddsls %s >%s',p_f,output_p_f);
%unix_command = sprintf('/home/scb/fgm/bin/ddsls %s >%s',p_f,output_p_f);
%unix(unix_command);

mext = mexext;
if strcmp(mext,'mexglx') % running on x86
	unix_command = sprintf('/home/scb/fgm/bin86/ddsls %s >%s',p_f,output_p_f);
elseif strcmp(mext,'mexsol') % running on Solaris/SPARC
	unix_command = sprintf('/home/scb/fgm/bin/ddsls %s >%s',p_f,output_p_f);
else
	error('Cannot determine operating system/platform.')
end
unix(unix_command);

file_path_and_name = output_p_f;
