function output_file = create_file(time,cl_id,mode)
%CREATE_FILE  create a header file for FGM
%
% output_file = create_file(time,cl_id,mode)
%
%	Create a header file fot the day of "time" and saves it the disk
%
% Input:
%	time -time in epoch
%	mode -b/n, means burst or normal
%

%Written by Robert Isaksson in the summer of -03

narginchk(3,3)

if ~(isnumeric(time) && time>0), error('TIME must be epoch'), end
if ~(isnumeric(cl_id) && any(cl_id==(1:4))), error('CL_ID must be 1..4'), end
if ~(ischar(mode) && (strcmp(mode,'b') || strcmp(mode,'n')))
  error('MODE must be n or b')
end

time_u = fromepoch(time);

%should be in the form ex:020302*fb.01
y = datestr(time_u,11);
m = datestr(time_u,5);
d = datestr(time_u,7);

path = sprintf('%s/DDS/',c_ctl(0,'data_path'));

%the filename
filename = sprintf('%s%s%s*f%s.0%d',y,m,d,mode,cl_id);
p_f = sprintf('%s%s',path,filename);

output_file = sprintf('t_%s%s%sf%s.0%d',y,m,d,mode,cl_id);

mext = mexext;
if strcmp(mext,'mexglx') % running on x86
  unix_command = sprintf('/home/scb/fgm/bin/ddsls %s >%s',p_f,output_file);
elseif strcmp(mext,'mexa64')
  unix_command = sprintf('/home/scb/fgm/bin64/ddsls %s >%s',p_f,output_file);
elseif strcmp(mext,'mexmaci64')
  unix_command = sprintf('/Users/huishan/Software/matlab/irfu-matlab/fgm/clfgm/bin64/ddsls %s >%s',p_f,output_file);
elseif strcmp(mext,'mexsol') % running on Solaris/SPARC
  error('SPARC is not supported at the moment.')
else
  error('Cannot determine operating system/platform.')
end
[s,h] = unix(unix_command);
if s~=0, warning('IRFU:unix','output from %s:\n%s', unix_command, h), end
