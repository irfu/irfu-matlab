function B=c_ri_get_B(from,to,cl_nr,mode,path_output)
% function B=c_ri_get_B(from,to,cl_nr,mode,path_output)
%
%Input:
% from -time in epoch
% to	-time in epoch
% cl_nr	-cluster number
% mode - 'b' or 'n', burst or normal
%
%Output:
% if no output argument then
%   saved to file:
%   path_output
%   filename ex: Ba_020302_F030100_T070301_b.01, (where the a stands for ascii)
% if output argument then return column vector B [isdat_time Bx By Bz]
%
%Descrition of the function:
% Download B-data from /data/cluster/DDS/filename1
% and saves the data in /share/robert/B_data/filename2
% (ex filename1: 020302*fY.XX  where Y = (b/n) XX =(01/02/03/04))
% (ex filename1: Ba_020302_F030100_T070403_Y.XX  where Y = (b/n) XX =(01/02/03/04))
%
%Using:
% ddscut
% fgmtel
% fgmcal
% fgmhrt
% fgmvec
% 
%Work method:
% Download B-data with "ddscut" from /data/cluster/DDS/filename
% (ex filename: 020302*fY.XX  where Y = (b/n) XX =(01/02/03/04)
% and saves it in a temporary file. Then the data is unpacked, calibrate and save in a format
% which matlab can read using "fmgtel - fgmcal -fgmhrt-fgmvec"
%
%Error:
%
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if nargin == 4
path_output = [pwd '/'];
end

from_U =fromepoch(from);
to_U = fromepoch(to);

y = datestr(from_U,11);
m = datestr(from_U,5);
d = datestr(from_U,7);
d_s = sprintf('%s%s%s',y,m,d);
fhh_mm_ss = datestr(from_U,13);
thh_mm_ss = datestr(to_U,13);

fhhmmss = sprintf('F%s%s%s', fhh_mm_ss(1:2) ,fhh_mm_ss(4:5) ,fhh_mm_ss(7:8));
thhmmss = sprintf('T%s%s%s', thh_mm_ss(1:2) ,thh_mm_ss(4:5) ,thh_mm_ss(7:8));

d_path = sprintf('/data/cluster/DDS/');

d_source = sprintf('%s%s*f%s.0%d',d_path,d_s,mode,cl_nr);
tfn = sprintf('tB_%s_%s_%s_%s.0%d',d_s,fhhmmss,thhmmss,mode,cl_nr');
to_file = sprintf('%s%s',path_output,tfn);

%cuts out the time intervall and creates a temporary file
%unix_command = sprintf('/home/scb/fgm/bin86/ddscut -b %d -e %d %s > %s',from,to,d_source,to_file);
unix_command = sprintf('/home/scb/fgm/bin/ddscut -b %d -e %d %s > %s',from,to,d_source,to_file);
unix(unix_command);

d_source = to_file;
fn = sprintf('Ba_%s_%s_%s_%s.0%d',d_s,fhhmmss,thhmmss,mode,cl_nr');
to_file = sprintf('%s%s',path_output,fn);
%get_fgm = sprintf('/home/scb/fgm/bin86/fgmtel %s | /home/scb/fgm/bin86/fgmcal | /home/scb/fgm/bin86/fgmhrt -a %s%s*ga.0%d |/home/scb/fgm/bin86/fgmvec > %s ',d_source,d_path,d_s,cl_nr,to_file);
get_fgm = sprintf('/home/scb/fgm/bin/fgmtel %s | /home/scb/fgm/bin/fgmcal | /home/scb/fgm/bin/fgmhrt -a %s%s*ga.0%d |/home/scb/fgm/bin/fgmvec > %s ',d_source,d_path,d_s,cl_nr,to_file);

if nargout,  % return B
  to_file='/tmp/sckmvnskjaqwedasdawd';
  get_fgm = sprintf('/home/scb/fgm/bin/fgmtel %s | /home/scb/fgm/bin/fgmcal | /home/scb/fgm/bin/fgmhrt -a %s%s*ga.0%d |/home/scb/fgm/bin/fgmvec > %s ',d_source,d_path,d_s,cl_nr,to_file);
  unix(get_fgm);
  fp = fopen(to_file);
  [S,count] = fscanf(fp,'%4d%1s%2d%1s%2d%1s%2d%1s%2d%1s%6f%2s%f%f%f',[15,inf]);
  fclose(fp);
  unix(['rm ' to_file]);
  [r_s, c_s] = size(S);
  if r_s == 15,
    temp_B = S([1 3 5 7 9 11 13 14 15],:)';
    B = [toepoch(temp_B(:,1:6)) temp_B(:,7:9)];
  else,
    B=[-1 0 0 0];      % no data
  end
else
  %download an unpack the downloaded data
  unix_command = sprintf('%s',get_fgm);
  unix(unix_command);
  disp(['saved to file ' to_file])

  %removes the temporary file
  unix_command =sprintf('rm %s',d_source);
  unix(unix_command);
end
