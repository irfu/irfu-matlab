function B=c_ri_get_B(from,to,cl_nr,mode,path_output)
% function B=c_ri_get_B(from,to,cl_nr,mode,path_output)
%
%Input:
% from  - time in epoch
% to	  - time in epoch
% cl_nr	- cluster number
% mode  - 'b' or 'n', burst or normal
%
%Output:
% if no output argument then
%   saved to file:
%   path_output
%   filename ex: Ba_020302_F030100_T070301_b.01, (where the a stands for ascii)
% if output argument then return column vector B [isdat_time Bx By Bz]
%
%Description of the function:
% Download B-data from /data/cluster/DDS/filename1
% (ex filename1: 020302*fY.XX  where Y = (b/n) XX =(01/02/03/04))
% (ex filename1: Ba_020302_F030100_T070403_Y.XX  where Y = (b/n) XX =(01/02/03/04))
%
%Work method:
% Download B-data with "ddscut" from /data/cluster/DDS/filename
% (ex filename: 020302*fY.XX  where Y = (b/n) XX =(01/02/03/04)
% and saves it in a temporary file. Then the data is unpacked, calibrate and save in a format
% which matlab can read using "fmgtel - fgmcal -fgmhrt-fgmvec"
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

if nargin == 4
path_output = [pwd '/'];
end

if and(~strcmp(mode,'b'),~strcmp(mode,'n'))
	error('MODE must be n or b')
end

% test the operating system we are running on
mext = mexext;
if strcmp(mext,'mexglx') % running on x86
	CMDPATH = '/home/scb/fgm/bin86/';
elseif strcmp(mext,'mexsol') % running on Solaris/SPARC
	CMDPATH = '/home/scb/fgm/bin/';
else
	error('Cannot determine operating system/platform.')
end

ddscut = [CMDPATH 'ddscut'];
fgmtel = [CMDPATH 'fgmtel'];
fgmcal = [CMDPATH 'fgmcal'];
fgmhrt = [CMDPATH 'fgmhrt'];
fgmvec = [CMDPATH 'fgmvec'];
d_path = sprintf('/data/cluster/DDS/');

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


d_source = sprintf('%s%s*f%s.0%d',d_path,d_s,mode,cl_nr);
tfn = sprintf('tB_%s_%s_%s_%s.0%d',d_s,fhhmmss,thhmmss,mode,cl_nr');
to_file = sprintf('%s%s',path_output,tfn); 


%cuts out the time intervall and creates a temporary file
unix_command = sprintf('%s -b %d -e %d %s > %s',ddscut,from,to,d_source,to_file);
unix(unix_command);

d_source = to_file;
fn = sprintf('Ba_%s_%s_%s_%s.0%d',d_s,fhhmmss,thhmmss,mode,cl_nr');
to_file = sprintf('%s%s',path_output,fn);
disp(['Reading FGM. ' d_s ' ' fhhmmss '-' thhmmss ', s/c' num2str(cl_nr) ]);

FGMPATH = '/share/fgm_cal';
[s,h] = unix('hostname');
if ~strcmp(h(1:end-1),'sanna.irfu.se')
	FGMPATH = ['/net/sanna/export' FGMPATH];
end

if nargout,  % return B
  to_file=tempname;
  unix_command = ['export FGMPATH; FGMPATH=' FGMPATH '; ' fgmtel ' ' d_source ' | ' fgmcal ' | ' fgmhrt ' -a ' d_path d_s '*ga.0' num2str(cl_nr) ' > ' to_file];
  unix(['/bin/sh -c ''' unix_command '''']);
  to_file_attr=dir(to_file);
  if to_file_attr.bytes>0,
    fvs = fgmvec_stream(to_file);tavail(fvs,[])
    dat = get(fvs, 'data', 'b', ['T00:00:00Z' 'T24:00:00Z']);
    B=[rem(dat.time,1)*3600*24+toepoch(fromepoch(from).*[1 1 1 0 0 0]) dat.b];
    close(fvs);
  else
    disp('Zero size FGM data file!');
    B=[-1 NaN NaN NaN];
  end
  unix(['rm ' to_file]);
else
  %download an unpack the downloaded data
  unix_command = ['export FGMPATH; FGMPATH=' FGMPATH '; ' fgmtel ' ' d_source ' | ' fgmcal ' | ' fgmhrt ' -a ' d_path d_s '*ga.0' num2str(cl_nr) ' | ' fgmvec ' > ' to_file];
  unix(['/bin/sh -c ''' unix_command '''']);
  disp(['saved to file ' to_file])

  %removes the temporary file
  unix_command =sprintf('rm %s',d_source);
  unix(unix_command);
end
unix(['rm ' FGMPATH '/tmp_att']);
