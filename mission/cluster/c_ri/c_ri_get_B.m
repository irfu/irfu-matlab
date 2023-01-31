function B=c_ri_get_B(from,to,cl_id,mode,path_output)
%C_RI_GET_B  read FGM data from DDS files
%
% B = c_ri_get_B(from,to,cl_id,mode,path_output)
%
%	Download B-data from /data/cluster/DDS/filename1
%
%	ex filename1:
%		020302*fY.XX  where Y = (b/n) XX =(01/02/03/04))
%	ex filename1:
%		Ba_020302_F030100_T070403_Y.XX  where Y = (b/n) XX =(01/02/03/04))
%
%
% Input:
%	from,to  - time in epoch
%	cl_id	- cluster number
%	mode  - 'b' or 'n', burst or normal
%
% Output:
%	if no output argument then
%		saved to file: path_output
%		filename ex: Ba_020302_F030100_T070301_b.01, (a stands for ascii)
%	if output argument then return column vector B [isdat_time Bx By Bz]
%
% Work method:
%	Download B-data with "ddscut" from /data/cluster/DDS/filename
%	(ex filename: 020302*fY.XX  where Y = (b/n) XX =(01/02/03/04)
%	and saves it in a temporary file. Then the data is unpacked,
%	calibrate and save in a format which matlab can read
%	using "fmgtel - fgmcal -fgmhrt-fgmvec"
%

%Written by Robert Isaksson in the summer of -03

narginchk(4,5)

%--------------------- the beginning --------------------------
B=[]; % return nothing if no data available

if nargin == 4, path_output = [pwd '/']; end

if ~(isnumeric(from) && from >0), error('FROM must be epoch'), end
if ~(isnumeric(to) && to>=from), error('TO must be larger or equal to FROM'), end
if ~(isnumeric(cl_id) && any(cl_id==(1:4))), error('CL_ID must be 1..4'), end
if ~(ischar(mode) && (strcmp(mode,'b') || strcmp(mode,'n')))
  error('MODE must be n or b')
end

% test the operating system we are running on
mext = mexext;
if strcmp(mext,'mexglx') % running on x86
  CMDPATH = '/home/scb/fgm/bin/';
elseif strcmp(mext,'mexa64') % running on amd64
  CMDPATH = '/home/scb/fgm/bin64/';
elseif strcmp(mext,'mexmaci64')
  CMDPATH = '/Users/huishan/Software/matlab/irfu-matlab/fgm/clfgm/bin64/';
else
  error('Cannot find CMDPATH for the current operating system/platform.')
end

ddscut = [CMDPATH 'ddscut'];
fgmtel = [CMDPATH 'fgmtel'];
fgmcal = [CMDPATH 'fgmcal'];
fgmhrt = [CMDPATH 'fgmhrt'];
fgmvec = [CMDPATH 'fgmvec'];
d_path = sprintf('%s/DDS/',c_ctl(0,'data_path'));

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


d_source = sprintf('%s%s*f%s.0%d',d_path,d_s,mode,cl_id);
tfn = sprintf('tB_%s_%s_%s_%s.0%d',d_s,fhhmmss,thhmmss,mode,cl_id');
to_file = sprintf('%s%s',path_output,tfn);

fromstr = sprintf('%04d-%02d-%02dT%02d:%02d:%02d.%03dZ', from_U(1),...
  from_U(2), from_U(3), from_U(4), from_U(5), fix(from_U(6)),...
  fix((from_U(6)-fix(from_U(6)))*100));
tostr = sprintf('%04d-%02d-%02dT%02d:%02d:%02d.%03dZ', to_U(1),...
  to_U(2), to_U(3), to_U(4), to_U(5), fix(to_U(6)),...
  fix((to_U(6)-fix(to_U(6)))*100));
% cuts out the time intervall and creates a temporary file
unix_command = sprintf('%s -b %s -e %s %s > %s',ddscut,fromstr,tostr,d_source,to_file);
[s,h] = unix(unix_command);
if s~=0, warning('IRFU:unix','output from %s:\n%s', unix_command, h), end

d_source = to_file;
fn = sprintf('Ba_%s_%s_%s_%s.0%d',d_s,fhhmmss,thhmmss,mode,cl_id');
to_file = sprintf('%s%s',path_output,fn);
irf_log('dsrc',['Reading FGM ' d_s ' ' fhhmmss '-' thhmmss ', s/c' num2str(cl_id) ]);

FGMPATH = sprintf('%s/cal/fgm/',c_ctl(0,'data_path'));
if ~exist(FGMPATH,'dir'), error('FGMPATH does not exist'), end

if nargout  % return B
  to_file=tempname;
  unix_command = ['export FGMPATH; FGMPATH=' FGMPATH '; ' fgmtel ' ' d_source ' | ' fgmcal ' | ' fgmhrt ' -a ' d_path d_s '*ga.0' num2str(cl_id) ' > ' to_file];
  [stat, res] = unix(['/bin/sh -c ''' unix_command '''']);
  if stat
    warning('IRFU:unix','error when running DP pipe:\n%s', unix_command);
    if ~isempty(res), warning('IRFU:unix','output:\n%s', unix_command, res); end
    return;
  end
  to_file_attr=dir(to_file);
  if to_file_attr.bytes == 0, return;end
  fvs = fgmvec_stream(to_file);
  ta=tavail(fvs); % debb=tavail(fvs,[]);disp(debb);clear debb;
  if min(ta)>0
    try
      dat = get(fvs, 'data', 'b', ['T00:00:00Z' 'T24:00:00Z']);
      B=[rem(dat.time,1)*3600*24+toepoch(fromepoch(from).*[1 1 1 0 0 0])...
        dat.b];
    catch
      % get() can produce error
      irf_log('dsrc','FGM Error')
    end
  end
  close(fvs);
  [s,h] = unix(['rm ' to_file]);
  if s~=0, warning('IRFU:unix','output from %s:\n%s', ['rm ' to_file], h); end
else
  % download an unpack the downloaded data
  unix_command = ['export FGMPATH; FGMPATH=' FGMPATH '; ' fgmtel ' ' d_source ' | ' fgmcal ' | ' fgmhrt ' -a ' d_path d_s '*ga.0' num2str(cl_id) ' | ' fgmvec ' > ' to_file];
  [s,h] = unix(['/bin/sh -c ''' unix_command '''']);
  if s~=0, warning('IRFU:unix','output from %s:\n%s', unix_command, h), end
  disp(['saved to file ' to_file])
end

% remove the temporary files
[s,h] = unix(['rm ' d_source]);
if s~=0, warning('IRFU:unix','output from %s:\n%s', ['rm ' d_source], h), end

