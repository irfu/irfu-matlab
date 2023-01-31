function caa_sp_batch(fname)
%CAA_SP_BATCH  produce summary plots
%
% CAA_SP_BATCH(FNAME)
%       Produce summary plots for dirs listed in FNAME
%
% See also: CAA_PL_SUMMARY_L1
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

old_pwd = pwd;
BASE_DIR = '/data/caa/l1';
dirs = textread(fname,'%s');
errors=0;
errlog=[fname '.failed'];

if isempty(dirs), disp('NO DIRS'), cd(old_pwd), return, end

for d=1:length(dirs)
  s = dirs{d};
  if length(s)>12 && strcmp(BASE_DIR,s(1:12)), s = s(14:end); end
  try
    caa_pl_summary_l1(...
      sprintf('%s-%s-%sT%s:00:00Z',s(1:4),s(10:11),s(12:13),s(15:16)),...
      3*3600,	[BASE_DIR '/' s],'saveps')
  catch ME
    cd(old_pwd);
    errors=errors+1;
    fprintf(1,'Error encountered for %s:\n',s);
    disp(ME.identifier);
    disp(ME.message);
    fid=fopen(errlog,'a');
    fprintf(fid,'Error encountered for %s:\n',s);
    fprintf(fid,'%s\n',ME.identifier);
    fprintf(fid,'%s\n',ME.message);
    fprintf(fid,'\n');
    fclose(fid);
  end
end

cd(old_pwd);
fprintf(1,'%i errors encountered. See %s.\n',errors,errlog);

