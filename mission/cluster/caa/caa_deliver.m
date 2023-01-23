function caa_deliver(sp)
%CAA_DELIVER  deliver CEF files to the CAA
%
% caa_deliver([sp])
%   sp - storage path
%	compress and deliver CEF files to the CAA
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<1, sp=pwd; end

old_pwd = pwd;
cd(sp);
dname = pwd;
ii = find(dname=='/');

fname = [dname(ii(end)+1:end) '.tgz'];
fmask = 'C[1-4]_CP_EFW_L[1-3]_*_V*.cef';

irf_log('proc',['compressing ' fname])
s = unix(['tar zcvf ' fname ' ' fmask],'-echo');
if s~=0, irf_log('save',['problem compressing ' fname]),cd(old_pwd),return, end

irf_log('proc',['uploading ' fname ' to the CAA'])
s = unix(...
  ['echo "put ' fname '"| sftp -b - efw@caa1.estec.esa.int:/c/data-20/EFW/'],...
  '-echo');
if s~=0, irf_log('save',['problem uploading ' fname]),cd(old_pwd),return, end

s = unix(...
  ['echo "' fname ' ' epoch2iso(date2epoch(now),1) '">>/data/caa/log/delivered.log'],...
  '-echo');
if s~=0, irf_log('save','problem updating log '),cd(old_pwd),return, end

irf_log('proc',['moving ' fname ' to storage'])
s = unix(['scp ' fname ' yuri@amanda:/data/caa/delivered && rm ' fname],'-echo');

if s~=0, irf_log('save',['problem moving ' fname]),cd(old_pwd),return, end

irf_log('proc','cleaning up')
unix(['rm ' fmask],'-echo');

irf_log('proc','moving summary plots to storage')
s = unix('mv *.pdf *.png /data/caa/splots/','-echo');
if s~=0, irf_log('save','problem moving summary plots'),cd(old_pwd),return, end

cd(old_pwd)
