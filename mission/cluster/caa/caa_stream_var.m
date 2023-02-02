function res = caa_stream_var(start_time,dt,dsetName,varName)
%CAA_STREAM_VAR  stream a variable from the CAA/CSA
%
% res = CAA_STREAM_VAR(start_time,dt,dsetName,varName)
%
% See also: CAA_DOWNLOAD

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

currentDir = pwd; tempDir = tempname;
try
  mkdir(tempDir);
  cd(tempDir);
  caa_download(start_time + [0 dt],dsetName,'stream')
  cd(['CAA' filesep dsetName]);
  d=dir('*.cef.gz');
  if isempty(d)
    res = [];
  else
    res = c_caa_cef_var_get(varName,d.name);
  end
catch
  irf_log('dsrc',['Error downloading from CAA: ' dsetName '/' varName])
end
cd(currentDir);
if exist(tempDir,'dir'), rmdir(tempDir,'s'); end