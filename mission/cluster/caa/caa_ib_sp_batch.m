function caa_ib_sp_batch(iso_t,dt,cl_id,int2)
%CAA_IB_SP_BATCH  produce summary plots for internal burst
%
% CAA_IB_SP_BATCH(ISO_T,DT,CL_ID,[INT2])
%
% INT2 - plot interval half-width
%
% List of internal bursts is obtained from listing isdat cluster/burst dir.
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

DP = '/data/cluster';
DB = 'db:9';

vars0 = {'tmode','fdm','efwt','ibias','p','e','a','sax','r','v','bfgm'};
vars1 = {'whip','sweep','bdump','probesa','p','ps' 'dies','die','pburst','dieburst'};

B_DT = 300;
B_DELTA = 60;

if nargin < 4, int2 = 30; end

[st,dt] = irf_stdt(iso_t,dt);

t = fromepoch(st);
t0 = toepoch([t(1) t(2) t(3) 0 0 0]);
t = fromepoch(st+dt+ 86400 -1);
t1 = toepoch([t(1) t(2) t(3) 0 0 0]);

for day=t0:86400:t1
  fname = irf_fname(day,1);
  files = dir([DP '/burst/' fname(1:6) '*we.0' num2str(cl_id)]);
  if isempty(files), continue, end

  nfiles = length(files);
  irf_log('load',...
    ['20' fname(1:6) ' C' num2str(cl_id) ' : ' num2str(nfiles) ' bursts(s)'])
  for f=1:nfiles
    s = files(f).name;
    bstart = iso2epoch(['20' s(1:2) '-' s(3:4) '-' s(5:6) 'T' ...
      s(7:8) ':' s(9:10) ':' s(11:12) 'Z']);
    sp = [pwd() '/' irf_fname(bstart)];
    if ~exist(sp,'dir'), mkdir(sp), end
    cdb = ClusterDB(DB,DP,sp);
    data = getData(cdb,bstart-B_DELTA, B_DT, cl_id,'pburst');
    if isempty(data)
      data = getData(cdb,bstart-B_DELTA, B_DT, cl_id,'eburst');
    end
    if isempty(data), irf_log('load','NO DATA'), continue, end
    st_int = data{2}(1,1) -int2; clear data

    no_data = 0;
    for v=1:length(vars0)
      data = getData(cdb,st_int, 2*int2, cl_id,vars0{v});
      if isempty(data) && (strcmp(vars0{v},'tmode') || strcmp(vars0{v},'fdm'))
        irf_log('load','No EFW data')
        no_data = 1;
        break
      end
    end
    if no_data, continue, end
    cp = ClusterProc(sp);
    for v=1:length(vars1)
      getData(cp,cl_id,vars1{v});
    end
    figure(61), clf
    summaryPlot(cp,cl_id,'fullb','ib','st',st_int,'dt',2*int2)
    orient tall
    print('-dpdf',['C' num2str(cl_id),'_EFW_IB_' irf_fname(bstart)])
  end
end