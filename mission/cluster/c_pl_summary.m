function c_pl_summary(varargin)
%C_PL_SUMMARY   Make EFW summary plots
%
% c_pl_summary([sp],[sc_list],[options])
% Input:
%   sp - storage path. [default '.']
%   sc_list - list of SC [default 1:4]
%   options :
%      'noproc' - do not recalculate anything
%      'sp' - storage directory; // default: '.'
%      'sc_list' - list of SC;   // default: 1:4
%      'ang_limit' - angle limit for E.B=0 (help @ClusterProc/getData)
%      ++ options to summaryPlot (see help summary plot)
%
% It is necessary to obtain the Sunward offset and amplitude factor
% before running this routine (use c_cal_gui)
%
% Examples:
%   c_pl_summary
%   % make plots using data in the current directory,
%   % for SC 1 to 4, with reprocessing of data.
%
%   c_pl_summary('noproc')
%   % same as above, but without reprocessing
%
%   c_pl_summary('noproc','sc_list',[2 4],'fullb')
%   % make plots for SC 2 and 4 without reprocessing.
%   % 'fullb' is passed to summaryPlot
%
%   sc_list = 2:4;
%   sp = '/home/yuri/caa-data/20020304';
%   st = toepoch([2002 03 04 10 00 00]);
%   dt = 30*60;
%   %load the necessary data
%   c_get_batch(st,dt,sc_list,'sp',sp);
%   %correct the S-offset and amplitude
%   c_cal_gui
%   %do plots
%   c_pl_summary('sp',sp,'sc_list',sc_list,'ang_limit',10)
%   % zoom in (1 min after ST, 20 sec interval) with full B FGM
%   c_pl_summary('noproc','sp',sp,'sc_list',sc_list,'st',st+60,'dt',20,'fullb')
%
% See also ClusterProc/summaryPlot, ClusterProc/getData
%
% $Revision$  $Date$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

do_proc = 1;
sp = '.';
sc_list = 1:4;
ang_limit = 15;

if nargin>0, have_options = 1; args = varargin;
else, have_options = 0; args = '';
end

while have_options
  l = 2;
  if ~isempty(args)
    switch(args{1})
      case 'sp'
        if ischar(args{2}), sp = args{2};
        else, irf_log('fcal','wrongArgType : sp must be string')
        end
      case 'sc_list'
        if isnumeric(args{2}), sc_list = args{2};
        else, irf_log('fcal','wrongArgType : sc_list must be numeric')
        end
      case 'ang_limit'
        if isnumeric(args{2}), ang_limit = args{2};
        else, irf_log('fcal','wrongArgType : ang_limit must be numeric')
        end
      case 'noproc'
        do_proc = 0; l=1;
      otherwise
        break
    end
    if length(args) >= l
      args = args(l+1:end);
      if isempty(args), break, end
    else, break
    end
  end
end

cp=ClusterProc(sp);

for cl_id=sc_list
  if do_proc
    getData(cp,cl_id,'edbs','ang_blank','ang_limit',ang_limit);
    getData(cp,cl_id,'edb','ang_blank','ang_limit',ang_limit);
    getData(cp,cl_id,'vedbs');
  end

  figure(cl_id), clf, orient tall
  splot_options = [{cp} {cl_id} args];
  summaryPlot(splot_options{:})
end
