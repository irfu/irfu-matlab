function c_export_ps(st,varargin)
%C_EXPORT_PS export figures 1:4 into PS and PDF
%
% function c_export_ps(st,[sc_list],[option,value])
% Input:
%	st - isdat epoch, defines filename YYYYMMDD_HHMM
%	if st is omitted, we try to guess it from userdata of figure 1
%	sc_list - list of sc [optional], default 1:4
%	option [optional]:
%	sp  - storage path
%	suf - suffix to add // file name becomes YY...HMM_SUF
%	pdf - export PDF only
%	ps  - export PS only
%
% Example:
%	c_export_ps(toepoch([2002 03 04 10 00 00]),[2 4],'suf','zoom1','pdf')
%
% See also C_EXPORT_PNG, IRF_FNAME
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


if nargin>2, have_options = 1; args = varargin;
else, have_options = 0;
end

sp = '.';
suf = '';
ex_ps = 1;
ex_pdf = 1;
sc_list = 1:4;

while have_options
  l = 2;
  switch(args{1})
    case 'sc_list'
      if length(args) < 2, error('SC_LIST requires a value'), end
      if isnumeric(args{2}), sc_list = args{2};
      else, irf_log('fcal','wrong ArgType : suf must be numeric')
      end
    case 'sp'
      if length(args) < 2, error('SP requires a value'), end
      if ischar(args{2}), sp = args{2};
      else, irf_log('fcal','wrong ArgType : sp must be string')
      end
    case 'suf'
      if length(args) < 2, error('SUF requires a value'), end
      if ischar(args{2}), suf = ['_' args{2}];
      else, irf_log('fcal','wrong ArgType : suf must be string')
      end
    case {'pdf','nops'}
      ex_ps = 0;
      l = 1;
    case {'nopdf','ps'}
      ex_pdf = 0;
      l = 1;
    otherwise
      irf_log('fcal',['Option ''' args{1} '''not recognized'])
  end
  if length(args) > l, args = args(l+1:end);
  else, break
  end
end

old_pwd=pwd;
cd(sp)
if nargin<1
  figs=get(0,'children');if length(figs)<4,sc_list=figs;end
  ud=get(figs(1),'userdata');
  if isfield(ud,'t_start_epoch'), st=ud.t_start_epoch;
  else, error('please specify ts')
  end
else
  if ischar(st), st = iso2epoch(st); end
end

for cl_id=sc_list(:)'
  figure(cl_id)
  orient tall
  fn = sprintf('EFW_C%d_%s%s',cl_id,irf_fname(st),suf);
  irf_log('save',['saving ' fn])
  if ex_ps
    print( gcf, '-dpsc2', fn)
  end
  if ex_pdf
    print( gcf, '-dpdf', fn)
  end
end

cd(old_pwd)
