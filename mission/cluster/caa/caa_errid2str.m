function [res_l, res_s] = caa_errid2str(id)
%CAA_ERRID2STR  convert ns_ops error ID to string
%
% [res_l, res_s] = caa_errid2str(id)
%
% Numbers in this file must be in sync with
% /data/cluster/caa-control/text.xsl
%
% See also CAA_STR2ERRID

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

switch id
  case 0
    r.s = 'no_tm';
    r.l = 'No TM';
  case 1
    r.s = 'bad_tm';
    r.l = 'Bad TM';
  case 2
    r.s = 'bad_data';
    r.l = 'Bad data';
  case 11
    r.s = 'no_p1';
    r.l = 'No probe 1';
  case 12
    r.s = 'no_p2';
    r.l = 'No probe 2';
  case 13
    r.s = 'no_p3';
    r.l = 'No probe 3';
  case 14
    r.s = 'no_p4';
    r.l = 'No probe 4';
  case 15
    r.s = 'hxonly';
    r.l = 'HXONLY sampling mode';
  case 16
    r.s = 'bad_bias';
    r.l = 'Bad bias';
  case 17
    r.s = 'bad_hx';
    r.l = 'Bad HX';
  case 18
    r.s = 'bad_lx';
    r.l = 'Bad LX';
  case 19
    r.s = 'high_bias';
    r.l = 'High bias saturation';
  case 112
    r.s = 'no_p12';
    r.l = 'No probe 12';
  case 132
    r.s = 'no_p32';
    r.l = 'No probe 32';
  case 134
    r.s = 'no_p34';
    r.l = 'No probe 34';
  case 255
    r.s = 'no_10Hz_filt';
    r.l = 'No 10Hz filter';
  case 260
    r.s = 'no_spin_fits';
    r.l = 'No on-board spin fits';
  case 261
    r.s = 'spec_bias';
    r.l = 'Special bias settings';
  case 999
    r.s = 'TBD';
    r.l = 'Problem to be defined';
  case 9999
    r.s = 'unknown';
    r.l = 'Unknown problem';
  otherwise
    error('Unknown problem ID')
end

if nargout==1
  res_l = r.l;
elseif nargout>1
  res_l = r.l;
  res_s = r.s;
else
  disp(r.l)
end
