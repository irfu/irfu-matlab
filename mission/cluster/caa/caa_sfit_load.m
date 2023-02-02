function [res,msg] = caa_sfit_load(cl_id)
%CAA_SFIT_LOAD  Load EFW spinfit data
%
% res = caa_sfit_load(cl_id)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

res = [];

[probePair,flag_lx] = caa_sfit_probe(cl_id);
if flag_lx
  vs = sprintf('diELXs%dp%d',cl_id,probePair/10);
else
  vs = sprintf('diEs%dp%d',cl_id,probePair);
end

[ok,diEs,msg] = c_load(vs);
if ~ok || isempty(vs), return, end

res = struct('diEs',diEs,'probePair',probePair,'flagLX',flag_lx);
