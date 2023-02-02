function [st_out,dt_out] = caa_ns_ops_int(st,dt,ns_ops,errlist)
%CAA_NS_OPS_INT  split/truncate interval according to EFW NS_OPS
%
% [st_out,dt_out] = caa_ns_ops_int(st,dt,ns_ops)
%
% See also: C_CTL, CAA_GET_NS_OPS_INT
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if isempty(ns_ops), error('Empty NS_OPS'), end

if nargin<4, errlist = []; end

% Remove records which cover permanent problems (as loss of
% probes, filters, etc.) as these must be programmed separately
ns_ops(ns_ops(:,2)==-1,:) = [];

% Problem covers the whole interval
ii = find( ns_ops(:,1)<=st & ns_ops(:,1)+ns_ops(:,2)>=st+dt );
if ~isempty(ii)
  for j=1:length(ii)
    if match_err(ns_ops(ii(j),4),errlist)
      % no/bad data - remove the interval
      irf_log('proc',prob_s(ns_ops(ii(j),:)))
      irf_log('proc',	'blanking the whole interval')
      st_out = [];
      dt_out = [];
      return
    else, irf_log('proc',prob_s(ns_ops(ii(j),:),1))
    end
  end
end
ns_ops(ii,:) = [];

% Problems starts inside the interval and ends after the interval
while 1
  ii = find( ns_ops(:,1)<st+dt & ns_ops(:,1)>st & ns_ops(:,1)+ns_ops(:,2)>=st+dt);
  if isempty(ii), break, end
  if match_err(ns_ops(ii(1),4),errlist)
    % no/bad data - truncate the interval
    irf_log('proc',prob_s(ns_ops(ii(1),:)))
    dt = ns_ops(ii(1),1) - st;
    irf_log('proc',	['truncating interval: setting DT to ' num2str(dt)])
  else, irf_log('proc',prob_s(ns_ops(ii(1),:),1))
  end
  % clear already processed records
  ns_ops(ii(1),:) = [];
end

% Problems starts before the interval and ends inside the interval
while 1
  ii = find( ns_ops(:,1)<=st & ns_ops(:,1)+ns_ops(:,2)>st & ns_ops(:,1)+ns_ops(:,2)<=st+dt);
  if isempty(ii), break, end
  if match_err(ns_ops(ii(1),4),errlist)
    % no/bad data - truncate the interval
    irf_log('proc',prob_s(ns_ops(ii(1),:)))
    et = st + dt;
    st = ns_ops(ii(1),1) + ns_ops(ii(1),2);
    dt = et - st;
    irf_log('proc',	['truncating interval: setting START_TIME to ' epoch2iso(st,1)])
  else, irf_log('proc',prob_s(ns_ops(ii(1),:),1))
  end
  % clear already processed records
  ns_ops(ii(1),:) = [];
end

st_out = st; dt_out = dt;
% Problem is inside the interval
found = 1;
while found
  found = 0;
  for in=1:length(st_out)
    ii = find( ns_ops(:,1)>st_out(in) & ns_ops(:,1)+ns_ops(:,2)<st_out(in)+dt_out(in));
    if ~isempty(ii)
      if match_err(ns_ops(ii(1),4),errlist)
        % no/bad data - truncate the interval
        irf_log('proc',prob_s(ns_ops(ii(1),:)))
        st = st_out(in);
        dt = dt_out(in);
        st_out(in+1:end+1) = st_out(in:end);
        dt_out(in+1:end+1) = dt_out(in:end);
        dt_out(in) = ns_ops(ii(1),1) - st;
        st_out(in+1) = ns_ops(ii(1),1) + ns_ops(ii(1),2);
        dt_out(in+1) = st + dt - st_out(in+1);
        irf_log('proc',	['splitting: ' ...
          epoch2iso(st,1) ' -- ' epoch2iso(st+dt,1)])
        irf_log('proc',	['to       : ' ...
          epoch2iso(st_out(in),1) ' -- ' epoch2iso(st_out(in)+dt_out(in),1)])
        irf_log('proc',	['and      : ' ...
          epoch2iso(st_out(in+1),1) ' -- ' epoch2iso(st_out(in+1)+dt_out(in+1),1)])
        % clear already processed records
        ns_ops(ii(1),:) = [];
        found = 1;
        break
      else
        irf_log('proc',prob_s(ns_ops(ii(1),:),1))
        % clear already processed records
        ns_ops(ii(1),:) = [];
      end
      
    end
  end
  if ~found, break, end
end

function res = match_err(opcode,errlist)
% See if OPCODE matches error condition
res = (opcode<10 && opcode>0) || (~isempty(errlist) && any(opcode==errlist));

function ss = prob_s(ns_ops_rec,warn)
if nargin<2, warn=0; end
if warn, s = 'WARNING: ';
else, s = 'PROBLEM: ';
end
ss = [s caa_errid2str(ns_ops_rec(4)) ' ' epoch2iso(ns_ops_rec(1),1)...
  ' -- ' epoch2iso(ns_ops_rec(1)+ns_ops_rec(2),1)];
