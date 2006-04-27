function ints_out = caa_get_ns_ops_int(st,dt,ns_ops,ids)
%CAA_GET_NS_OPS_INT  get interval conerging particular problem in EFW NS_OPS
%
% ints_out = caa_get_ns_ops_int(st,dt,ns_ops,ids)
%
% See also: C_CTL
%
% $Id$

% Copyright 2006, Yuri Khotyaintsev

if isempty(ns_ops), error('Empty NS_OPS'), end

% Remove records which cover permanent problems (as loss of 
% probes, filters, etc.) as these must be programmed separately
ii = find(ns_ops(:,2)==-1);
ns_ops(ii,:) = [];

ints_out = [];

id = caa_str2errid(ids);
ns_ops = ns_ops(ns_ops(:,4)==id,:);

if isempty(ns_ops), return, end

% Problem covers the whole interval
ii = find( ns_ops(:,1)<=st & ns_ops(:,1)+ns_ops(:,2)>=st+dt );
if ~isempty(ii)
	ints_out = [st st+dt];
	return
end

% Problems starts inside the interval and ends after the interval
while 1
	ii = find( ns_ops(:,1)<st+dt & ns_ops(:,1)>st & ns_ops(:,1)+ns_ops(:,2)>=st+dt);
	if isempty(ii), break, end
	
	ints_out = [ints_out; ns_ops(ii(1),1) st+dt];
	
	% clear already processed records
	ns_ops(ii(1),:) = [];
end

% Problems starts before the interval and ends inside the interval
while 1
	ii = find( ns_ops(:,1)<=st & ns_ops(:,1)+ns_ops(:,2)>st & ns_ops(:,1)+ns_ops(:,2)<=st+dt);
	if isempty(ii), break, end
	
	ints_out = [ints_out; st ns_ops(ii(1),1)+ns_ops(ii(1),2)];
	
	% clear already processed records
	ns_ops(ii(1),:) = [];
end

% Problem is inside the interval
ii = find( ns_ops(:,1)>st & ns_ops(:,1)+ns_ops(:,2)<st+dt );
if isempty(ii), return, end

ints_out = [ints_out; ns_ops(ii(:),1) ns_ops(ii(:),1)+ns_ops(ii(:),2)];
