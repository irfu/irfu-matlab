function caa_reg_skel(year,month,outdir)
%CAA_REG_SKEL  make skeleton XML file from coverage data
%
% caa_reg_skel(year,month,[outdir])
%

% Copyright 2005 Yuri Khotyaintsev

DT = 20*60; % 30 minutes

if nargin<3, outdir=pwd; end

if month==12
  nexty = year + 1;
  nextm = 1;
else
  nexty = year;
  nextm = month + 1;
end

ts = toepoch([year month 1 00 00 00]);
te = toepoch([nexty nextm 1 00 00 00]);

[onoff,cover] = caa_read_coverage('/data/cluster/caa');
[plan, plan_ind] = caa_pl_coverage(onoff,cover, [ts te],0);

if isempty(plan) && isempty(plan_ind), disp('plans are empty. nothing to do'), return, end

if ~isempty(plan_ind) && ~isempty(plan)
  plan = sortrows([plan; plan_ind(:,1:2)],1);
end

while 1
  ii = find(plan(2:end,1)-plan(1:end-1,2)<DT);
  if ~isempty(ii)
    plan(ii(end),2) = plan(ii(end)+1,2);
    plan(ii(end)+1,:) = [];
  else, break
  end
end

if plan(1,1)<ts+DT, plan(1,1) = ts; end
if plan(end,2)>te-DT, plan(end,2) = te; end

if isempty(plan), disp('plan is empty'), return, end

plan =  sortrows(reshape(plan,length(plan(:,1))*2,1),1);
m_s = num2str(month);
if length(m_s)<2, m_s =['0' m_s]; end
tit = [num2str(year) '-' m_s];
fname = ['regions-' tit '-skel.xml'];

old_pwd = pwd;
cd(outdir)

fid = fopen(fname,'w');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<!DOCTYPE regions SYSTEM "https://www.space.irfu.se/regions.dtd">\n');
fprintf(fid,'\n<!--\nList of regions crossed by Cluster for %s\n-->\n<regions>\n',tit);
fprintf(fid,'<title>%s</title>\n',tit);
if ts~=plan(1)
  fprintf(fid,'<region start="%s" desc="nodata"/>\n',epoch2iso(ts,1));
end
desc={'undef','nodata'};
for j=1:length(plan)
  fprintf(fid,'<region start="%s" desc="%s"/>\n',epoch2iso(plan(j),1),...
    desc{(j/2==round(j/2))+1});
end
fprintf(fid,'</regions>\n',tit);
fclose(fid);

cd(old_pwd)

