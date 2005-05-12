function caa_make_jobs(year,month,region_s,outdir)
%CAA_MAKE_JOBS  make job files
%
% caa_make_jobs(year,month,region_s,[outdir])
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

ddir = '/data/cluster/caa/';

if nargin<4, outdir = pwd; end

if strcmp(region_s,'mp'), rid = 1;
elseif strcmp(region_s,'sh'), rid = 2;
elseif strcmp(region_s,'bs'), rid = 3;
elseif strcmp(region_s,'sw'), rid = 4;
elseif strcmp(region_s,'psh'), rid = 5;
elseif strcmp(region_s,'psp'), rid = 6;
elseif strcmp(region_s,'lo'), rid = 7;
elseif strcmp(region_s,'cusp'), rid = 8;
elseif strcmp(region_s,'nodata'), rid = 0;
elseif strcmp(region_s,'undef'), rid = 999;
else, error('unkrovn region')
end

if month==12
	nexty = year + 1;
	nextm = 1;
else
	nexty = year;
	nextm = month + 1;
end

ts = toepoch([year month 1 00 00 00]);
te = toepoch([nexty nextm 1 00 00 00]);

[onoff,cover] = caa_read_coverage(ddir);
[plan, plan_ind] = caa_pl_coverage(onoff,cover, [ts te],0);

m_s = num2str(month);
if length(m_s)<2, m_s =['0' m_s]; end
tit = [num2str(year) '-' m_s];
fname = ['regions-' tit '.xml'];

t_file = tempname;
%disp(['/usr/local/bin/xsltproc ' ddir 'region-text.xsl ' ddir fname '>' t_file]);
%return
[s,w] = unix(['/usr/local/bin/xsltproc ' ddir 'region-text.xsl ' ddir fname '>' t_file]);

if s>0, disp(w), end

[iso_t,region]=textread(t_file,'%s %d',-1);
[s,w] = unix(['rm -f ' t_file]);

t = iso2epoch(cell2mat(iso_t));
ii = find(region==rid);

if ii(end)==length(region), regions = [t(ii(1:end-1)) t(ii(1:end-1)+1);t(end) te];
else, regions = [t(ii) t(ii+1)];
end

old_pwd = pwd;
cd(outdir)

for part=1:6
	job = [];
	ts_tmp = toepoch([year month (part-1)*5+1 00 00 00]);
	if part==6, te_tmp = te;
	else, te_tmp = toepoch([year month (part)*5+1 00 00 00]);
	end
	ints = plan(find( plan(:,2)>=ts_tmp & plan(:,1)<=te_tmp),:);
	if isempty(ints), continue, end
	regs = plan(find( regions(:,2)>=ts_tmp & regions(:,1)<=te_tmp),:);
	if isempty(ints), continue, end
	
	for j=1:length(regs(:,1))
		ii = find( (ints(:,1)>=regs(j,1) &  ints(:,1)<=regs(j,2)) | ...
			(ints(:,2)>=regs(j,1) & ints(:,2)<=regs(j,2)));
		for i=ii
			if ints(i,1)<regs(j,1), ints(i,1) = regs(j,1); end
			if ints(i,2)>regs(j,2), ints(i,2) = regs(j,2); end
			job = [job; ints(i,1) ints(i,2)-ints(i,1)];
		end
	end
	if isempty(job), continue, end
	
	fname = ['job-' tit '-p' num2str(part) '.dat'];
	irf_log('save',['writing ' fname])
	fid = fopen(fname,'w');
	for j=1:length(job(:,1))
		fprintf(fid,'%s %d\n',epoch2iso(job(j,1)),job(j,2));
	end
	fclose(fid);
end

cd(old_pwd)
