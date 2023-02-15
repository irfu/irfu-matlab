function caa_make_jobs(year,month,region_s,outdir)
%CAA_MAKE_JOBS  make job files
%
% caa_make_jobs(year,month,region_s,[outdir])
%

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
elseif strcmp(region_s,'az'), rid = 9;
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
  regs = regions(find( regions(:,2)>=ts_tmp & regions(:,1)<=te_tmp),:);
  if isempty(ints), continue, end
  
  for j=1:length(regs(:,1))
    ii = find( (ints(:,1)>=regs(j,1) & ints(:,1)<regs(j,2)) | ...
      (ints(:,2)>regs(j,1) & ints(:,2)<=regs(j,2)));
    if ~isempty(ii)
      ii = torow(ii(:));
      for i=ii
        ts_t = ints(i,1);
        te_t = ints(i,2);
        if ints(i,1)<regs(j,1), ts_t = regs(j,1); end
        if ints(i,2)>regs(j,2), te_t = regs(j,2); end
        job = [job; ts_t te_t-ts_t];
        disp(['1 adding ' epoch2iso(job(end,1),1) '-' ...
          epoch2iso(job(end,1)+job(end,2),1)])
      end
    end
    % intervals which cover the whole region
    ii = find(ints(:,1)<regs(j,1) & ints(:,2)>regs(j,2));
    if ~isempty(ii)
      job = [job; regs(j,1) regs(j,2)-regs(j,1)];
      disp(['2 adding ' epoch2iso(job(end,1),1) '-' ...
        epoch2iso(job(end,1)+job(end,2),1)])
    end
  end
  if ~isempty(job)
    job = sortrows(job,1);
    fname = ['job-' tit '-p' num2str(part) '-' region_s];
    irf_log('save',['writing ' fname])
    fid = fopen([fname '.dat'],'w');
    for j=1:length(job(:,1))
      fprintf(fid,'%s %d\n',epoch2iso(job(j,1),1),job(j,2));
    end
    fclose(fid);
    
    pl_cov(onoff,cover,[ts_tmp te_tmp])
    hold on
    for j=1:length(job(:,1))
      irf_plot([ job(j,1)+[0 job(j,2)]; 4*[1 1]+.2]','k-x')
    end
    hold off
    title(['EFW data coverage and CAA planning (' fname ')'])
    print('-dpng',fname)
    print( gcf, '-dpsc2', fname)
    unix(['/usr/local/bin/ps2pdf12 ' fname '.ps; rm -f ' fname '.ps']);
  end
end

cd(old_pwd)

function pl_cov(onoff,cover,tint)
clf
axes;
set(gca,'YLim',[0 5]);
hold on
for cl_id=1:4
  oo = onoff{cl_id};
  co = cover{cl_id};
  io = find(oo(:,1)>=tint(1) & oo(:,1)<=tint(2));
  if isempty(io)
    io = find(oo(:,1)<tint(1));
    if isempty(io)
      io = find(oo(:,1)>tint(2));
      io = io(1);
    else, io = io(end);
    end
  end
  if oo(io(1),2)==0
    if io(1)>1, io=[io(1)-1; io]; end
  end
  if oo(io(end),2)==1
    if io(end)<length(oo(:,1)), io=[io; io(end)+1]; end
  end
  
  ic = find(co(:,2)>=tint(1) & co(:,1)<=tint(2));
  
  for j=1:2:length(io)-1
    irf_plot([[oo(io(j),1) oo(io(j+1),1)]; (5-cl_id)*[1 1]]','g-x')
  end
  for j=1:length(ic)
    if co(ic(j),3)==1
      irf_plot([[co(ic(j),1) co(ic(j),2)]; (5-cl_id)*[1 1]-.2]','r-x')
    else
      irf_plot([[co(ic(j),1) co(ic(j),2)]; (5-cl_id)*[1 1]-.2]')
    end
  end
  
end
hold off
set(gca,'YTick',[1 2 3 4],'YTickLabel',['C4'; 'C3'; 'C2'; 'C1'])
irf_zoom(tint,'x',gca)
