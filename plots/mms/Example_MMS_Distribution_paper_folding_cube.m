% Initialize database
mms.db_init('local_file_db','/Volumes/mms');

% Load data
ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z');
iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
iPDist_counts = iPDist; iPDist_counts.data = (iPDist.data./iPDistErr.data).^2;
 
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);', ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)

L_gse = [1 0 0];
M_gse = [0 1 -0.2];
N_gse = cross(L_gse,M_gse);

c_eval('mvaB? = gseB?*[L_gse; M_gse; N_gse]'';',ic)
c_eval('mvaE? = gseE?*[L_gse; M_gse; N_gse]'';',ic)
c_eval('mvaVExB? = gseVExB?*[L_gse; M_gse; N_gse]''; mvaVExB?.name = ''E LMN'';',ic)

c_eval('tsLgse? = irf.ts_vec_xyz(iPDist.time,repmat(L_gse,iPDist.length,1));',ic)
c_eval('tsMgse? = irf.ts_vec_xyz(iPDist.time,repmat(M_gse,iPDist.length,1));',ic)
c_eval('tsNgse? = irf.ts_vec_xyz(iPDist.time,repmat(N_gse,iPDist.length,1));',ic)

c_eval('tsLdsl? = mms_dsl2gse(tsLgse?,defatt?,-1);',ic)
c_eval('tsMdsl? = mms_dsl2gse(tsMgse?,defatt?,-1);',ic)
c_eval('tsNdsl? = mms_dsl2gse(tsNgse?,defatt?,-1);',ic)

%% Plot distribution
% Structure for folding
%   o
% o o o
%   o
%   o

% Prepare distribution for plotting
ic = 3;
time = irf_time('2017-07-11T22:34:01.000Z','utc>EpochTT');
%time = irf_time('2017-07-11T22:33:50.000Z','utc>EpochTT');
tint_dist = time + 1*0.5*0.150*[-1 1]; % includes one distribution

% Taking a N-point moving mean and removing data below the one-count level
nMovMean = 5;
pdist_movmean = iPDist.movmean(nMovMean,'RemoveOneCounts',iPDist_counts);

vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];
elim = [0 Inf];

pdist = pdist_movmean.elim(elim).tlim(tint_dist);

% Set up coordinate system
t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
Ldsl = Tdsl(1,:);
Mdsl = Tdsl(2,:);
Ndsl = Tdsl(3,:);


c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   


% Set up figure
fig = figure('Position',[1     1   550   696]);

nrows = 4;
ncols = 3;
h(1) = subplot('position',[0.3883    0.7212    0.2583    0.2037]);
h(2) = subplot('position',[0.1300    0.5175    0.2583    0.2037]);
h(3) = subplot('position',[0.3883    0.5175    0.2583    0.2037]);
h(4) = subplot('position',[0.6466    0.5175    0.2583    0.2037]);
h(5) = subplot('position',[0.3883    0.3138    0.2583    0.2037]);
h(6) = subplot('position',[0.3883    0.1101    0.2583    0.2037]);

%compact_panels(h,0.0,0.0)

% Plot distributions
nSmooth = 0;
doPlotB = 1;
doPlotExB = 1;
vlim = 2500;

if doPlotB
  B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
  B_ = mean(B__.data,1);
  b = B_/norm(B_);
  b = b*vlim;
end

isub = 1;
if 1 % f(L,M)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
  vdf.plot_plane(hca,'smooth',nSmooth)
  axis(hca,'square')
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_M (km/s)';
  if doPlotB % plot B direction
    xlim = hca.XLim; % not needed, because done later
    ylim = hca.YLim;
    hold(hca,'on')
    quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    irf_legend(hca,sprintf('B = %.1f nT',sum(sqrt(B_.^2))),[0.98 0.98],'fontsize',7,'color','k')
  end     
  if 1 % plot ExB
    hold(hca,'on')
    hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
    hold(hca,'off')    
  end
end

isub = 3;
if 1 % f(L,N)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',Ldsl,Ndsl,'vint',vint_M);
  vdf.plot_plane(hca,'smooth',nSmooth)
  axis(hca,'square')
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_N (km/s)';

  if doPlotB % plot B direction
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    quiver(-b(1),-b(3),2*b(1),2*b(3),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end     
  if 1 % plot ExB
    hold(hca,'on')
    hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
    hold(hca,'off')    
  end
end

isub = 4;
if 1 % f(M,N)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',Mdsl,Ndsl,'vint',vint_L);
  vdf.plot_plane(hca,'smooth',nSmooth)
  axis(hca,'square')
  hca.XLabel.String = 'v_M (km/s)';
  hca.YLabel.String = 'v_N (km/s)';

  if doPlotB % plot B direction
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    quiver(-b(2),-b(3),2*b(2),2*b(3),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end     
  if 1 % plot ExB
    hold(hca,'on')
    hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
    hold(hca,'off')    
  end
end

isub = 2;
if 1 % f(M,N)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',Mdsl,Ndsl,'vint',vint_L);
  vdf.plot_plane(hca,'smooth',nSmooth)
  axis(hca,'square')
  hca.XLabel.String = 'v_M (km/s)';
  hca.YLabel.String = 'v_N (km/s)';

  if doPlotB % plot B direction
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    quiver(-b(2),-b(3),2*b(2),2*b(3),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end     
  if 1 % plot ExB
    hold(hca,'on')
    hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
    hold(hca,'off')    
  end
  hca.XDir = 'reverse';
end


isub = 5;
if 1 % f(L,M)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
  vdf.plot_plane(hca,'smooth',nSmooth)
  axis(hca,'square')
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_M (km/s)';

  if doPlotB % plot B direction
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end     
  if 1 % plot ExB
    hold(hca,'on')
    hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
    hold(hca,'off')    
  end
  hca.YDir = 'reverse';
end  


isub = 6;
if 1 % f(L,N)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',Ldsl,Ndsl,'vint',vint_M);
  vdf.plot_plane(hca,'smooth',nSmooth)
  axis(hca,'square')
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_N (km/s)';

  if doPlotB % plot B direction
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    quiver(-b(1),-b(3),2*b(1),2*b(3),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end     
  if 1 % plot ExB
    hold(hca,'on')
    hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
    hold(hca,'off')    
  end
  hca.YDir = 'reverse';
end

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1); delete(hb)
c_eval('h(?).FontSize = 5;',1:numel(h))
c_eval('h(?).XAxisLocation = "origin";',1:numel(h))
c_eval('h(?).YAxisLocation = "origin";',1:numel(h))
c_eval('h(?).Layer = "top";',1:numel(h))
c_eval('h(?).XTick = sort([0:-1000:-vlim 1000:1000:vlim]);',1:numel(h))
c_eval('h(?).YTick = sort([0:-1000:-vlim 1000:1000:vlim]);',1:numel(h))
c_eval('h(?).XTickLabelRotation = 90;',1:numel(h))
c_eval('h(?).YLabel.Rotation = 90;',1:numel(h))
c_eval('h(?).YLabel.VerticalAlignment = "top";',1:numel(h))
c_eval('h(?).YLabel.VerticalAlignment = "bottom";',[4])
c_eval('h(?).YLabel.HorizontalAlignment = "right";',1:numel(h))
c_eval('h(?).XLim = vlim*[-1 1];',1:numel(h))
c_eval('h(?).YLim = vlim*[-1 1];',1:numel(h))

linkprop(h,{'CLim'});
colormap(irf_colormap('thermal'))

irf_legend(h(1),{t_dist_center.utc('yyyy:mm:dd'),t_dist_center.utc('HH:MM:SS.mmm')}',[0.02 0.98],'color','k','fontsize',7)
hcb = colorbar(h(1),'location','south');
hcb.Position(3) = 0.1;
hcb.Position(4) = 0.01;
hcb.YLabel.String = 'f_i (s/m^4)';

if 1 % Make lines for scissors
  %%
  pos = h(1).Position;
  x = [pos(1) pos(1)+pos(3)*0.15 pos(1)+pos(3)*0.85 pos(1)+pos(3)*1]; 
  y = [pos(2)+pos(4) pos(2)+pos(4)*1.15 pos(2)+pos(4)*1.15 pos(2)+pos(4)*1.0]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':')

  pos = h(2).Position;
  x = [pos(1) pos(1)+pos(3)*0.15 pos(1)+pos(3)*0.85 pos(1)+pos(3)*1]; 
  y = [pos(2)+pos(4) pos(2)+pos(4)*1.15 pos(2)+pos(4)*1.15 pos(2)+pos(4)*1.0]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':')

  pos = h(4).Position;
  x = [pos(1) pos(1)+pos(3)*0.15 pos(1)+pos(3)*0.85 pos(1)+pos(3)*1]; 
  y = [pos(2)+pos(4) pos(2)+pos(4)*1.15 pos(2)+pos(4)*1.15 pos(2)+pos(4)*1.0]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':')

  pos = h(2).Position;
  x = [pos(1) pos(1)+pos(3)*0.15 pos(1)+pos(3)*0.85 pos(1)+pos(3)*1]; 
  y = [pos(2) pos(2)-pos(4)*0.15 pos(2)-pos(4)*0.15 pos(2)]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':')

  pos = h(4).Position;
  x = [pos(1) pos(1)+pos(3)*0.15 pos(1)+pos(3)*0.85 pos(1)+pos(3)*1]; 
  y = [pos(2) pos(2)-pos(4)*0.15 pos(2)-pos(4)*0.15 pos(2)]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':')

  pos = h(2).Position;
  x = [pos(1) pos(1)-pos(3)*0.15 pos(1)-pos(3)*0.15 pos(1)]; 
  y = [pos(2) pos(2)+pos(4)*0.15 pos(2)+pos(4)*0.85 pos(2)+pos(4)]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':')

  pos = h(4).Position;
  x = [pos(1)+pos(3) pos(1)+pos(3)*1.15 pos(1)+pos(3)*1.15 pos(1)+pos(3)]; 
  y = [pos(2) pos(2)+pos(4)*0.15 pos(2)+pos(4)*0.85 pos(2)+pos(4)]; 
  annotation('line',x(1:2),y(1:2),'linestyle',':')
  annotation('line',x(2:3),y(2:3),'linestyle',':')
  annotation('line',x(3:4),y(3:4),'linestyle',':') 
end
