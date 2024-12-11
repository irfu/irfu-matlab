% Example of how to make a combined plot with time series panels on the
% left and subplot panels on the right. For example, one can use the time
% series as an overview and plot detailed particle distributions for
% several different times, and mark them in the time series panels.

%% Example 1 (empty): Initialize plot
[h1,h2] = initialize_combined_plot(5,2,2,3,'horizontal');
%% Example 1 (empty): Plot timeseries
% Plot your irf_plot panels here using axes handles h1
%% Example 1 (empty): Plot particle distributions
% Define which time(s) you want to plot the detailed distrubtions for
tint = irf.tint('2015-10-16T10:33:30.20Z/2015-10-16T10:33:30.50Z'); %#ok<NASGU>
c_eval('times = desDist?.tlim(tint).time;',ic)

% Make one plot for each time or time interval
for ii = 1:times.length
  tint = times(ii);

  % Mark the time in the time panels
  if exist('hTimeMark','var'); delete(hTimeMark); end % take away mark if plotted before
  hTimeMark = irf_pl_mark(h1,tint.epochUnix','black');

  % Plot same plot for the four spacecraft
  for ic = 1:4
    % Plot pitchangles/projection/skymap using handles h2
  end

  pause(1)
  % Print ...
end

%% Example 2: Initialize plot
[h1,h2] = initialize_combined_plot(5,2,2,3,'vertical');
%[h1,h2] = initialize_combined_plot(5,2,2,0.4,'vertical');
%% Example 2: Plot timeseries
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
ic = 1;

hca = irf_panel('B');
c_eval('irf_plot(hca,{gseB?.tlim(tint).x,gseB?.tlim(tint).y,gseB?.tlim(tint).z,gseB?.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B_{?,GSE}',ic),'(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.99 0.9]);

hca = irf_panel('n');
c_eval('irf_plot(hca,{ne?.tlim(tint),ni?.tlim(tint)},''comp'');',ic)
hca.YLabel.String = {irf_ssub('n_{?}',ic),'(cm^{-3})'};
irf_legend(hca,{'n_e','n_i'},[0.99 0.9]);

hca = irf_panel('vi');
c_eval('irf_plot(hca,{Vi?.tlim(tint).x,Vi?.tlim(tint).y,Vi?.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {irf_ssub('v_i_{?}',ic),'(km/s)'};
irf_legend(hca,{'v_i_x','v_i_y','v_i_z'},[0.99 0.9]);

hca = irf_panel('ve');
c_eval('irf_plot(hca,{Ve?.tlim(tint).x,Ve?.tlim(tint).y,Ve?.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {irf_ssub('v_e_{?}',ic),'(km/s)'};
irf_legend(hca,{'v_e_x','v_e_y','v_e_z'},[0.99 0.9]);

hca = irf_panel('E');
c_eval('irf_plot(hca,{gseE?.tlim(tint).x,gseE?.tlim(tint).y,gseE?.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {irf_ssub('E_{?,GSE}',ic),'(nT)'};
irf_legend(hca,{'E_x','E_y','E_z'},[0.99 0.9]);

irf_zoom(h1(:),'x',tint)
irf_plot_axis_align
%% Example 2: Plot particle distributions
ic = 1;
tintDist = irf.tint('2015-10-16T10:33:30.20Z/2015-10-16T10:33:30.50Z'); %#ok<NASGU>
c_eval('times = desDist?.tlim(tintDist).time;',ic)


for ii = 1:times.length
  tint = times(ii)+[-1 1]*0.5*0.03*0.5;
  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,tint.epochUnix','green');
  vlim = 15*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [0 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];

  for ic = 1:4
    c_eval('dist = desDist?;',ic)

    c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic);
    hatVe0 = double(irf_norm(Ve0));

    % Get mean magnetic field direction
    c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(tint).time).data);',ic);
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic);
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);

    vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e'};%0;hatVe0,'V_e'};

    % Projection coordinate system
    if 1
      x = hatB0;
      y = hatExB0;
      z = cross(x,y);
    else
      x = [1 0 0]; %#ok<UNRCH>
      y = [0 1 0];
      z = [0 0 1];
    end



    isub = ic;

    if 0 % Plot psd 0 90 180
      hca = h2(isub); isub = isub + 1; %#ok<UNRCH>
      psdtint = tint;%+[-1 1]*0.03;
      c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
      TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
      %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      %hca.Title.String = irf_ssub('C?',ic);
      hca.Title.String = TitleStr;
    end

    if 0 % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;       %#ok<UNRCH>
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};

      % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};
    end

    % Plot project ion onto a plane
    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end

    hca = h2(isub); %isub = isub + 1;
    mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
    hca.Title.String = titleStr;
    colormap(hca,strCMap)

    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end
  end
  pause(1)
  %cn.print([irf_ssub('BvnP_1proj_mms1234_gseish',1) irf_time(times(ii),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end

%% Example 3 (Ohm's law, 4 sc): Transform the data into LMN coordinates
% Set up coordinate system
coordSystem = 4;
switch coordSystem
  case 1 % N: minimum variance of B
    [out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
    L = v(1,:); M = v(2,:); N = v(3,:);
    coordLabels = {'L','M','N'};
    lmn = [N;-M;L];
  case 2 % N: minimum variance of J
    [out,l,v] = irf_minvar(gseJ1.tlim(irf.tint('2015-10-16T10:33:23.000Z/2015-10-16T10:33:32.000Z')));
    L = -v(2,:); M = -v(1,:); N = v(3,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
  case 3 % N: maximum variance of E
    tint = irf.tint('2015-10-16T10:33:30.100Z/2015-10-16T10:33:30.400Z');
    [out,l,v] = irf_minvar(gseE3.tlim(tint));
    L = v(2,:); M = -v(3,:); N = -v(1,:);
    coordLabels = {'N','-M','L'};
    lmn = [N;-M;L];
  case 4 % N: magnetosheath side normal derived from mms1 and mms4
    gseVec14 = gseR4-gseR1; gseVec14 = gseVec14.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec14.time,M);
    gseNorm14 = gseVec14.cross(gseM);
    gseNormalMSH = gseNorm14/gseNorm14.abs;

    N = -gseNormalMSH.data;
    M = M;
    L = cross(mshM,mshN);
  case 5 % N: magnetosphere side normal derived from mms1 and mms4
    gseVec34 = gseR4-gseR3; gseVec34 = gseVec34.resample(tint.start);
    gseM = irf.ts_vec_xyz(gseVec34.time,M);
    gseNorm34 = gseVec34.cross(gseM);
    gseNormalMSP = gseNorm34/gseNorm34.abs;

    N = -gseNormalMSP.data;
    M = M;
    L = cross(mspM,mspN);
end
% Rotate data
c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]);')
c_eval('mvaB? = irf.ts_vec_xyz(gseB?.time,[gseB?.dot(L).data gseB?.dot(M).data gseB?.dot(N).data]);')
c_eval('mvaE? = irf.ts_vec_xyz(gseE?.time,[gseE?.dot(L).data gseE?.dot(M).data gseE?.dot(N).data]);')
c_eval('mvaVe? = irf.ts_vec_xyz(gseVe?.time,[gseVe?.dot(L).data gseVe?.dot(M).data gseVe?.dot(N).data]);')
c_eval('mvaVi? = irf.ts_vec_xyz(gseVi?.time,[gseVi?.dot(L).data gseVi?.dot(M).data gseVi?.dot(N).data]);')
c_eval('mvaJ? = irf.ts_vec_xyz(gseJ?.time,[gseJ?.dot(L).data gseJ?.dot(M).data gseJ?.dot(N).data]); mvaJ?.units = gseJ?.units;')
mvaJcurl = irf.ts_vec_xyz(Jcurl.time,[Jcurl.dot(L).data Jcurl.dot(M).data Jcurl.dot(N).data]); mvaJcurl.units = Jcurl.units;
%gradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,Pe1,Pe2,Pe3,Pe4);
c_eval('mvaPe? = mms.rotate_tensor(gsePe?,''rot'',L,M,N); mvaPe? = irf.ts_tensor_xyz(mvaPe?.time,mvaPe?.data); mvaPe?.units = Pe?.units;',ic)
gradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4);
mvaGradPe = irf.ts_vec_xyz(gradPe.time,[gradPe.dot(L).data gradPe.dot(M).data gradPe.dot(N).data]);
c_eval('mvaVexB? = mvaVe?.resample(mvaB?.time).cross(mvaB?)*1e-3; mvaVexB?.units = ''mV/m''')
c_eval('mvaVixB? = mvaVi?.resample(mvaB?.time).cross(mvaB?)*1e-3; mvaVixB?.units = ''mV/m''')
c_eval('mvaJxB? = mvaJ?.resample(mvaB?.time).cross(mvaB?);')
mvaAvE = (mvaE1+mvaE2.resample(mvaE1.time)+mvaE3.resample(mvaE1.time)+mvaE4.resample(mvaE1.time))/4;
mvaAvVe = (mvaVe1+mvaVe2.resample(mvaVe1.time)+mvaVe3.resample(mvaVe1.time)+mvaVe4.resample(mvaVe1.time))/4;
mvaAvVi = (mvaVi1+mvaVi2.resample(mvaVi1.time)+mvaVi3.resample(mvaVi1.time)+mvaVi4.resample(mvaVi1.time))/4;
mvaAvB = (mvaB1+mvaB2.resample(mvaB1.time)+mvaB3.resample(mvaB1.time)+mvaB4.resample(mvaB1.time))/4;
mvaAvJ = (mvaJ1+mvaJ2.resample(mvaJ1.time)+mvaJ3.resample(mvaJ1.time)+mvaJ4.resample(mvaJ1.time))/4; mvaAvJ.units = mvaJ1.units;
mvaAvVexB = (mvaVexB1+mvaVexB2.resample(mvaVexB1.time)+mvaVexB3.resample(mvaVexB1.time)+mvaVexB4.resample(mvaVexB1.time))/4;
mvaAvVixB = (mvaVixB1+mvaVixB2.resample(mvaVixB1.time)+mvaVixB3.resample(mvaVixB1.time)+mvaVixB4.resample(mvaVixB1.time))/4;
avNe = (ne1.resample(ne1.time)+ne2.resample(ne1.time)+ne3.resample(ne1.time)+ne4.resample(ne1.time))/4;

% Diamagnetic drift
avVDe = -1*mvaGradPe.cross(mvaAvB.resample(mvaGradPe.time))/(-1*e)/avNe/mvaAvB.abs/mvaAvB.abs*1e-6*1e9*1e-9*1e-3;
avVDi = -1*mvaGradPe.cross(mvaAvB.resample(mvaGradPe.time))/(+1*e)/avNe/mvaAvB.abs/mvaAvB.abs*1e-6*1e9*1e-9*1e-3;
%% Example 3 (Ohm's law, 4 sc): Ohm's law terms
units = irf_units;
e = units.e; % add the minus below
mvaOhmGradPe = mvaGradPe.z/avNe.resample(mvaGradPe.time)/e*1e-9*1e-6; mvaOhmGradPe.units = 'mV/m';
mvaOhmVexB = mvaAvVexB; mvaOhmVexB.units = 'mV/m';
mvaOhmVixB = mvaAvVixB; mvaOhmVixB.units = 'mV/m';
mvaOhmJxB_a = mvaAvJ.resample(mvaAvB.time).cross(mvaAvB)/avNe/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_a.units = 'mV/m';
mvaOhmJxB_b = mvaJcurl.resample(mvaAvB.time).cross(mvaAvB)/avNe/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_b.units = 'mV/m';
mvaOhmJxB_c = (mvaJxB1/ne1+mvaJxB2.resample(mvaJxB1.time)/ne2+mvaJxB3.resample(mvaJxB1.time)/ne3+mvaJxB4.resample(mvaJxB1.time)/ne4)/4/e*1e-9*1e-9*1e-6*1e3; mvaOhmJxB_c.units = 'mV/m';
mvaOhmJxB = mvaOhmJxB_c;
% irf_plot({mvaOhmJxB_a,mvaOhmJxB_b,mvaOhmJxB_c},'comp'); legend('<J_{fpi}>x<B>','J_{curl}x<B>','<J_{fpi}xB>')
%% Example 3 (Ohm's law, 4 sc): Initialize plot
[h1,h2] = initialize_combined_plot(9,3,2,3,'horizontal');
%[h1,h2] = initialize_combined_plot(5,2,2,0.4,'vertical');
%% Example 3 (Ohm's law, 4 sc): Plot timeseries with Ohm's terms
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
%tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_plot(hca,{mvaAvE.z,-1*mvaOhmVexB.z,-1*mvaOhmGradPe.z,-1*mvaOhmVexB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % Ohms law ions: N
  hca = irf_panel('Ohms law ions: N');
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_plot(hca,{mvaAvE.z,mvaOhmJxB.z,-1*mvaOhmVixB.z,-1*mvaOhmGradPe.z,mvaOhmJxB.z.resample(mvaAvE.time)-1*mvaOhmVixB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_legend(hca,{'E','jxB/ne','-v_{i}xB','-\nabla \cdot P_e/ne','jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % Ohms law 2: N
  hca = irf_panel('Ohms law 2: N');
  set(hca,'ColorOrder',mms_colors('xyzba'))
  irf_plot(hca,{mvaAvE.z+1*mvaOhmVexB.resample(mvaAvE.time).z,-1*mvaOhmGradPe.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))
  irf_legend(hca,{'E+v_{e}xB','-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.98 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 0 % JN
  hca = irf_panel('JL'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h1(:),'x',tintZoom)
irf_zoom(h1(:),'y')
irf_plot_axis_align
%% Example 3 (Ohm's law, 4 sc): Initialize plot
[h1,h2] = initialize_combined_plot(10,3,2,3,'horizontal');
%[h1,h2] = initialize_combined_plot(5,2,2,0.4,'vertical');
%% Example 3 (Ohm's law, 4 sc): Plot timeseries with other terms
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
%tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % Pe N
  hca = irf_panel('Pe N');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{mvaPe1.tlim(tint).zz,mvaPe2.tlim(tint).zz,mvaPe3.tlim(tint).zz,mvaPe4.tlim(tint).zz},'comp');
  %hca.YLabel.String = {irf_ssub('P_{eN}',ic),'(nPa)'};
  irf_plot(hca,{mvaPe1.tlim(tint).trace/3,mvaPe2.tlim(tint).trace/3,mvaPe3.tlim(tint).trace/3,mvaPe4.tlim(tint).trace/3},'comp');
  hca.YLabel.String = {irf_ssub('P_{e}',ic),'(nPa)'};
end
if 1 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
    mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
    mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
    mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp');
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_plot(hca,{mvaAvE.z,-1*mvaOhmVexB.z,-1*mvaOhmGradPe.z,-1*mvaOhmVexB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.98],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.02 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % Ohms law ions: N
  hca = irf_panel('Ohms law ions: N');
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_plot(hca,{mvaAvE.z,mvaOhmJxB.z,-1*mvaOhmVixB.z,-1*mvaOhmGradPe.z,mvaOhmJxB.z.resample(mvaAvE.time)-1*mvaOhmVixB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_legend(hca,{'E','jxB/ne','-v_{i}xB','-\nabla \cdot P_e/ne','jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[0.98 0.98],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.02 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 0 % JN
  hca = irf_panel('JL'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h1(:),'x',tintZoom)
irf_zoom(h1(:),'y')
irf_plot_axis_align
%% Example 3 (Ohm's law, 4 sc): Initialize plot
[h1,h2] = initialize_combined_plot(10,3,2,3,'horizontal');
%[h1,h2] = initialize_combined_plot(5,2,2,0.4,'vertical');
%% Example 3 (Ohm's law, 4 sc): Plot timeseries with other terms
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
%tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % Pe N
  hca = irf_panel('Pe N');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{mvaPe1.tlim(tint).zz,mvaPe2.tlim(tint).zz,mvaPe3.tlim(tint).zz,mvaPe4.tlim(tint).zz},'comp');
  %hca.YLabel.String = {irf_ssub('P_{eN}',ic),'(nPa)'};
  irf_plot(hca,{mvaPe1.tlim(tint).trace/3,mvaPe2.tlim(tint).trace/3,mvaPe3.tlim(tint).trace/3,mvaPe4.tlim(tint).trace/3},'comp');
  hca.YLabel.String = {irf_ssub('P_{e}',ic),'(nPa)'};
end
if 1 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'v_e_L','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
    mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
    mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
    mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp');
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_plot(hca,{mvaAvE.z,-1*mvaOhmVexB.z,-1*mvaOhmGradPe.z,-1*mvaOhmVexB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.98],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.02 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % Ohms law ions: N
  hca = irf_panel('Ohms law ions: N');
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_plot(hca,{mvaAvE.z,mvaOhmJxB.z,-1*mvaOhmVixB.z,-1*mvaOhmGradPe.z,mvaOhmJxB.z.resample(mvaAvE.time)-1*mvaOhmVixB.z.resample(mvaAvE.time)-mvaOhmGradPe.z.resample(mvaAvE.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_legend(hca,{'E','jxB/ne','-v_{i}xB','-\nabla \cdot P_e/ne','jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[0.98 0.98],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc av',ic),[0.02 0.95],'color','k','fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 0 % JN
  hca = irf_panel('JL'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h1(:),'x',tintZoom)
irf_zoom(h1(:),'y')
irf_plot_axis_align
%% Example 3 (Ohm's law, 4 sc): Plot 2x sc positions + 4x projection, electrons

tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:34:00.00Z');
%tt = irf.tint('2015-10-16T10:33:26.350Z',0.1);
%tt = irf.tint('2015-10-16T10:33:30.440Z',0.1)+(-0.18);
%tt = irf.tint('2015-10-16T10:33:27.920Z',0.1);
%tt=tt(1);

times = desDist1.tlim(tintZoom).time;


for it = 119;1:times.length  % 644 - ring distributions
  time = times(it);
  tint = time;
  tt = time;


  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,tt.epochUnix','green');

  vectors = {'B',{mvaB1.resample(tt(1)),mvaB2.resample(tt(1)),mvaB3.resample(tt(1)),mvaB4.resample(tt(1))},2;...
    'Ve',{mvaVe1.resample(tt(1)),mvaVe2.resample(tt(1)),mvaVe3.resample(tt(1)),mvaVe4.resample(tt(1))},70};

  isub = 1;
  hca = h2(isub); isub = isub + 1;
  lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 -1 0])
  axis(hca,'square')

  hca = h2(isub); isub = isub + 1;
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 0 1])
  axis(hca,'square')

  %for ic = 1:4;
  %hca = h2(isub); isub = isub + 1;
  %mms.plot_projection(hca,desDist,'tint',times(ii),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  vlim = 15*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [0 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];


  for ic = 1:4
    c_eval('dist = desDist?;',ic)

    c_eval('Ve0 = gseVe?.resample(time).data;',ic);
    hatVe0 = double(irf_norm(Ve0));

    % Get mean magnetic field direction
    c_eval('B0 = gseB?.resample(time).data;',ic);
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = gseE?.resample(time).data;',ic);
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);

    vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

    % Projection coordinate system
    if 1
      x = hatB0;
      y = hatExB0;
      z = cross(x,y);
    else
      x = [1 0 0]; %#ok<UNRCH>
      y = [0 1 0];
      z = [0 0 1];
    end



    %isub = ic;

    if 0 % Plot psd 0 90 180
      hca = h2(isub); isub = isub + 1; %#ok<UNRCH>
      psdtint = tint;%+[-1 1]*0.03;
      c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
      TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
      %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      %hca.Title.String = irf_ssub('C?',ic);
      hca.Title.String = TitleStr;
    end

    if 0 % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;      %#ok<UNRCH>
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};

      % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};
    end

    % Plot project ion onto a plane
    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end

    hca = h2(isub); isub = isub + 1;
    mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
    hca.Title.String = titleStr;
    colormap(hca,strCMap)

    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end
  end
  pause(1)
  %cn.print([irf_ssub('ohm_sc_pos_4sc',1) irf_time(time,'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end

if 0
  hca = h2(3); %#ok<UNRCH>
  tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
  mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
  c_eval('tts? = double((mvaB?.tlim(tintQuiver).time.ttns-mvaB?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
  c_eval('R? = mvaR?-mvaR0;')
  velocity = [0 0 30];
  c_eval('xyz? = irf.ts_vec_xyz(mvaB?.tlim(tintQuiver).time,tts?*velocity);')
  c_eval('quivers? = mvaB?.tlim(tintQuiver);')
  c_eval('positions? = xyz?+R?.resample(xyz?);')
  c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?); hold(hca,''off'');')
  hca.XLabel.String = 'N';
  hca.YLabel.String = '-M';
  hca.ZLabel.String = 'L';

  hca = h2(4);
  tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
  mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
  c_eval('tts? = double((mvaVe?.tlim(tintQuiver).time.ttns-mvaVe?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
  c_eval('R? = mvaR?-mvaR0;')
  velocity = [0 0 30];
  c_eval('xyz? = irf.ts_vec_xyz(mvaVe?.tlim(tintQuiver).time,tts?*velocity);')
  c_eval('quivers? = mvaVe?.tlim(tintQuiver);')
  c_eval('positions? = xyz?+R?.resample(xyz?);')
  c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?,mms_colors(''?'')); hold(hca,''off'');')
  hca.XLabel.String = 'N';
  hca.YLabel.String = '-M';
  hca.ZLabel.String = 'L';
end
%% Example 3 (Ohm's law, 4 sc): Plot 2x sc positions + 4x projection, ions

tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:35.00Z');
%tt = irf.tint('2015-10-16T10:33:26.350Z',0.1);
%tt = irf.tint('2015-10-16T10:33:30.440Z',0.1)+(-0.18);
%tt = irf.tint('2015-10-16T10:33:27.920Z',0.1);
%tt=tt(1);

times = disDist1.tlim(tintZoom).time;


for it = 1:times.length
  time = times(it);
  tint = time;
  tt = time;


  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,tt.epochUnix','green');

  vectors = {'B',{mvaB1.resample(tt(1)),mvaB2.resample(tt(1)),mvaB3.resample(tt(1)),mvaB4.resample(tt(1))},2;...
    'Vi',{mvaVi1.resample(tt(1)),mvaVi2.resample(tt(1)),mvaVi3.resample(tt(1)),mvaVi4.resample(tt(1))},70};

  isub = 1;
  hca = h2(isub); isub = isub + 1;
  lim = 12; xlims = lim*[-1 1]; ylims = lim*[-1 1]; zlims = lim*[-1 1];
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 -1 0])
  axis(hca,'square')

  hca = h2(isub); isub = isub + 1;
  plot_lmn3D(hca,mvaR1,mvaR2,mvaR3,mvaR4,[0 0 1;0 -1 0;1 0 0],{'N','-M','L'},'vectors',vectors)
  hca.XLim = xlims; hca.YLim = ylims; hca.ZLim = zlims;
  view(hca,[0 0 1])
  axis(hca,'square')

  %for ic = 1:4;
  %hca = h2(isub); isub = isub + 1;
  %mms.plot_projection(hca,desDist,'tint',times(ii),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  vlim = 1*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [3 9];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];


  for ic = 1:4
    c_eval('dist = disDist?;',ic)

    c_eval('Ve0 = gseVe?.resample(time).data;',ic);
    hatVe0 = double(irf_norm(Ve0));
    c_eval('Vi0 = gseVi?.resample(time).data;',ic);
    hatVi0 = double(irf_norm(Vi0));

    % Get mean magnetic field direction
    c_eval('B0 = gseB?.resample(time).data;',ic);
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = gseE?.resample(time).data;',ic);
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);

    vectors = {hatB0,'B';hatE0,'E';hatVi0,'V_i';L,'L';M,'M';N,'N'};%0;hatVe0,'V_e'};

    % Projection coordinate system
    if 0
      x = hatB0; %#ok<UNRCH>
      y = hatExB0;
      z = cross(x,y);
    elseif 0 %#ok<IFCDUP>
      x = [1 0 0]; %#ok<UNRCH>
      y = [0 1 0];
      z = [0 0 1];
    else
      x = -N;
      y = L;
      z = -M;
    end



    %isub = ic;

    if 0 % Plot psd 0 90 180
      hca = h2(isub); isub = isub + 1; %#ok<UNRCH>
      psdtint = tint;%+[-1 1]*0.03;
      c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
      TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
      %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      %hca.Title.String = irf_ssub('C?',ic);
      hca.Title.String = TitleStr;
    end

    if 0 % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;       %#ok<UNRCH>
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};

      % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};
    end

    % Plot project ion onto a plane
    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end

    hca = h2(isub); isub = isub + 1;

    vlabel = {'-N','L','M'};
    mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim,'vlabel',vlabel);
    titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
    hca.Title.String = titleStr;
    colormap(hca,strCMap)

    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end
  end
  pause(1)
  cn.print([irf_ssub('ohm_sc_pos_4sc_ions',1) irf_time(time,'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end

if 0
  hca = h2(3); %#ok<UNRCH>
  tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
  mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
  c_eval('tts? = double((mvaB?.tlim(tintQuiver).time.ttns-mvaB?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
  c_eval('R? = mvaR?-mvaR0;')
  velocity = [0 0 30];
  c_eval('xyz? = irf.ts_vec_xyz(mvaB?.tlim(tintQuiver).time,tts?*velocity);')
  c_eval('quivers? = mvaB?.tlim(tintQuiver);')
  c_eval('positions? = xyz?+R?.resample(xyz?);')
  c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?); hold(hca,''off'');')
  hca.XLabel.String = 'N';
  hca.YLabel.String = '-M';
  hca.ZLabel.String = 'L';

  hca = h2(4);
  tintQuiver = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:31.00Z');
  mvaR0 = (mvaR1+mvaR2.resample(mvaR1)+mvaR3.resample(mvaR1)+mvaR4.resample(mvaR1))/4;
  c_eval('tts? = double((mvaVe?.tlim(tintQuiver).time.ttns-mvaVe?.tlim(tintQuiver).time.ttns(1)))*1e-9;')
  c_eval('R? = mvaR?-mvaR0;')
  velocity = [0 0 30];
  c_eval('xyz? = irf.ts_vec_xyz(mvaVe?.tlim(tintQuiver).time,tts?*velocity);')
  c_eval('quivers? = mvaVe?.tlim(tintQuiver);')
  c_eval('positions? = xyz?+R?.resample(xyz?);')
  c_eval('hold(hca,''on''); plot_quivers(hca,quivers?,positions?,mms_colors(''?'')); hold(hca,''off'');')
  hca.XLabel.String = 'N';
  hca.YLabel.String = '-M';
  hca.ZLabel.String = 'L';
end
%%
ic = 1;
tintDist = irf.tint('2015-10-16T10:33:23.20Z/2015-10-16T10:33:34.50Z'); %#ok<NASGU>
c_eval('times = desDist?.tlim(tintDist).time;',ic)


for ii = 100;1:times.length;
  tint = times(ii)+[-1 1]*0.5*0.03*0.5;
  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,tint.epochUnix','green');
  vlim = 15*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [0 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];

  for ic = 1:4
    c_eval('dist = disDist?;',ic)

    c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic);
    hatVe0 = double(irf_norm(Ve0));

    % Get mean magnetic field direction
    c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(tint).time).data);',ic);
    hatB0 = double(irf_norm(B0));
    c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic);
    hatE0 = double(irf_norm(E0));
    hatExB0 = cross(hatE0,hatB0);

    vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e'};%0;hatVe0,'V_e'};

    % Projection coordinate system
    if 1
      x = hatB0;
      y = hatExB0;
      z = cross(x,y);
    else
      x = [1 0 0]; %#ok<UNRCH>
      y = [0 1 0];
      z = [0 0 1];
    end



    isub = ic;

    if 0 % Plot psd 0 90 180
      hca = h2(isub); isub = isub + 1; %#ok<UNRCH>
      psdtint = tint;%+[-1 1]*0.03;
      c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
      TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
      %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
      %hca.Title.String = irf_ssub('C?',ic);
      hca.Title.String = TitleStr;
    end

    if 0 % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;       %#ok<UNRCH>
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};

      % Plot skymap for a given energy
      hca = h2(isub); isub = isub + 1;
      c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
      %hca.Title.String = hca.Title.String{2};
    end

    % Plot project ion onto a plane
    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end

    hca = h2(isub); %isub = isub + 1;
    mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
    hca.Title.String = titleStr;
    colormap(hca,strCMap)

    if 0
      hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
      mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
      hca.Title.String = '';
      colormap(hca,strCMap)
    end
  end
  pause(1)
  %cn.print([irf_ssub('BvnP_1proj_mms1234_gseish',1) irf_time(times(ii),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end

%% Example 5 (Overview plot, 1 sc): Construct all needed data
ic = 1;
% Pitch angle distribution
% load /Users/Cecilia/Data/MMS/2015Oct16/desDistPA1.mat
% load /Users/Cecilia/Data/MMS/2015Oct16/disDistPA1.mat

% c_eval('[desDistPA?,eThetaPA?,eEnergyPA?,eTintPA?] = mms.get_pitchangledist(desDist?,ephi?,etheta?,estepTable?,eenergy0?,eenergy1?,dmpaB?);',ic)
% c_eval('[disDistPA?,iThetaPA?,iEnergyPA?,iTintPA?] = mms.get_pitchangledist(disDist?,iphi?,itheta?,istepTable?,ienergy0?,ienergy1?,dmpaB?);',ic)
%c_eval('eAnisParPerp?.t = desDistPA?.epochUnix; eAnisParPerp?.f = eEnergyPA?; eAnisParPerp?.p = log10(squeeze(desDistPA1.data(:,:,1)'')./squeeze(desDistPA1.data(:,:,end)''));',ic)

ePApap.t = desDistPA1.time.epochUnix;
ePApap.f = mean(iEnergyPA1(1:2,:));
ePApap.p = desDistPA1.data(:,:,1)./desDistPA1.data(:,:,end);
ePApap.f_label = 'E_i (eV)';
ePApap.p_label = 'f_{par}/f_{apar}';

ePApep.t = desDistPA1.time.epochUnix;
ePApep.f = mean(iEnergyPA1(1:2,:));
ePApep.p = (desDistPA1.data(:,:,1)+desDistPA1.data(:,:,end))./mean(desDistPA1.data(:,:,6:7),3)/2;
ePApep.f_label = 'E_e (eV)';
ePApep.p_label = 'f_{par/apar}/f_{perp}';

iPApap.t = disDistPA1.time.epochUnix;
iPApap.f = mean(iEnergyPA1(1:2,:));
iPApap.p = disDistPA1.data(:,:,1)./disDistPA1.data(:,:,end);
iPApap.f_label = 'E_i (eV)';
iPApap.p_label = 'f_{par}/f_{apar}';

iPApep.t = disDistPA1.time.epochUnix;
iPApep.f = mean(iEnergyPA1(1:2,:));
iPApep.p = (disDistPA1.data(:,:,1)+disDistPA1.data(:,:,end))./mean(disDistPA1.data(:,:,6:7),3)/2;
iPApep.f_label = 'E_i (eV)';
iPApep.p_label = 'f_{par/apar}/f_{perp}';
%% Example 5 (Overview plot, 1 sc): Initialize plot
[h1,h2] = initialize_combined_plot(7,2,2,3,'vertical');
%[h1,h2] = initialize_combined_plot(5,2,2,0.4,'vertical');
%% Example 5 (Overview plot, 1 sc): Plot timeseries
tint = irf.tint('2015-10-16T10:32:50.00Z/2015-10-16T10:34:10.00Z');
%tintZoom = irf.tint('2015-10-16T10:33:26.00Z/2015-10-16T10:33:31.00Z');
ic = 1;

if 1 % B, LMN
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
end

if 1 % Electrons, omni
  hca = irf_panel('Electrons Omni');
  c_eval('irf_plot(hca,eDEFomni?);',ic)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 1 % Electrons, anis
  hca = irf_panel('Electron anisotropy a/par/perp');
  c_eval('irf_plot(hca,ePApep);',ic)
  hca.CLim = 2*[-1 1];
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  colormap(hca,cn.cmap('bluered3'))
end
if 1 % Electrons, anis
  hca = irf_panel('Electron anisotropy par/apar');
  c_eval('irf_plot(hca,ePApap);',ic)
  hca.CLim = 2*[-1 1];
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  colormap(hca,cn.cmap('bluered3'))
end
if 1 % Ions, omni
  hca = irf_panel('Ions Omni');
  c_eval('irf_plot(hca,iDEFomni?);',ic)
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 1 % Ions, anis
  hca = irf_panel('Ion anisotropy a/par/perp');
  c_eval('irf_plot(hca,iPApep);',ic)
  hca.CLim = 2*[-1 1];
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  colormap(hca,cn.cmap('bluered3'))
end
if 1 % Ions, anis
  hca = irf_panel('Ion anisotropy par/apar');
  c_eval('irf_plot(hca,iPApap);',ic)
  hca.CLim = 2*[-1 1];
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  colormap(hca,cn.cmap('bluered3'))
end
if 0 % Anisotropy, electrons
  hca = irf_panel('B'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'B_L','B_M','B_N'},[0.98 0.9],'fontsize',12);
end

if 0 % ne
  hca = irf_panel('ne'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors(irf_ssub('?',ic)))
  c_eval('irf_plot(hca,ne?);',ic)
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % Ve LMN
  hca = irf_panel('Ve'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaVe?.x,mvaVe?.y,mvaVe?.z},''comp'');',ic)
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'v_L','v_M','v_N'},[0.98 0.9],'fontsize',12);
end
if 0 % E
  hca = irf_panel('E'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'E_L','E_M','E_N'},[0.98 0.9],'fontsize',12);
end
if 0 % J, comparing with curl and 4 sc av, 3 panels
  hca = irf_panel('JL'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('1234ab'))
  c_eval('lines = irf_plot(hca,{mvaJ?.x,mvaAvJ.x,mvaJcurl.x},''comp'');',ic)
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  irf_legend(hca,{irf_ssub('J_?',ic),'J_{av,fpi}','J_{curl}'},[0.02 0.9],'fontsize',12);

  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  c_eval('lines = irf_plot(hca,{mvaJ?.y,mvaAvJ.y,mvaJcurl.y},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  irf_legend(hca,{irf_ssub('J_?',ic),'J_{av,fpi}','J_{curl}'},[0.02 0.9],'fontsize',12);

  hca = irf_panel('JN');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  c_eval('lines = irf_plot(hca,{mvaJ?.z,mvaAvJ.z,mvaJcurl.z},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  irf_legend(hca,{irf_ssub('J_?',ic),'J_{av,fpi}','J_{curl}'},[0.02 0.9],'fontsize',12);
end
if 0 % J
  hca = irf_panel('J'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('lines = irf_plot(hca,mvaJ?);',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_L','J_M','J_N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1200 2000];
end
if 0 % P
  hca = irf_panel('Pe'); %#ok<UNRCH>
  %set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('lines = irf_plot(hca,{mvaPe?.xx,mvaPe?.yy,mvaPe?.zz,mvaPe?.xy,mvaPe?.xz,mvaPe?.yz},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/km^2)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'P_{LL}','P_{MM}','P_{NN}','P_{LM}','P_{LN}','P_{MN}'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyzb'))
  c_eval('irf_plot(hca,{mvaE?.z,-1*mvaOhmVexB?.z,-1*mvaOhmGradPe?.z,-1*mvaOhmVexB?.z.resample(mvaE?.time)-mvaOhmGradPe?.z.resample(mvaE?.time)},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzb'))
  irf_legend(hca,{'E','-v_{e}xB','-\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,['v = ' num2str(velocity,'%.0f') ' km/s'],[0.98 0.95],'color',[0 0 0],'fontsize',12);
end
if 0 % Ohms law ions: N
  hca = irf_panel('Ohms law ions: N'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyazb'))
  c_eval('irf_plot(hca,{mvaE?.z,mvaOhmJxB?.z,-1*mvaOhmVixB?.z,-1*mvaOhmGradPe?.z,mvaOhmJxB?.z.resample(mvaE?.time)-1*mvaOhmVixB?.z.resample(mvaE?.time)-mvaOhmGradPe?.z.resample(mvaE?.time)},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyazb'))
  irf_legend(hca,{'E','jxB/ne','-v_{i}xB','-\nabla \cdot P_e/ne','jxB/ne-v_{i}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,['v = ' num2str(velocity,'%.0f') ' km/s'],[0.98 0.95],'color',[0 0 0],'fontsize',12);
end
if 0 % Ohms law, E+vexB: N
  hca = irf_panel('Ohms law 2: N'); %#ok<UNRCH>
  set(hca,'ColorOrder',mms_colors('xyzba'))
  c_eval('irf_plot(hca,{mvaE?.z+1*mvaOhmVexB?.resample(mvaE?.time).z,-1*mvaOhmGradPe?.z},''comp'');',ic)
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))
  irf_legend(hca,{'E+v_{e}xB','-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,['v = ' num2str(velocity,'%.0f') ' km/s'],[0.98 0.95],'color',[0 0 0],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end

h1(1).Title.String = irf_ssub('MMS ?',ic);
irf_zoom(h1(:),'x',tint)
%irf_zoom(h1(:),'y')
irf_plot_axis_align
%% Example 5 (Overview plot, 1 sc): Plot particle distributions
ic = 1;
tintDist = irf.tint('2015-10-16T10:33:23.20Z/2015-10-16T10:33:34.50Z');
c_eval('times = desDist?.tlim(tintDist).time;',ic)
tint = irf.tint('2015-10-16T10:33:29.50Z',0.3)+2.4*[1 1];



for ii = 1%:tint.length;
  %tint = times(ii)+[-1 1]*0.5*0.03*0.5;
  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,tint.epochUnix','green');

  isub = 1;
  hca = h2(isub); isub = isub + 1;
  c_eval('PA = squeeze(mean(disDistPA?.tlim(tint).data,1))*1e30;',ic)
  c_eval('energy = mean(iEnergyPA?(1:2,:),1);',ic)
  loglog(hca,energy,mean(PA(:,1:2),2),energy,mean(PA(:,6:7),2),energy,mean(PA(:,11:12),2))

  if 1 % Plot psd 0 90 180
    hca = h2(isub); isub = isub + 1;
    psdtint = tint;
    c_eval('mms.plot_pitchangles(hca,disDist?,dmpaB?,''tint'',psdtint,''xlim'',h2(1).XLim);',ic)
    TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
    %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];
    %hca.Title.String = irf_ssub('C?',ic);
    hca.Title.String = TitleStr;
  end
  if 0 % Plot skymap for a given energy
    hca = h2(isub); isub = isub + 1;       %#ok<UNRCH>
    c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
    %hca.Title.String = hca.Title.String{2};

    % Plot skymap for a given energy
    hca = h2(isub); isub = isub + 1;
    c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
    %hca.Title.String = hca.Title.String{2};
  end
  if 0 % Plot project ion onto a plane
    hca = h2(isub); isub = isub + 1;  %#ok<UNRCH>
    mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
    hca.Title.String = '';
    colormap(hca,strCMap)
  end


  titleStr = {irf_ssub('MMS ?',ic),[irf_time(tint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};
  hca.Title.String = titleStr;
  colormap(hca,strCMap)

end
pause(1)
%cn.print([irf_ssub('BvnP_1proj_mms1234_gseish',1) irf_time(times(ii),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
