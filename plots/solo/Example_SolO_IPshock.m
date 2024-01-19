% Script to plot a few data from interplanetary shocks observed by MMS.
% The user is asked to select which event and which reference frame to plot
% the data in.
%
% The script shows how to convert shock related plasma data to
% the normal incidence frame. The same procedure can be used at the bow
% shock.
%
% /AJ

refFrameNum = irf_ask('Which frame? (1: s/c frame, 2: NI frame) [%] > ','refFrameNum',2);
switch refFrameNum
  case 1
    refFrame = 'scf';
  case 2
    refFrame = 'nif';
end

tint = irf.tint('2021-10-30T21:59:00/2021-10-30T22:07:00');
% time for upstream
tintu = irf.tint('2021-10-30T22:01:00/2021-10-30T22:01:30');
% time for downstream
tintd = irf.tint('2021-10-30T22:02:21.375000Z/2021-10-30T22:04:08.625000Z');

% Solar wind energy table for this event. In the NI frame, a new energy
% table is needed.
%Eg = logspace(.3,3.93,32);
%inpE = 1;

% add more events here as they come in
% tintu/tintd are not needed for s/c frame plot

% velocity limit in spacecraft frame
vnlim = [0,1000]; % km/s
vxlim = [-1000,0]; % km/s

% define velocity grid
vn1D = linspace(vnlim(1),vnlim(2),100);
vx1D = linspace(vxlim(1),vxlim(2),100);


% Number of Monte Carlo iterations per bin. Decrease to improve
% performance, increase to improve plot.
nMC = 2e2;


%% Get data (needs a database initialized)
% get ion distribution (mms.get_data is somewhat slow)
% also get errors

iPDist = solo.get_data('pas_vdf',tint);
iPDist = iPDist.(cell2mat(fieldnames(iPDist)));
%iPDistErr = solo.get_data('PDERRi_fpi_brst_l2',tint,ic);

% ignore psd where count is 1
%iPDist.data(iPDist.data<1.1*iPDistErr.data) = 0;

%
B = solo.get_data('B_srf_brst',tint);

%
Ni = solo.get_data('Ni',tint);
Ne = solo.get_data('Nescpot',tint);
Vi = solo.get_data('Vi_srf',tint);
Ti = solo.get_data('Ti',tint);

% sc position (if irfu NAS24 is mounted)
%R = mms.get_data('R_gse',tint);
%c_eval('R? = irf.ts_vec_xyz(R.time,R.gseR?(:,1:3));')

% otherwise replace with
%c_eval('R? = mms.db_get_ts(''mms?_mec_srvy_l2_epht89q'',''mms?_mec_r_gse'',tint);')





%% Set up- and downstream parameters
if strcmp(refFrame,'nif')
  % structure of plasma parameters
  plp = [];
  plp.Bu = double(mean(B.tlim(tintu).data));
  plp.nu = double(mean(Ni.tlim(tintu).data));
  plp.Vu = double(mean(Vi.tlim(tintu).data));
  plp.Bd = double(mean(B.tlim(tintd).data));
  plp.nd = double(mean(Ni.tlim(tintd).data));
  plp.Vd = double(mean(Vi.tlim(tintd).data));

  % get normal and shock speed
  nst = irf_shock_normal(plp);

  % let's use mixed mode 3 as normal
  nvec = nst.n.mx3;

  % and the SB for shock speed
  Vsh = nst.Vsh.mf.mx3;

  % then add info to plp
  plp.nvec = nvec;
  plp.Vsh = Vsh;
  plp.ref_sys = 'nif';

  % set coordinate system
  t2vec = cross(nvec,plp.Bu)/norm(cross(nvec,plp.Bu));
  t1vec = cross(t2vec,nvec);
  % rotation matrix to shock-aligned coordinate system
  Rmat = [nvec;t1vec;t2vec];

  % Get more shock parameters
  shp = irf_shock_parameters(plp);

  % velocity of the NI frame in the SC frame
  Vnif = shp.Vnifu;
end


%% reduce ion distribution
switch refFrame
  case 'nif'
    tic
    % reduced distribution along normal vector
    f1Dn = iPDist.reduce('1D',nvec,'vg',vn1D,'nMC',nMC);
    toc
  case 'scf'
    tic
    % reduced distribution along x
    f1Dx = iPDist.reduce('1D',[1,0,0],'vg',vx1D,'nMC',nMC);
    toc
end


%% convert ion dist to NIF
if strcmp(refFrame,'nif')
  % can be done in other directions too
  f1DnNIF = f1Dn;
  f1DnNIF.depend{1} = f1Dn.depend{1}-dot(Vnif,nvec);
  fOmniNIF = iPDist.omni(Vnif);
end

%% plot the data



% size of figure in cm
figSize = [12,14];


% axis handle array
h = irf_plot(6,'newfigure');

% figure handle
fig = h(1).Parent;


% Set parameters
fig.PaperUnits = 'centimeters';
xSize = figSize(1); ySize = figSize(2);
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
fig.PaperPosition = [xLeft yTop xSize ySize];
fig.Position = [10 10 xSize*50 ySize*50];
fig.PaperPositionMode = 'auto';

switch refFrame
  case 'nif'
    % magnetic field magnitude panel
    hca = irf_panel(h,'Babs');
    irf_plot(hca,B.abs,'linewidth',1.5) % plot magnitude
    ylabel(hca,'$B$ [nT]','interpreter','latex')

    % magnetic field vector panel
    hca = irf_panel(h,'Brtn');
    irf_plot(hca,B*Rmat','linewidth',1.5) % plot vector components
    hold(hca,'on')
    hl = irf_plot(hca,[tint.epochUnix,[0;0]],'--','color',[1,1,1]*.5,'linewidth',1.8);
    uistack(hl,'bottom')
    ylabel(hca,'$\mathbf{B}$ [nT]','interpreter','latex')
    irf_legend(hca,{'$B_n$';'$B_{t1}$';'$B_{t2}$'},[1.02,0.98],'fontsize',18,...
      'interpreter','latex')

    % density panel
    hca = irf_panel(h,'n');
    irf_plot(hca,Ni,'linewidth',1.5)
    hca.ColorOrder = [0,0,0;1,0,0];
    hold(hca,'on')
    irf_plot(hca,Ne-10,'linewidth',1.5)
    ylabel(hca,'$n$ [cm$^{-3}$]','interpreter','latex')
    irf_legend(hca,{'$n_i$';'$n_e$'},[1.02,0.98],'fontsize',18,...
      'interpreter','latex')
    hca.YLim(1) = 0;


    % ion velocity vector panel
    hca = irf_panel(h,'V');
    irf_plot(hca,(Vi-Vnif)*Rmat','linewidth',1.5)
    hold(hca,'on')
    hl = irf_plot(hca,[tint.epochUnix,[0;0]],'--','color',[1,1,1]*.5,'linewidth',1.8);
    uistack(hl,'bottom')
    ylabel(hca,'$\mathbf{V}_i$ [km/s]','interpreter','latex')
    irf_legend(hca,{'$V_n$';'$V_{t1}$';'$V_{t2}$'},[1.02,0.98],'fontsize',18,...
      'interpreter','latex')

    % Plot omni distribution
    hca = irf_panel(h,'omni');
    irf_spectrogram(hca,fOmniNIF.convertto('s^3/m^6').specrec,'donotshowcolorbar');
    hold(hca,'on')
    hcb1 = colorbar(hca);
    hca.YScale = 'log';
    ylabel(hcb1,{'$\log_{10} f_i$','[s$^{3}$\,m$^{-6}$]'},'interpreter','latex','fontsize',15)
    ylabel(hca,'$E_i$ [eV]','interpreter','latex')
    %irf_colormap(hca,'waterfall')
    % hack
    %hca.CLim = [-15.5,-8.5];
    hcb1.Ticks = -100:3:100;

    % Plot reduced distribution
    hca = irf_panel(h,'pdist');
    irf_spectrogram(hca,f1DnNIF.specrec('1D_velocity'),'donotshowcolorbar');
    hold(hca,'on')
    hcb2 = colorbar(hca);
    ylabel(hcb2,{'$\log_{10} F_i$',' [s m$^{-4}$]'},'interpreter','latex','fontsize',15)
    ylabel(hca,'$v_{n}$ [km/s]','interpreter','latex')
    % don't use irf_colormap - it is ridicolously slow
    %irf_colormap(hca,'waterfall')


    %        title(h(1),['\bf{MMS',num2str(ic),' - IP shock in the NI frame}'],...
    %           'interpreter','latex')



  case 'scf'
    % magnetic field magnitude panel
    hca = irf_panel(h,'Babs');
    irf_plot(hca,B.abs,'linewidth',1.5) % plot magnitude
    ylabel(hca,'$B$ [nT]','interpreter','latex')

    % magnetic field vector panel
    hca = irf_panel(h,'Brtn');
    irf_plot(hca,B,'linewidth',1.5) % plot vector components
    hold(hca,'on')
    hl = irf_plot(hca,[tint.epochUnix,[0;0]],'--','color',[1,1,1]*.5,'linewidth',1.8);
    uistack(hl,'bottom')
    ylabel(hca,'$\mathbf{B}$ [nT]','interpreter','latex')
    irf_legend(hca,{'$B_x$';'$B_{y}$';'$B_{z}$'},[1.02,0.98],'fontsize',18,...
      'interpreter','latex')

    % density panel
    hca = irf_panel(h,'n');
    irf_plot(hca,Ni,'linewidth',1.5)
    hca.ColorOrder = [0,0,0;1,0,0];
    hold(hca,'on')
    irf_plot(hca,Ne-10,'linewidth',1.5)
    ylabel(hca,'$n$ [cm$^{-3}$]','interpreter','latex')
    irf_legend(hca,{'$n_i$';'$n_e$'},[1.02,0.98],'fontsize',18,...
      'interpreter','latex')
    hca.YLim(1) = 0;


    % ion velocity vector panel
    hca = irf_panel(h,'V');
    irf_plot(hca,Vi,'linewidth',1.5)
    hold(hca,'on')
    hl = irf_plot(hca,[tint.epochUnix,[0;0]],'--','color',[1,1,1]*.5,'linewidth',1.8);
    uistack(hl,'bottom')
    ylabel(hca,'$\mathbf{V}_i$ [km/s]','interpreter','latex')
    irf_legend(hca,{'$V_x$';'$V_{y}$';'$V_{z}$'},[1.02,0.98],'fontsize',18,...
      'interpreter','latex')

    % Plot omni distribution
    hca = irf_panel(h,'omni');
    irf_spectrogram(hca,iPDist.omni.convertto('s^3/m^6').specrec,'donotshowcolorbar');
    hold(hca,'on')
    hcb1 = colorbar(hca);
    hca.YScale = 'log';
    ylabel(hcb1,{'$\log_{10} f_i$','[s$^{3}$\,m$^{-6}$]'},'interpreter','latex','fontsize',15)
    ylabel(hca,'$E_i$ [eV]','interpreter','latex')
    %irf_colormap(hca,'waterfall')
    % hack
    hcb1.Ticks = -100:3:100;

    % Plot reduced distribution
    hca = irf_panel(h,'pdist');
    irf_spectrogram(hca,f1Dx.specrec('1D_velocity'),'donotshowcolorbar');
    hold(hca,'on')
    hcb2 = colorbar(hca);
    ylabel(hcb2,{'$\log_{10} F_i$',' [s m$^{-4}$]'},'interpreter','latex','fontsize',15)
    ylabel(hca,'$v_{n}$ [km/s]','interpreter','latex')
    % irf_colormap(hca,'waterfall')


    title(h(1),['\bf{MMS',num2str(ic),' - IP shock in the s/c frame}'],...
      'interpreter','latex')

end


% make figure look nice

irf_zoom(h,'x',tint)

irf_plot_axis_align(h);
irf_plot_ylabels_align(h)

for ii = 1:length(h)
  h(ii).FontSize = 18;
  grid(h(ii),'off')

  h(ii).LineWidth = 1.8;
  h(ii).Layer = 'top';

  irf_zoom(h(ii),'y',h(ii).YLim)

  h(ii).YLabel.Units = 'normalized';
  h(ii).YLabel.Position(1) = -.14;

  if ii ~= length(h)
    h(ii).XTickLabel = '';
  else
    h(ii).XTickLabelRotation = 0;
  end
end

% fix colorbars
hcb1.Position(1) = .82;
hcb2.Position(1) = .82;

hcb1.FontSize = 18;
hcb2.FontSize = 18;

hcb1.LineWidth = 1.8;
hcb2.LineWidth = 1.8;

% a bit ugly code
hcb1.Position([2,4]) = h(5).Position([2,4])+[1,-2]*.005;
hcb2.Position([2,4]) = h(6).Position([2,4])+[1,-2]*.005;