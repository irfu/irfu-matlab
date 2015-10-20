% MMS.SKYMAP_EXAMPLE
%   Loads and plots a skymap of ion and electron phase space density as 
%   well as B and electron and ion moments directions, for a given time
%   interval.

%% Load data
% Find files which include the start of this time interval and load them
tint = irf.tint('2015-08-15T12:59:40',1);
sc = 3;

% Magnetic field survey data
c_eval('[B?fg,dobjB] = mms.cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint(1));',sc)

% Electric field burst data
c_eval('[dslE?,dobjE] = mms.cn_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint(1));',sc)

% Ion skymap
c_eval('[disDist?,dobjDisDist] = mms.cn_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint(1));',sc)

% Electron skymap
c_eval('[desDist?,dobjDesDist] = mms.cn_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint(1));',sc)

% Electron moments
c_eval('[~,dobj] = mms.cn_get_ts(''mms?_fpi_brst_l1b_des-moms'',[],tint(1));',sc)
dobjDesMoms = dobj;

% Electron density
c_eval('ne? = mms.variable2ts(get_variable(dobj,''mms?_des_numberDensity''));',sc)

% Electron pressure
c_eval('peXX = mms.variable2ts(get_variable(dobj,''mms?_des_PresXX'')); ',sc)
c_eval('peYY = mms.variable2ts(get_variable(dobj,''mms?_des_PresYY'')); ',sc)
c_eval('peZZ = mms.variable2ts(get_variable(dobj,''mms?_des_PresZZ'')); ',sc)
c_eval('pe? = irf.ts_scalar(peXX.time,(peXX.data + peYY.data + peZZ.data)/3);',sc)
c_eval('pe?.units = peXX.units;',sc)
c_eval('pe?.userData = peXX.userData;',sc)

% Electron temperature
c_eval('TeXX = mms.variable2ts(get_variable(dobj,''mms?_des_TempXX''));',sc)
c_eval('TeYY = mms.variable2ts(get_variable(dobj,''mms?_des_TempYY''));',sc)
c_eval('TeZZ = mms.variable2ts(get_variable(dobj,''mms?_des_TempZZ''));',sc)
c_eval('Te? = irf.ts_scalar(TeXX.time,(TeXX.data + TeYY.data + TeZZ.data)/3);',sc)
c_eval('Te?.units = TeXX.units;',sc)
c_eval('Te?.userData = TeXX.userData;',sc)

% Electron velocity
c_eval('veX = mms.variable2ts(get_variable(dobj,''mms?_des_bulkX''));',sc)
c_eval('veY = mms.variable2ts(get_variable(dobj,''mms?_des_bulkY''));',sc)
c_eval('veZ = mms.variable2ts(get_variable(dobj,''mms?_des_bulkZ''));',sc)
c_eval('ve? = irf.ts_vec_xyz(veX.time,[veX.data veY.data veZ.data]);',sc)
c_eval('ve?.units = veX.units;',sc)
c_eval('ve?.userData = veX.userData;',sc)

% Downsample electron moments
c_eval('fs = 1/(ne?.time(2)-ne?.time(1));',sc)
fny = fs/2;
c_eval('pe?_lowres = irf_filt(pe?,0,fny/2,fs,5);',sc)
c_eval('Te?_lowres = irf_filt(Te?,0,fny/2,fs,5);',sc)
c_eval('ne?_lowres = irf_filt(ne?,0,fny/2,fs,5);',sc)

% Ion moments
c_eval('[~,dobj] = mms.cn_get_ts(''mms?_fpi_brst_l1b_dis-moms'',[],tint(1));',sc)
dobjDisMoms = dobj;

% Ion density
c_eval('ni? = mms.variable2ts(get_variable(dobj,''mms?_dis_numberDensity''));',sc)

% Ion pressure 
c_eval('piXX = mms.variable2ts(get_variable(dobj,''mms?_dis_PresXX'')); ',sc)
c_eval('piYY = mms.variable2ts(get_variable(dobj,''mms?_dis_PresYY'')); ',sc)
c_eval('piZZ = mms.variable2ts(get_variable(dobj,''mms?_dis_PresZZ'')); ',sc)
c_eval('pi? = irf.ts_scalar(piXX.time,(piXX.data + piYY.data + piZZ.data)/3);',sc)
c_eval('pi?.units = piXX.units;',sc)
c_eval('pi?.userData = piXX.userData;',sc)

% Ion temperature    
c_eval('TiXX = mms.variable2ts(get_variable(dobj,''mms?_dis_TempXX''));',sc)
c_eval('TiYY = mms.variable2ts(get_variable(dobj,''mms?_dis_TempYY''));',sc)
c_eval('TiZZ = mms.variable2ts(get_variable(dobj,''mms?_dis_TempZZ''));',sc)
c_eval('Ti? = irf.ts_scalar(TiXX.time,(TiXX.data + TiYY.data + TiZZ.data)/3);',sc)
c_eval('Ti?.units = TiXX.units;',sc)
c_eval('Ti?.userData = TiXX.userData;',sc)

% Ion velocity
c_eval('viX = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkX''));',sc)
c_eval('viY = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkY''));',sc)
c_eval('viZ = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkZ''));',sc)
c_eval('vi? = irf.ts_vec_xyz(viX.time,[viX.data viY.data viZ.data]);',sc)
c_eval('vi?.units = viX.units;',sc)
c_eval('vi?.userData = viX.userData;',sc)

% Downsample ion moments
c_eval('fs = 1/(ni?.time(2)-ni?.time(1));',sc)
fny = fs/2;
c_eval('pi?_lowres = irf_filt(pi?,0,fny/2,fs,5);',sc)
c_eval('Ti?_lowres = irf_filt(Ti?,0,fny/2,fs,5);',sc)
c_eval('ni?_lowres = irf_filt(ni?,0,fny/2,fs,5);',sc)

%% Plot a skymap fxor one time and several energy levels with vph and B0 inserted
% For this part vi?, ve?, B?fg, dslE?, desDist?, disDist? are needed
% Choose different time interval
% tint = irf.tint('2015-08-15T12:55:20.5',0.6);
% tint = irf.tint('2015-08-15T12:57:16.6',0.3);

% Get vectors for a small time interval
% I didn't find a mean function for TSeries, like 
% B?fg.tlim(tint).data.mean or B?fg.tlim(tint).mean
c_eval('B0 = mean(B?fg.tlim(tint).data);',sc); hatB0 = double(irf_norm(B0));
c_eval('E0 = mean(dslE?.tlim(tint).data);',sc); hatE0 = double(irf_norm(E0));
c_eval('vi = mean(vi?.tlim(tint).data);',sc); hatVi = double(irf_norm(vi));
c_eval('ve = mean(ve?.tlim(tint).data);',sc); hatVe = double(irf_norm(ve));
c_eval('vExB? = cross(dslE?,B?fg.resample(dslE?)); vExB = mean(vExB?.tlim(tint).data);',sc);  hatExB = double(irf_norm(vExB));

% Set up coordinates for skymap plot
c_eval('desDist = desDist?;',sc);
c_eval('disDist = disDist?;',sc);

% Set up coordinate system separately for ions and electrons if they for
% some reason would be different from each other.
r = 1; % radius of sphere
% Electrons
phi_edges = linspace(0,2*pi,size(desDist.data,3)+1); % azimuthal angle bin edges?
theta_edges = linspace(0,pi,size(desDist.data,4)+1); % polar angle bin edges?
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
eX = -r*sin(THETA).*cos(PHI); % '-' because the data shows which direction the particles where coming from
eY = -r*sin(THETA).*sin(PHI);
eZ = -r*cos(THETA);
% Ions
phi_edges = linspace(0,2*pi,size(disDist.data,3)+1); % azimuthal angle bin edges?
theta_edges = linspace(0,pi,size(disDist.data,4)+1); % polar angle bin edges?
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
iX = -r*sin(THETA).*cos(PHI); % '-' because the data shows which direction the particles where coming from
iY = -r*sin(THETA).*sin(PHI);
iZ = -r*cos(THETA);

% Plot a few energy levels on a skymap for both ion and electron data
nrows = 3;
ncols = 4;
nSubPlots = nrows*ncols;
for k = 1:nSubPlots; h(k) = subplot(nrows,ncols,k); end

% Choose which energy levels to plot
ionEnergyLevels = [10 15 20 25]; % between 1:32
electronEnergyLevels = [5 9 12 15 ]; % between 1:32

[etId,~] = desDist.time.tlim(tint); 
[itId,~] = disDist.time.tlim(tint); 

for k = 1:ncols
    % Plot ion and electron moments
    hca = subplot(nrows,ncols,1:4);
    axes(hca);
    veline = irf_plot(hca,ve3); hold(hca,'on'); 
    viline = irf_plot(hca,vi3); hold(hca,'off');
    %irf_legend(hca,{'v_{ex}','v_{ey}','v_{ez}','v_{ix}','v_{iy}','v_{iz}'},[0.9 0.9])
    legend([veline(:); viline(:)],{'v_{ex}','v_{ey}','v_{ez}','v_{ix}','v_{iy}','v_{iz}'},'location','EastOutside')
    irf_pl_mark(hca,tint.epochUnix','red')
    irf_zoom(hca,'x',[vi3.time.start.epochUnix vi3.time.stop.epochUnix])
    title(hca,'Ion and electron velocity moments')
    
    % Plot electron skymap data in second row
    hca = h(ncols*1+k);    
    axes(hca)
    C = squeeze(nanmean(desDist.data(etId,electronEnergyLevels(k),:,:),1))';
    hs = surf(hca,eX,eY,eZ,C);
    hc = colorbar('peer',hca);
    axis(hca,'square')
    axis(hca,'equal')
    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    %titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(electronEnergyLevels(k))]};
    titleString = {[tint(1).utc ' + ' num2str(tint.stop-tint.start) ' s'],['Energy level = ' num2str(electronEnergyLevels(k))]};
    hca.Title.String = titleString;   
    hc.YLabel.String = 'Electron phase spase density';
    shading flat;
    
    % Plot vectors
    hold(hca,'on');
    hold on
    scale = 1.5;
    v_ExB_hat = double(irf_norm(vExB));
    quiver3(hca,0,0,0,hatExB(1),hatExB(2),hatExB(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVi(1),hatVi(2),hatVi(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVe(1),hatVe(2),hatVe(3),scale,'linewidth',2)
    quiver3(hca,-scale*hatB0(1),-scale*hatB0(2),-scale*hatB0(3),hatB0(1),hatB0(2),hatB0(3),2*scale,'linewidth',2)%quiver3(hca,0,0,0,hatB0(1),hatB0(2),hatB0(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatE0(1),hatE0(2),hatE0(3),scale,'linewidth',2)
    
    % Label vectors
    scale = 1.7;
    text(scale*hatExB(1),scale*hatExB(2),scale*hatExB(3),'v_{ExB}','fontsize',14)
    text(scale*hatVi(1),scale*hatVi(2),scale*hatVi(3),'v_{i}','fontsize',14)
    text(scale*hatVe(1),scale*hatVe(2),scale*hatVe(3),'v_{e}','fontsize',14)
    text(scale*hatB0(1),scale*hatB0(2),scale*hatB0(3),'B_{0}','fontsize',14)
    text(scale*hatE0(1),scale*hatE0(2),scale*hatE0(3),'E_{0}','fontsize',14)
    hold(hca,'off');
    hold off
    
    % Plot ion skymap data in bottom row
    hca = h(ncols*2+k);    
    axes(hca)
    C = squeeze(nanmean(disDist.data(itId,ionEnergyLevels(k),:,:),1))';
    hs = surf(hca,iX,iY,iZ,C);
    hc = colorbar('peer',hca);
    axis(hca,'square')
    axis(hca,'equal')
    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    %titleString = {tint(1).utc,tint(2).utc,['Energy level = ' num2str(ionEnergyLevels(k))]};
    titleString = {[tint(1).utc ' + ' num2str(tint.stop-tint.start) ' s'],['Energy level = ' num2str(ionEnergyLevels(k))]};
    hca.Title.String = titleString;   
    hc.YLabel.String = 'Ion phase spase density';
    shading flat;
    
    % Plot vectors
    hold(hca,'on');
    hold on
    scale = 1.5;
    v_ExB_hat = double(irf_norm(vExB));
    quiver3(hca,0,0,0,hatExB(1),hatExB(2),hatExB(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVi(1),hatVi(2),hatVi(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatVe(1),hatVe(2),hatVe(3),scale,'linewidth',2)        
    quiver3(hca,-scale*hatB0(1),-scale*hatB0(2),-scale*hatB0(3),hatB0(1),hatB0(2),hatB0(3),2*scale,'linewidth',2) %quiver3(hca,0,0,0,hatB0(1),hatB0(2),hatB0(3),scale,'linewidth',2)
    quiver3(hca,0,0,0,hatE0(1),hatE0(2),hatE0(3),scale,'linewidth',2)
    
    % Label vectors
    scale = 1.7;
    text(scale*hatExB(1),scale*hatExB(2),scale*hatExB(3),'v_{ExB}','fontsize',14)
    text(scale*hatVi(1),scale*hatVi(2),scale*hatVi(3),'v_{i}','fontsize',14)
    text(scale*hatVe(1),scale*hatVe(2),scale*hatVe(3),'v_{e}','fontsize',14)
    text(scale*hatB0(1),scale*hatB0(2),scale*hatB0(3),'B_{0}','fontsize',14)
    text(scale*hatE0(1),scale*hatE0(2),scale*hatE0(3),'E_{0}','fontsize',14)
    hold(hca,'off');
    hold off
end
% Set  viewing angles
for ii = 9:12; view(h(ii),hatB0); end
for ii = 5:8; view(h(ii),hatB0); end


