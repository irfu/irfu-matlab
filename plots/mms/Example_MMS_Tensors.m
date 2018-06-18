% Shows implementation of TSeries tensor.

%% Create tensor of artifical data
utcT1 = '2002-03-04T09:30:00Z'; 
time = EpochTT(utcT1):1:(EpochTT(utcT1)+10);
data = repmat(magic(3),1,1,time.length); data = permute(data,[3 1 2]);
T = irf.ts_tensor_xyz(time,data);
T.name = 'tensor';

%% Perform operations and take out components
T.xy
T.zz
T.trace
T.tlim(tint);

%% Make tensor of real data using mms.psd_moments
tint = irf.tint('2015-10-16T10:33:10.00Z/2015-10-16T10:34:10.00Z'); %#ok<NASGU>
db_info = datastore('mms_db');
ic = 1;

% Load electron distributions with depend data
c_eval('tmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l1b/des-dist/2015/10/16/mms?_fpi_brst_l1b_des-dist_20151016103254_v1.1.0.cdf'']);',ic);
c_eval('pdist? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_brstSkyMap_dist''));',ic);
c_eval('energy0? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_energy0'');',ic);
c_eval('energy1? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_energy1'');',ic);
c_eval('phi? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_brstSkyMap_phi''));',ic);
c_eval('theta? = get_variable(tmpDataObj?,''mms?_des_brstSkyMap_theta'');',ic);
c_eval('stepTable? = mms.variable2ts(get_variable(tmpDataObj?,''mms?_des_stepTable_parity''));',ic);

% Load spacecraft potential
c_eval('P?brst=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',tint);',ic);

% Calculate moments
c_eval('emoments? = mms.psd_moments(pdist?,phi?,theta?,stepTable?,energy0?,energy1?,P?brst,''electron'');',ic);

% Pressure
c_eval('Pet? = emoments?.P_psd.tlim(tint);',ic);
c_eval('Pdata? = nan(Pet?.length,3,3);',ic)
c_eval(['Pdata?(:,1,1) = Pet?.data(:,1);',...
        'Pdata?(:,1,2) = Pet?.data(:,2);',...
        'Pdata?(:,1,3) = Pet?.data(:,3);',...
        'Pdata?(:,2,1) = Pet?.data(:,2);',...
        'Pdata?(:,2,2) = Pet?.data(:,4);',...
        'Pdata?(:,2,3) = Pet?.data(:,5);',...
        'Pdata?(:,3,1) = Pet?.data(:,3);',...
        'Pdata?(:,3,2) = Pet?.data(:,5);',...
        'Pdata?(:,3,3) = Pet?.data(:,6);',...
        ],ic);
c_eval('P? = irf.ts_tensor_xyz(Pet?.time,Pdata?);',ic)

% Temperature
c_eval('Tet? = emoments?.T_psd.tlim(tint);',ic);
c_eval('Tdata? = nan(Tet?.length,3,3);',ic)
c_eval(['Tdata?(:,1,1) = Tet?.data(:,1);',...
        'Tdata?(:,1,2) = Tet?.data(:,2);',...
        'Tdata?(:,1,3) = Tet?.data(:,3);',...
        'Tdata?(:,2,1) = Tet?.data(:,2);',...
        'Tdata?(:,2,2) = Tet?.data(:,4);',...
        'Tdata?(:,2,3) = Tet?.data(:,5);',...
        'Tdata?(:,3,1) = Tet?.data(:,3);',...
        'Tdata?(:,3,2) = Tet?.data(:,5);',...
        'Tdata?(:,3,3) = Tet?.data(:,6);',...
        ],ic);
c_eval('T? = irf.ts_tensor_xyz(Tet?.time,Tdata?);',ic)

%% Make tensor of real data using mms.get_data
tint = irf.tint('2015-10-16T10:33:10.00Z/2015-10-16T10:34:10.00Z');
brstPe1 = mms.get_data('Pe_fpi_brst',tint,1); % c_eval('brstPe? = mms.get_data(''Pe_fpi_brst'',tint,?);')
brstTe1 = mms.get_data('Te_fpi_brst',tint,1);

%% Rotate tensors
% need to load gseB?
ic = 1;
c_eval('facTe? = mms.rotate_tensor(brstTe?,''fac'',gseB?); facTe? = irf.ts_tensor_xyz(facTe?.time,facTe?.data);',ic)
c_eval('hca = irf_plot({facTe?.xx,facTe?.yy,facTe?.zz,facTe?.xy,facTe?.xz,facTe?.yz},''comp'');',ic)
irf_legend(hca,{'||','\perp_1','\perp_2','xy','xz','yz'},[0.98 0.95])

%% Plot pressure
% irf_plot can not yet handle TSeries with tensorOrder = 2, so its
% necessary to pick out the components and use 'comp'
% c_eval('P? = P?fpi;') % pick out if you want to plot fpi or irfu moments
% c_eval('P? = brstPe?;') % pick out if you want to plot fpi or irfu moments

h = irf_plot(6);
isub = 1;
hca = h(isub); isub = isub + 1;
irf_plot(hca,{P1.trace/3,P2.trace/3,P3.trace/3,P4.trace/3},'comp')
ylabel(hca,{'trace(P)/3','(nPa)'})
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.95])

hca = h(isub); isub = isub + 1;
irf_plot(hca,{P1.xx,P2.xx,P3.xx,P4.xx},'comp')
ylabel(hca,{'P_{xx}','(nPa)'})
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.95])

for ic = [1 2 3 4] % Pxx Pyy Pzz
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{P?.xx,P?.yy,P?.zz,P?.xy,P?.xz,P?.yz},''comp'')',ic)
  irf_legend(hca,{'xx','yy','zz','xy','xz','yz'},[0.98 0.95])
  ylabel(hca,{irf_ssub('P?',ic),'(nPa)'})
end

tintZoom = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:34:00.00Z');
irf_zoom(h,'x',tintZoom)

%% Plot temperature
% irf_plot can not yet handle TSeries with tensorOrder = 2, so its
% necessary to pick out the components
% c_eval('T? = T?fpi;') % pick out if you want to plot fpi or irfu moments

h = irf_plot(6);
isub = 1;
hca = h(isub); isub = isub + 1;
irf_plot(hca,{T1.trace/3,T2.trace/3,T3.trace/3,T4.trace/3},'comp')
ylabel(hca,{'trace(T)/3','(eV)'})
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.95])

hca = h(isub); isub = isub + 1;
irf_plot(hca,{T1.xx,T2.xx,T3.xx,T4.xx},'comp')
ylabel(hca,{'T_{xx}','(eV)'})
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.95])

for ic = [1 2 3 4] % Txx Tyy Tzz
  hca = h(isub); isub = isub + 1;
  c_eval('irf_plot(hca,{T?.xx,T?.yy,T?.zz},''comp'')',ic)
  ylabel(hca,{irf_ssub('T?',ic),'(eV)'})
end

tintZoom = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:34:00.00Z');
irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')

