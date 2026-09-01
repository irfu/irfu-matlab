%% Test function of PDist.reduce

%% Load data
mms.db_init('local_file_db','/Volumes/mms');

tint = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:34:30.00Z');

ic = 3;
c_eval('dmpaB = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
ne = mms.get_data('Ne_fpi_brst_l2',tint,ic);
c_eval('dbcsVe = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('scPot = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
ePDist = mms.get_data('PDe_fpi_brst_l2',tint,ic);
iPDist = mms.get_data('PDi_fpi_brst_l2',tint,ic);
iPDistErr = mms.get_data('PDERRi_fpi_brst_l2',tint,ic);
iPDist_counts = iPDist; iPDist_counts.data = (iPDist_counts.data./iPDistErr.data).^2;


%% Compare reduced distribution to model Maxwellian PDist
units = irf_units;

n = 1; % cc
Ti = 5000; % eV
Te = 1000; % eV
vd = [1 0 0]; % km/s, breaks for exactly zero
B = [1 0.01 0];
m = units.mp;
debug_mult_factor = 1;%pi*1e12;

w = @(T,m) sqrt(2*units.eV*T/(m)); % m/s
f1D_max = @(v_kms, T_eV, n_cc, vd_kms, m) n_cc*1e6/((pi)^(1/2)*w(T_eV,m).^1)*exp(-(v_kms*1e3-vd_kms*1e3).^2./w(T_eV,m)./w(T_eV,m));
f2D_max = @(v_kms, T_eV, n_cc, vd_kms, m) n_cc*1e6/((pi)^(2/2)*w(T_eV,m).^2)*exp(-(v_kms*1e3-vd_kms*1e3).^2./w(T_eV,m)./w(T_eV,m));
f3D_max = @(v_kms, T_eV, n_cc, vd_kms, m) n_cc*1e6/((pi)^(3/2)*w(T_eV,m).^3)*exp(-(v_kms*1e3-vd_kms*1e3).^2./w(T_eV,m)./w(T_eV,m));


nMC = 1000;

it = 1;
pdist = iPDist(it);
time = pdist.time;
tsB = irf.ts_vec_xyz(time,B);
tsVsc = irf.ts_scalar(time,0);
tsN = irf.ts_scalar(time,n);
tsV = irf.ts_vec_xyz(time,vd);
tsT = irf.ts_tensor_xyz(time,[1 0 0; 0 1 0; 0 0 1]*(Ti));

i_model_dist = mms.make_model_dist(pdist,tsB,tsVsc,tsN,tsV,tsT);

f2D = i_model_dist.reduce('2D',[1 0 0],[0 1 0],'nMC',nMC);
f1D = i_model_dist.reduce('1D',[1 0 0],'nMC',nMC);

nRows = 3; nCols = 1;
ip = 0; h = gobjects(0);
for iRow = 1:nRows
  for iCol = 1:nCols
    ip = ip + 1;
    h(ip) = subplot(nRows,nCols,ip);
  end
end

isub = 1;

hca = h(isub); isub = isub + 1;
f2D.plot_plane(hca)
hca.Title.String = sprintf('n_{MC} = %g', nMC);

hca = h(isub); isub = isub + 1;
plot(hca, f1D.depend{1}, f1D_max(f1D.depend{1},Ti,n,vd(1),units.mp)*debug_mult_factor)
hold(hca,'on')
plot(hca, f1D.depend{1}, f1D.data)
hold(hca,'off')
hca.XLabel.String = 'v';
hca.YLabel.String = 'f (s/m^4)';
legend(hca,{'Reduced PDist','Maxwellian'})


hca = h(isub); isub = isub + 1;
plotyy(hca, f1D.depend{1}, f1D_max(f1D.depend{1},Ti,n,vd(1),units.mp)*debug_mult_factor, f1D.depend{1}, f1D.data)
hca.XLabel.String = 'v';
hca.YLabel.String = 'f (s/m^4)';
legend(hca,{'Reduced PDist','Maxwellian'})
hca.Title.String = 'Double Y-scale for debugging';

h = findobj(gcf,'type','axes');
c_eval('h(?).FontSize = 20;',1:numel(h))
c_eval('h(?).LineWidth = 2;',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))


%% Effect of initialization of particles
units = irf_units;

E1 = 1000; % eV
E2 = 3000; % eV

f_v = @(E,m) sqrt(2*units.eV*E/m);
f_E = @(v,m) m*v.^2/2/units.eV;

v1 = f_v(E1,units.mp);
v2 = f_v(E2,units.mp);

nMP = 100000;

v = v1 + (v2-v1)*rand(nMP,1);
E = E1 + (E2-E1)*rand(nMP,1);

v_from_E = f_v(E,units.mp);
E_from_v = f_E(v,units.mp);

nRows = 2; nCols = 1;
ip = 0; h = gobjects(0);
for iRow = 1:nRows
  for iCol = 1:nCols
    ip = ip + 1;
    h(ip) = subplot(nRows,nCols,ip);
  end
end

isub = 1;

hca = h(isub); isub = isub + 1;
histogram(hca,v*1e-3,20)
hold(hca,'on')
histogram(hca,v_from_E*1e-3,20)
hold(hca,'off')
hca.XLabel.String = 'v (km/s)';
legend(hca,{'Randomized in v (how int\_sph\_dist does it)','Randomized in E'},'Box','off')
irf_legend(hca,{sprintf('<v> = %.0f km/s',mean(v*1e-3)), sprintf('<v> = %.0f km/s',mean(v_from_E*1e-3))},[0.02 0.98])


hca = h(isub); isub = isub + 1;
histogram(hca,E_from_v,20)
hold(hca,'on')
histogram(hca,E,20)
hold(hca,'off')
hca.XLabel.String = 'E (eV)';
legend(hca,{'Randomized in v (how int\_sph\_dist does it)','Randomized in E'},'Box','off')
irf_legend(hca,{sprintf('<E> = %.0f eV',mean(E)), sprintf('<E> = %.0f eV',mean(E_from_v))},[0.02 0.98])

h(1).Title.String = {sprintf('Initialization of macroparticles, N_{MP} = %g', nMP),'Illustration of spread of macroparticles within a given instrument bin'};
linkprop(h,{'YLim'});


%% Effect of initialization of particles, for FPI energy grid
units = irf_units;

energy = iPDist.depend{1}(1,:);
de_minus = iPDist.ancillary.delta_energy_minus;
de_plus = iPDist.ancillary.delta_energy_plus;
energy_edges = [energy(:,1) - de_minus(:,1) energy + de_plus];
nE = numel(energy);

mass = units.mp;
nMP = 100000;
f_v = @(E,m) sqrt(2*units.eV*E/m);
f_E = @(v,m) m*v.^2/2/units.eV;

v_from_v_mean = zeros(nE,1);
v_from_E_mean = zeros(nE,1);
E_from_E_mean = zeros(nE,1);
E_from_v_mean = zeros(nE,1);

for iE = 1:nE
  E1 = energy(iE) - de_minus(iE); % eV
  E2 = energy(iE) + de_plus(iE); % eV

  v1 = f_v(E1,units.mp);
  v2 = f_v(E2,units.mp);

  % Randomize particle in either v or E
  v_from_v = v1 + (v2-v1)*rand(nMP,1);
  E_from_E = E1 + (E2-E1)*rand(nMP,1);

  % Convert v to E and vice versa
  v_from_E = f_v(E_from_E,mass);
  E_from_v = f_E(v_from_v,mass);

  % Collect averages
  v_from_v_mean(iE) = mean(v_from_v);
  v_from_E_mean(iE) = mean(v_from_E);
  E_from_v_mean(iE) = mean(E_from_v);
  E_from_E_mean(iE) = mean(E_from_E);

end

% Plot results as a function of energy
nRows = 2; nCols = 2;
ip = 0; h = gobjects(0);
for iRow = 1:nRows
  for iCol = 1:nCols
    ip = ip + 1;
    h(ip) = subplot(nRows,nCols,ip);
  end
end

isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,energy,v_from_v_mean,energy,v_from_E_mean)
hca.YLabel.String = '<v> (km/s)';
hca.XLabel.String = 'E (eV)';
legend(hca,{'Randomized in v (how int\_sph\_dist does it)','Randomized in E'},'Box','off')

hca = h(isub); isub = isub + 1;
plot(hca,energy,E_from_v_mean,energy,E_from_E_mean)
hca.YLabel.String = '<E> (eV)';
hca.XLabel.String = 'E (eV)';
legend(hca,{'Randomized in v (how int\_sph\_dist does it)','Randomized in E'},'Box','off')

hca = h(isub); isub = isub + 1;
plot(hca,energy,v_from_v_mean-v_from_E_mean)
hca.YLabel.String = '<v>-<v> (km/s)';
hca.XLabel.String = 'E (eV)';
legend(hca,{'Randomized in v (how int\_sph\_dist does it)','Randomized in E'},'Box','off')

hca = h(isub); isub = isub + 1;
semilogy(hca,energy,E_from_v_mean-E_from_E_mean)
hca.YLabel.String = '<E>-<E> (eV)';
hca.XLabel.String = 'E (eV)';
legend(hca,{'Randomized in v (how int\_sph\_dist does it)','Randomized in E'},'Box','off')


h(1).Title.String = {sprintf('Initialization of macroparticles, N_{MP} = %g', nMP)};

