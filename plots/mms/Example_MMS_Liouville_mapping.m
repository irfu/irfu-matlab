% A routine to compute and plot reduced electron distributions from FPI
% Then apply Liouville mapping, assuming parallel potential solely
% responsible for electron acceleration.
%
% Compare with: Wilder, F. D., et al. (2016), GRL, 43, 5909?5917,
% doi:10.1002/2016GL069473.
%
% The Example is fairly slow. Approx 4 min.
%
% Originally 'Example_MMS_reduced_ele_dist' written by A. Johlander
% Modified to include Liouville mapping by J. D. White


%% Set parameters and get data
% time interval
% tint = irf.tint('2015-09-19T10:08:14/2015-09-19T10:08:18');
tint = irf.tint('2018-05-05T17:12:18/2018-05-05T17:12:26');

% times to make lines
% t1 = irf.time_array('2015-09-19T10:08:16.275');
% t2 = irf.time_array('2015-09-19T10:08:16.335');
t1 = irf.time_array('2018-05-05T17:12:22.080');
t2 = irf.time_array('2018-05-05T17:12:22.180');

% sc number
ic = 3;

% color/y-limit
clim = 10.^[-4.5,-0.5]; % s m^-4

% define velocity grid
vg = linspace(-45e3,45e3,200); % km/s

% Number of Monte Carlo iterations per bin. Decrease to improve
% performance time, increase to improve plot.
nMC = 1e3;

% velocity limit in plot
vlim = 50e3; % km/s

% define physical constants
units = irf_units;
qe = units.e;
me = units.me;

% get distribution function
c_eval('ePDist = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint));',ic)
ePDist = ePDist.tlim(tint);

% get magnetic field in DMPA (since ePDist is in DMPA)
c_eval('B = mms.get_data(''B_dmpa_fgm_brst_l2'',tint,?);',ic)

% Get the spacecraft potential (in this example, it is not so important)
c_eval('scPot = mms.get_data(''V_edp_brst_l2'',tint,?);',ic)

% remove flux from bottom two energy levels to make it more like in the
% paper
ePDist.data(:,1:2,:,:) = 0;


%% Reduce distribution
tic
% reduced distribution along B
f1D = ePDist.reduce('1D',B,'vg',vg,'nMC',nMC,'scpot',scPot);
toc

%% Plot reduced distribution as a time series
% make figure
h = irf_plot(2,'newfigure');

% plot magnetic field
hca = irf_panel(h,'Bxyz');
irf_plot(hca,B)
hca.YLabel.Interpreter = 'tex';hline = findobj(gcf, 'type', 'line');
ylabel(hca,'B_{dmpa} [nT]')
irf_legend(hca,{'Bx';'By';'Bz'},[1.02,0.9])

% Plot reduced distribution
hca = irf_panel(h,'pdist');
[~,hcb] = irf_spectrogram(hca,f1D.specrec('1D_velocity'));
hcb.Label.String = 'log_{10}F_e [s m^{-4}]';
ylabel(hca,'V_{||} [km/s]')
colormap('jet')
irf_zoom(hca,'y',[min(vg),max(vg)])
irf_plot_axis_align(h)
irf_zoom(h,'x',tint)

% plot yellow over ROI, enclosed by dashed lines
irf_pl_mark(h,[t1 t2],'y')
irf_pl_mark(h,t1,'k')
irf_pl_mark(h,t2,'k')

h(1).Title.String = ['MMS ',num2str(ic)];

%% plot distribution as lines for the two selected lines

% get indices for the times
it1 = interp1(ePDist.time.epochUnix,1:length(ePDist),t1.epochUnix,'nearest');
it2 = interp1(ePDist.time.epochUnix,1:length(ePDist),t2.epochUnix,'nearest');

% matlab colours
col = [0.8500    0.3250    0.0980;...
       0.0000    0.4470    0.7410;...
       0.0000    0.8000    0.0000];

% get f_e and v_par
c_eval('f_e? = f1D(it?).data;',[1,2])
c_eval('v_par? = f1D(it?).depend{1};',[1,2])

% Threshold to ignore spacecraft effects, necessary for next steps
v_par_thresh = 2000; % km/s, a bit arbitrary. Could perhaps automate.
c_eval('v_par?_thr = v_par? .* (v_par? > v_par_thresh);',[1,2])
c_eval('f_e?_thr = f_e? .* (v_par? > v_par_thresh);',[1,2])

% Find out which distribution (1 or 2) is more energised (?)
% Using subscript b for energised beam distribution
% Using subscript c for colder distribution
c_eval('mean_fv? = mean(f_e? .* v_par?_thr);',[1,2])
warm_dist = 1; cold_dist = 2;
if mean_fv2 > mean_fv1
  warm_dist = 2; cold_dist = 1;
end
c_eval(['f_eb = f_e?; f_eb_thr = f_e?_thr;'...
        'v_parb = v_par?; v_parb_thr = v_par?_thr;'],warm_dist)
c_eval(['f_ec = f_e?; f_ec_thr = f_e?_thr;' ...
        'v_parc = v_par?; v_parc_thr = v_par?_thr;'],cold_dist)

% Method 1 of finding beam
% Find gradient sign change to find beam, used to define start of ROI.
% Using warm_dist because beam is generally clearer than cold_dist
grad_fvb = gradient(f_eb,v_parb_thr);
grad_fvb(~isfinite(grad_fvb)) = 0; % remove ±inf/NaN
%%% this line could easily throw an empty array if no positive grad found
%%% could also find a false beam imposter... need to fix this
%%% I think clearly this method is not a great idea
idxb = find(grad_fvb > 0.000001,1,'last'); % find beam index

% Method 2 of finding beam: after removing spacecraft effects, find index
% of maximum f_eb. This removes the same beam imposter error from before,
% and will always find an index (whereas gradient method may not)
[~,idxb] = max(f_eb_thr);

% Definitions: 
% Start of ROI: idx = beam index + 3 (to get away from maximum)
% End of ROI: from roi_a to end / scale_fact
% ROI array: from start to end of ROI, interpolated from warm_dist
scale_fact = 2.5;
% num_pts = 30;
% xbq = linspace(v_parb(idx1),v_parb(round(idx1+(end-idx1)/scale_fact)),num_pts);
yb_pts = f_eb(idxb:round(idxb+(end-idxb)/scale_fact));
xb_pts = v_parb(idxb:round(idxb+(end-idxb)/scale_fact));
% roi_array = interp1(xb_pts,yb_pts,xbq);

% find indices (idx2) of closest points in cold_dist
% idx2 = [];
% for k = 1:length(yb_pts)
%   [~,idxk] = min(abs(f_ec_thr - yb_pts(k)));
%   idx2 = [idx2 idxk]; % maybe find more efficient way to do this
% end

% find first/last indices (idx2a/b) of closest points in cold_dist to beam
[~,idxc_start] = min(abs(f_ec_thr - yb_pts(1)));
[~,idxc_end] = min(abs(f_ec_thr - yb_pts(end)));

% find nearest x,y coords of points in cold_dist
yc_pts_nearest = f_ec(idxc_start:idxc_end);
xc_pts_nearest = v_parc(idxc_start:idxc_end);

% interpolate to find matching y coords (Liouville's theorem)
% need to extrapolate for final xc_pt in array (else it's NaN)
xc_pts = interp1(yc_pts_nearest,xc_pts_nearest,yb_pts,'linear','extrap');
yc_pts = yb_pts;

% calculate array of phi using conservation of energy, assuming parallel
% potential causes electron acceleration, by mapping warm onto cold dist
% for phi_array in volts, need velocities in m/s
phi_array = me*((xb_pts.*1E3).^2 - (xc_pts.*1E3).^2)/(2*qe);
phi_mean = mean(phi_array);
phi_err = std(phi_array); % is this a suitable error in phi?

% map cold onto warm distribution using phi_mean, this will be plotted
% to show the validity of the assumptions made for these two dists
% for xb_liouv in km/s, need to divide phi_mean by 1E6
xb_liouv = sqrt(2*qe*phi_mean/(me*1E6) + xc_pts.^2);

% initiate figure
hca = irf_plot(1,'newfigure');
hold(hca,'on')
% plot warm and cold dists
c_eval(['plot(hca,v_par?,f_e?,''Color'',col(?,:),' ...
        '''linewidth'',2)'],[1,2])
% plot points selected to map from
c_eval('plot(hca,xb_pts,yb_pts,''x'',''Color'',col(?,:))',warm_dist)
c_eval(['plot(hca,xc_pts,yc_pts,''o'',''Color'',col(?,:),' ...
        '''MarkerFaceColor'',col(?,:))'],cold_dist)
% plot Lioville-mapped points
plot(hca,xb_liouv,yb_pts,'o','Color',col(3,:),'MarkerFaceColor', ...
     col(3,:),LineStyle='-')
hold(hca,'off')
hca.YScale = 'log';
hca.YLim = clim;
% legends show time centers
hca.ColorOrder = col; % set colour order for legends
irf_legend(hca,{ePDist(it1).time.toUtc; ...
                ePDist(it2).time.toUtc; ...
                ['Liouville Mapping, \Delta\phi_{||} = ' ...
                num2str(phi_mean,3) '±' num2str(phi_err,3) ...
                'V']},[0.98,0.98])
% hca.XLim = [min(vg),max(vg)];
%labels
xlabel(hca,'V_{||} [km/s]')
ylabel(hca,'f_{e||} [s/m^4]')
hca.Title.String = ['MMS ',num2str(ic)];