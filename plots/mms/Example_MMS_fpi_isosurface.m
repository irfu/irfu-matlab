% Script to plot a 3D isosurface of a FPI particle distribution. The
% distribution is converted from a spherical grid to a cartesian one to
% facilitate plotting. So far, the example only includes an electron
% distribution. Ions should be added.
%

%% set input

% time interval (short)
% butterfly distribution from (Zhang et al., 2021, 10.1029/2021GL096056)
tint = irf.tint('2015-10-14T06:15:25.75/2015-10-14T06:15:26.25');
% time of distribution (nearest)
t = irf.time_array('2015-10-14T06:15:26.0',0);

% spacecraft number
ic = 1;

%% get data
ePDist = mms.get_data('PDe_fpi_brst_l2',tint,ic);
Bdmpa = mms.get_data('B_dmpa_fgm_brst_l2',tint,ic);


%% call function to convert distribution to cartesian grid
%
u = irf_units;

% get index of measurement nearest t
it = interp1(ePDist.time.epochUnix,1:ePDist.length,t.epochUnix,'nearest');

M = u.me;
%M = u.mp;

% set velocity grid (same in all directions)
vmax = 3e7; % m/s
nvg = 100;
vg = linspace(-vmax,vmax,nvg);

% call function
tic
[VXG,VYG,VZG,Fg] = get3Ddist(ePDist.convertto('s^3/m^6'),it,vg,M);
toc


%% plot isosurface

B0 = Bdmpa.resample(ePDist.time(it)).data;
b0 = B0/norm(B0);

figure;

fsurf = 2e-16;

% isosurface does not allow for axis input, this might be an issue if you
% have many panels. It seems better to let isosurface initiate the axis
isosurface(VXG*1e-4,VYG*1e-4,VZG*1e-4,Fg,fsurf); axis equal

% then set axis handle hca
hca = gca;

hold(hca,'on')
Blen = 1500; % "length" of B arrow

xlabel(hca,'$v_x$ ','interpreter','latex')
ylabel(hca,'$v_y$ ','interpreter','latex')
zlabel(hca,'$v_z$ [10$^4$~m\,s$^{-1}$]','interpreter','latex')

plot3(hca,[-1,1]*b0(1)*Blen,[-1,1]*b0(2)*Blen,[-1,1]*b0(3)*Blen,'k','linewidth',4)

hca.Box = 'off';
grid on
hca.LineWidth = 2;
hca.FontSize = 18;

% ugly title
title(hca,['f_{iso} = ',num2str(fsurf),' s^3 m^{-6}'],'interpreter','tex')


%% function
function [VXG,VYG,VZG,Fg] = get3Ddist(dist,it,vg,M)
% GET3DDIST get 3D distribution in cartesian grid (DMPA coordinates)
%  [VXG,VYG,VZG,F] = get3Ddist(dist,vg,M) get velocity mesh grids VXG, VYG,
%  VZG, given PDist dist, velocity array vg, and mass M.

% elementary charge
qe = 1.6022e-19;

emat = double(dist.depend{1}); % in eV
energy = emat(it,:);
v = sqrt(2*energy*qe/M); % m/s

% azimuthal angle
phi = double(dist.depend{2}(it,:)); % in degrees
%phi = phi+180;
%phi(phi>360) = phi(phi>360)-360;
phi = phi-180;
phi = phi*pi/180; % in radians

% elevation angle
th = double(dist.depend{3}); % polar angle in degrees
th = th-90; % elevation angle in degrees
th = th*pi/180; % in radians


nAz = length(phi);
nEle = length(th);
nV = length(v);

%
% 3D matrices for instrumental bin centers
TH = repmat(th,nV,1,nAz);       % [phi,th,v]
TH = permute(TH,[1,3,2]);       % [v,phi,th]
PHI = repmat(phi,nV,1,nEle);    % [v,phi,th]
VEL = repmat(v,nAz,1,nEle);     % [phi,v,th]
VEL = permute(VEL,[2,1,3]);     % [v,phi,th]

F = double(squeeze(dist.data(it,:,:,:)));

% instrument grid in DMPA
[VX,VY,VZ] = sph2cart(PHI,TH,VEL);

% make grid

% allow for different grids in future ?
vxg = vg; vyg = vg; vzg = vg;

[VXG,VYG,VZG] = meshgrid(vxg,vyg,vzg);

% get the cartesian distribution, this is the slow part
Fg = griddata(VX,VY,VZ,F,VXG,VYG,VZG,'linear');

end
