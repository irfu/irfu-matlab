% Examples of Partial moments calculations using mms.psd_moments.

%% Example 1. Compute moments of ions with positive and negative V_x individually 

% Load ion skymap data...

x = zeros(length(dist.time),32,16);
for ii = 1:length(dist.time)
	x(ii,:,:) = -cosd(phi.data(ii,:)')*sind(theta.data);
end
x = repmat(x,1,1,1,32);
x = squeeze(permute(x,[1 4 2 3]));
xpositive = zeros(size(x)); xnegative = zeros(size(x));
xpositive(x > 0) = 1; xnegative(x < 0) = 1;

imoms = mms.psd_moments(dist,phi,theta,stepTable,energy0,energy1,SCpot,'ion');
imomspos = mms.psd_moments(dist,phi,theta,stepTable,energy0,energy1,SCpot,'ion','partialmoms',xpositive);
imomsneg = mms.psd_moments(dist,phi,theta,stepTable,energy0,energy1,SCpot,'ion','partialmoms',xnegative);

%% Example 2. Compute electron moments for pitch-angles < 90 deg and > 90 deg.

% Load skymap data and B field

Bxyz = Bxyz.resample(dist);
Bvec = Bxyz/Bxyz.abs;
Bvecx = repmat(Bvec.data(:,1),1,32,32,16);
Bvecy = repmat(Bvec.data(:,2),1,32,32,16);
Bvecz = repmat(Bvec.data(:,3),1,32,32,16);
x = zeros(length(dist.time),32,16);
y = zeros(length(dist.time),32,16);
z = zeros(length(dist.time),32,16);

for ii = 1:length(dist.time)
	x(ii,:,:) = -cosd(phi.data(ii,:)')*sind(theta.data);
	y(ii,:,:) = -sind(phi.data(ii,:)')*sind(theta.data);
	z(ii,:,:) = -ones(32,1)*cosd(theta.data);
end
x = repmat(x,1,1,1,32);
x = squeeze(permute(x,[1 4 2 3]));
y = repmat(y,1,1,1,32);
y = squeeze(permute(y,[1 4 2 3]));
z = repmat(z,1,1,1,32);
z = squeeze(permute(z,[1 4 2 3]));

thetab = acosd(x.*Bvecx+y.*Bvecy+z.*Bvecz);

pitchl90 = zeros(size(x)); pitchg90 = zeros(size(x));
pitchl90(thetab < 90) = 1; pitchg90(thetab > 90) = 1;

emoms = mms.psd_moments(dist,phi,theta,stepTable,energy0,energy1,SCpot,'electron');
emomsl90 = mms.psd_moments(dist,phi,theta,stepTable,energy0,energy1,SCpot,'electron','partialmoms',pitchl90);
emomsg90 = mms.psd_moments(dist,phi,theta,stepTable,energy0,energy1,SCpot,'electron','partialmoms',pitchg90);

ne = TSeries(emoms.n_psd.time,[emoms.n_psd.data emomsl90.n_psd.data emomsg90.n_psd.data]);
Te = irf.ts_scalar(emoms.T_psd.time,(emoms.T_psd.data(:,1)+emoms.T_psd.data(:,4)+emoms.T_psd.data(:,3))/3);
Tel90 = irf.ts_scalar(emomsl90.T_psd.time,(emomsl90.T_psd.data(:,1)+emomsl90.T_psd.data(:,4)+emomsl90.T_psd.data(:,3))/3);
Teg90 = irf.ts_scalar(emomsg90.T_psd.time,(emomsg90.T_psd.data(:,1)+emomsg90.T_psd.data(:,4)+emomsg90.T_psd.data(:,3))/3);
Te = TSeries(Te.time,[Te.data Tel90.data Teg90.data]);