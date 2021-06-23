% Examples of Partial moments calculations using mms.psd_moments.

%% Example 1. Compute moments of ions with positive and negative V_x individually
% Load ion data
ic = 1; % Spacecraft number
Tint = irf.tint('2015-11-04T04:57:50.00Z/2015-11-04T04:58:30.00Z');
c_eval('iPDist = mms.get_data(''PDi_fpi_brst_l2'',Tint,?);',ic);
c_eval('SCpot = mms.get_data(''V_edp_brst_l2'',Tint,?);',ic);

x = zeros(length(iPDist.time),32,16);
for ii = 1:length(iPDist.time)
  x(ii,:,:) = -cosd(iPDist.depend{1,2}(ii,:)')*sind(iPDist.depend{1,3});
end
x = repmat(x,1,1,1,32);
x = squeeze(permute(x,[1 4 2 3]));
xpositive = zeros(size(x)); xnegative = zeros(size(x));
xpositive(x > 0) = 1; xnegative(x < 0) = 1;

imoms = mms.psd_moments(iPDist,SCpot);
imomspos = mms.psd_moments(iPDist,SCpot,'partialmoms',xpositive);
imomsneg = mms.psd_moments(iPDist,SCpot,'partialmoms',xnegative);

%% Example 2. Compute electron moments for pitch-angles < 90 deg and > 90 deg.

ic = 1; % Spacecraft number
Tint = irf.tint('2015-12-06T23:38:27.00Z/2015-12-06T23:38:37.00Z');
c_eval('ePDist = mms.get_data(''PDe_fpi_brst_l2'',Tint,?);',ic);
c_eval('SCpot = mms.get_data(''V_edp_brst_l2'',Tint,?);',ic);
c_eval('Bxyz = mms.get_data(''B_dmpa_brst_l2'',Tint,?);',ic);

Bxyz = Bxyz.resample(ePDist);
Bvec = Bxyz/Bxyz.abs;
Bvecx = repmat(Bvec.data(:,1),1,32,32,16);
Bvecy = repmat(Bvec.data(:,2),1,32,32,16);
Bvecz = repmat(Bvec.data(:,3),1,32,32,16);
x = zeros(length(ePDist.time),32,16);
y = zeros(length(ePDist.time),32,16);
z = zeros(length(ePDist.time),32,16);

for ii = 1:length(ePDist.time)
  x(ii,:,:) = -cosd(ePDist.depend{1,2}(ii,:)')*sind(ePDist.depend{1,3});
  y(ii,:,:) = -sind(ePDist.depend{1,2}(ii,:)')*sind(ePDist.depend{1,3});
  z(ii,:,:) = -ones(32,1)*cosd(ePDist.depend{1,3});
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

emoms = mms.psd_moments(ePDist,SCpot);
emomsl90 = mms.psd_moments(ePDist,SCpot,'partialmoms',pitchl90);
emomsg90 = mms.psd_moments(ePDist,SCpot,'partialmoms',pitchg90);

neall = TSeries(emoms.n_psd.time,[emoms.n_psd.data emomsl90.n_psd.data emomsg90.n_psd.data]);
Te = emoms.T_psd.trace/3;
Tel90 = emomsl90.T_psd.trace/3;
Teg90 = emomsg90.T_psd.trace/3;
Te = TSeries(Te.time,[Te.data Tel90.data Teg90.data]);