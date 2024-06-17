function [paddist,theta,energy,tint] = get_pitchangledist(varargin)
% GET_PITCHANGLEDIST compute pitch angle distributions from burst mode data
%
% Examples:
% paddist = mms.get_pitchangledist(pdist,B,'angles',24)
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,B,tint,'angles',[0 45 90 135 180])
%
% Inputs:
%       pdist - electron or ion particle distribution in PDist format (must be skymap)
%       B     - magnetic field in TSeries (must be in the same coordinate system as pdist)
%       tint  - optional time or time interval (TSeries format). If a single element paddist is an
%       array of the pitch-angle distribution at the time closest to tint.
%       If tint has two elements paddist is a Tseries/PDist over tint.
% Output:
%       paddist - Particle pitch angle distribution (PDist format)
%
% Options:
%       'angles' - flag to set number of pitch angles used (single value Default is
%       12) or pitch angle edges (array). When counts are high (plasma is dense) 24 pitch angles for
%       electrons may be used.
%       'meanorsum' - 'mean' (Default) or 'sum'; + @ 2018-03-20;
%
% Written by D. B. Graham
%

anglevec = 15:15:180; % Default pitch angles. 15 degree angle widths
dangle = median(diff(anglevec))*ones(1,length(anglevec));
tmpnargin = nargin;
meanorsum = 'mean';
have_char = 0;

if isa(varargin{tmpnargin-1},'char'), have_char = 1; end

while have_char
  if isa(varargin{tmpnargin-1},'char')
    if strcmp(varargin{tmpnargin-1},'angles')
      if (length(varargin{tmpnargin}) == 1)
        nangles = varargin{tmpnargin};
        nangles = floor(nangles); % Make sure input is integer
        dangle = 180/nangles;
        anglevec = dangle:dangle:180;
        dangle = dangle*ones(1,length(anglevec));
        irf.log('notice','User defined number of pitch angles.')
      elseif (length(varargin{tmpnargin}) > 1)
        anglevec = varargin{tmpnargin};
        dangle = diff(anglevec);
        anglevec = anglevec(2:end);
        irf.log('notice','User defined pitch angle limits.')
      else
        irf.log('warning','angles parameter not understood. ')
      end
    elseif strcmp(varargin{tmpnargin-1},'meanorsum')
      meanorsum = varargin{tmpnargin};
    else
      irf.log('warning', 'parameters not understood. ')
    end
  end
  tmpnargin = tmpnargin-2;
  have_char = 0;
  if isa(varargin{tmpnargin-1},'char'), have_char = 1; end
end

pitcha = anglevec-dangle/2;

%anglevec
% Input check
findnearest = 0;
if isa(varargin{1},'PDist') && strcmp('skymap',varargin{1}.type)
  pdist = varargin{1};
  B = varargin{2};
  if (tmpnargin == 3)
    tint = varargin{3};
    if(length(tint) > 2)
      irf.log('critical','Format of tint is wrong.');
      return;
    end
    findnearest = 1;
    if (length(tint) == 2)
      findnearest = 0;
      pdist = pdist.tlim(tint);
      B = B.tlim(tint);
    end
  end
  tmpPhi = pdist.depend{1,2};
  if size(tmpPhi,1) == 1 % only one value for one time or for all times
    phi = TSeries(pdist.time,repmat(pdist.depend{1,2},pdist.length,1));
  else
    phi = TSeries(pdist.time,pdist.depend{1,2});
  end
  theta = pdist.depend{1,3};
else
  irf.log('critical','PDist must be skymap.');
  return;
end

energy = pdist.depend{1,1};

% Check size of energy
numechannels = length(energy(1,:));
lengthphi = length(phi.data(1,:));
lengththeta = length(theta);


B = B.resample(pdist);
Bvec = B/B.abs;
Bvecx = repmat(Bvec.data(:,1),1,numechannels,lengthphi,lengththeta);
Bvecy = repmat(Bvec.data(:,2),1,numechannels,lengthphi,lengththeta);
Bvecz = repmat(Bvec.data(:,3),1,numechannels,lengthphi,lengththeta);

if 0
if findnearest == 1
  [~,tintpos] = min(abs(pdist.time-tint));
  tint = pdist.time(tintpos);
  Bvecx = ones(numechannels,lengthphi,lengththeta)*Bvec.data(tintpos,1);
  Bvecy = ones(numechannels,lengthphi,lengththeta)*Bvec.data(tintpos,2);
  Bvecz = ones(numechannels,lengthphi,lengththeta)*Bvec.data(tintpos,3);
  pdist = TSeries(pdist.time(tintpos),pdist.data(tintpos,:,:,:));
  energy = energy(tintpos,:);
end

x = zeros(length(pdist.time),lengthphi,lengththeta);
y = zeros(length(pdist.time),lengthphi,lengththeta);
z = zeros(length(pdist.time),lengthphi,lengththeta);

for ii = 1:length(pdist.time)
  x(ii,:,:) = -cosd(phi.data(ii,:)')*sind(theta);
  y(ii,:,:) = -sind(phi.data(ii,:)')*sind(theta);
  z(ii,:,:) = -ones(lengthphi,1)*cosd(theta);
end

xt = repmat(x,1,1,1,numechannels);
xt = squeeze(permute(xt,[1 4 2 3]));
yt = repmat(y,1,1,1,numechannels);
yt = squeeze(permute(yt,[1 4 2 3]));
zt = repmat(z,1,1,1,numechannels);
zt = squeeze(permute(zt,[1 4 2 3]));

%thetab = acosd(xt.*Bvecx+yt.*Bvecy+zt.*Bvecz);
thetab = acosd(xt.*squeeze(Bvecx)+yt.*squeeze(Bvecy)+zt.*squeeze(Bvecz));


c_eval('dist? = pdist.data;',1:length(anglevec));
%dist1(thetab > anglevec(1)) = NaN;
for jj = 1:(length(anglevec))
  c_eval('dist?(thetab < (anglevec(?)-dangle(?))) = NaN;',jj);
  c_eval('dist?(thetab > anglevec(?)) = NaN;',jj);
end
%c_eval('dist?(thetab < (anglevec(end)-dangle(?))) = NaN;',length(anglevec));
if strcmp(meanorsum, 'mean')
  c_eval('dist? =  squeeze(irf.nanmean(irf.nanmean(dist?,4),3));',1:length(anglevec));
elseif strcmp(meanorsum, 'sum')
  c_eval('dist? =  squeeze(irf.nansum(irf.nansum(dist?,4),3));',1:length(anglevec));
elseif strcmp(meanorsum, 'sum_weighted')
  c_eval('sr? = pdist.solidangle.data;',1:length(anglevec));
  c_eval('sr?(isnan(dist?)) = NaN;',1:length(anglevec));
  c_eval('dist? =  squeeze(irf.nansum(irf.nansum(dist?.*sr?,4),3));',1:length(anglevec));
  c_eval('sumsr? =  squeeze(  irf.nansum(irf.nansum(sr?,4),3));',1:length(anglevec));
  c_eval('dist? = dist?./sumsr?;',1:length(anglevec));
end
%paddistarr = cat(3,dist15,dist30,dist45,dist60,dist75,dist90,dist105,dist120,dist135,dist150,dist165,dist180);
paddistarr = dist1;
for ii = 2:length(anglevec)
  c_eval('paddistarr = cat(3,paddistarr,dist?);',ii);
end

theta = pitcha;
else
  
  [VX,VY,VZ] = pdist.v;
  VABS = sqrt(VX.^2 + VY.^2 + VZ.^2);
  VXnorm = VX./VABS;
  VYnorm = VY./VABS;
  VZnorm = VZ./VABS;
  
  angle_b = acosd(VXnorm.*Bvecx + VYnorm.*Bvecy + VZnorm.*Bvecz);
  
  % Bin pitchangles
  pitch_angle_bin_edges = 0:15:180;
  pitch_angle_bin_edges = [0 anglevec];
  dangle = diff(pitch_angle_bin_edges);
  theta = pitch_angle_bin_edges(1:end-1) + 0.5*dangle;
  nPitchangles = numel(pitch_angle_bin_edges)-1;
  bins = discretize(angle_b,pitch_angle_bin_edges);

  % Accumulate all the data into proper grid 
  energy_bin_edges = [pdist.ancillary.energy(1,:) - pdist.ancillary.delta_energy_minus(1,:) pdist.ancillary.energy(1,end) + pdist.ancillary.delta_energy_plus(1,end)];
  nEnergies = numel(pdist.depend{1}(1,:));  
  nTimes = pdist.length;
  
  datasize = size(pdist.data);
  
  bins_pa = discretize(angle_b(:),pitch_angle_bin_edges);
  bins_en = repmat(discretize(pdist.depend{1}(:),energy_bin_edges),[prod(datasize(3:4)),1]); %repmat((1:nEnergies)',[nTimes 1]);  
  bins_time = repmat((1:nTimes)',[prod(datasize(2:end)),1]);
  %bins_en = discretize(angle_b,pitch_angle_bin_edges);
  %for it = 1:nTimes
%     for iE = 1:nEnergies
      %bins = discretize(squeeze(angle_b(:,iE,:)),pitch_angle_bin_edges);
      %data = squeeze(pdist.data(:,iE,:));
      subs = [bins_time,bins_en,bins_pa];
      vals = pdist.data(:);
      paddistarr = accumarray(subs,vals,[nTimes,nEnergies,nPitchangles],@nanmean);
      
%       A(:,iE,:) = accumarray(subs,vals,[nTimes,nPitchangles],@nanmean);
%       N(:,iE,:) = accumarray([(1:nTimes)' bins],(pdist.data(:,iE,:,:)>0),[nTimes,nPitchangles],@sum);
%     end
  
  1;
  %end
end
if findnearest == 0
  paddist = PDist(pdist.time,paddistarr,'pitchangle',energy,theta);
  paddist.units = varargin{1}.units;
  paddist.species = varargin{1}.species;
  paddist.ancillary = varargin{1}.ancillary;
  paddist.ancillary.meanorsum = meanorsum;
  paddist.ancillary.delta_pitchangle_minus = dangle(1:end)*0.5;
  paddist.ancillary.delta_pitchangle_plus = dangle(1:end)*0.5;
else
  paddist = squeeze(paddistarr);
end
end
