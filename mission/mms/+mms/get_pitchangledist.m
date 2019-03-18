function [paddist,theta,energy,tint] = get_pitchangledist(varargin) 
% GET_PITCHANGLEDIST compute pitch angle distributions from burst mode data
%
% Examples:
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,B,[tint]) - For v0.2.0 data
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,phi,theta,stepTable,energy0,energy1,B,[tint]) - For v1.0.0 or higher data
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,B,[tint]) - For PDist structure
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,B,[tint],'angles',24) - For PDist structure
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,B,[tint],'angles',[0 45 90 135 180]) - For PDist structure
%
% Computes the pitch angle distributions from l1b brst particle data. 
%
% Inputs:
%       pdist - electron or ion particle distribution in TSeries or PDist
%       format
%       B     - magnetic field in TSeries
%       phi   - TSeries of all phi angles of distribution.
%       theta - 1D array or structure of theta angles. 
%       stepTable - TSeries of stepping table between energies.
%       energy0 - 1D array or structure of energy table 0 (burst)
%       energy1 - 1D array or structure of energy table 1 (burst)
%       tint - optional time or time interval (TSeries format). If a single element paddist is an
%       array of the pitch-angle distribution at the time closest to tint.
%       If tint has two elements paddist is a Tseries over tint. 
% Output: 
%       paddist - Particle pitch angle distribution
%       theta - Pitch angles
%       energy - energy as a vector or array corresponding to paddist
%       tint - interval or time of the returned distibution(s)
%
% Options: 
%       'angles' - flag to set number of pitch angles used (single value Default is
%       12) or pitch angle edges (array). When counts are high (plasma is dense) 24 pitch angles for
%       electrons may be used. 
%       'meanorsum' - 'mean' (Default) or 'sum'; + @ 2018-03-20;
%
% Support for v0.2.0 should be removed when it is not needed
% 
% Written by D. B. Graham
%

tic;

numechannels = 32;
lengthphi = 32;
lengththeta = 16;
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

% Input check
rtrnTS = 1;
if isa(varargin{1},'PDist')
    if strcmp('skymap',varargin{1}.type)
        irf.log('warning','PDist is skymap format; computing pitch angle distribution.');
        pdist = varargin{1};
        B = varargin{2};
        tmpPhi = pdist.depend{1,2};
        if size(tmpPhi,1) == 1 % only one value for one time or for all times
          phi = TSeries(pdist.time,repmat(pdist.depend{1,2},pdist.length,1));
        else
          phi = TSeries(pdist.time,pdist.depend{1,2});
        end
          
        theta = pdist.depend{1,3};
        if isfield(pdist.ancillary, 'esteptable')
            stepTable = TSeries(pdist.time,pdist.ancillary.esteptable);
        else
            stepTable = TSeries(pdist.time, zeros(length(pdist.time), 1));
        end
        if and(isfield(pdist.ancillary, 'energy0'), isfield(pdist.ancillary, 'energy1'))
            energy0 = pdist.ancillary.energy0;
            energy1 = pdist.ancillary.energy1;            
        else    
            if isfield(pdist.ancillary, 'energy')
                energy0 = pdist.ancillary.energy(1, :);
                energy1 = pdist.ancillary.energy(2, :);
            else
                irf.log('warning', 'no data for energy0 & energy1.');
            end
        end
        noangles = 0;        
        numechannels = size(pdist.depend{1},2);
        if (tmpnargin == 3)
            tint = varargin{3};
            if(length(tint) > 2)
                irf.log('critical','Format of tint is wrong.');
                return; 
            end
            rtrnTS = 0;
            if (length(tint) == 2)
                rtrnTS = 1;
                pdist = pdist.tlim(tint);
                B = B.tlim(tint);
                phi = phi.tlim(tint);
                stepTable = stepTable.tlim(tint);
            end
        end
    else
        irf.log('critical','PDist must be skymap.');
        return;        
    end
elseif (tmpnargin == 2 || tmpnargin==3)
    pdist = varargin{1};
    B = varargin{2};
    noangles = 1;
    if (tmpnargin == 3)
        tint = varargin{3};
        if(length(tint) > 2)
            irf.log('critical','Format of tint is wrong.');
            return; 
        end
        rtrnTS = 0;
        if(length(tint) == 2)
            rtrnTS = 1;
            pdist = pdist.tlim(tint);
            B = B.tlim(tint);
        end
    end
    irf.log('warning','No angles passed. Default values used.');
elseif (tmpnargin==7 || tmpnargin==8)
    pdist = varargin{1};
    phi = varargin{2};
    theta = varargin{3};
    stepTable = varargin{4};
    energy0 = varargin{5};
    energy1 = varargin{6};
    if isstruct(theta)
        theta = theta.data;
    end
    if isstruct(energy0)
        energy0 = energy0.data;
    end
    if isstruct(energy1)
        energy1 = energy1.data;
    end
    B = varargin{7};
    noangles = 0;
    if (tmpnargin == 8)
        tint = varargin{8};
        if(length(tint) > 2)
            irf.log('critical','Format of tint is wrong.');
            return; 
        end
        rtrnTS = 0;
        if(length(tint) == 2)
            rtrnTS = 1;
            pdist = pdist.tlim(tint);
            phi = phi.tlim(tint);
            stepTable = stepTable.tlim(tint);
            B = B.tlim(tint);
        end
    end
    irf.log('notice','Angles passed; using these to calculate PAD.');
else
    irf.log('critical','Input not recognized.');
    return;
end

% Check size of energy
if noangles == 0
    numechannels = length(energy0);
    lengthphi = length(phi.data(1,:));
    lengththeta = length(theta);
end

B = B.resample(pdist);
Bvec = B/B.abs;
Bvecx = repmat(Bvec.data(:,1),1,numechannels,lengthphi,lengththeta);
Bvecy = repmat(Bvec.data(:,2),1,numechannels,lengthphi,lengththeta);
Bvecz = repmat(Bvec.data(:,3),1,numechannels,lengthphi,lengththeta);


if(rtrnTS == 0)
    [~,tintpos] = min(abs(pdist.time-tint));
    tint = pdist.time(tintpos);
    Bvecx = ones(numechannels,lengthphi,lengththeta)*Bvec.data(tintpos,1);
    Bvecy = ones(numechannels,lengthphi,lengththeta)*Bvec.data(tintpos,2);
    Bvecz = ones(numechannels,lengthphi,lengththeta)*Bvec.data(tintpos,3);
    pdist = TSeries(pdist.time(tintpos),pdist.data(tintpos,:,:,:));
    if(noangles == 0)
        stepTable = TSeries(stepTable.time(tintpos),stepTable.data(tintpos));
    end
end

if noangles
    % Define angles
    dangle = 180/16;
    phi = dangle*(0:31)+dangle/2;
    theta = dangle*(0:15)+dangle/2;
    x = -cosd(phi')*sind(theta);
    y = -sind(phi')*sind(theta);
    z = -ones(lengthphi,1)*cosd(theta);
    x = repmat(x,1,1,length(pdist.time));
    y = repmat(y,1,1,length(pdist.time));
    z = repmat(z,1,1,length(pdist.time));
    x = permute(x,[3 1 2]);
    y = permute(y,[3 1 2]);
    z = permute(z,[3 1 2]);
    [~,energy] = hist([log10(10),log10(30e3)],32);
    energy = 10.^energy;
    energy = ones(length(pdist.time),1)*energy;
else
    x = zeros(length(pdist.time),lengthphi,lengththeta);
    y = zeros(length(pdist.time),lengthphi,lengththeta);
    z = zeros(length(pdist.time),lengthphi,lengththeta);

    for ii = 1:length(pdist.time)
        x(ii,:,:) = -cosd(phi.data(ii,:)')*sind(theta);
        y(ii,:,:) = -sind(phi.data(ii,:)')*sind(theta);
        z(ii,:,:) = -ones(lengthphi,1)*cosd(theta);
    end

    if isa(pdist,'PDist') && size(pdist.depend{1},1) == pdist.length
      energy = pdist.depend{1};
    else % make energy table from ancillary data
      energy = ones(length(pdist.time),1)*energy0;
      for ii = 1:length(pdist.time)
          if stepTable.data(ii)
           energy(ii,:) = energy1;
          end
      end
    end
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
    c_eval('sumsr? =  squeeze(irf.nansum(irf.nansum(sr?,4),3));',1:length(anglevec));
    c_eval('dist? = dist?./sumsr?;',1:length(anglevec));
end    
%paddistarr = cat(3,dist15,dist30,dist45,dist60,dist75,dist90,dist105,dist120,dist135,dist150,dist165,dist180);
paddistarr = dist1;
for ii = 2:length(anglevec)
    c_eval('paddistarr = cat(3,paddistarr,dist?);',ii);
end

if rtrnTS
    paddist = TSeries(pdist.time,paddistarr);
    tint = irf.tint(paddist.time.start.utc,paddist.time.stop.utc);
else
    paddist = squeeze(paddistarr);
end

theta = pitcha;
toc;
if (isa(varargin{1},'PDist') && length(tint)==2)
    paddist = PDist(pdist.time,paddistarr,'pitchangle',energy,theta);
    paddist.units = pdist.units;
    paddist.species = pdist.species;
    paddist.ancillary = varargin{1}.ancillary;
    paddist.ancillary.meanorsum = meanorsum;
    paddist.ancillary.delta_pitchangle_minus = dangle(1:end)*0.5;
    paddist.ancillary.delta_pitchangle_plus = dangle(1:end)*0.5;
end

end
