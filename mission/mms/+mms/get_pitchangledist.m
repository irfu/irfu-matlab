function [paddist,theta,energy,tint] = get_pitchangledist(varargin) 
% GET_PITCHANGLEDIST compute pitch angle distributions from burst mode data
%
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,B,[tint]) - For v0.2.0 data
% [paddist,theta,energy,tint] = mms.get_pitchangledist(pdist,phi,theta,stepTable,energy0,energy1,B,[tint]) - For v1.0.0 or higher data
%
%
% Computes the pitch angle distributions from l1b brst particle data. 
%
% Inputs:
%       pdist - electron or ion particle distribution in TSeries
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
% Support for v0.2.0 should be removed when it is not needed
% 
% Written by D. B. Graham
%

tic;

anglevec = [15 30 45 60 75 90 105 120 135 150 165 180];
pitcha = anglevec-7.5;
numechannels = 32;
lengthphi = 32;
lengththeta = 16;

% Input check
rtrnTS = 1;
if (nargin == 2 || nargin==3),
    pdist = varargin{1};
    B = varargin{2};
    noangles = 1;
    if (nargin == 3),
        tint = varargin{3};
        if(length(tint) > 2),
            irf_log('proc','Format of tint is wrong.')
            return; 
        end
        rtrnTS = 0;
        if(length(tint) == 2),
            rtrnTS = 1;
            pdist = pdist.tlim(tint);
            B = B.tlim(tint);
        end
    end
    irf_log('proc','No angles passed. Default values used.')
elseif (nargin==7 || nargin==8)
    pdist = varargin{1};
    phi = varargin{2};
    theta = varargin{3};
    stepTable = varargin{4};
    energy0 = varargin{5};
    energy1 = varargin{6};
    if isstruct(theta),
        theta = theta.data;
    end
    if isstruct(energy0),
        energy0 = energy0.data;
    end
    if isstruct(energy1),
        energy1 = energy1.data;
    end
    B = varargin{7};
    noangles = 0;
    if (nargin == 8),
        tint = varargin{8};
        if(length(tint) > 2),
            irf_log('proc','Format of tint is wrong.')
            return; 
        end
        rtrnTS = 0;
        if(length(tint) == 2),
            rtrnTS = 1;
            pdist = pdist.tlim(tint);
            stepTable = stepTable.tlim(tint);
            B = B.tlim(tint);
        end
    end
    irf_log('proc','Angles passed; using these to calculate PAD.')
else
    irf_log('proc','Input not recognized.')
    return;
end

B = B.resample(pdist);
Bvec = B/B.abs;
Bvecx = repmat(Bvec.data(:,1),1,numechannels,lengthphi,lengththeta);
Bvecy = repmat(Bvec.data(:,2),1,numechannels,lengthphi,lengththeta);
Bvecz = repmat(Bvec.data(:,3),1,numechannels,lengthphi,lengththeta);


if(rtrnTS == 0),
    [~,tintpos] = min(abs(pdist.time-tint));
    tint = pdist.time(tintpos);
    Bvecx = ones(32,32,16)*Bvec.data(tintpos,1);
    Bvecy = ones(32,32,16)*Bvec.data(tintpos,2);
    Bvecz = ones(32,32,16)*Bvec.data(tintpos,3);
    pdist = TSeries(pdist.time(tintpos),pdist.data(tintpos,:,:,:));
    if(noangles == 0),
        stepTable = TSeries(stepTable.time(tintpos),stepTable.data(tintpos));
    end
end

if noangles,
    % Define angles
    dangle = 180/16;
    phi = dangle*[0:31]+dangle/2;
    theta = dangle*[0:15]+dangle/2;
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

    for ii = 1:length(pdist.time);
        x(ii,:,:) = -cosd(phi.data(ii,:)')*sind(theta);
        y(ii,:,:) = -sind(phi.data(ii,:)')*sind(theta);
        z(ii,:,:) = -ones(lengthphi,1)*cosd(theta);
    end
    energy = ones(length(pdist.time),1)*energy0;
    for ii = 1:length(pdist.time);
        if stepTable.data(ii),
        	energy(ii,:) = energy1;
        end
    end
end

xt = repmat(x,1,1,1,32);
xt = squeeze(permute(xt,[1 4 2 3]));
yt = repmat(y,1,1,1,32);
yt = squeeze(permute(yt,[1 4 2 3]));
zt = repmat(z,1,1,1,32);
zt = squeeze(permute(zt,[1 4 2 3]));

thetab = acosd(xt.*Bvecx+yt.*Bvecy+zt.*Bvecz);

c_eval('dist? = pdist.data;',anglevec);
dist15(thetab > 15) = NaN;
for jj = 2:(length(anglevec)-1)
	c_eval('dist?(thetab < (?-15)) = NaN;',anglevec(jj));
	c_eval('dist?(thetab > ?) = NaN;',anglevec(jj));
end 
dist180(thetab < 165) = NaN;
c_eval('dist? =  squeeze(irf.nanmean(irf.nanmean(dist?,4),3));',anglevec);

paddistarr = cat(3,dist15,dist30,dist45,dist60,dist75,dist90,dist105,dist120,dist135,dist150,dist165,dist180);

if rtrnTS,
    paddist = TSeries(pdist.time,paddistarr);
    tint = irf.tint(paddist.time.start.utc,paddist.time.stop.utc);
else
    paddist = squeeze(paddistarr);
end

theta = pitcha;
toc;

end