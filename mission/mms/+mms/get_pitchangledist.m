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

if(rtrnTS == 0),
    [~,tintpos] = min(abs(pdist.time-tint));
    tint = pdist.time(tintpos);
    Bvec = TSeries(Bvec.time(tintpos),Bvec.data(tintpos,:));
    pdist = TSeries(pdist.time(tintpos),pdist.data(tintpos,:,:,:));
    if(noangles == 0),
        stepTable = TSeries(stepTable.time(tintpos),stepTable.data(tintpos));
    end
end

numechannels = 32;

if noangles,
    % Define angles
    dangle = 180/16;
    phi = dangle*[0:31]+dangle/2;
    theta = dangle*[0:15]+dangle/2;
    x = -cosd(phi')*sind(theta);
    y = -sind(phi')*sind(theta);
    [~,energy] = hist([log10(10),log10(30e3)],32);
    energy = 10.^energy;
    energy = ones(length(pdist.time),1)*energy;
    lengthphi = 32;
else
    lengthphi = length(phi.data(1,:));
    xt = zeros(length(pdist.time),lengthphi,length(theta));
    yt = zeros(length(pdist.time),lengthphi,length(theta));
    for ii = 1:length(pdist.time);
        xt(ii,:,:) = -cosd(phi.data(ii,:)')*sind(theta);
        yt(ii,:,:) = -sind(phi.data(ii,:)')*sind(theta);
    end
    energy = ones(length(pdist.time),1)*energy0;
    for ii = 1:length(pdist.time);
        if stepTable.data(ii),
        	energy(ii,:) = energy1;
        end
    end
end

z = -ones(lengthphi,1)*cosd(theta);

if 0, % Not including solid angles. Leads to spurious results for Bz around zero. 
z2 = ones(lengthphi,1)*sind(theta);
dangle = pi/16;
solida = dangle*dangle*z2;

for jj=1:numechannels;
    allsolid(jj,:,:) = solida;
end
for ii = 1:length(pdist.time);
        allsolid2(ii,:,:,:) = allsolid;
end
end

%pdist.data = pdist.data.*allsolid2;

anglevec = [15 30 45 60 75 90 105 120 135 150 165 180];
pitcha = anglevec-7.5;

paddistarr = zeros(length(pdist.time),32,12);

for ii = 1:length(pdist.time);
    dist1 = squeeze(pdist.data(ii,:,:,:));
    if(noangles == 0),
        x = squeeze(xt(ii,:,:));
        y = squeeze(yt(ii,:,:));
    end
    thetab = acosd(x*Bvec.data(ii,1)+y*Bvec.data(ii,2)+z*Bvec.data(ii,3));
    c_eval('pos? = ones(lengthphi,length(theta));',anglevec); 
    %c_eval('sol? = solida;',anglevec);    
    pos15(thetab > 15) = NaN;
    %sol15(thetab > 15) = NaN;
    for jj = 2:(length(anglevec)-1)
        c_eval('pos?(thetab < (?-15)) = NaN;',anglevec(jj));
        c_eval('pos?(thetab > ?) = NaN;',anglevec(jj));
        %c_eval('sol?(thetab < (?-15)) = NaN;',anglevec(jj));
        %c_eval('sol?(thetab > ?) = NaN;',anglevec(jj));
    end
    pos180(thetab < 165) = NaN;
    %sol180(thetab < 165) = NaN;
    
    c_eval('dist? = dist1;',anglevec);
    
    for kk = 1:numechannels;
        c_eval('dist?(kk,:,:)  = squeeze(dist1(kk,:,:)).*pos?;',anglevec);
    end 
    %c_eval('dist? =  squeeze(irf.nanmean(irf.nanmean(dist?,3),2))/irf.nanmean(irf.nanmean(sol?));',anglevec);
    c_eval('dist? =  squeeze(irf.nanmean(irf.nanmean(dist?,3),2));',anglevec);
    paddistarr(ii,:,:) = [dist15 dist30 dist45 dist60 dist75 dist90 dist105 dist120 dist135 dist150 dist165 dist180];
end

if rtrnTS,
    paddist = TSeries(pdist.time,paddistarr);
    tint = irf.tint(paddist.time.start.utc,paddist.time.stop.utc);
else
    paddist = squeeze(paddistarr);
end

theta = pitcha;

end