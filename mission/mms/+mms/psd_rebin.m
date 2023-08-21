function [pdistr,phir,energyr] = psd_rebin(pdist,phi,energy0,energy1,stepTable)
%PSD_REBIN Convert burst mode distribution into 64 energy channel distribution
%
%  [pdistr,phir,energyr] = mms.psd_rebin(pdist,phi,energy0,energy1,stepTable)
%
%  Functions takes the burst mode distribution sampled in two energy tables
%  and converts to a single energy table with 64 energy channels. Time
%  resolution is halved and phi angles are averaged over adjacent times.
%
%  Input:
%       pdist - particle distribution in TSeries format
%       phi - phi angles in TSeries
%       energy0 - energy table 0 (structure or array)
%       energy1 - energy table 1 (structure or array)
%       stepTable - TSeries of stepping table between energies (burst)
%
%  Output:
%       pdistr - rebinned particle distribution
%       phir - recalculated phi angle TSeries
%       energyr - revised energy table
%
%  Notes: I'm assuming no gaps in the burst data interval. If there is a
%  gap use tlim before running. To be updated later.
%
% Written by D. B. Graham

tic;
if isstruct(energy0)
  energy0 = energy0.data;
end
if isstruct(energy1)
  energy1 = energy1.data;
end

stepTable = stepTable.data;

% Sort energy levels
energyr = [energy0 energy1];
energyr = sort(energyr);

% Define new times
deltat = median(diff(pdist.time.epochUnix))/2;
newtimes = pdist.time(1:2:end-1)+deltat;
pdistr = zeros(length(newtimes),64,32,16);
phir = zeros(length(newtimes),32);
newelnum = 1;

phis = circshift(phi.data,1,2);
phis(:,1) = phis(:,1)-360;

for ii=1:2:length(pdist.time)-1
  if phi.data(ii,1) > phi.data(ii+1,1)
    phir(newelnum,:) = (phi.data(ii,:)+phis(ii+1,:))/2;
    pdisttemp = circshift(squeeze(pdist.data(ii+1,:,:,:)),1,2);

    if stepTable(ii)
      pdistr(newelnum, 2:2:64,:,:) = pdist.data(ii,:,:,:);
      pdistr(newelnum, 1:2:63,:,:) = pdisttemp;
    else
      pdistr(newelnum, 1:2:63,:,:) = pdist.data(ii,:,:,:);
      pdistr(newelnum, 2:2:64,:,:) = pdisttemp;
    end
  else
    phir(newelnum,:) = (phi.data(ii,:)+phi.data(ii+1,:))/2;

    if stepTable(ii)
      pdistr(newelnum,2:2:64,:,:) = pdist.data(ii,:,:,:);
      pdistr(newelnum,1:2:63,:,:) = pdist.data(ii+1,:,:,:);
    else
      pdistr(newelnum,1:2:63,:,:) = pdist.data(ii,:,:,:);
      pdistr(newelnum,2:2:64,:,:) = pdist.data(ii+1,:,:,:);
    end
  end

  newelnum = newelnum+1;
end

phir = TSeries(newtimes,phir);
pdistr = TSeries(newtimes,pdistr);
toc;

end

