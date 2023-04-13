function [phiEB,vbest,dirbest,thetas,corrs] = lhwaveanalysis(varargin)
% LHWAVEANALYSIS calculate lower-hybrid wave properties from MMS data
%
% [phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(Tint,Exyz,Bscm,Bxyz,ne,option,optionvalue)
%
% Input:
%   Tint - time interval
%   Exyz - Electric field (TSeries)
%   Bscm - Search coil magnetic field (TSeries)
%   Bxyz - Background magnetic field (TSeries)
%   ne - number density (TSeries)
%
%   Options:
%       'lhfilt' - Filter for LH fluctuations. For one element it is the minimum
%       frequency in the highpass filter. For two elements the fields are
%       bandpassed between the frequencies.
%       'blpass' - set maximum frequency for low-pass filter of background
%       magnetic field (FGM)
%       'plot' - set to 1 to plot figure.
%       'vmax' - maximum speed (km/s), if vbest = vmax, increase interval
%
% Example:
%   Tintl = irf.tint('2015-12-14T01:17:39.00Z/2015-12-14T01:17:43.00Z');
%   Bxyz=mms.db_get_ts('mms2_fgm_brst_l2','mms2_fgm_b_gse_brst_l2',Tintl);
%   Exyz=mms.db_get_ts('mms2_edp_brst_l2_dce','mms2_edp_dce_gse_brst_l2',Tintl);
%   Bscm=mms.db_get_ts('mms2_scm_brst_l2_scb','mms2_scm_acb_gse_scb_brst_l2',Tintl);
%   ne = mms.db_get_ts('mms2_fpi_brst_l2_des-moms','mms2_des_numberdensity_brst',Tintl);
%   Tint = irf.tint('2015-12-14T01:17:40.20Z/2015-12-14T01:17:41.50Z');
%   [phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(Tint,Exyz,Bscm,Bxyz,ne,'lhfilt',[5 100],'blpass',5,'plot',1);
%
% Written by D. B. Graham. Based in part on irf_match_phibe_dir,
% irf_match_phibe_v, and irf_match_phibe_vis.

% Check input
if (nargin < 5)
  help lhwaveanalysis;
  return;
end

tints = varargin{1};
Exyz = varargin{2};
Bscm = varargin{3};
Bxyz = varargin{4};
ne = varargin{5};

ic = Bscm.name(4);

args=varargin(6:end);
if numel(args)>0
  options=1;
else
  options=0;
end

% Default bandpasses
minfreq = 10;
maxfreq = 0; % No low-pass filter
lowpassBxyz = 2;
plotfigure = 0;
frange = 0;
vmax = 2e3; % 2000 km/s

while options
  l = 2;
  switch(lower(args{1}))
    case 'lhfilt'
      if numel(args)>1 && isnumeric(args{2})
        if length(args{2}) == 1
          minfreq = args{2};
        elseif length(args{2}) == 2
          minfreq = args{2}(1);
          maxfreq = args{2}(2);
          frange = 1;
        else
          irf.log('critical','lfbandpass not recognized.');
          return;
        end
      end
    case 'blpass'
      if numel(args)>1 && isnumeric(args{2})
        lowpassBxyz = args{2};
      end
    case 'plot'
      if numel(args)>1 && isnumeric(args{2})
        if args{2}
          plotfigure = 1;
        end
      end
    case 'vmax'
      vmax = args{2};
    otherwise
      irf.log('warning',['Unknown flag: ' args{1}])
      l=1;
      break
  end
  args = args(l+1:end);
  if isempty(args), options=0; end
end

% Bandpass filter data
fBxyz = 1/median(diff(Bxyz.time.epochUnix));
Bxyz = irf_filt(Bxyz,0,lowpassBxyz,fBxyz,5);
Exyz = Exyz.resample(Bscm);
ne = ne.resample(Bscm);
Bxyz = Bxyz.resample(Bscm);
Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);

fBscm = 1/median(diff(Bscm.time.epochUnix));
Bscmfac = irf_filt(Bscmfac,minfreq,maxfreq,fBscm,5);
Exyz = irf_filt(Exyz,minfreq,maxfreq,fBscm,5);

Units = irf_units; % Use IAU and CODATA values for fundamental constants.
qe = Units.e;
mu = Units.mu0;
Bmag = Bxyz.abs.data;
phiB = (Bscmfac.data(:,3)).*Bmag*1e-18./(ne.data*qe*mu*1e6);
phiB = irf.ts_scalar(Bscmfac.time,phiB);

% short buffer so phi_E does not begin at zero.
tint = tints+[-0.2 0.2];
Exyz = Exyz.tlim(tint);
phiBs = phiB.tlim(tints);

%Rotate Exyz into field-aligned coordinates
Bxyzs = Bxyz.tlim(tints);
Bmean = mean(Bxyzs.data,1);
Bvec = Bmean/norm(Bmean);
Rtemp = [1 0 0];
R2 = cross(Bvec,Rtemp);
R2 = R2/sqrt(R2(1)^2+R2(2)^2+R2(3)^2);
R1 = cross(R2,Bvec);
ER1 = Exyz.data(:,1)*R1(1)+Exyz.data(:,2)*R1(2)+Exyz.data(:,3)*R1(3);
ER2 = Exyz.data(:,1)*R2(1)+Exyz.data(:,2)*R2(2)+Exyz.data(:,3)*R2(3);
ER3 = Exyz.data(:,1)*Bvec(1)+Exyz.data(:,2)*Bvec(2)+Exyz.data(:,3)*Bvec(3);
Efac = irf.ts_vec_xyz(Exyz.time,[ER1 ER2 ER3]);

% Find best direction
thetas = 0:1:360;
corrs = zeros(1,length(thetas));

for ii = 1:length(thetas)
  Etemp = cosd(thetas(ii))*Efac.data(:,1)+sind(thetas(ii))*Efac.data(:,2);
  Etemp = TSeries(Exyz.time,Etemp);
  phitemp = irf_integrate(Etemp);
  phitemp = phitemp.tlim(tints);
  phitemp.data = phitemp.data-mean(phitemp.data);
  corrs(ii) = xcorr(phiBs.data,phitemp.data,0,'coeff');
end

[~,corrpos] = max(corrs);
Ebest = cosd(thetas(corrpos))*Efac.data(:,1)+sind(thetas(corrpos))*Efac.data(:,2);
Ebest = TSeries(Exyz.time,Ebest);
phibest = irf_integrate(Ebest);
phibest = phibest.tlim(tints);
phibest.data = phibest.data-mean(phibest.data);
thetabest = thetas(corrpos);
dirbest = R1*cosd(thetabest)+R2*sind(thetabest);
dirbestround = round(dirbest,2);

%Find best speed
vphvec = 1e1:1e0:vmax; % Maximum velocity may need to be increased in rare cases
corrv = zeros(1,length(vphvec));

for ii=1:length(vphvec)
  phiEtemp = phibest.data*vphvec(ii);
  corrv(ii)=sum(abs(phiEtemp-phiBs.data).^2);
end

[~,corrvpos] = min(corrv);
if corrvpos == length(vphvec)
  irf.log('warning',sprintf('Wave speed > vmax = %g km/s. Increase search interval using input option ''vmax''.',vmax));
end
phiEbest = phibest.data*vphvec(corrvpos);
phiEbest = TSeries(phiBs.time,phiEbest);
vbest = vphvec(corrvpos);

phiEB = TSeries(phiBs.time,[phiEbest.data phiBs.data]);

phiEmax = max(abs(phiEB.data(:,1)));
phiBmax = max(abs(phiEB.data(:,2)));

if plotfigure
  fn=figure;
  set(fn,'Position',[10 10 600 600])
  h(1)=axes('position',[0.1 0.58 0.8 0.4]);
  h(2)=axes('position',[0.1 0.07 0.8 0.4]);
  ud=get(fn,'userdata');
  ud.subplot_handles=h;
  set(fn,'userdata',ud);
  set(fn,'defaultLineLineWidth',2);
  
  h(1)=irf_panel('phi');
  irf_plot(h(1),phiEB);
  ylabel(h(1),'\phi (V)','Interpreter','tex','fontsize',14);
  irf_legend(h(1),{'\phi_{E}','\phi_{B}'},[0.1 0.12],'fontsize',14)
  if frange
    freqlab = [num2str(minfreq) ' Hz < f < ' num2str(maxfreq) ' Hz'];
  else
    freqlab = ['f > ' num2str(minfreq) ' Hz'];
  end
  vlab = ['v = ' num2str(vbest) ' km s^{-1}'];
  xdir = num2str(dirbestround(1));
  ydir = num2str(dirbestround(2));
  zdir = num2str(dirbestround(3));
  dirlab = ['dir: [' xdir ',' ydir ',' zdir ']'];
  irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',14)
  c_eval('irf_legend(h(1),''MMS ?'',[0.01 0.98],''color'',''k'',''fontsize'',14)',ic);
  irf_legend(h(1),freqlab,[0.95 0.2],'color','k','fontsize',14)
  irf_legend(h(1),vlab,[0.95 0.13],'color','k','fontsize',14)
  irf_legend(h(1),dirlab,[0.95 0.06],'color','k','fontsize',14)
  irf_legend(h(1),['|\phi_E|_{max} = ' num2str(round(phiEmax,1))],[0.90 0.98],'color','k','fontsize',14)
  irf_legend(h(1),['|\phi_B|_{max} = ' num2str(round(phiBmax,1))],[0.90 0.90],'color','k','fontsize',14)
  irf_zoom(h(1),'x',tints);
  
  plot(h(2),thetas,corrs);
  [maxcorr,ind] = max(corrs);
  hold(h(2),'on')
  plot(h(2),thetas(ind),maxcorr,'ro');
  hold(h(2),'off')
  xlabel(h(2),'\theta (deg)','Interpreter','tex','fontsize',14);
  ylabel(h(2),'C_{\phi}','Interpreter','tex','fontsize',14);
  axis(h(2),[0 360 -1 1]);
  irf_legend(h(2),['C_{\phi max} = ' num2str(round(maxcorr,2))],[0.95 0.06],'color','k','fontsize',14)
  irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',14)
  
  set(h(1:2),'fontsize',14);
  set(gcf,'color','w');
  xtickangle(h,0)

end

end
