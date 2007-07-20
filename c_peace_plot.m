function h= c_peace_plot(peace_spec,peace_pa)
%C_PEACE_PLOT  plot standard PEACE spectrograms (parallel, perp, antiparallel, ..)
%
% h = c_peace_plot(peace_spec,peace_pa)
%
% Input: 
%     peace_spec: PEACE spectrograms
%     peace_pa: PEACE pitch angles
%
% Output: 
%      hout: handle of all panels
%
%    See also C_PEACE_SPECTRA
%
% $Id$

if exist('peace_pa','var'),
    nsubplots=length(peace_spec.p)+1;
else
    nsubplots=length(peace_spec.p);
end

for jj=1:nsubplots,
    h(jj)=irf_subplot(nsubplots,1,-jj);
end

hspec=caa_spectrogram(h(1:length(peace_spec.p)),peace_spec);
set(hspec,'Yscale','log');
axes(hspec(1));
hc=colorbar;cax=caxis;ylabel(hc,peace_spec.p_label{1});
for jj=2:length(hspec),
    axes(hspec(jj));
    caxis(cax);hc=colorbar;ylabel(hc,peace_spec.p_label{jj});
end

if exist('peace_pa','var'),
    axes(h(end)); % plot angles
    if jj==3, % plot angles only when plotting par,perp,antipar (lazy programming)
        colorbar;hh=get(gca,'position');colorbar off;set(gca,'position',hh);
        htmp=irf_plot({[peace_spec.t peace_pa{1}(:,1)],[peace_spec.t peace_pa{2}(:,1)],[peace_spec.t peace_pa{3}(:,1)]},'comp','linestyle','.');
        ylabel('pitch angle [deg]'); 
        set(htmp,'ylim',[0 180],'ytick',[0 30 60 90 120 150 180])
    end
end

irf_figmenu

