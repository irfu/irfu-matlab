function hout = c_peace_plot(hin,peace_spec,peace_pa,i_spectra)
%C_PEACE_PLOT  plot standard PEACE spectrograms (parallel, perp, antiparallel, ..)
%
% h = c_peace_plot(peace_spec,peace_pa)
%
% Input: 
%     peace_spec: PEACE spectrograms
%     peace_pa: PEACE pitch angles
%     i_spectra: which spectra to plot, if 0 plot angles
%
% Output: 
%      hout: handle of all panels
%
%    See also C_PEACE_SPECTRA
%
% $Id$

error(nargchk(1,4,nargin))
if nargin==0,
    help c_peace_plot;
elseif nargin==1, % assume only peace spectrogram as input
    peace_spec=hin;hin=[];
elseif nargin==2,
elseif nargin==3,
    i_spectra=1:length(peace_spec);
elseif nargin==4,
else
    help c_peace_plot;return
end

nsubplots=length(i_spectra);

if length(hin==nsubplots),
    h=hin;
else
    for jj=1:nsubplots,
        h(jj)=irf_subplot(nsubplots,1,-jj);
    end
end

for jj=1:nsubplots,
    if i_spectra(jj)==0, % plot angles 
        disp('assumes that peace spectra have par/perp/antipar');
        axes(h(jj));cla;
        htmp=irf_plot({[peace_spec.t peace_pa{1}(:,1)],[peace_spec.t peace_pa{2}(:,1)],[peace_spec.t peace_pa{3}(:,1)]},'comp','linestyle','.');
        ylabel('pitch angle [deg]');
        set(htmp,'ylim',[0 180],'ytick',[0 30 60 90 120 150 180])
    else,
        peace_spec_comp=peace_spec;
        peace_psec_comp.p=peace_spec.p(i_spectra(jj));
        peace_psec_comp.p_label=peace_spec.p_label(i_spectra(jj));
        hspec=caa_spectrogram(h(jj),peace_spec_comp);
        set(hspec,'xtick',[]);
        set(hspec,'Yscale','log');
        if jj==1,
            hc=colorbar;cax=caxis; % all panels have common colorbar axis
        else
            caxis(cax);hc=colorbar;
        end
        ylabel(hc,peace_spec_comp.p_label);
    end
end

irf_figmenu

