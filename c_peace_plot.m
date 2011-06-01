function h = c_peace_plot(hin,peace_spec,i_spectra)
%C_PEACE_PLOT  plot standard PEACE spectrograms (parallel, perp, antiparallel, ..)
%
% h = c_peace_plot(peace_spec)
% h = c_peace_plot(peace_spec,i_spectra)
% h = c_peace_plot(hin,peace_spec,i_spectra)
%
% Input: 
%        hin: handle of panels in which to plot
%     peace_spec: PEACE spectrograms
%     i_spectra: the number of spectra to plot
%               if 0 line-plot angles, if -1 plot all 
%
% Output: 
%      h: handle of all panels
%
%    See also C_PEACE_SPECTRA
%
% $Id$

error(nargchk(1,3,nargin))
if nargin==0,
    help c_peace_plot;
elseif nargin==1, % assume c_peace_plot(peace_spec)
    peace_spec=hin;hin=[];i_spectra=1;
elseif nargin==2, % assume c_peace_plot(peace_spec,i_spectra)
    i_spectra=peace_spec;peace_spec=hin;hin=[];
elseif nargin==3, 
else
    help c_peace_plot;return
end

if i_spectra==-1,
    i_spectra=length(peace_spec.p):-1:1;
end

nsubplots=length(i_spectra);

if length(hin)==nsubplots,
    h=hin;
else
    h=1:nsubplots; % preallocate h
    for jj=1:nsubplots,
        h(jj)=irf_subplot(nsubplots,1,-jj);
    end
end

for jj=1:nsubplots,
    if i_spectra(jj)==0, % plot angles 
        disp('assumes that peace spectra have par/perp/antipar');
        axes(h(jj));cla; %#ok<LAXES>
        htmp=irf_plot({[peace_spec.t peace_spec.pa{1}(:,1)],[peace_spec.t peace_spec.pa{2}(:,1)],[peace_spec.t peace_spec.pa{3}(:,1)]},'comp','linestyle','.');
        ylabel('pitch angle [deg]');
        axis tight
        set(htmp,'ylim',[0 180],'ytick',[0 30 60 90 120 150 180])
    else
        peace_spec_comp=peace_spec;
        peace_spec_comp.p=peace_spec.p(i_spectra(jj));
        peace_spec_comp.p_label=peace_spec.p_label(i_spectra(jj));
        hspec=irf_spectrogram(h(jj),peace_spec_comp);
        if strcmp(peace_spec_comp.f_unit,'eV'),
            set(hspec,'xticklabel',[]);
            set(hspec,'Yscale','log');
            set(gca,'ytick',[10 100 1000 10000])
            grid on;
        else
           set(hspec,'xticklabel',[]);
            set(hspec,'Yscale','lin');
            set(gca,'ytick',[30 60 90 120 150 180 210 240 270 300 330])
            grid on;
        end
        if jj==1,
            hc=colorbar('peer',h(jj));
            caxis_save = caxis(h(jj)); % all panels have common colorbar axis
        else
            caxis(h(jj),caxis_save);
            hc=colorbar('peer',h(jj));
        end
        ylabel(hc,peace_spec_comp.p_label);
    end
end

irf_figmenu
irf_timeaxis

