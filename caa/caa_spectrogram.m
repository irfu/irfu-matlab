function caa_spectrogram(h,Pxx,F)
%CAA_SPECTROGRAM  plot fft
%
% caa_spectrogram(h,Pxx,F)
%
% See also CAA_SPECTROGRAM
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

error(nargchk(2,3,nargin))

ndata = size(Pxx,1);
ncomp = size(Pxx,2);
nf = size(Pxx,3);
C = reshape(Pxx(:,2,:),ndata,nf);
t = Pxx(:,1,1); t = t(:);
mm = min(min(C));

for jj=1:ndata
	C(jj,find(isnan(C(jj,:)))) = mm;
end
ud=get(gcf,'userdata');
ii = find(~isnan(t));
if isfield(ud,'t_start_epoch'), 
	t_start_epoch = ud.t_start_epoch;
elseif t(ii(1))> 1e8, % set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
	t_start_epoch = t(ii(1));
	ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
	irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
else
	t_start_epoch = 0;
end

	
pcolor(t-t_start_epoch,F,log10(C'))
load caa/cmap.mat
colormap(cmap)
shading flat
colorbar('vert')
add_timeaxis
set(gca,'TickDir','out')
ylabel('frequency [Hz]')
