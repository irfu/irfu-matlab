function hout = caa_spectrogram(h,specrec,Pxx,F)
%CAA_SPECTROGRAM  plot power spectrum
%
% [h] = caa_spectrogram([h],specrec)
% [h] = caa_spectrogram([h],t,Pxx,F)
%
% See also CAA_SPECTROGRAM
%
% $Id$

% Copyright 2005,2006 Yuri Khotyaintsev

error(nargchk(1,4,nargin))

if nargin==1, specrec = h; h = [];
elseif nargin==3, specrec.p = specrec; specrec.f = Pxx; specrec.t = h; h = [];
elseif nargin==4, specrec.t = specrec; specrec.p = Pxx; specrec.f = F; 
end

ndata = length(specrec.t);
if ndata<1, if nargout>0, hout=h; end, return, end
ncomp = length(specrec.p);

load caa/cmap.mat

if isempty(h), clf, for comp=1:ncomp, h(comp) = irf_subplot(ncomp,1,-comp); end, end

% If H is specified, but is shorter than NCOMP, we plot just first 
% length(H) spectra
for comp=1:min(length(h),ncomp)
	
	for jj=1:ndata
		specrec.p{comp}(jj,find(isnan(specrec.p{comp}(jj,:)))) = 1e-15;
	end
	
    ud=get(gcf,'userdata');
	ii = find(~isnan(specrec.t));
	if isfield(ud,'t_start_epoch'), 
		t_start_epoch = ud.t_start_epoch;
	elseif specrec.t(ii(1))> 1e8, 
		% Set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
		t_start_epoch = specrec.t(ii(1));
		ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
		irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
	else
		t_start_epoch = 0;
	end

	axes(h(comp));
	
	% Special case when we have only one spectrum
	% We duplicate it
	if ndata==1
		dt = .5/specrec.f(2);
		specrec.t = [specrec.t-dt; specrec.t+dt];
		specrec.p(comp) = {[specrec.p{comp}; specrec.p{comp}]};
	end
	pcolor(specrec.t-t_start_epoch,specrec.f,log10(specrec.p{comp}'))
	
	colormap(cmap)
	shading flat
%	colorbar('vert')
%	set(gca,'TickDir','out','YScale','log')
	set(gca,'TickDir','out')
	ylabel('frequency [Hz]')
	if comp==min(length(h),ncomp), add_timeaxis;
	else set(gca,'XTicklabel','')
	end
end

if nargout>0, hout=h; end
