function hout = caa_spectrogram(h,t,Pxx,F)
%CAA_SPECTROGRAM  plot power spectrum
%
% [h] = caa_spectrogram([h],specrec)
% [h] = caa_spectrogram([h],t,Pxx,F)
%
% See also CAA_POWERFFT
%
% $Id$

% Copyright 2005-2007 Yuri Khotyaintsev

error(nargchk(1,4,nargin))

if nargin==1, specrec = h; h = [];
elseif nargin==2
	specrec = t;
elseif nargin==3
	if size(t,2) == length(h), t = t'; end
	if iscell(t), specrec.p = t;
	else specrec.p = {t};
	end
	specrec.f = Pxx; 
	specrec.t = h; 
	h = [];
elseif	nargin==4
	specrec.t = t;
	if size(Pxx,2) == length(t), Pxx = Pxx'; end
	if iscell(Pxx),specrec.p = Pxx;
	else specrec.p = {Pxx};
	end
	specrec.f = F;
end

specrec.t = double(specrec.t);
specrec.f = double(specrec.f);

ndata = length(specrec.t);
if ndata<1, if nargout>0, hout=h; end, return, end
ncomp = length(specrec.p);

load caa/cmap.mat

if isempty(h), clf, for comp=1:ncomp, h(comp) = irf_subplot(ncomp,1,-comp); end, end

% If H is specified, but is shorter than NCOMP, we plot just first 
% length(H) spectra
for comp=1:min(length(h),ncomp)
	
	for jj=1:ndata
		specrec.p{comp}(jj,isnan(specrec.p{comp}(jj,:))) = 1e-15;
	end
	
    ud = get(gcf,'userdata');
	ii = find(~isnan(specrec.t));
	if isfield(ud,'t_start_epoch'), 
		t_start_epoch = double(ud.t_start_epoch);
	elseif specrec.t(ii(1))> 1e8, 
		% Set start_epoch if time is in isdat epoch
		% Warn about changing t_start_epoch
		t_start_epoch = double(specrec.t(ii(1)));
		ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
		irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
	else
		t_start_epoch = double(0);
	end

	axes(h(comp));
	
	% Special case when we have only one spectrum
	% We duplicate it
	if ndata==1
		dt = double(.5/specrec.f(2));
		specrec.t = [specrec.t-dt; specrec.t+dt];
		specrec.p(comp) = {[specrec.p{comp}; specrec.p{comp}]};
    end
    if ~isfield(specrec,'f_unit'), % if not specified assume units are Hz
        if max(specrec.f) > 2000, % check whether to use kHz 
            specrec.f=specrec.f*double(1e-3);
            specrec.f_unit='kHz';
        else,
            specrec.f_unit='Hz';
        end
        if ~isfield(specrec,'f_label')
            specrec.f_label=['frequency [' specrec.f_unit ']'];
        end
    end
    
	pcolor(double(specrec.t-t_start_epoch),double(specrec.f),...
		double(log10(specrec.p{comp}')))
	
	colormap(cmap)
    shading flat
    %	colorbar('vert')
    %	set(gca,'TickDir','out','YScale','log')
    set(gca,'TickDir','out')
    %check ylabel
    if ~isfield(specrec,'f_label')
        if ~isfield(specrec,'f_unit'),
            specrec.f_unit='a.u.';
        end
        specrec.f_label=['[' specrec.f_unit ']'];
    end
    ylabel(specrec.f_label)

    if comp==min(length(h),ncomp), add_timeaxis;
    else set(gca,'XTicklabel','')
    end
end

if nargout>0, hout=h; end
