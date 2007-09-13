function hout = caa_spectrogram(h,t,Pxx,F,dt,dF)
%CAA_SPECTROGRAM  plot power spectrum in logarithimic scale
%
% [h] = caa_spectrogram([h],specrec)
% [h] = caa_spectrogram([h],t,Pxx,[f],[dt],[df])
%
% Input:
%    specrec - structure including spectra
%              specrec.t  - time vector
%              specrec.f - frequency vector 
%              specrec.p - spectral density matrix (size(t)xsize(f))
%              specrec.dt - vector of dt interval for every t point (can be ommitted)
%              specrec.df - vector of dF interval for every frequency f point (can be ommitted)
%
% See also CAA_POWERFFT
%
% $Id$

% Copyright 2005-2007 Yuri Khotyaintsev

error(nargchk(1,6,nargin))

if nargin==1, 
    specrec = h; h = [];
    if ~isfield(specrec,'dt'), specrec.dt=[];end
    if ~isfield(specrec,'df'), specrec.df=[];end
elseif nargin==2 % caa_spectrogram(h,t)
	specrec = t;
    if ~isfield(specrec,'dt'), specrec.dt=[];end
    if ~isfield(specrec,'df'), specrec.df=[];end
elseif nargin==3 % caa_spectrogram(t,Pxx,F)
	if size(t,2) == length(h), t = t'; end
	if iscell(t), specrec.p = t;
	else specrec.p = {t};
	end
	specrec.f = Pxx; 
	specrec.t = h; 
    if length(specrec.t)>1 % assume equidistant times 
        specrec.dt=(specrec.t(2)-specrec.t(1))/2;
    else
        specrec.dt=[]; % will be calculated later
    end
    specrec.df=[];
	h = [];
elseif	nargin==4 % caa_spectrogram(h,t,Pxx,F)
	specrec.t = t;
	if (size(Pxx,1) ~= length(t)) && (size(Pxx,2) == length(t)), Pxx = Pxx'; end
	if iscell(Pxx),specrec.p = Pxx;
	else specrec.p = {Pxx};
	end
	specrec.f = F;
    if length(specrec.t)>1 % assume equidistant times 
        specrec.dt=(specrec.t(2)-specrec.t(1))/2;
    else
        specrec.dt=[]; % will be calculated later
    end
    specrec.df=[];
elseif	nargin==5 % caa_spectrogram(h,t,Pxx,F,dt)
	specrec.t = t;
	if (size(Pxx,1) ~= length(t)) && (size(Pxx,2) == length(t)), Pxx = Pxx'; end
	if iscell(Pxx),specrec.p = Pxx;
	else specrec.p = {Pxx};
	end
	specrec.f = F;
	specrec.dt = dt;
    specrec.df=[];
elseif	nargin==6 % caa_spectrogram(h,t,Pxx,F,dt,df)
	specrec.t = t;
	if (size(Pxx,1) ~= length(t)) && (size(Pxx,2) == length(t)), Pxx = Pxx'; end
	if iscell(Pxx),specrec.p = Pxx;
	else specrec.p = {Pxx};
	end
	specrec.f = F;
	specrec.dt = dt;
	specrec.df = dF;
end

specrec.t = double(specrec.t);
specrec.f = double(specrec.f);
specrec.dt = double(specrec.dt);
specrec.df = double(specrec.df);

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
        specrec.dt = double(.5/specrec.f(2));
        %		specrec.t = [specrec.t-dt; specrec.t+dt];
        %		specrec.p(comp) = {[specrec.p{comp}; specrec.p{comp}]};
    end
    if ~isfield(specrec,'f_unit'), % if not specified assume units are Hz
        if max(specrec.f) > 2000, % check whether to use kHz
            specrec.f=specrec.f*double(1e-3);
            specrec.f_unit='kHz';
        else
            specrec.f_unit='Hz';
        end
        if ~isfield(specrec,'f_label')
            specrec.f_label=['frequency [' specrec.f_unit ']'];
        end
    end
    
    if min(size(specrec.f))==1, ff=double(specrec.f(:))';end % if f vector make it row vector
    tt=double(specrec.t(:));
    pp=specrec.p{comp};
    if ~isempty(specrec.df) % if frequency steps are given
        if min(size(specrec.df))==1, df=double(specrec.df(:))';end % if df vector make it row vector
        fnew=[ff ff];
        jj=1:length(ff);
        fnew(jj*2-1)=ff-df;
        fnew(jj*2)=ff+df;
        ff=fnew;
        ppnew=[pp pp];
        ppnew(:,jj*2-1)=pp;
        ppnew(:,jj*2)=NaN;
        pp=ppnew;
    end
    if ~isempty(specrec.dt) % if time steps are not given
        dt=double(specrec.dt(:));
        ttnew=[tt; tt];
        jj=1:length(tt);
        ttnew(jj*2-1)=tt-dt;
        ttnew(jj*2)=tt+dt;
        tt=ttnew;
        ppnew=[pp;pp];
        ppnew(jj*2-1,:)=pp;
        ppnew(jj*2,:)=NaN;
        pp=ppnew;
    end
    pcolor(double(tt-t_start_epoch),ff,double(log10(pp')))
    
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
