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
%              specrec.dt - vector of dt interval for every t point (can be omitted)
%              specrec.df - vector of dF interval for every frequency f point (can be omitted)
%                           df can be structure with two vectors df.plus and df.minus
%         specrec.f_label - label of f axis
%         specrec.p_label - label of colorbar
%
% See also CAA_POWERFFT
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

disp('');
disp('WARNING!!!!!')
disp('');
disp('caa_spectrogram is replaced with irf_spectrogram');
disp('caa_spectrogram will be removed in the near future');
disp('');

narginchk(1,6)

if nargin==1
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
  else, specrec.p = {t};
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
  else, specrec.p = {Pxx};
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
  else, specrec.p = {Pxx};
  end
  specrec.f = F;
  specrec.dt = dt;
  specrec.df=[];
elseif	nargin==6 % caa_spectrogram(h,t,Pxx,F,dt,df)
  specrec.t = t;
  if (size(Pxx,1) ~= length(t)) && (size(Pxx,2) == length(t)), Pxx = Pxx'; end
  if iscell(Pxx),specrec.p = Pxx;
  else, specrec.p = {Pxx};
  end
  specrec.f = F;
  specrec.dt = dt;
  specrec.df = dF;
end

specrec.t = double(specrec.t);
specrec.f = double(specrec.f);
%specrec.dt = double(specrec.dt);
%specrec.df = double(specrec.df);
if iscell(specrec.p)
  ncomp=length(specrec.p);
elseif isnumeric(specrec.p)
  ncomp=1;
  specrec.p={specrec.p};
else
  disp('WARNING: cannot intepret input parameters in caa_spectrogram, returning.')
  return
end

ndata = length(specrec.t);
if ndata<1, if nargout>0, hout=h; end, return, end

load caa/cmap.mat

if isempty(h), clf, for comp=1:ncomp, h(comp) = irf_subplot(ncomp,1,-comp); end, end

% If H is specified, but is shorter than NCOMP, we plot just first
% length(H) spectra
for comp=1:min(length(h),ncomp)
  
  for jj=1:ndata
    %		specrec.p{comp}(jj,isnan(specrec.p{comp}(jj,:))) = 1e-15;
    specrec.p{comp}(jj,isnan(specrec.p{comp}(jj,:))) = NaN;
  end
  
  ud = get(gcf,'userdata');
  ii = find(~isnan(specrec.t));
  if isfield(ud,'t_start_epoch')
    t_start_epoch = double(ud.t_start_epoch);
  elseif specrec.t(ii(1))> 1e8
    % Set start_epoch if time is in isdat epoch
    % Warn about changing t_start_epoch
    t_start_epoch = double(specrec.t(ii(1)));
    ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
    irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
  else
    t_start_epoch = double(0);
  end
  
  % Special case when we have only one spectrum
  % We duplicate it
  if ndata==1
    specrec.dt = double(.5/specrec.f(2));
    %		specrec.t = [specrec.t-dt; specrec.t+dt];
    %		specrec.p(comp) = {[specrec.p{comp}; specrec.p{comp}]};
  end
  if ~isfield(specrec,'f_unit') % if not specified assume units are Hz
    if max(specrec.f) > 2000 % check whether to use kHz
      specrec.f=specrec.f*double(1e-3);
      specrec.f_unit='kHz';
    else
      specrec.f_unit='Hz';
    end
    if ~isfield(specrec,'f_label')
      specrec.f_label=['frequency [' specrec.f_unit ']'];
    end
  end
  
  if min(size(specrec.f))==1, ff=double(specrec.f(:))';
  else
    ff=double(specrec.f);
  end % if f vector make it row vector
  tt=double(specrec.t(:));
  pp=specrec.p{comp};
  if isempty(specrec.df) % if frequency steps are not given
    fnew=[ff ff];
    fnew(1)=ff(1)-0.5*(ff(2)-ff(1));
    fnew(end)=ff(end)+0.5*(ff(end)-ff(end-1));
    fnew(2:2:end-1)=0.5*(ff(1:end-1)+ff(2:end));
    fnew(3:2:end-1)=0.5*(ff(1:end-1)+ff(2:end));
    ff=fnew;
  else                   % if frequency steps are given
    if isstruct(specrec.df) % df.plus and df.minus should be specified
      dfplus=torow(double(specrec.df.plus(:))); % if df vector make it row vector
      dfminus=torow(double(specrec.df.minus(:))); % if df vector make it row vector
    else
      dfplus=torow(double(specrec.df(:))); % if df vector make it row vector
      dfminus=dfplus;
    end
    fnew=[ff ff];
    jj=1:size(ff,2);
    fnew(:,jj*2-1)=ff-dfminus;
    fnew(:,jj*2)=ff+dfplus;
    ff=fnew;
  end
  jj=1:size(pp,2);
  ppnew=[pp pp];
  ppnew(:,jj*2-1)=pp;
  ppnew(:,jj*2)=NaN;
  pp=ppnew;
  if ~isempty(specrec.dt) % if time steps are not given
    if isstruct(specrec.dt) % df.plus and df.minus should be specified
      dtplus=torow(double(specrec.dt.plus(:))); % if df vector make it row vector
      dtminus=torow(double(specrec.dt.minus(:))); % if df vector make it row vector
    else
      dtplus=tocolumn(double(specrec.dt(:))); % if df vector make it row vector
      dtminus=dtplus;
    end
    ttnew=[tt; tt];
    jj=1:length(tt);
    ttnew(jj*2-1)=tt-dtminus;
    ttnew(jj*2)=tt+dtplus;
    tt=ttnew;
    ppnew=[pp;pp];
    ppnew(jj*2-1,:)=pp;
    ppnew(jj*2,:)=NaN;
    pp=ppnew;
  end
  
  if min(size(ff))==1 % frequency is vector
    if any(min(pp)<0) % spectra include negative values linear spectrogram
      pcolor(h(comp),double(tt-t_start_epoch),ff,double(pp'))
    else
      pcolor(h(comp),double(tt-t_start_epoch),ff,log10(double(pp')))
    end
  else % frequency is matrix
    ttt = repmat(tt,1,size(ff,2));
    if any(min(pp)<0) % spectra include negative values linear spectrogram
      pcolor(h(comp),double(ttt-t_start_epoch),ff,double(pp))
    else
      pcolor(h(comp),double(ttt-t_start_epoch),ff,log10(double(pp)))
    end
  end
  %	colormap(cmap)
  shading(h(comp),'flat')
  %	colorbar('vert')
  %	set(gca,'TickDir','out','YScale','log')
  set(h(comp),'TickDir','out')
  %check ylabel
  if ~isfield(specrec,'f_label')
    if ~isfield(specrec,'f_unit')
      specrec.f_unit='a.u.';
    end
    specrec.f_label=['[' specrec.f_unit ']'];
  end
  ylabel(h(comp),specrec.f_label)
  
  if isfield(specrec,'p_label')
    if isa(h(comp),'handle'), hcb = colorbar(h(comp)); % HG2
    else, hcb = colorbar('peer',h(comp));
    end
    ylabel(hcb,specrec.p_label);
    irf_colorbar_fit_label_height(hcb);
  end
  if comp==min(length(h),ncomp), irf_timeaxis;
  else, set(h(comp),'XTicklabel','')
  end
end

if nargout>0, hout=h; end
