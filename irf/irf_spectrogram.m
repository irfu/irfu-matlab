function [hout,hcb,hlpe] = irf_spectrogram(varargin)
%function [hout,hcb] = irf_spectrogram(h,t,Pxx,F,dt,dF)
%IRF_SPECTROGRAM  plot spectrogram
%
% [h, hcb] = irf_spectrogram(h,specrec,'option1','option2',..)
% [h, hcb] = irf_spectrogram(h,t,Pxx,[f],[dt],[df])
%
% Input:
%          h   - axis handle
%          hcb - colorbar handle
%    specrec - structure including spectra
%              specrec.t  - time vector (epoch unix)
%              specrec.f  - frequency vector (can be also matrix the size specrec.p)
%              specrec.p  - spectral density matrix (size(t)xsize(f))
%              specrec.dt - Description of the ~half-width of each separate spectrum in the plots. (Can be omitted)
%                           Specified in seconds. Multiple allowed forms:
%                           (1) Numeric value(s). Half-width of each spectrum.
%                               (a) scalar, or
%                               (b) 1D vector, one value for every spectrum.
%                           (2) Struct with fields .plus & .minus. Width before and after timestamp. The fields can be
%                               (a) scalar, or
%                               (b) 1D vector, one value for every spectrum.
%              specrec.df - vector of dF interval for every frequency f point (can be omitted)
%                           df can be structure with two vectors df.plus and df.minus
%                           (can be also matrix the size of specrec.p)
%              specrec.f_label   - label of f axis
%              specrec.p_label   - label of colorbar
%              specrec.plot_type - 'lin' or 'log'
%    options - 'lin' - plot spectrogram values in linear scale
%              'log' - (default) plot spectrogram values in log10 scale
%              'donotfitcolorbarlabel' - do not shrink colorbar label fonts to fit axes size
%              'donotshowcolorbar' - do not create a colorbar in plot
%              'donotassumetimeaxis' - do not assume that x axis is time, e.g. plotting pitch angle spectrograms
%
% See also IRF_POWERFFT

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

[h,args,nargs] = axescheck(varargin{:});
if isempty(h)
  fig = get(groot,'CurrentFigure'); h = get(fig,'Children');
end

%% Defaults
flagLog = true;            % want log10(data) dy default
f_multiplier = 1;          % default value using Hz units when units not specified, can be overwritten later if kHz makes labels more reasonable
fitColorbarLabel = true;   % fit font size of colorbar label to fit into axes size
showColorbar = true;
flagReducePlot = false;    % by default do not reduce spectrogram size to fit window pixels
useTStartEpoch = true;

%% Check input
if nargs==1 || ischar(args{2})    % irf_spectrogram(specrec,[options])
  specrec = args{1};
  if ~isfield(specrec,'dt'), specrec.dt=[];end
  if ~isfield(specrec,'df'), specrec.df=[];end
  for iArgs = 2:numel(args)
    flagValue = args{iArgs};
    if ischar(flagValue)
      switch lower(flagValue)
        case 'donotfitcolorbarlabel'
          fitColorbarLabel = false;
        case 'log'
          flagLog = true;
        case 'lin'
          flagLog = false;
        case 'reduce' % number of spectra equal to half axes size in pixels
          flagReducePlot = true;
        case 'donotshowcolorbar'
          showColorbar = false;
        case 'donotassumetimeaxis'
          useTStartEpoch = false;
        otherwise
          errStr= ['irf_spectrogram(), unknown flag:' flagValue];
          irf.log('critical',errStr);
          error('irf_spectrogram:unknown_flag',errStr);
      end
    end
  end
elseif nargs==3 % irf_spectrogram(t,Pxx,F)
  t=args{1};Pxx=args{2};F=args{3};
  if size(Pxx,2) == length(t), Pxx = Pxx'; end
  if iscell(Pxx), specrec.p = Pxx;
  else, specrec.p = {Pxx};
  end
  specrec.f = F;
  specrec.t = t;
  specrec.dt=[]; % will be calculated later
  specrec.df=[];
elseif	nargs==4 % irf_spectrogram(t,Pxx,F,dt)
  t=args{1};Pxx=args{2};F=args{3};dt=args{4};
  specrec.t = t;
  if (size(Pxx,1) ~= length(t)) && (size(Pxx,2) == length(t)), Pxx = Pxx'; end
  if iscell(Pxx),specrec.p = Pxx;
  else, specrec.p = {Pxx};
  end
  specrec.f = F;
  specrec.dt = dt;
  specrec.df=[];
elseif	nargin==5 % irf_spectrogram(t,Pxx,F,dt,df)
  t=args{1};Pxx=args{2};F=args{3};dt=args{4};dF=args{5};
  specrec.t = t;
  if (size(Pxx,1) ~= length(t)) && (size(Pxx,2) == length(t)), Pxx = Pxx'; end
  if iscell(Pxx),specrec.p = Pxx;
  else, specrec.p = {Pxx};
  end
  specrec.f = F;
  specrec.dt = dt;
  specrec.df = dF;
end

if isfield(specrec,'plot_type') && ...
    strcmpi(specrec.plot_type,'lin')
  flagLog = 0;
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
  disp('WARNING: cannot interpret input parameters in irf_spectrogram, returning.')
  return
end

ndata = length(specrec.t);
if ndata<1, if nargout>0, hout=h; end, return, end

%% Plot spectrogram
%load caa/cmap.mat

% Initiate figure if handles not given
if isempty(h)
  h=irf_plot(ncomp,'newfigure');
end

% If H is specified, but is shorter than NCOMP, we plot just first
% length(H) spectra
for comp=1:min(length(h),ncomp)
  
  specrec.p{comp}(isnan(specrec.p{comp})) = NaN; % WHY is this done? NaN = NaN already.
  
  ud = get(gcf,'userdata');
  ii = find(~isnan(specrec.t));
  if isfield(ud,'t_start_epoch')
    t_start_epoch = double(ud.t_start_epoch);
  elseif specrec.t(ii(1))> 1e8
    % Set start_epoch if time is in isdat epoch
    % Warn about changing t_start_epoch
    t_start_epoch = double(specrec.t(ii(1)));
    ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
    irf.log('notice',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
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
  if ~isfield(specrec,'f_unit') && ~isfield(specrec,'f_label') % if not specified assume units are Hz
    if max(specrec.f) > 2000 % check whether to use kHz
      specrec.f=specrec.f*double(1e-3);
      f_multiplier=1e-3;
      specrec.f_unit='kHz';
    else
      specrec.f_unit='Hz';
    end
    if ~isfield(specrec,'f_label')
      specrec.f_label=['f [' specrec.f_unit ']'];
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
    fnew(:,1)=ff(:,1)-0.5*(ff(:,2)-ff(:,1));
    fnew(:,end)=ff(:,end)+0.5*(ff(:,end)-ff(:,end-1));
    fnew(:,2:2:end-1)=0.5*(ff(:,1:end-1)+ff(:,2:end));
    fnew(:,3:2:end-1)=0.5*(ff(:,1:end-1)+ff(:,2:end));
    ff=fnew;
  else                   % if frequency steps are given
    if isstruct(specrec.df)                % if df is structure df.plus and df.minus should be specified
      if numel(specrec.df.plus)==1           % if df.plus is scalar
        dfplus=double(specrec.df.plus);   %    assign it
      elseif min(size(specrec.df.plus))==1   % if df.plus is vector
        dfplus=double(specrec.df.plus(:))'; %    make df.plus row vector
        dfplus=repmat(dfplus,size(ff,1),1); %    replicate to size of ff
      else                                    % if df.plus is matrix
        dfplus=double(specrec.df.plus);     %    assign it
      end
      if numel(specrec.df.minus)==1          % if df.minus is scalar
        dfminus=double(specrec.df.minus); %    assign it
      elseif min(size(specrec.df.minus))==1  % if df.minus is vector
        dfminus=double(specrec.df.minus(:))';%    make df.minus row vector
        dfminus=repmat(dfminus,size(ff,1),1);%    replicate to size of ff
      else                                    % if df.minus is matrix
        dfminus=double(specrec.df.minus);   %     assign it
      end
    else
      if min(size(specrec.df))==1       % if df is vector or scalar
        dfplus=double(specrec.df(:))'; %    make df row vector
      else                               % if df is matrix
        dfplus=double(specrec.df);     %    assign it
      end
      dfminus=dfplus;
    end
    dfplus=dfplus*f_multiplier;
    dfminus=dfminus*f_multiplier;
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
  if ~isempty(specrec.dt) % if time steps are given
    if isstruct(specrec.dt) % dt.plus and dt.minus should be specified
      dtplus=double(specrec.dt.plus(:)); % if dt vector make it column vector
      dtminus=double(specrec.dt.minus(:)); % if dt vector make it column vector
    else
      dtplus=double(specrec.dt(:)); % if dt vector make it column vector
      dtminus=dtplus;
    end
    ttnew=zeros(numel(tt),1);
    jj=1:length(tt);
    ttnew(jj*2-1)=tt-dtminus;
    ttnew(jj*2)=tt+dtplus;
    tt=ttnew;
  else
    ttnew=[tt;tt];
    ttnew(1)=tt(1)-0.5*(tt(2)-tt(1));
    ttnew(end)=tt(end)+0.5*(tt(end)-tt(end-1));
    ttnew(2:2:end-1)=0.5*(tt(1:end-1)+tt(2:end));
    ttnew(3:2:end-1)=0.5*(tt(1:end-1)+tt(2:end));
    tt=ttnew;
    
  end
  ppnew=[pp;pp];
  ppnew(1:2:end,:)=pp;
  ppnew(2:2:end,:)=NaN;
  pp=ppnew;
  if min(size(ff))~= 1 % ff is matrix
    ffnew=zeros(size(ff).*[2 1]);
    ffnew(1:2:end,:)=ff;
    ffnew(2:2:end,:)=ff;
    ff=ffnew;
  end
  
  tag=get(h(comp),'tag'); % keep tag during plotting
  ud=get(h(comp),'userdata'); % keep tag during plotting
  cData = double(pp); % plot double
  cData = cData'; % cData should be size length(Y)xlength(X)
  ff = ff';
  if flagLog
    cData = log10(cData);
  end
  if useTStartEpoch
    xData = double(tt-t_start_epoch)';
  else 
    xData = double(tt)';
  end
  if flagReducePlot
    axes(h(comp));
    hlpe = SpecPlotReducer(xData,ff,cData);
  else
    pcolor(h(comp),xData,ff,cData)
  end
  set(h(comp),'tag',tag);
  set(h(comp),'userdata',ud);
  zoom_in_if_necessary(h(comp)); %
  
  shading(h(comp),'flat')
  set(h(comp),'TickDir','out')
  %check ylabel
  if ~isfield(specrec,'f_label')
    if ~isfield(specrec,'f_unit')
      specrec.f_unit='a.u.';
    end
    specrec.f_label=['[' specrec.f_unit ']'];
  end
  ylabel(h(comp),specrec.f_label)
  
  if showColorbar
    if isfield(specrec,'p_label')
      if isa(h(comp),'handle'), hcb = colorbar(h(comp)); % HG2
      else, hcb = colorbar('peer',h(comp));
      end
      drawnow
      posCb = get(hcb,'Position');
      posAx = get(h(comp),'Position');
      drawnow
      set(hcb,'TickDir','out','Position',...
        [posCb(1) posCb(2)+posCb(4)*0.05 posCb(3)*.75 posCb(4)*0.9])
      set(h(comp),'Position',[posAx(1) posAx(2) (posCb(1)-posAx(1))*0.97 posAx(4)])
      ylabel(hcb,specrec.p_label);
      if fitColorbarLabel
        irf_colorbar_fit_label_height(hcb);
      end
    else
      hcb{comp}= [];
    end
  else
    hcb{comp}= [];
  end
  if comp==min(length(h),ncomp)
    irf_timeaxis;
  else
    set(h(comp),'XTicklabel','')
  end
end

if nargout>0, hout=h; end

function zoom_in_if_necessary(h)
  ud=get(h,'userdata');
  if isfield(ud,'zoom_x')
    disp('zooming in the updated plot')
    irf_zoom(h,'x',ud.zoom_x);
    if ud.zoom_x(1) > 1e8 && ud.zoom_x(1) < 1e10 % isdat epoch
      irf_timeaxis(h,'nolabel');
    end
end
