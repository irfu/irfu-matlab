function irf_timeaxis( h, t_start_epoch, xlabels, xlabeltitle )
%IRF_TIMEAXIS  add time axis
%
% function irf_timeaxis( h, t_start_epoch, xlabels, xlabeltitle );
% function irf_timeaxis( h, xlabels, xlabeltitle );
% function irf_timeaxis( h, t_start_epoch );
% function irf_timeaxis( h, 'usefig' );      % to use t_start_epoch from the figure
% function irf_timeaxis( h, 'date' );        % to add xlabel with date
% function irf_timeaxis( h, 'nodate' );      % do not add xlabel with date
% function irf_timeaxis( h, 'nolabels' );    % do not add any labels (only ticks)
%
% Adds time axis in hh:mm:ss format to axis given by handles h.
% If t_start_epoch is given, then adds so many seconds to the time axis.
% For example if time on axis is in seconds from the beginning of 0300 UT
% 01-Jan-2001, then use irf_timeaxis(h, toepoch([2001 01 01 03 00 00])).
%
% if figures user_data.t_start_epoch is defined use that as t_start_epoch
%
% If xlabels are defined then adds x-tra labels in addition to time.
% xlabels format is column vector [time lab1 lab2 ...] where lab1 are
% numerical values and time is in isdat_epoch. Program then interpolates
% to the time of labels xlabeltitle = {'LAB1' 'LAB2' ..}; is the str
% for labels.

addLabelsToLastHandle = true;  % can be changed by 'nolabels' argument
addDateLabel          = true;  % default is to add date label to last handle
addXtraXlabels        = false; % default add only time labels on x axis
useFiguresTEpochStart = false;

if nargin == 0
  h = gca;
end
if nargin > 2
  addXtraXlabels = true;
  if nargin == 3
    xlabeltitle = xlabels; % TODO bad style, redo all input with varargin
    xlabels = t_start_epoch;
    t_start_epoch =[];
    if isa(xlabels,'TSeries')
      xlabels = [xlabels.time.epochUnix xlabels.data];
    end
    useFiguresTEpochStart = true;
  end
  for j=1:numel(h)
    ud=get(h,'UserData');
    if ~isstruct(ud), clear ud; ud=struct; end
    ud.xlabels=xlabels;
    ud.xlabeltitle=xlabeltitle;
    set(h(j),'UserData',ud);
  end
end

hh = reshape( h, 1, numel(h) );
clear h;
h = hh;

if (nargin >= 2) && (ischar(t_start_epoch))
  if strcmp(t_start_epoch,'date')
    addDateLabel = true;
  elseif strcmp(t_start_epoch,'nodate')
    addDateLabel = false;
  elseif strcmp(t_start_epoch,'nolabels')
    addLabelsToLastHandle = false;
    addDateLabel = false;
    remove_extra_xlabel_handles(h);
    for j=1:numel(h) % clean extra xlabels if present
      ud=get(h(j),'UserData');
      if isfield(ud,'xlabels'), ud=rmfield(ud,'xlabels');end
      if isfield(ud,'xlabeltitle'), ud=rmfield(ud,'xlabeltitle');end
      set(h(j),'UserData',ud);
    end
  elseif strcmp(t_start_epoch,'usefig')
    useFiguresTEpochStart = 1;
  end
end

if ~exist('t_start_epoch','var') || ischar(t_start_epoch) || useFiguresTEpochStart
  user_data = get(gcf,'userdata');
  if isfield(user_data,'t_start_epoch')
    t_start_epoch = double(user_data.t_start_epoch);
  else
    t_start_epoch = double(0);
  end
end

for j=1:numel(h)
  %    xlabel(h(j),'');
  if addDateLabel==0, xlabel(h(j),''); end
  tint = get(h(j),'xlim') + t_start_epoch;
  res  = timeaxis(tint);
  set( h(j), 'XTick', res{1} - t_start_epoch );
  if j == numel(h)
    if ~addXtraXlabels
      ud=get(h(j),'UserData');
      if isfield(ud,'xlabels') % add extra labels
        addXtraXlabels = true;
        xlabels=ud.xlabels;
        xlabeltitle=ud.xlabeltitle;
      else                                 % add only time labels
        set( h(j), 'XTickLabel', res{2} );
      end
    end
  else
    set( h(j), 'XTickLabel','');
  end

  if addXtraXlabels  % xlabels should be added
    set( h(j), 'XTickLabel','');
    lab    = res{2};
    xcoord = res{1};
    h_xlabels=gobjects(size(res{1},2),size(xlabels,2)); % allocated handles for text labels
    remove_extra_xlabel_handles(h(j))
    for ii = 1:size(res{1},2)
      if ~strcmp(lab(ii),' ')
        ax = double(axis(h(j)));
        if isa(xlabels,'TSeries')
          xlabels = [xlabels.time.epochUnix xlabels.data];
        end
        mm = irf_resamp( xlabels, xcoord(ii));
        for jj = 1:length(mm)
          if jj==1 % the first line is time
            str = lab(ii);
          else % other lines are xlabels
            str = [repmat(' \newline',1,jj-1) num2str(mm(jj),3)];
          end
          h_xlabels(ii,jj) = text(xcoord(ii)-t_start_epoch, ax(3), str);
          set(h_xlabels(ii,jj), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top','parent',h(j),'units','normalized');
        end
      end
    end
    ud=get(h(j),'UserData');
    ud.h_xlabels=h_xlabels;
    set(h(j),'UserData',ud);

    % Add titles
    h_xlabeltitle=1:size(xlabeltitle,2); % allocated handles for text labels
    str = 'UT      ';
    for jj = 0:size(xlabeltitle,2)
      if jj>0
        addDateLabel=0; % if more than one line in xlabels, remove date
        str      = [repmat(' \newline',1,jj) xlabeltitle{jj} '     '];
      end
      h_xlabeltitle(jj+1) = text( ax(1),ax(3), str, 'parent',h(j));
      set(h_xlabeltitle(jj+1), 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'parent',h(j),'units','normalized');
    end
    ud=get(h(j),'UserData');
    ud.h_xlabeltitle=h_xlabeltitle;
    set(h(j),'UserData',ud);
  end
end

% determine the time stamp corrsponding to first labeled XTick, so that
% xlabel will have the date corresponding to that tick.
xTickLast=get(h(end),'XTick');  xTickLabelLast=get(h(end),'XTickLabel');
tFirstTick = [];
if ~isempty(xTickLabelLast)
  for i=1:length(xTickLast)
    if strcmp(xTickLabelLast{i},' '), continue, end
    tFirstTick = xTickLast(i); break
  end
end
xlimlast=get(h(end),'XLim');
if isempty(tFirstTick), tFirstTick = xlimlast(1); end
time_label = irf_time(tFirstTick + t_start_epoch,'epoch>utc_yyyy-mm-dd UTC');
if addDateLabel == 1 && addXtraXlabels ~= 1 && diff(xlimlast)<=3600*24*100
  xlabel(h(end),time_label);  % add data only if no extra xlabels
end
if addLabelsToLastHandle == 0
  set(h(end),'XTickLabel',' ');
end
return;

function remove_extra_xlabel_handles(h)
for j=1:numel(h) % clean extra xlabels if present
  hca=h(j);
  ud=get(hca,'UserData');
  if isempty(ud)
    hca.XLabel.String = '';
    hca.XTickLabel = '';
  else
    if isfield(ud,'h_xlabels') % remove old handles
      for jj=1:numel(ud.h_xlabels)
        if ishandle(ud.h_xlabels(jj)) && ud.h_xlabels(jj)~=0, delete(ud.h_xlabels(jj));end
      end
      ud=rmfield(ud,'h_xlabels');
    end
    if isfield(ud,'h_xlabeltitle') % remove old handles
      for jj=1:numel(ud.h_xlabeltitle)
        if ishandle(ud.h_xlabeltitle(jj)) && ud.h_xlabeltitle(jj)~=0, delete(ud.h_xlabeltitle(jj));end
      end
      ud=rmfield(ud,'h_xlabeltitle');
    end
    set(hca,'UserData',ud);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
