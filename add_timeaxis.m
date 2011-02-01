function add_timeaxis( h, t_start_epoch, xlabels, xlabeltitle )
%ADD_TIMEAXIS  add time axis
%
% function add_timeaxis( h, t_start_epoch, xlabels, xlabeltitle );
% function add_timeaxis( h, t_start_epoch );
% function add_timeaxis( h, 'usefig' );      % to use t_start_epoch from the figure
% function add_timeaxis( h, 'date' );        % to add xlabel with date
% function add_timeaxis( h, 'nodate' );      % do not add xlabel with date
% function add_timeaxis( h, 'nolabels' );    % do not add any labels (only ticks)
%
% Adds time axis in hh:mm:ss format to axis given by handles h.
% If t_start_epoch is given, then adds so many seconds to the time axis.
% For example if time on axis is in seconds from the beginning of 0300 UT
% 01-Jan-2001, then use add_timeaxis(h, toepoch([2001 01 01 03 00 00])).
%
% if figures user_data.t_start_epoch is defined use that as t_start_epoch
%
% If xlabels are defined then adds x-tra labels in addition to time.
% xlabels format is column vector [time lab1 lab2 ...] where lab1 are
% numerical values and time is in isdat_epoch. Program then interpolates
% to the time of labels xlabeltitle = {'LAB1' 'LAB2' ..}; is the str
% for labels.
%
% $Id$

flag_labels=1; % default is to add labels to the last axis handle, can be changed by 'nolabels' argument
flag_date=1;   % default is to add date labels
  if nargin == 0
     h = gca;
  end

  hh = reshape( h, 1, numel(h) );
  clear h;
  h = hh;

  flag_usefig = 0;

  if (nargin >= 2) && (ischar(t_start_epoch))
     if strcmp(t_start_epoch,'date')
        flag_date = 1;
     elseif strcmp(t_start_epoch,'nodate')
        flag_date = 0;
     elseif strcmp(t_start_epoch,'nolabels')
        flag_labels = 0;
        flag_date = 0;
	 elseif strcmp(t_start_epoch,'usefig')
        flag_usefig = 1;
     end
  end

  if ~exist('t_start_epoch','var') || ischar(t_start_epoch) || flag_usefig
    user_data = get(gcf,'userdata');
    if isfield(user_data,'t_start_epoch')
      t_start_epoch = double(user_data.t_start_epoch);
    else
      t_start_epoch = double(0);
    end
  end

  for j=1:length(h)
      xlabel(h(j),'');
      ax = axis;axis(axis);
      tint = ax(1:2) + t_start_epoch;
      res  = timeaxis(tint);
      set( h(j), 'XTick', res{1} - t_start_epoch );
      if j == length(h)
         set( h(j), 'XTickLabel', res{2} );
      else
         set( h(j), 'XTickLabel','');
      end

      if nargin > 2  % xlabels should be added
          set( h(j), 'XTickLabel','');
          lab    = res{2};
          xcoord = res{1};
          for ii = 1:size(res{1},2)
              if ~strcmp(lab(ii),' ')
                  ax = axis;
                  mm = irf_resamp( xlabels, xcoord(ii));
                  for jj = 1:length(mm)
                      if jj==1, % the first line is time
                          str = lab(ii);
                      else % other lines are xlabels
                          str = [repmat(' \newline',1,jj-1) num2str(mm(jj),3)];
                      end
                      outhandle   = text(xcoord(ii)-t_start_epoch, ...
					  	ax(3)-abs(ax(3)-ax(4))/100, str);
                      set( outhandle, 'HorizontalAlignment', 'center', ...
                          'VerticalAlignment', 'top', 'FontSize', 10);
                  end
              end
          end

          % Add titles
          str = 'UT';
          for jj = 0:size(xlabeltitle,2),
              if jj>0,
                  flag_date=0; % if more than one line in xlabels, remove date
                  str      = [repmat(' \newline',1,jj) xlabeltitle{jj}];
              end
              outhandle = text( ax(1)-abs(ax(2)-ax(1))/20, ...
			  	ax(3)-abs(ax(3)-ax(4))/100, str );
              set( outhandle, 'HorizontalAlignment', 'right', ...
                  'VerticalAlignment', 'top', 'FontSize', 10);
          end
      end
  end

  start_time = fromepoch( ax(1) + t_start_epoch );
  time_label = datestr( datenum(start_time),1 );
  if flag_date == 1, xlabel(time_label);  end
  if flag_labels == 0, set(gca,'XTickLabel',' '); end
  return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
