function add_timeaxis( h, t_start_epoch, xlabels, xlabeltitle );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function add_timeaxis( h, t_start_epoch, xlabels, xlabeltitle );
% function add_timeaxis( h, t_start_epoch );
% function add_timeaxis( h, 'date' );        % to add xlabel with date
% function add_timeaxis( h, 'nodate' );      % do not add xlabel with date
% function add_timeaxis( h, 'nolabels' );    % do not add any labels (only ticks)
%
% Adds time axis in hh:mm:ss format to axis given by handles h.
% If t_start_epoch is given, then adds so many seconds to the time axis.
% For example if time on axis is in seconds from the beginning of 0300 UT
% 01-Jan-2001, then use add_timeaxis(h, toepoch([2001 01 01 03 00 00])).
%
% If xlabels are defined then adds x-tra labels in addition to time.
% xlabels format is column vector [time lab1 lab2 ...] where lab1 are
% numerical values and time is in isdat_epoch. Program then interpolates
% to the time of labels xlabeltitle = {'LAB1' 'LAB2' ..}; is the string
% for labels.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_labels=1; % default is to add labels to the last axis handle, can be changed by 'nolabels' argument
flag_date=1;   % default is to add date labels
  if nargin == 0
     h = gca;
  end

  hh = reshape( h, 1, prod(size(h)) );
  clear h;
  h = hh;


%  if size(h,2) == 1
%     flag_date = 0;
%  else
%     flag_date = 1;
%  end

  if nargin < 2
     t_start_epoch=0;
  end

  if (nargin >= 2) & (ischar(t_start_epoch))
     if strcmp(t_start_epoch,'date')
        flag_date = 1;
     elseif strcmp(t_start_epoch,'nodate')
        flag_date = 0;
     elseif strcmp(t_start_epoch,'nolabels')
        flag_labels = 0;
     end
     t_start_epoch = 0;
  end

  for j=1:length(h)
      axes( h(j) );xlabel('');
      ax   = axis;
      tint = ax(1:2) + t_start_epoch;
      res  = timeaxis(tint);
      set( h(j), 'XTick', res{1} - t_start_epoch );
      if j == length(h)
         set( h(j), 'XTickLabel', res{2} );
      else
         set( h(j), 'XTickLabel','');
      end

      if nargin > 2
         axes( h(j) );
         ax   = axis;
         tint = ax(1:2) + t_start_epoch;
         res  = timeaxis(tint);
         set( h(j), 'XTick', res{1} - t_start_epoch );
         if j == length(h)
            set( h(j), 'XTickLabel', res{2} );
         else
            set( h(j), 'XTickLabel','');
         end

         lab    = res{2};
         xcoord = res{1};
         for ii = 1:size(res{1},2)
             if ~strcmp(lab(ii),' ')
                ax = axis;
                mm = av_interp( xlabels, xcoord(ii) + t_start_epoch );
                string_init = ' \newline';
                for jj = 2:length(mm)
                    string      = [string_init num2str(mm(jj),3)];
                    string_init = [string_init ' \newline'];
                    outhandle   = text( xcoord(ii), ax(3), string );
                    set( outhandle, 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'top', 'FontSize', 10);
                end
             end
         end

         % Add titles
         string_init = ' \newline';
         outhandle = text( ax(1), ax(3), 'UT' );
         set( outhandle, 'HorizontalAlignment', 'right', ...
              'VerticalAlignment', 'top', 'FontSize', 10);
         for jj = 1:size(xlabeltitle,2)
             string      = [string_init xlabeltitle{jj}];
             string_init = [string_init ' \newline'];
             outhandle = text( ax(1), ax(3), string );
             set( outhandle, 'HorizontalAlignment', 'right', ...
                  'VerticalAlignment', 'top', 'FontSize', 10);
         end
     end
  end

  start_time = fromepoch( ax(1) + t_start_epoch );
  time_label = datestr( datenum(start_time),1 );
  %disp(time_label);
  if flag_date == 1
     xlabel(time_label);
  end
  if flag_labels == 0,
    set(gca,'xticklabel',' ');
  end
  return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
