function [hax,hl]=irf_plot_new(varargin)
% [hax,hl]=irf_plot_new(varargin)

%% Check input
[ax,args,nargs] = irf.axescheck(varargin{:});
x=args{1}; args=args(2:end); original_args=args;
if isempty(x), eS='nothing to plot'; irf.log('critical','eS'),error(eS),end

hl = [];
% Check if single number argument, then use syntax IRF_PLOT(number)
if isnumeric(x) && numel(x)==1, init_figure()
  if nargout>0, hax=cTmp; end, return
end

% Default values that can be overriden by options
dt = 0;
flag_yy = 0;
scaleyy = 1;
plot_type = '';
marker = '-';
flag_plot_all_data=1;
flag_colorbar=1;
check_input_options()
if ~iscell(x), x = {x}; end

switch plot_type
  case '', plot_in_single_panel()
  case 'comp', plot_components_in_separate_panels()
  case 'subplot', plot_variables_in_separate_panels()
  otherwise, error('should not be here')
end

return % return from main function

  function plot_variables_in_separate_panels()

  end

  function  plot_components_in_separate_panels()

  end

  function plot_in_single_panel()
    if isempty(ax) % if empty axis use current axis GCA
      if isempty(get(0,'CurrentFigure')) % there is no figure open
        ax = irf_panel(randStr);
      else, ax = gca;
      end
    end
    hca = ax(1);
    ts = irf_plot_start_epoch(x{1}.time);
    tag=get(hca,'tag'); ud=get(hca,'userdata'); % keep tag/userdata during plotting

    %if length(x)==2 && strcmp(get_units(x{1}),get_units(x{2}))
    flagHold = ishold(hca);
    if flagHold, flagHolding = true; else, flagHolding = false; end
    for iVar = 1:length(x)
      h = plot(hca, x{iVar}.time.epochUnix-ts-dt, x{iVar}.data, marker, args{:});
      hl=[hl; h];
      if iVar == 1 && ~flagHold && length(x)>1
        hold(hca,'on'), flagHolding = true;
      end
    end
    if ~flagHold && flagHolding, hold(hca,'off'), end

    grid(hca,'on');
    set(hca,'tag',tag); set(hca,'userdata',ud); % restore
    %zoom_in_if_necessary(hca);
    % Put YLimits so that no labels are at the end (disturbing in multipanel plots)
    if ~ishold(hca), irf_zoom(hca,'y'); end % automatic zoom only if hold is not on
    %ylabel(hca,get_label());
    irf_timeaxis(hca)
    hax=hca;

    function st = randStr
      symbols = ['a':'z' 'A':'Z' '0':'9'];
      MAX_ST_LENGTH = 10;
      stLength = randi(MAX_ST_LENGTH);
      nums = randi(numel(symbols),[1 stLength]);
      st = symbols (nums);
    end
  end

  function check_input_options()
    have_options = 0;
    if nargs > 1, have_options = 1; end
    while have_options
      l = 1;
      switch(lower(args{1}))
        case 'newfigure'
          hax=initialize_figure(x);
        case 'subplot'
          plot_type = 'subplot';
        case 'comp'
          plot_type = 'comp';
        case 'dt'
          if nargs>1
            if isnumeric(args{2})
              dt = args{2};
              l = 2;
            else
              irf.log('critical','wrongArgType : dt must be numeric')
              error('dt must be numeric');
            end
          else, irf.log('critical','wrongArgType : dt value is missing')
            error('dt value missing');
          end
        case 'tint'
          if nargs>1 && isnumeric(args{2})
            tint = args{2};
            l = 2;
            flag_plot_all_data=0;
            original_args(find(strcmpi(original_args,'tint'))+ [0 1])=[]; % remove tint argument
          else % TODO implement string tint
            irf.log('critical','wrongArgType : tint must be numeric')
            error('tint must be numeric')
          end
        case 'yy'
          if nargs>1
            if isnumeric(args{2})
              flag_yy = 1;
              scaleyy = args{2};
              l = 2;
            else
              irf.log('critical','wrongArgType : yy must be numeric')
              error('yy must be numeric');
            end
          else
            irf.log('critical','wrongArgType : yy value is missing')
            error('yy value missing');
          end
        case 'linestyle'
          marker = args{2};
          l = 2;
        case 'nocolorbar'
          flag_colorbar=0;
          l = 1;
        otherwise
          marker = args{1};
          args = args(2:end);
          break
      end
      args = args(l+1:end);
      if isempty(args), break, end
    end
  end
  function init_figure
    if x>=1 && x<=20
      % check if there is 'newfigure' argument
      if numel(args)>=2 && ischar(args{2}) && strcmpi(args{2},'newfigure')
        hax=initialize_figure(x,'newfigure');
      elseif numel(args)>=2 && ischar(args{2}) && strcmpi(args{2},'reset')
        hax=initialize_figure(x,'reset');
      else, hax=initialize_figure(x);
      end
    else, irf.log('warning','Max 20 subplots supported ;)');
    end
  end

end
