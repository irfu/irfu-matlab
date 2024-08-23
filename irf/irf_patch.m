function hp_out = irf_patch(varargin)
% irf_patch - plots timeseries patches
%
%   tint = irf.tint('2000-01-01T00:00:00.00Z',60);
%   time_vector = (tint.start:1:tint.stop)-tint.start;
%   tsScalar0 = irf.ts_scalar(tint.start:1:tint.stop,cos(2*pi*time_vector/60));
%   tsScalar1 = irf.ts_scalar(tint.start:1:tint.stop,cos(2*pi*time_vector/60)+0.3*abs(randn(size(time_vector))));
%   tsScalar2 = irf.ts_scalar(tint.start:1:tint.stop,cos(2*pi*time_vector/60)-0.3*abs(randn(size(time_vector))));
%   tsScalar3 = irf.ts_scalar(tint.start:1:tint.stop,0.5*cos(2*pi*time_vector/60)+0.1*abs(randn(size(time_vector))));
%   randMat = abs(rand(length(time_vector),100)); edgesRand = 0:0.1:1;
%   randDist = histc(randMat,edgesRand,2);
%   cumSumRandDist = cumsum(randDist(:,1:10),2);
%   tsScalarMat = irf.ts_scalar(tint.start:1:tint.stop,cumSumRandDist);
%
%   % Example 1:
%     hp = irf_patch({tsScalar2,tsScalar1}); hold on;
%     hl = irf_plot(tsScalar0); hold off;
%
%   % Example 2:
%     divider = 0.1;
%     hp = irf_patch({tsScalar0,divider});
%
%   % Example 3:
%     divider = 0.1;
%     hp = irf_patch({tsScalar0,divider},'larger'); hold on; % !NB: bugs around the overlap points, use with care!
%     hp = irf_patch({tsScalar0,divider},'smaller'); hold on;
%     hl = irf_plot(tsScalar0); hold off; hl.Color = [0 0 0];
%
%   % Example 4:
%     hp = irf_patch({tsScalar1,tsScalar0}); hold on;
%     hp = irf_patch({tsScalar2,tsScalar0}); hold on;
%     hl = irf_plot(tsScalar0); hold off; hl.Color = [0 0 0];
%
%   % Example 5:
%     hp = irf_patch({tsScalar1,tsScalar3},'larger'); hold on;
%     hp = irf_patch({tsScalar1,tsScalar3},'smaller'); hold on;
%     hl = irf_plot(tsScalar1); hl.Color = [0 0 0];
%     hl = irf_plot(tsScalar3); hl.Color = [0 0 0]; hold off;
%
%   Example 6:
%     irf_patch(tsScalarMat);
%     labels = arrayfun(@(x,y) {[num2str(x) ' > n > ' num2str(y)]}, edgesRand(1:1:end-1),edgesRand(2:1:end));
%     hpatches=findall(gca,'Type','patch');
%     legend(hpatches(end:-1:1),labels,'location','eastoutside')

%

[ax,args,nargs] = axescheck(varargin{:});

if isempty(ax)
  ax = gca;
end

axis(ax);

plotOption = 'difference';
if numel(args) == 2 && ischar(args{2})
  plotOption = args{2};
end

if numel(args) == 3 && ischar(args{3})
  doStack = args{3};
end


if iscell(args{1})
  patch1 = args{1}{1};
  if isa(args{1}{2},'TSeries')
    patch2 = args{1}{2};
  elseif isnumeric(args{1}{2}) && size(args{1}{2},1)==1
    if size(args{1}{2},2) == 1
      patch2 = patch1.clone(patch1.time,repmat(args{1}{2},patch1.length,size(patch1.data,2)));
    elseif size(args{1}{2},2) == size(patch1.data,2)
      patch2 = patch1.clone(patch1.time,repmat(args{1}{2},patch1.length,1));
    end
  end

  switch plotOption
    case 'difference'
    case 'larger'
      indT = find(patch1.data<patch2.data);
      patch2.data(indT,:) = patch1.data(indT,:);
    case 'smaller'
      indT = find(patch1.data>patch2.data);
      patch2.data(indT) = patch1.data(indT,:);
  end
  %dtStart =
  %ts = ;
  tPatch1 = patch1.time.epochUnix-t_start_epoch(patch1.time.epochUnix);
  tPatch2 = patch2.time.epochUnix-t_start_epoch(patch2.time.epochUnix);

  allPatch = irf.ts_scalar(patch1.time,[patch1.data patch2.data]);
  hl = irf_plot(ax,allPatch); for ii = 1:numel(hl); hl(ii).Visible = 'off'; end; hold(ax,'off')

  %hl = irf_plot(ax,patch1); for ii = 1:numel(hl); hl(ii).Visible = 'off'; end; hold(ax,'on')
  %hl = irf_plot(ax,patch2); for ii = 1:numel(hl); hl(ii).Visible = 'off'; end; hold(ax,'off')

  for iComp = 1:size(patch1.data,2)
    xPatch = [tPatch1; NaN; tPatch2(end:-1:1)];
    yPatch = [patch1.data(:,iComp); NaN; patch2.data(end:-1:1,iComp)];

    isNan = isnan(yPatch);
    xPatch(isNan) = [];
    yPatch(isNan) = [];

    hp = patch(xPatch,yPatch,'k','Parent', ax);
    hp.FaceAlpha = 0.2;
    hp.EdgeColor = hl(iComp).Color;
    hp.FaceColor = hl(iComp).Color;
    hp_out(iComp) = hp;
  end

elseif isa(args{1},'TSeries')
  tsData = args{1};
  tmpData = tsData.data;
  tsData = irf.ts_scalar(tsData.time,[zeros(tsData.length,1) tmpData]); % add zero level
  nDim2 = size(tsData.data,2);
  cmap = colormap('jet'); cmap = flipdim(cmap(fix(linspace(1,64,nDim2-1)),:),1);
  tPatch = tsData.time.epochUnix-t_start_epoch(tsData.time.epochUnix);

  hold(ax,'on')
  colors = get(ax,'ColorOrder');
  % plot white line of max value so that irf_zoom works
  maxData = max(tsData.data,[],2);
  aa = irf_plot(ax,irf.ts_scalar(tsData.time,maxData*1.1),'w'); aa.Visible = 'off';
  aa = irf_plot(ax,irf.ts_scalar(tsData.time,maxData*0),'w');  aa.Visible = 'off';
  % plot patches
  for iP = 1:nDim2-1
    xPatch = [tPatch; NaN; tPatch(end:-1:1)];
    yPatch = [tsData.data(:,iP); NaN; tsData.data(end:-1:1,iP+1)];


    isNan = isnan(yPatch);
    xPatch(isNan) = [];
    yPatch(isNan) = [];

    hp = patch(xPatch,yPatch,'k','Parent', ax);
    hp.FaceAlpha = 0.9;
    hp.EdgeColor = 'none';
    hp.FaceColor = cmap(iP,:);

    if all(tsData.data(:,iP) == tsData.data(:,iP+1))
      hp.Visible = 'off';
    end

    %pause
    hp_out(iP) = hp;
  end

  hold(ax,'off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function t_st_e = t_start_epoch(t)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gives back the value of t_start_epoch of the figure
    % if not  set, sets t_start_epoch of the figure
    ud = get(gcf,'userdata');
    ii = find(~isnan(t));
    if ~isempty(ii), valid_time_stamp = t(ii(1)); else, valid_time_stamp = []; end

    if isfield(ud,'t_start_epoch')
      t_st_e = double(ud.t_start_epoch);
    elseif ~isempty(valid_time_stamp)
      if valid_time_stamp > 1e8
        % Set start_epoch if time is in isdat epoch
        % Warn about changing t_start_epoch
        t_st_e = double(valid_time_stamp);
        ud.t_start_epoch = t_st_e;
        set(gcf,'userdata',ud);
        irf.log('notice',['user_data.t_start_epoch is set to ' ...
          epoch2iso(t_st_e,1)]);
      else
        t_st_e = double(0);
      end
    else
      t_st_e = double(0);
    end

  end
end
