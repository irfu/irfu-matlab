function t_st_e = irf_plot_start_epoch(t)
%IRF_PLOT_START_EPOCH  get/set t_start_epoch of the figure
%
% t_st_e = irf_plot_start_epoch(t)
%
% Gives back the value of t_start_epoch (epochUnix) of the figure
% if not  set, sets t_start_epoch of the figure

ud = get(gcf,'userdata');
if isfield(ud,'t_start_epoch')
  t_st_e = double(ud.t_start_epoch);
  return;
end

if isa(t,'GenericTimeArray'), epochUnix = t.epochUnix; 
elseif isa(t,'int64'), epochUnix = EpochUnix(t).epochUnix;
elseif isa(t,'double'), epochUnix = t;
else
  errS = 'Unrecognized type of T'; irf.log('critical',errS), error(errS)
end

ii = find(~isnan(epochUnix));
if ~isempty(ii), valid_time_stamp = epochUnix(ii(1)); 
else, valid_time_stamp = EpochUnix('2010-01-01T00:00:00Z').epochUnix;
end

t_st_e = double(valid_time_stamp);
ud.t_start_epoch = t_st_e;
set(gcf,'userdata',ud);
irf.log('notice',['user_data.t_start_epoch is set to ' ...
  EpochUnix(t_st_e).utc]);
  