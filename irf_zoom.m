function irf_zoom(interval,c,axis_handles,t_ref)
%IRF_ZOOM   Zoom in
%
%  irf_zoom([min max],c,axis_handles,t_ref)
%  irf_zoom([min max],c,axis_handles)
%  irf_zoom({[yyyy mm dd hh mm ss] [yyyy mm dd hh mm ss]},c,axis_handles)
%  irf_zoom({[yyyy mm dd hh mm ss] [yyyy mm dd hh mm ss]},c,axis_handles,t_ref)
%  left side of the date vectors can be skipped, then one uses the 
%  values from axis
%  c='x' for x-axis, 'y' for y-axis
%  t_ref is isdat_epoch of time=0 point
%
% $Id$

flag_use_t_start_epoch=0; % if 1 use userdata.t_start_epoch as tref

if nargin==0, help irf_zoom, return; end
if nargin<4,
  % Try to read the reference time from figures user_data variable
  user_data=get(gcf,'userdata');
  if isfield(user_data,'t_start_epoch')
    t_ref=user_data.t_start_epoch;
    flag_use_t_start_epoch=1;
  else
    t_ref=0;
  end
end

if nargin < 3, axis_handles=gca;end
if nargin == 3, 
  ax = reshape(axis_handles,1,numel(axis_handles));
  clear axis_handles;
  axis_handles = ax;
end
 
if size(interval,2) ~= 2, disp('zooming interval in wrong format');return;end

if strcmp(c,'x'), 
  if iscell(interval),  % Simplified time zooming
    ax=axis;
    if ( ax(1)+t_ref>1e8 && ax(1)+t_ref<1e10 ),
      int_min=fromepoch(ax(1)+t_ref);
      int_max=fromepoch(ax(2)+t_ref);
      int_min(7-size(interval{1},2):6)=interval{1};
      int_max(7-size(interval{2},2):6)=interval{2};
      clear interval;
      interval=[toepoch(int_min) toepoch(int_max)]-t_ref;
    end
  else % Interval must be vector with two values
    if flag_use_t_start_epoch, % Account for reference time from userdata.t_start_epoch
      interval=interval-t_ref;
    end
  end
end

% Make interval finite if it has only one point
if diff(interval)==0, interval(2)=interval(1)+1e-10; end 

% Remove XTickLabel and XLabel from all panels but the last one
if c=='x' && length(axis_handles)>1
	p = cell2mat(get(axis_handles,'Position'));
	pymin = min(p(:,2));
end

for h=axis_handles
	if c=='x'
		set(h,'XLim',interval);
		if ax(1)+t_ref>1e8 && ax(1)+t_ref<1e10
			if flag_use_t_start_epoch % Read t_ref from userdata.t_start_epoch
				p = get(h,'position');
				if length(axis_handles)>1 && p(2)==pymin, add_timeaxis(h);
				else add_timeaxis(h,'nolabels');
				end
			else
				add_timeaxis(h,t_ref);
				if length(axis_handles)>1
					p = get(h,'position');
					if p(2)>pymin, xlabel(h,''), set(h,'XTickLabel',''), end
				end
			end
		end
	end
	if c=='y'
		set(h,'Ylim',interval);
	end
end
