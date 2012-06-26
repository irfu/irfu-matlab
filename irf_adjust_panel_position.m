


% ================================================================
% used to adjust the panel postion, which is produced by irf_plot. 
% It should be used behind all subplots, and before irf_zoom
% see also irf_zoom
% 
% created by Huishan  (huishanf@gmail.com)
% ================================================================



for ii=1:size(h,2)
    pos_panel = get(h(ii), 'Position');
    panel_dx(ii) = pos_panel (3);
end

back_space=[0 0 max(panel_dx)-min(panel_dx) 0];

%Adjust the position
for ii=1:size(h,2)
    pos = get(h(ii), 'Position');
    if pos(3)>(min(panel_dx)+0.2*back_space(3))
        set(h(ii), 'Position',pos-back_space);
    end
end


