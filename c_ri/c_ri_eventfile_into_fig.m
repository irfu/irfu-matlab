function h = c_ri_eventfile_into_fig(time_interval,path_events,panels)
%function h = c_ri_eventfile_into_fig(time_interval,path_events,panels)
%
%Input:
% time_interval - isdat_epoch [start_time end_time]
% path_events - path where event files are located, ex './'
% panels - structure with list of panels to plot, ex. {'Bx','By','B','B1','Ex','Vps'}
%
%Output:
%  h - handle to figures

global AV_DEBUG; if isempty(AV_DEBUG), debug=0;else, debug=AV_DEBUG;end
n_panels=size(panels,1);
i_fig=1;
plot_command=struct(...
  'Bx','c_pl_tx(B1,B2,B3,B4,2);ylabel(''B_X [nT] GSE'');', ...
  'By','c_pl_tx(B1,B2,B3,B4,3);ylabel(''B_Y [nT] GSE'');', ...
  'Bz','c_pl_tx(B1,B2,B3,B4,4);ylabel(''B_Z [nT] GSE'');', ...
  'B' ,'c_pl_tx(av_abs(B1),av_abs(B2),av_abs(B3),av_abs(B4),2);ylabel(''B [nT] GSE'');', ...
  'B1' ,'av_tplot(av_abs(B1));ylabel(''B [nT] GSE, sc1'');', ...
  'test','test' ...
)

file_list=dir([path_events '*F*t*T*t*.mat']);
for i_file=1:size(file_list,1),
  if c_ri_timestr_within_tint(file_list(i_file).name,time_interval),
     load([path_events file_list(i_file)]);
     figure(i_fig);i_panel=1;
     for i_panel=1:size(panels,1),
        h(i_fig,i_panel)=av_subplot(n_panels,1,-i_panel);i_panel=i_panel+1;
        eval(['plot_command.' panels{i_panel}]);
     end
     add_timeaxis(h(i_fig));
     i_fig=i_fig+1;
  end
end


