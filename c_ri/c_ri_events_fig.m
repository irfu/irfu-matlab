function c_ri_events_fig(time_interval,p_data,p_MP,p_Bp,p_E,p_Out, period, fig_type)
% function c_ri_events_fig(time_interval,p_data,p_MP,p_Bp,p_E,p_Out, period, fig_type)
%
%Input:
%
%Output:
%
%Descrition of the function:
%
%Using:
% 


event_files=dir([p_E 'E_' '*.mat']);

for i_event_file=1:size(event_files,1),
  file_name = event_files(i_event_file).name;
  if c_ri_timestr_within_intervall_E(file_name,fromepoch(time_interval(1)),fromepoch(time_interval(2))) == 1,
     c_ri_events_into_fig(p_E,file_name,period,p_Out,p_Bp,p_MP,p_data);
  end
end

