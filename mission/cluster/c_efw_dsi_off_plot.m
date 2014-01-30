%c_efw_dsi_off_plot

t=toepoch([2001 02 01 00 00 0]):86400:toepoch([2010 01 01 00 00 0]); 
t = t' + 600;



%% ms offsets
for cl_id = 2:4
  off=zeros(length(t),2);
  off(:,1) = t;
  for idx=1:length(t)
    off(idx,2) = c_efw_dsi_off(t(idx),cl_id,[t(idx) -30]);
  end
  c_eval('msoff?=off;',cl_id);
  clear off
end

c_pl_tx msoff?
  
  