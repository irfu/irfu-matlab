%c_efw_dsi_off_plot

t=toepoch([2001 02 01 00 00 0]):86400:toepoch([2010 01 01 00 00 0]); 
t = t' + 600;



%% ms offsets
for cl_id = 1:4
  off=zeros(length(t),2);
  off(:,1) = t;
  for idx=1:length(t)
    off(idx,2) = c_efw_dsi_off(t(idx),cl_id,[t(idx) -30]);
  end
  c_eval('msoff?=off;',cl_id);
  clear off
end

figure
c_pl_tx msoff?
ylabel('Ex offset [mV/m]')
xlabel('Time [year]')
legend('C1','C2','C3','C4')
print -dpdf ms_offsets

%% sh offsets
for cl_id = 1:4
  off=zeros(length(t),2);
  off(:,1) = t;
  for idx=1:length(t)
    off(idx,2) = c_efw_dsi_off(t(idx),cl_id);
  end
  c_eval('shoff?=off;',cl_id);
  clear off
end
figure
c_pl_tx shoff?
ylabel('Ex offset [mV/m]')
xlabel('Time [year]')
legend('C1','C2','C3','C4')
print -dpdf sh_offsets
  
  