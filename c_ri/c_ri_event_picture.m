function a = c_ri_event_picture(time_of_events,dt,angles,ampl,path_out)
%function a = c_ri_event_picture(time_of_events,dt,angles,ampl,path_out)
%
% fast solution to plot B overview for the events 
% time_of_events,ampl are obtained from c_ri_angles_and_ampl

per=dt;
f_events=time_of_events;
f_count=size(f_events,1);
%-------------------------------------
%writes an ascii file for the reduced number of points
p_m = [fromepoch(f_events(:,1)) f_events(:,2:5)];
p_and_f = sprintf('%sEvent_lists.txt',path_out);
disp(['saved: ' p_and_f]);
fp = fopen(p_and_f,'w+');

pr_r = sprintf('Time                 angle ampl1 ampl2 mode\n');
fwrite(fp,pr_r);

for k = 1:f_count
  dstr = R_datestring2(p_m(k,1:6));
  pr_r = sprintf('%s %5.0f %4.0f %4.0f    %s\n',dstr,p_m(k,7:10));
  fwrite(fp,pr_r);
end

fclose(fp);

for g =1:f_count
  t = f_events(g,1);
  s_t = fromepoch(t);
  e_t = fromepoch(t); 
  
  [B1,B2,B3,B4]=c_get_bfgm(t+ [-per per]);whos B1 B2 B3 B4;
  c_eval('B?=av_abs(B?);');
  fg = figure;
  %plots B1
  h(1)=av_subplot(6,1,-1);
  av_tplot(B1);hold on;
  set(gca,'xlim',t+[-per per]);
  plot(t,0,'xk')
  ylabel('B1, nT')
  add_timeaxis
  legend('Bx','By','Bz')
  hold off
  
  h(2)=av_subplot(6,1,-2);
  c_pl_tx(B1,B2,B3,B4,5);hold on;
  set(gca,'xlim',t+[-per per]);
  plot(t,0,'xk')
  ylabel('|B|, nT')
  add_timeaxis
  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
  
  %plots Bx for all clusters
  h(3)=av_subplot(6,1,-3);
  c_pl_tx(B1,B2,B3,B4,2);hold on;
  set(gca,'xlim',t+[-per per]);
  plot(t,0,'xk')
  ylabel('Bx, nT')
  add_timeaxis
  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
  hold off
  
  %plots By for all clusters
  subplot(6,1,4)
  hold on
  c_pl_tx(B1,B2,B3,B4,3);
  set(gca,'xlim',t+[-per per]);
  plot(t,0,'xk')
  ylabel('By, nT')
  add_timeaxis
  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
  hold off
  
  %plots Bz for all clusters
  subplot(6,1,5)
  hold on
  c_pl_tx(B1,B2,B3,B4,4);
  set(gca,'xlim',t+[-per per]);
  plot(t,0,'xk')
  ylabel('Bz, nT')
  add_timeaxis
  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
  hold off
  
  
  %plotting data
  load '.c_ri_parameters.mat'
  
  
  subplot(6,1,6)
  
  %right side 
  text(2.5,3,['angles: ' int2str(angles(g,2:7)) ' degrees']);
  text(2.5,4.5,['|B|: ' int2str(ampl(g,:)) ' nT']);
  axis([1,4,1,6])
  axis off
  
  p_and_f_picture = sprintf('%sF_%s',p_up,c_ri_datestring_file(fromepoch(t)));
  print( fg, '-djpeg', p_and_f_picture);
  disp(['saving: ' p_and_f_picture]);
  close(fg)
end

