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
  disp([num2str(g) '.event. ' datestr(fromepoch(t),31)]);
  
  [B1,B2,B3,B4]=c_get_bfgm(t+ [-per per]);
  c_eval('try B?=av_abs(B?);catch B?=[NaN NaN NaN NaN NaN];end');
  %c_eval('[ts,te,tm?]=createEFWModeTableFDM(''disco:10'',s_t,10,?,''tm'');');
  %title_str=['sc mode (0-normal, 1-burst): ' num2str([tm1 tm2 tm3 tm4])];
  title_str=[' angles: ' int2str(angles(g,2:7)) ' degrees'];
  title_str=[title_str '|B|: ' int2str(ampl(g,:)) ' nT'];

  fg = figure('visible','off');
  set(gcf,'Units','centimeters')
  set(gcf,'Position',[1 1 15 20])
  %plots B1
  h(1)=av_subplot(6,1,-1);
  av_tplot(B1);
  ylabel('B1, nT')
  title(title_str);
  %  legend('Bx','By','Bz')
  
  h(2)=av_subplot(6,1,-2);
  c_pl_tx(B1,B2,B3,B4,5);hold on;
  set(gca,'xlim',t+[-per per]);
  ylabel('|B|, nT')
%  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')

  %plots Bx for all clusters
  h(3)=av_subplot(6,1,-3);
  c_pl_tx(B1,B2,B3,B4,2);
  ylabel('Bx, nT')
%  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
  
  %plots Bx for all clusters
  h(4)=av_subplot(6,1,-4);
  c_pl_tx(B1,B2,B3,B4,3);
  ylabel('By, nT')
%  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
  
  %plots Bx for all clusters
  h(5)=av_subplot(6,1,-5);
  c_pl_tx(B1,B2,B3,B4,4);
  ylabel('Bz, nT')
%  legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')

  av_pl_mark(t + [-.1 .1],h,'y');
  av_zoom(t+[-per per],'x',h);
  add_timeaxis(h);
%  legend
  
  %plotting data
  load '.c_ri_parameters.mat'
  
  
  subplot(6,1,6)
  %right side 
 
  axis off;  

  p_and_f_picture = sprintf('%sF_%s',path_out,c_ri_datestring_file(fromepoch(t)));
  eps2png = '/usr/local/bin/eps2png';

  print('-depsc2',[p_and_f_picture '.eps'])
  unix(sprintf('(%s %s.eps)',eps2png,p_and_f_picture));

  disp(['saving: ' p_and_f_picture '.png']);
  close(fg)
end

