function a = c_ri_events_into_pictures(path_eve, file_events, per, p_up, path_Bp,p_MP)
%
%Input:
% path_eve, file_events - path and name of file to be loaded
%                    this file should contain:
%                    [time in epoch | angle | amplitude | amplitude | mode]
% per -the length of the filter if two events are within per, then they are classed as one event
% p_up -path to output
%
%Output:
%
%Descrition of the function:
%
%Using:
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
p_and_f =sprintf('%s%s',path_eve,file_events);
load(p_and_f)

nr_events = length(time_of_events(:,1));

%builds the data points into events
dt_ev = diff(time_of_events(:,1));

[i_end,c] = size(dt_ev);

%-----------------------------------
% this section reduces the number of events
%first point is always an event
f_events = time_of_events(1,:);
f_count = 1;
for i = 1:i_end

  %new event
  if dt_ev(i) > per
    f_count = f_count +1;
    f_events(f_count,:) = time_of_events(i+1,:);

    % if there is a new point with greater angle
  elseif f_events(f_count,2) < time_of_events(i+1,2)
    f_events(f_count,:) = time_of_events(i+1,:);
  end
end
a=f_events;
%-------------------------------------
%writes an ascii file for the reduced number of points
p_m = [fromepoch(f_events(:,1)) f_events(:,2:5)];
p_and_f = sprintf('%sR%st.txt',p_up,file_events(2:length(file_events)-4));
disp(['saved: ' p_and_f]);
fp = fopen(p_and_f,'w+');

pr_r = sprintf('Time                 angle ampl1 ampl2 mode\n');
fwrite(fp,pr_r);

for k = 1:f_count
  dstr = R_datestring2(p_m(k,1:6));
  pr_r = sprintf('%s %5.0f %4.0f %4.0f    %s\n',dstr,p_m(k,7:9),char(p_m(k,10)));
  fwrite(fp,pr_r);
end

fclose(fp);

%---------------------------------------
% ls all MP_*.* files

% ls all Bp_*.* files
f_prefix = 'MP_';
s_cat = pwd;
cd(p_MP);

%writing the ls-command to a file, which can be opened by matlab
ls_out = sprintf('%sls_MP.txt',p_up);
unix_command = sprintf('ls %s*.mat >%s' ,f_prefix, ls_out );
unix(unix_command);

for g =1:f_count
  t = f_events(g,1);
  s_t = fromepoch(t);
  e_t = fromepoch(t);

  fp = fopen(ls_out, 'r');

  % continue until end of file
  while feof(fp) == 0
    file_name = fgetl(fp);

    if c_ss_timestr_within_intervall_MP(file_name,s_t,e_t) == 1
      %load passing_MP,dist_t,dist2MP,p_solarwind
      load(file_name);

      %plotting data
      variable = who('dist_t');
      if isempty(variable)
        dist_t = 0;
      end
      if dist_t(1,1) ~= 0
        t_dMP = time2row(t,dist_t);
        dMP(g) = dist_t(t_dMP,2);
      else
        dMP =0;
      end

    end
  end
  fclose(fp);
end

cd(s_cat);

%---------------------------------------
% ls all Bp_*.* files
f_prefix = 'Bp_';
s_cat = pwd;
cd(path_Bp);

%writing the ls-command to a file, which can be opened by matlab
ls_out = sprintf('%sls_BP.txt',p_up);
unix_command = sprintf('ls %s*.mat >%s' ,f_prefix, ls_out );
unix(unix_command);

for g =1:f_count
  t = f_events(g,1);
  s_t = fromepoch(t);
  e_t = fromepoch(t);

  fp = fopen(ls_out, 'r');

  % continue until end of file
  while feof(fp) == 0
    file_name = fgetl(fp);

    if c_ri_timestr_within_intervall(file_name,s_t,e_t) == 1
      %load B1,B2,B3,B4
      load(file_name);

      mode = file_name(length(file_name)-4);

      m1 = time2row(t-5,B1);
      m2 = time2row(t+5,B1);

      B1_eve = B1(m1:m2,:);
      B2_eve = B2(m1:m2,:);
      B3_eve = B3(m1:m2,:);
      B4_eve = B4(m1:m2,:);

      % fills datagaps, this can be removed if the process get to slow
      B1_eve = c_ri_fill_datagaps(B1_eve);
      B2_eve = c_ri_fill_datagaps(B2_eve);
      B3_eve = c_ri_fill_datagaps(B3_eve);
      B4_eve = c_ri_fill_datagaps(B4_eve);

      [angles, ampl] = c_ri_angles_and_ampl(B1_eve,B2_eve,B3_eve,B4_eve);


      fg = figure;
      %plots B1
      subplot(6,1,1);
      hold on
      irf_plot(B1_eve);
      set(gca,'xlim',t+[-5 5]);
      plot(t,0,'xk')
      ylabel('B1, nT')
      irf_timeaxis
      legend('Bx','By','Bz')
      hold off

      %plots the amplitude of all 4 B
      ampl1 = [B1_eve(:,1), ampl(:,1)];
      ampl2 = [B1_eve(:,1), ampl(:,2)];
      ampl3 = [B1_eve(:,1), ampl(:,3)];
      ampl4 = [B1_eve(:,1), ampl(:,4)];

      subplot(6,1,2)
      hold on
      c_pl_tx(ampl1,ampl2,ampl3,ampl4,2);
      set(gca,'xlim',t+[-5 5]);
      plot(t,0,'xk')
      ylabel('|B|, nT')
      irf_timeaxis
      legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
      hold off

      %plots Bx for all clusters
      subplot(6,1,3)
      hold on
      c_pl_tx(B1_eve,B2_eve,B3_eve,B4_eve,2);
      set(gca,'xlim',t+[-5 5]);
      plot(t,0,'xk')
      ylabel('Bx, nT')
      irf_timeaxis
      legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
      hold off

      %plots By for all clusters
      subplot(6,1,4)
      hold on
      c_pl_tx(B1_eve,B2_eve,B3_eve,B4_eve,3);
      set(gca,'xlim',t+[-5 5]);
      plot(t,0,'xk')
      ylabel('By, nT')
      irf_timeaxis
      legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
      hold off

      %plots Bz for all clusters
      subplot(6,1,5)
      hold on
      c_pl_tx(B1_eve,B2_eve,B3_eve,B4_eve,4);
      set(gca,'xlim',t+[-5 5]);
      plot(t,0,'xk')
      ylabel('Bz, nT')
      irf_timeaxis
      legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
      hold off


      %plotting data
      variable = who('min_ampl');
      if isempty(variable)
        min_ampl = 0;
      end

      variable = who('min_angle');
      if isempty(variable)
        min_angle = 0;
      end

      variable = who('dist2MP');
      if isempty(variable)
        dist2MP = 3;
      end

      variable = who('p_solarwind');
      if isempty(variable)
        p_solarwind = 2;
      end

      if dMP(1,1) ~= 0
        distance = dMP(g);
      end
      variable = who('distance');
      if isempty(variable)
        distance = 0;
      end

      r= time2row(t,angles);

      [amaxim, amaximp] = max(angles(r,2:7));
      [apr1, apr2] =  ind2nr(amaximp);

      distan_str = sprintf('%2.1f',distance);
      subplot(6,1,6)
      %left side
      text(1,4.5,R_datestring2(fromepoch(t)));
      text(1,3,['angle: ' int2str(f_events(g,2)) ' degrees cl ' int2str(apr1) ' to cl ' int2str(apr2) ]);
      text(1,1.5,['min angle: ' int2str(min_angle) ' degrees']);
      text(1,0,['distance to MP: ' distan_str ' RE']);
      text(1,-1.5, ['burst or normal, (b/n) ' mode]);

      %right side
      text(2.5,3,['angles: ' int2str(angles(r,2:7)) ' degrees']);
      text(2.5,4.5,['|B|: ' int2str(ampl(r,:)) ' nT']);
      text(2.5,1.5,['minimum |B|: ' int2str(min_ampl) ' nT']);
      text(2.5,0,['distance to class position as within MP: +-' int2str(dist2MP)  '  RE']);
      text(2.5,-1.5,['solarwind preasure: ' int2str(p_solarwind) ' nT']);
      axis([1,4,1,6])
      axis off

      p_and_f_picture = sprintf('%sF_%s',p_up,c_ri_datestring_file(fromepoch(t)));
      print( fg, '-djpeg', p_and_f_picture);
      disp(['saving: ' p_and_f_picture]);

      close(fg)
    end

  end
  fclose(fp);
end
cd(s_cat);

%------------------------------------------------------

