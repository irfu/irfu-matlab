function a = c_ri_events_into_fig(time_interval,path_eve, file_events, per, p_up, path_Bp,p_MP,path_data)
%function a = c_ri_events_into_fig(time_interval,path_eve, file_events, per, p_up, path_Bp,p_MP,path_data)
%
%Input:
% time_interval - isdat_epoch [start_time end_time]
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
global AV_DEBUG
if isempty(AV_DEBUG), debug=0;else, debug=AV_DEBUG;end


%--------------------- the beginning --------------------------
load([path_eve file_events]);
time_of_events=av_t_lim(time_of_events,time_interval);if debug, whos time_of_events;end
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
  elseif  f_events(f_count,2) < time_of_events(i+1,2)
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
  pr_r = sprintf('%s %5.0f %4.0f %4.0f    %s\n',dstr,p_m(k,7:10));
  fwrite(fp,pr_r);
end

fclose(fp);

%---------------------------------------
% ls all MP_*.* files

MP_files=dir([p_MP 'MP_*.mat']);

for g =1:f_count
  t = f_events(g,1);
  s_t = fromepoch(t);
  e_t = fromepoch(t);

  for i_MP_file=1:size(MP_files,1)
    file_name = MP_files(i_MP_file).name;
    if c_ri_timestr_within_intervall_M(file_name,s_t,e_t) == 1
      %load passing_MP,dist_t,dist2MP,p_solarwind
      load([p_MP file_name]);
      if debug, disp(['MP file: ' p_up file_name]);end
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

  %---------------------------------------
  % read data files - E field, P

  data_files=dir([path_data 'F*T*.mat']);
  for i_data_file=1:size(data_files,1)
    file_name = data_files(i_data_file).name;
    if c_ri_timestr_within_tint(file_name,[t t]) == 1
      load([path_data file_name]);
      if debug, disp(['Data file: ' path_data file_name]);end
    end
  end

  %---------------------------------------
  % ls all Bp_*.* files
  Bp_files=dir([path_Bp 'Bp_*.mat']);
  for i_Bp_file=1:size(Bp_files,1)
    file_name = Bp_files(i_Bp_file).name;

    if c_ri_timestr_within_intervall(file_name,s_t,e_t) == 1
      %load B1,B2,B3,B4
      load([path_Bp file_name]);
      if debug, disp(['Bp file: ' path_Bp file_name]);end
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


      clf;clear h;
      np=8;ip=1;
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(B1_eve,B2_eve,B3_eve,B4_eve,2);
      hold on;
      plot(t,0,'xk');
      ylabel('Bx, nT')
      legend('cl 1', 'cl 2' , 'cl 3', 'cl 4')
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(B1_eve,B2_eve,B3_eve,B4_eve,3);
      hold on;
      plot(t,0,'xk');
      ylabel('By, nT')
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(B1_eve,B2_eve,B3_eve,B4_eve,4);
      hold on;
      plot(t,0,'xk');
      ylabel('Bz, nT')
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(av_abs(B1_eve),av_abs(B2_eve),av_abs(B3_eve),av_abs(B4_eve),5);
      hold on;
      plot(t,0,'xk');
      ylabel('|B|, nT')
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(wE1,wE2,wE3,wE4,3);
      hold on;
      plot(t,0,'xk');
      ylabel('p34, mV/m')
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(wE1,wE2,wE3,wE4,4);
      hold on;
      plot(t,0,'xk');
      ylabel('p12, mV/m')
      h(ip)=irf_subplot(np,1,-ip);ip=ip+1;
      c_pl_tx(P1,P2,P3,P4,3);
      hold on;
      ud=get(gcf,'userdata');
      if isfield(ud,'t_start_epoch'), 	ts=ud.t_start_epoch;else, ts=0;end
      plot(t-ts,0,'xk');
      ylabel('Vps, [cc]')

      irf_zoom(t+[-5 5],'x',h);
      irf_timeaxis(h);
      set(h,'YLimMode','auto');

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

      h(ip)=subplot(np,1,ip);ip=ip+1;
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
      p_and_f_picture = [p_up 'BEP_' c_ri_datestring_file(fromepoch(t)) '.png'];
      print( gcf, '-dpng', p_and_f_picture);
      disp(['saving: ' p_and_f_picture]);
    end
  end
end

%------------------------------------------------------

