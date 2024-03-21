function c_ri_run_class_angle_as_event(s_tm, e_tm, min_angle, min_ampl, path_input, path_output)
%
% c_ri_run_class_angle_as_event(s_t, e_t, min_angle, min_ampl, path_input, path_output)
%
%Input:
% s_t - ex:[2002 03 02 3 0 0]. Can be in matrix form
% e_t - ex:[2002 03 02 3 0 0]. Can be in matrix form
% min_angle -the minimum angle to be called a event
% min_ampl -the minimum ampl to be called a event
% path_input, path_output - end with / ex: /share/ang/.important is that the filename ends with "/"
%                         (Can use default)
%
%Output:
% Saves to file: /path_output/E_F20020302t010303_T20020302t010303.mat
% time_of_events - [time in epoch | angle | max amplitude | mode(in acssi code)]
%
%Descrition of the function:
% Runs "calc_angel_as_event" for all the files the input cataloge and the result is saved "path_output"
%
%Using:
% c_ri_timestr_within_intervall
% class_angle_as_event
%
%Work method:
%
%Error:
% If no events are found, no file is saved
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
[nr_mp_cross,c ] = size(s_tm);
for mpc = 1:nr_mp_cross
  s_t = s_tm(mpc,:);
  e_t = e_tm(mpc,:);

  if nargin == 4
    path_input = [pwd '/'];
    path_output = [pwd '/'];
  end

  %writing the ls-command to a file, which can be opened by matlab
  file_prefix = 'Ap_';
  start_path = pwd;

  cd(path_input);

  ls_out = sprintf('%sls_B_data.txt',path_output);
  unix_command = sprintf('ls %s*.mat >%s' ,file_prefix, ls_out );
  unix(unix_command);

  fp = fopen(ls_out, 'r');
  % continue until end of file

  time_of_events = 0;

  while feof(fp) == 0
    file_name = fgetl(fp);
    f_length = length(file_name);

    %excluding bad filenames
    if strcmp(file_name(1:3),file_prefix) && strcmp(file_name(f_length-3:f_length), '.mat')

      if c_ri_timestr_within_intervall(file_name,s_t,e_t) == 1

        disp(['loading : ' file_name])
        load(file_name);

        new_events = 0;
        mode = file_name(f_length-4);
        new_events = class_angle_as_event(angles, ampl, min_angle, min_ampl,mode);
        time_of_events = add_A2M(time_of_events, new_events);
      end

    end

  end


  [r_t,c_t] = size(time_of_events);

  if r_t == 1 && c_t == 1
    disp('no events found, nothing saved to file')
  else
    f_date = c_ri_datestring_file(s_t);
    t_date = c_ri_datestring_file(e_t);
    p_and_f = sprintf('%sE_F%s_T%s',path_output, f_date, t_date);

    %save section
    disp(['saving: ' p_and_f])
    save(p_and_f, 'time_of_events', 'min_angle','min_ampl')

  end

end % nr of MP_cross

cd(start_path);
